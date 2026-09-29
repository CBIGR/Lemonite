#!/usr/bin/env bash
set -euo pipefail
# Run LemonTree clustering + regulator assignment for the phospho-clustering variant, using
# the pipeline's own JAR inside the pipeline .sif. We invoke LemonTree directly (not the
# nextflow main.nf) because main.nf re-runs its own preprocessing from raw data/ and has no
# phospho-as-primary preprocessing_type — see PIPELINE_ISSUES.md. This reproduces the same
# LemonTree steps the pipeline would run (ganesh -> tight_clusters -> regulators per layer).
#
# Usage: run_lemontree.sh <variant-dir> [n_clusters] [min_weight]
VARIANT_DIR=$1
NCLUST=${2:-10}          # Gibbs samplers (pipeline default 100; 10 for a tractable first pass)
MINW=${3:-0.25}

SIF=/home/borisvdm/repo/LemonIte/nextflow/lemontree-pipeline_v1.0.0.sif
JAR=/opt/lemontree/lemontree_v3.1.1.jar
REPO=/home/borisvdm/repo/LemonIte
PREP="$VARIANT_DIR/Preprocessing"
OUT="$VARIANT_DIR/Lemon_out"
mkdir -p "$OUT/Lemon_results"

sing() { singularity exec --cleanenv --bind "$VARIANT_DIR:$VARIANT_DIR" --bind "$REPO:$REPO" "$SIF" "$@"; }

EXPR="$PREP/LemonPreprocessed_expression.txt"
COMPLETE="$PREP/LemonPreprocessed_complete.txt"

MAXPAR=${4:-10}          # max concurrent Gibbs samplers (box has 12 cores)
echo "=== [1/3] Clustering: $NCLUST Gibbs samplers (<=$MAXPAR parallel) on $(wc -l < "$EXPR") rows ==="
for i in $(seq 1 "$NCLUST"); do
    sing java -cp "$JAR" lemontree.modulenetwork.RunCli \
        -task ganesh -data_file "$EXPR" -output_file "$OUT/Lemon_results/cluster_$i" \
        > "$OUT/cluster_$i.log" 2>&1 &
    # throttle to MAXPAR concurrent jobs
    while (( $(jobs -rp | wc -l) >= MAXPAR )); do sleep 2; done
done
wait
echo "clusters done: $(ls "$OUT"/Lemon_results/cluster_* 2>/dev/null | wc -l)"

echo "=== [2/3] Tight clusters (consensus) ==="
rm -f "$OUT/clusterfile"
for cf in "$OUT"/Lemon_results/cluster_*; do echo "$cf" >> "$OUT/clusterfile"; done
sing java -cp "$JAR" lemontree.modulenetwork.RunCli \
    -task tight_clusters -data_file "$EXPR" -cluster_file "$OUT/clusterfile" \
    -output_file "$OUT/tight_clusters.txt" -node_clustering false \
    -min_weight "$MINW" -min_clust_size 10 -min_clust_score 2 \
    > "$OUT/tight.log" 2>&1
echo "tight_clusters: $(grep -c . "$OUT/tight_clusters.txt" 2>/dev/null || echo 0) lines"

echo "=== [3/3] Regulator assignment per layer (parallel) ==="
# The 5 layers are independent (same tight_clusters, different reg_file/output) and each
# LemonTree JVM is single-threaded, so run them concurrently. Prefix:reg_list_file.
for pair in TFs:tfs Metabolites:metabolites Proteins:proteins Lipids:lipids KinaseActivity:kinaseactivity; do
    prefix="${pair%%:*}"; rf="$PREP/${pair##*:}.txt"
    [ -s "$rf" ] || { echo "  skip $prefix (no reg list)"; continue; }
    echo "  launching $prefix ($(wc -l < "$rf") regulators)..."
    sing java -cp "$JAR" lemontree.modulenetwork.RunCli \
        -task regulators -data_file "$COMPLETE" -reg_file "$rf" \
        -cluster_file "$OUT/tight_clusters.txt" -output_file "$OUT/$prefix" \
        > "$OUT/reg_$prefix.log" 2>&1 &
done
wait
for prefix in TFs Metabolites Proteins Lipids KinaseActivity; do
    [ -s "$OUT/$prefix.topreg.txt" ] && echo "  OK $prefix -> $prefix.topreg.txt" || echo "  [FAIL] $prefix (see reg_$prefix.log)"
done
echo "=== DONE $VARIANT_DIR ==="
