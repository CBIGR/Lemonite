#!/usr/bin/env bash
set -euo pipefail

# Stage a self-contained run directory for one phospho variant, combining the
# EXISTING GBM regulator outputs (TFs/Metabolites/Lipids, kept unchanged) with the
# newly-computed Phospho regulator outputs, so the nextflow downstream scripts can be
# run over the phospho-augmented network without re-clustering.
#
# Usage: stage_downstream.sh <variant>    (variant = all_phosphosites | top2000_variable)
#
# Produces  <RESULTS_ROOT>/<variant>/run/{Lemon_out,Preprocessing,data}  ready for the
# network / PKN / heatmap / enrichment / overview / summary scripts.

VARIANT=${1:?"usage: stage_downstream.sh <all_phosphosites|top2000_variable>"}

BASE=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum
RESULTS_ROOT=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum_phospho
DATA=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/data

VDIR="$RESULTS_ROOT/$VARIANT"
RUN="$VDIR/run"
PHOS="$VDIR/Lemon_out"     # where prepare_phospho_regulator.py + LemonTree wrote Phospho.*

mkdir -p "$RUN/Lemon_out" "$RUN/Preprocessing"

# --- Lemon_out: base regulator outputs (unchanged) + Phospho -----------------
# Base regulators that entered the original network: Lovering (TFs), Metabolites, Lipids.
for pfx in Lovering Metabolites Lipids; do
    for suf in allreg.txt randomreg.txt topreg.txt; do
        cp -f "$BASE/Lemon_out/${pfx}.${suf}" "$RUN/Lemon_out/" 2>/dev/null || true
    done
done
# Phospho outputs from this variant
for suf in allreg.txt randomreg.txt topreg.txt; do
    cp -f "$PHOS/Phospho.${suf}" "$RUN/Lemon_out/"
done
# Modules (kept as-is) + helper lists
cp -f "$BASE/Lemon_out/tight_clusters.txt"  "$RUN/Lemon_out/"
cp -f "$BASE/Lemon_out/clusters_list.txt"   "$RUN/Lemon_out/" 2>/dev/null || true
cp -f "$BASE/Lemon_out/DESeq_groups.txt"    "$RUN/Lemon_out/" 2>/dev/null || true
cp -f "$BASE/Lemon_out/ensemble_mapping.txt" "$RUN/Lemon_out/" 2>/dev/null || true

# --- Preprocessing: expression + complete (+phospho rows) + mappings ---------
cp -f "$BASE/Preprocessing/LemonPreprocessed_expression.txt" "$RUN/Preprocessing/"
cp -f "$BASE/Preprocessing/DESeq_groups.txt"                 "$RUN/Preprocessing/"
cp -f "$BASE/Preprocessing/name_map.csv"                     "$RUN/Preprocessing/" 2>/dev/null || true
cp -f "$BASE/Preprocessing/name_mapping.tsv"                 "$RUN/Preprocessing/" 2>/dev/null || true
cp -f "$BASE/Preprocessing/ensemble_mapping.txt"            "$RUN/Preprocessing/" 2>/dev/null || true
# Use the phospho-augmented complete matrix so phospho regulator names resolve
cp -f "$PHOS/LemonPreprocessed_complete_phospho.txt" "$RUN/Preprocessing/LemonPreprocessed_complete.txt"

# --- data: symlink case-study data dir (regulator lists, metabolite files) ---
rm -rf "$RUN/data"
ln -sfn "$DATA" "$RUN/data"

echo "[ok] staged $RUN"
echo "  Lemon_out regulator files:"
ls -1 "$RUN/Lemon_out/"*.allreg.txt | xargs -n1 basename | sed 's/^/    /'
