#!/usr/bin/env bash
set -euo pipefail

# Run ONLY the LemonTree regulator-assignment step for the phosphoproteomics layer,
# against the EXISTING modules (tight_clusters.txt) of the Wang GBM Lemonite run.
# No re-clustering, no re-generation of tight clusters — modules are kept as-is.
#
# LemonTree is executed inside the pipeline's Singularity image so the exact same
# jar (/opt/lemontree/lemontree_v3.1.1.jar) and Java runtime are used as in the
# nextflow pipeline.
#
# Usage:
#   run_phospho_regulators.sh <variant_Lemon_out_dir>
# where <variant_Lemon_out_dir> contains:
#   tight_clusters.txt                     (copied from the base run; the modules to keep)
#   LemonPreprocessed_complete_phospho.txt (base complete matrix + phospho rows)
#   phospho.txt                            (regulator list)
# Produces in the same dir:  Phospho, Phospho.topreg.txt, Phospho.xml.gz, ...

OUT_DIR=${1:?"usage: run_phospho_regulators.sh <variant_Lemon_out_dir>"}
SIF=${SIF:-/home/borisvdm/repo/LemonIte/nextflow/lemontree-pipeline_v1.0.0.sif}
JAR=${JAR:-/opt/lemontree/lemontree_v3.1.1.jar}
PREFIX=${PREFIX:-Phospho}

OUT_DIR=$(readlink -f "$OUT_DIR")

for f in tight_clusters.txt LemonPreprocessed_complete_phospho.txt phospho.txt; do
    if [ ! -s "$OUT_DIR/$f" ]; then
        echo "Error: required input '$OUT_DIR/$f' is missing or empty" >&2
        exit 1
    fi
done

echo "====================================================="
echo "Phospho regulator assignment (existing modules kept)"
echo "  Singularity image : $SIF"
echo "  LemonTree jar     : $JAR"
echo "  Working dir       : $OUT_DIR"
echo "  Output prefix     : $PREFIX"
echo "  N phosphosites    : $(wc -l < "$OUT_DIR/phospho.txt")"
echo "  N modules (tight) : $(grep -c . "$OUT_DIR/tight_clusters.txt" || true)"
echo "====================================================="

singularity exec --cleanenv \
    --bind "$OUT_DIR:$OUT_DIR" \
    "$SIF" \
    java -cp "$JAR" lemontree.modulenetwork.RunCli \
        -task regulators \
        -data_file "$OUT_DIR/LemonPreprocessed_complete_phospho.txt" \
        -reg_file  "$OUT_DIR/phospho.txt" \
        -cluster_file "$OUT_DIR/tight_clusters.txt" \
        -output_file "$OUT_DIR/$PREFIX"

echo ""
echo "[ok] Phospho regulator assignment complete. Outputs:"
ls -1 "$OUT_DIR/${PREFIX}"* 2>/dev/null || echo "  (no ${PREFIX}* files produced — check the log above)"
