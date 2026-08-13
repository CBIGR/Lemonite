#!/bin/bash

#############################################################################
# Lemonite Pipeline Reproducibility Test (Wang GBM Dataset)
#############################################################################
# Runs the pipeline twice on the Wang GBM test dataset with the same
# random seed and compares key output files to verify reproducibility.
#
# Dataset: /home/borisvdm/Documents/PhD/Test_data_Wang
# Regulators: TFs, Metabolites, Lipidomics, Proteomics
# Settings: 2000 HVGs, 10 Gibbs sampling clustering runs
#
# Author: Vandemoortele Boris, CBIGR lab @ Ghent University
# Date: April 16, 2026
#############################################################################

set -euo pipefail

#############################################################################
# Configuration
#############################################################################

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INPUT_DIR="/home/borisvdm/Documents/PhD/Test_data_Wang"
SINGULARITY_IMAGE="${PIPELINE_DIR}/lemontree-pipeline_v1.0.0.sif"
RESULTS_BASE="${INPUT_DIR}/reproducibility_test"
RUN1_DIR="${RESULTS_BASE}/run1"
RUN2_DIR="${RESULTS_BASE}/run2"
WORK1_DIR="${RESULTS_BASE}/.work_run1"
WORK2_DIR="${RESULTS_BASE}/.work_run2"

# Pipeline parameters
TOP_N_GENES=2000
N_CLUSTERS=10
RANDOM_SEED=42
MAX_CPUS=5
ENRICHMENT_METHOD="ENRICHR"   # Use GSEA for offline-safe execution

# Regulator configuration
REGULATOR_TYPES="TFs:lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:lipidomics.tsv,Proteins:Proteomics_TFs.txt"

#############################################################################
# Helper functions
#############################################################################

print_header() {
    echo ""
    echo "╔══════════════════════════════════════════════════════════════════════════╗"
    printf  "║  %-72s║\n" "$1"
    echo "╚══════════════════════════════════════════════════════════════════════════╝"
    echo ""
}

print_section() {
    echo ""
    echo "──────────────────────────────────────────────────────────────────────────"
    echo "  $1"
    echo "──────────────────────────────────────────────────────────────────────────"
}

compare_files() {
    # Compare two files; print MATCH or DIFF with context
    local label="$1"
    local f1="$2"
    local f2="$3"

    if [ ! -f "$f1" ] && [ ! -f "$f2" ]; then
        echo "    [SKIP]  $label — file absent in both runs"
        return
    elif [ ! -f "$f1" ]; then
        echo "    [SKIP]  $label — missing in run 1: $f1"
        DIFF_COUNT=$((DIFF_COUNT + 1))
        return
    elif [ ! -f "$f2" ]; then
        echo "    [SKIP]  $label — missing in run 2: $f2"
        DIFF_COUNT=$((DIFF_COUNT + 1))
        return
    fi

    local lines1 lines2
    lines1=$(wc -l < "$f1")
    lines2=$(wc -l < "$f2")

    if diff -q "$f1" "$f2" > /dev/null 2>&1; then
        echo "    [MATCH] $label  (${lines1} lines)"
        MATCH_COUNT=$((MATCH_COUNT + 1))
    else
        local changed
        changed=$(diff "$f1" "$f2" | grep -c "^[<>]" || true)
        echo "    [DIFF]  $label  (run1: ${lines1} lines | run2: ${lines2} lines | ${changed} changed lines)"
        DIFF_COUNT=$((DIFF_COUNT + 1))
    fi
}

compare_numeric_file() {
    # Compare two TSV/CSV files ignoring floating-point noise (threshold 1e-6)
    local label="$1"
    local f1="$2"
    local f2="$3"

    if [ ! -f "$f1" ] || [ ! -f "$f2" ]; then
        compare_files "$label" "$f1" "$f2"
        return
    fi

    local max_diff
    max_diff=$(python3 - "$f1" "$f2" <<'PYEOF'
import sys, csv, math
f1, f2 = sys.argv[1], sys.argv[2]
sep = '\t' if f1.endswith(('.tsv', '.txt')) else ','
max_d = 0.0
with open(f1) as a, open(f2) as b:
    for row_a, row_b in zip(csv.reader(a, delimiter=sep), csv.reader(b, delimiter=sep)):
        for ca, cb in zip(row_a, row_b):
            try:
                max_d = max(max_d, abs(float(ca) - float(cb)))
            except ValueError:
                if ca != cb:
                    max_d = float('inf')
print(f"{max_d:.2e}")
PYEOF
)
    if [ "$max_diff" = "0.00e+00" ]; then
        echo "    [MATCH] $label  (max numeric diff = 0)"
        MATCH_COUNT=$((MATCH_COUNT + 1))
    elif python3 -c "import sys; sys.exit(0 if float('$max_diff') < 1e-6 else 1)" 2>/dev/null; then
        echo "    [~MATCH] $label  (max numeric diff = ${max_diff}, within 1e-6)"
        MATCH_COUNT=$((MATCH_COUNT + 1))
    else
        echo "    [DIFF]  $label  (max numeric diff = ${max_diff})"
        DIFF_COUNT=$((DIFF_COUNT + 1))
    fi
}

#############################################################################
# Step 0: Preflight checks
#############################################################################

print_header "LEMONITE REPRODUCIBILITY TEST — Wang GBM Dataset"

echo "Configuration:"
echo "   Input data:       $INPUT_DIR"
echo "   Pipeline:         $PIPELINE_DIR"
echo "   Container:        $SINGULARITY_IMAGE"
echo "   Top N genes:      $TOP_N_GENES"
echo "   Clustering runs:  $N_CLUSTERS"
echo "   Random seed:      $RANDOM_SEED"
echo "   Enrichment:       $ENRICHMENT_METHOD"
echo "   Output base:      $RESULTS_BASE"
echo ""

print_section "Step 0: Preflight checks"

if ! command -v nextflow &> /dev/null; then
    echo "ERROR: nextflow not found in PATH. Install from https://www.nextflow.io/"
    exit 1
fi
echo "   Nextflow: $(nextflow -version 2>&1 | head -1 | tr -d '\n')"

if ! command -v singularity &> /dev/null; then
    echo "ERROR: singularity not found in PATH."
    exit 1
fi
echo "   Singularity: $(singularity --version)"

if [ ! -f "$SINGULARITY_IMAGE" ]; then
    echo "ERROR: Singularity image not found: $SINGULARITY_IMAGE"
    echo "   Run ./build-singularity.sh first."
    exit 1
fi
echo "   Container: $(basename "$SINGULARITY_IMAGE") ($(du -h "$SINGULARITY_IMAGE" | cut -f1))"

if [ ! -d "$INPUT_DIR/data" ]; then
    echo "ERROR: Data directory not found: $INPUT_DIR/data"
    exit 1
fi

for f in Counts.tsv Metabolomics.txt lipidomics.tsv Proteomics_TFs.txt lovering_TF_list.txt Metadata.txt; do
    if [ ! -f "$INPUT_DIR/data/$f" ]; then
        echo "ERROR: Required data file missing: $INPUT_DIR/data/$f"
        exit 1
    fi
done
echo "   All required data files present"

mkdir -p "$RUN1_DIR" "$RUN2_DIR" "$WORK1_DIR" "$WORK2_DIR"
echo "   Output directories created"

#############################################################################
# Shared nextflow arguments
#############################################################################

NF_ARGS=(
    "${PIPELINE_DIR}/main.nf"
    --input_dir      "$INPUT_DIR"
    --top_n_genes    "$TOP_N_GENES"
    --n_clusters     "$N_CLUSTERS"
    --random_seed    "$RANDOM_SEED"
    --regulator_types "$REGULATOR_TYPES"
    --enrichment_method "$ENRICHMENT_METHOD"
    --max_cpus       "$MAX_CPUS"
    -profile         singularity
)

#############################################################################
# Step 1: Run 1
#############################################################################

print_section "Step 1: Pipeline Run 1  (seed=${RANDOM_SEED})"
echo "   Output:    $RUN1_DIR"
echo "   Work dir:  $WORK1_DIR"
echo ""

START1=$SECONDS
nextflow run "${NF_ARGS[@]}" \
    --output_dir "$RUN1_DIR" \
    -w "$WORK1_DIR" \
    -with-report "${RUN1_DIR}/execution_report.html" \
    -with-trace  "${RUN1_DIR}/execution_trace.txt"
DURATION1=$((SECONDS - START1))

echo ""
echo "   Run 1 completed in $(( DURATION1 / 60 ))m $(( DURATION1 % 60 ))s"

#############################################################################
# Step 2: Run 2
#############################################################################

print_section "Step 2: Pipeline Run 2  (seed=${RANDOM_SEED}, independent)"
echo "   Output:    $RUN2_DIR"
echo "   Work dir:  $WORK2_DIR"
echo ""

START2=$SECONDS
nextflow run "${NF_ARGS[@]}" \
    --output_dir "$RUN2_DIR" \
    -w "$WORK2_DIR" \
    -with-report "${RUN2_DIR}/execution_report.html" \
    -with-trace  "${RUN2_DIR}/execution_trace.txt"
DURATION2=$((SECONDS - START2))

echo ""
echo "   Run 2 completed in $(( DURATION2 / 60 ))m $(( DURATION2 % 60 ))s"

#############################################################################
# Step 3: Compare results
#############################################################################

print_section "Step 3: Comparing results"

MATCH_COUNT=0
DIFF_COUNT=0

# Locate run subdirectories (named by HVG/coherence params)
RUN1_SUBDIR=$(find "$RUN1_DIR" -maxdepth 1 -mindepth 1 -type d | head -1)
RUN2_SUBDIR=$(find "$RUN2_DIR" -maxdepth 1 -mindepth 1 -type d | head -1)

if [ -z "$RUN1_SUBDIR" ] || [ -z "$RUN2_SUBDIR" ]; then
    echo "   WARNING: Could not locate run subdirectories for comparison."
    echo "   Run 1 output: $RUN1_DIR"
    echo "   Run 2 output: $RUN2_DIR"
    exit 1
fi

echo ""
echo "   Run 1 subdir: $(basename "$RUN1_SUBDIR")"
echo "   Run 2 subdir: $(basename "$RUN2_SUBDIR")"
echo ""

# ── Preprocessing outputs ─────────────────────────────────────────────────
echo "  Preprocessing:"
for fname in \
    "LemonTree/Preprocessing/LemonPreprocessed_expression.txt" \
    "LemonTree/Preprocessing/LemonPreprocessed_complete.txt" \
    "LemonTree/Preprocessing/DESeq_groups.txt"; do
    compare_numeric_file "$fname" \
        "${RUN1_SUBDIR}/${fname}" \
        "${RUN2_SUBDIR}/${fname}"
done

# ── Clustering (module assignments) ───────────────────────────────────────
echo ""
echo "  Clustering:"
for fname in \
    "LemonTree/ModuleViewer_files/clusters_list.txt" \
    "LemonTree/ModuleViewer_files/tight_clusters_gene_sets.gmt"; do
    compare_files "$fname" \
        "${RUN1_SUBDIR}/${fname}" \
        "${RUN2_SUBDIR}/${fname}"
done

# Compare any module coherence score file
for f in "${RUN1_SUBDIR}"/LemonTree/Lemon_out/Module_coherence_scores*.txt; do
    [ -f "$f" ] || continue
    bname="LemonTree/Lemon_out/$(basename "$f")"
    compare_numeric_file "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

# ── Networks ──────────────────────────────────────────────────────────────
echo ""
echo "  Regulatory Networks:"
for f in "${RUN1_SUBDIR}"/Networks/LemonNetwork_*.txt; do
    [ -f "$f" ] || continue
    bname="Networks/$(basename "$f")"
    compare_files "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

for f in "${RUN1_SUBDIR}"/Networks/*.tsv; do
    [ -f "$f" ] || continue
    bname="Networks/$(basename "$f")"
    compare_files "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

# ── Regulator selection lists ─────────────────────────────────────────────
echo ""
echo "  Regulator Selection Lists:"
for f in "${RUN1_SUBDIR}"/LemonTree/Preprocessing/*.selected_regs_list.txt; do
    [ -f "$f" ] || continue
    bname="LemonTree/Preprocessing/$(basename "$f")"
    compare_files "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

# ── Enrichment ────────────────────────────────────────────────────────────
echo ""
echo "  Enrichment:"
for f in "${RUN1_SUBDIR}"/enrichment/*.csv; do
    [ -f "$f" ] || continue
    bname="enrichment/$(basename "$f")"
    compare_numeric_file "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

# ── PKN Evaluation ────────────────────────────────────────────────────────
echo ""
echo "  PKN Evaluation:"
for f in "${RUN1_SUBDIR}"/Evaluation/*.txt "${RUN1_SUBDIR}"/Evaluation/*.tsv; do
    [ -f "$f" ] || continue
    bname="Evaluation/$(basename "$f")"
    compare_numeric_file "$bname" "$f" "${RUN2_SUBDIR}/${bname}"
done

#############################################################################
# Step 4: Summary
#############################################################################

print_section "Reproducibility Summary"

TOTAL=$((MATCH_COUNT + DIFF_COUNT))
echo "   Compared files:  $TOTAL"
echo "   Matching:        $MATCH_COUNT"
echo "   Different:       $DIFF_COUNT"
echo ""

if [ "$DIFF_COUNT" -eq 0 ] && [ "$TOTAL" -gt 0 ]; then
    echo "   RESULT: ✅  FULLY REPRODUCIBLE — all $MATCH_COUNT compared files are identical"
    echo "           The --random_seed fix is working correctly."
elif [ "$DIFF_COUNT" -eq 0 ] && [ "$TOTAL" -eq 0 ]; then
    echo "   RESULT: ⚠️   No files were found to compare. Check output directories:"
    echo "           $RUN1_SUBDIR"
    echo "           $RUN2_SUBDIR"
else
    echo "   RESULT: ❌  NOT FULLY REPRODUCIBLE — $DIFF_COUNT/$TOTAL files differ"
    echo "           Inspect diff with:"
    echo "           diff ${RUN1_SUBDIR}/Networks/ ${RUN2_SUBDIR}/Networks/"
fi

echo ""
echo "   Run 1 output: $RUN1_SUBDIR"
echo "   Run 2 output: $RUN2_SUBDIR"
echo ""
echo "   Run 1 duration: $(( DURATION1 / 60 ))m $(( DURATION1 % 60 ))s"
echo "   Run 2 duration: $(( DURATION2 / 60 ))m $(( DURATION2 % 60 ))s"
echo ""
