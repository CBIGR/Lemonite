#!/usr/bin/env bash
set -uo pipefail

# Run the nextflow post-regulator-assignment steps for one phospho variant,
# reproducing what main.nf does AFTER POST_CLUSTERING, over the phospho-augmented
# regulator set (existing TFs/Metabolites/Lipids + Phospho). Steps:
#   1. NETWORK_GENERATION      lemontree_to_network.py
#   2. SUBNETWORK_GRAPHS       create_subnetwork_graphs.py
#   3. PKN_EVALUATION          evaluate_against_PKN.py
#   4. MODULE_VIEWER_HEATMAPS  module_viewer.py
#   5. ENRICHMENT_ANALYSIS     enrichment_analysis.R
#   6. MODULE_OVERVIEW         module_overview_interactive.py   (MegaGO NOT recomputed:
#                              --n_clusters 1 disables it; the base run's module-level
#                              MegaGO/rrvgo labels are injected afterwards, since the
#                              modules are identical to the initial run)
#   7. SUMMARY_REPORT          generate_summary_report.py
#
# Everything runs inside the pipeline Singularity image so the environment matches
# the real pipeline. All scripts are read from nextflow/scripts/ (host copy).
#
# Usage: run_downstream.sh <all_phosphosites|top2000_variable> [first_step] [last_step]
#   step numbers 1..7; default runs 1..7.

VARIANT=${1:?"usage: run_downstream.sh <all_phosphosites|top2000_variable> [first_step] [last_step]"}
FIRST=${2:-1}
LAST=${3:-7}

PROJ=/home/borisvdm/repo/LemonIte/nextflow
SIF=$PROJ/lemontree-pipeline_v1.0.0.sif
SCR=$PROJ/scripts
RESULTS_ROOT=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum_phospho
BASE=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum
RUN=$RESULTS_ROOT/$VARIANT/run
LOGDIR=$RESULTS_ROOT/$VARIANT/downstream_logs
mkdir -p "$LOGDIR"

# ---- pipeline parameters (match the base GBM run) ---------------------------
RUN_ID="phospho_${VARIANT}"
REG_TYPES="Lovering:lovering_TFs.txt,Metabolites:metabolites.txt,Lipids:lipids.txt,Phospho:phospho.txt"
# The base GBM run kept 46 modules, which corresponds to coherence_threshold 0.5
# (the config default is 0.6 → 36 modules). Match the initial analysis with 0.5.
COHERENCE=0.5
SEL_METHOD=percentage
TOP_PCT=2.0
FOLD_CUTOFF=2.0
ORGANISM=human
GROUP_COL=diagnosis
ENRICH_METHOD=EnrichR
PKN=$PROJ/PKN/Lemonite_PKN.tsv
MET_INTERACT=$PROJ/PKN/metabolite_gene_PKN.tsv
# Subnetwork figures use a viz-only PKN = main PKN + standalone phospho edges, so phospho
# regulator->target edges render (orange). PKN eval keeps using the main PKN (phospho eval
# is handled separately by evaluate_phospho_against_pkn.py). Falls back to main PKN if absent.
PHOSPHO_DIR=/home/borisvdm/repo/LemonIte/Wang_GBM/Lemonite_phospho
PKN_VIZ=$PHOSPHO_DIR/Lemonite_PKN_plus_phospho.tsv
# The viz PKN is a large, gitignored, regenerable artifact — rebuild it if missing.
if [ ! -f "$PKN_VIZ" ] && [ -f "$PHOSPHO_DIR/phospho_pkn.tsv" ]; then
    python3 "$PHOSPHO_DIR/build_viz_pkn.py" || true
fi
[ -f "$PKN_VIZ" ] || PKN_VIZ=$PKN

# run scripts inside the container, bound to the run dir + project
sing() {
    singularity exec --cleanenv --bind "$RUN:$RUN" --bind "$PROJ:$PROJ" --bind "$BASE:$BASE" "$SIF" "$@"
}

run_step() { # $1=num $2=name ; body reads from stdin
    local n=$1 name=$2
    if [ "$n" -lt "$FIRST" ] || [ "$n" -gt "$LAST" ]; then return 0; fi
    echo "======================================================================"
    echo " STEP $n: $name   (variant=$VARIANT)"
    echo "======================================================================"
}

cd "$RUN"

# =============================================================================
run_step 1 "NETWORK_GENERATION"
if [ "$FIRST" -le 1 ] && [ "$LAST" -ge 1 ]; then
    sing python3 "$SCR/lemontree_to_network.py" \
        --input_dir "$RUN" --output_dir "$RUN" --run_id "$RUN_ID" \
        --coherence_threshold "$COHERENCE" \
        --regulator_selection_method "$SEL_METHOD" \
        --top_n_percent_regulators "$TOP_PCT" \
        --regulator_fold_cutoff "$FOLD_CUTOFF" \
        --regulator_types "$REG_TYPES" \
        > "$LOGDIR/01_network.log" 2>&1
    echo "  exit=$? -> $LOGDIR/01_network.log"
fi

# =============================================================================
run_step 2 "SUBNETWORK_GRAPHS"
if [ "$FIRST" -le 2 ] && [ "$LAST" -ge 2 ]; then
    # build Prefix:selected_regs_list.txt list from ModuleViewer_files
    REGFILES=""
    for pfx in Lovering Metabolites Lipids Phospho; do
        f=$(find "$RUN/ModuleViewer_files" -name "${pfx}.selected_regs_list.txt" | head -1)
        [ -z "$f" ] && f=$(find "$RUN/Networks" -name "${pfx}.selected_regs_list.txt" | head -1)
        [ -n "$f" ] && REGFILES="${REGFILES:+$REGFILES,}${pfx}:${f}"
    done
    METARG=""; [ -f "$RUN/Preprocessing/name_map.csv" ] && METARG="--metabolite_mapping $RUN/Preprocessing/name_map.csv"
    NMARG="";  [ -f "$RUN/Preprocessing/name_mapping.tsv" ] && NMARG="--name_mapping $RUN/Preprocessing/name_mapping.tsv"
    PHOSPHO_DIR=/home/borisvdm/repo/LemonIte/Wang_GBM/Lemonite_phospho
    singularity exec --cleanenv --bind "$RUN:$RUN" --bind "$PROJ:$PROJ" --bind "$BASE:$BASE" \
        --bind "$PHOSPHO_DIR:$PHOSPHO_DIR" "$SIF" \
        python3 "$SCR/create_subnetwork_graphs.py" \
        --regulator_files "$REGFILES" \
        --clusters "$RUN/ModuleViewer_files/clusters_list.txt" \
        --pkn "$PKN_VIZ" $METARG $NMARG \
        --output_dir "$RUN/Networks/subnetworks" \
        > "$LOGDIR/02_subnetworks.log" 2>&1
    echo "  exit=$? -> $LOGDIR/02_subnetworks.log"
fi

# =============================================================================
# Flatten files into the run root the way the nextflow processes stage them:
# evaluate_against_PKN.py chdir's to workdir and globs LemonNetwork_*.txt / *2targets_*
# in CWD; enrichment_analysis.R reads input_dir/LemonPreprocessed_expression.txt.
flatten_run_root() {
    # Clear stale root-level network files first so an old module-count network from a
    # previous run can't win the glob in evaluate_against_PKN.py.
    rm -f "$RUN"/LemonNetwork_*.txt "$RUN"/*2targets_*.txt 2>/dev/null || true
    cp -f "$RUN"/Networks/LemonNetwork_*.txt        "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Networks/*2targets_*.txt           "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Networks/specific_modules.txt      "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Networks/Module_coherence_scores.txt "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/ModuleViewer_files/*.selected_regs_list.txt        "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/ModuleViewer_files/clusters_list.txt              "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Preprocessing/LemonPreprocessed_expression.txt    "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Preprocessing/LemonPreprocessed_complete.txt      "$RUN"/ 2>/dev/null || true
    cp -f "$RUN"/Preprocessing/DESeq_groups.txt                    "$RUN"/ 2>/dev/null || true
}
if [ "$FIRST" -le 5 ]; then flatten_run_root; fi

run_step 3 "PKN_EVALUATION"
if [ "$FIRST" -le 3 ] && [ "$LAST" -ge 3 ]; then
    NETFILE=$(find "$RUN/Networks" -name "LemonNetwork_*.txt" | head -1)
    CLUSTERF=$(find "$RUN" -maxdepth 1 -name "clusters_list.txt" | head -1)
    [ -z "$CLUSTERF" ] && CLUSTERF=$(find "$RUN" -name "clusters_list.txt" | head -1)
    REGFILES=""
    for pfx in Lovering Metabolites Lipids Phospho; do
        f=$(find "$RUN/Networks" "$RUN/ModuleViewer_files" -name "${pfx}.selected_regs_list.txt" 2>/dev/null | head -1)
        [ -n "$f" ] && REGFILES="${REGFILES:+$REGFILES,}${pfx}:${f}"
    done
    METMAP="$RUN/Preprocessing/name_map.csv"; [ -f "$METMAP" ] || METMAP=/tmp/name_map.csv
    sing python3 "$SCR/evaluate_against_PKN.py" \
        --workdir "$RUN" \
        --network_file "$NETFILE" \
        --pkn_file "$PKN" \
        --metabolite_interactions_file "$MET_INTERACT" \
        --metabolite_mapping_file "$METMAP" \
        --regulator_files "$REGFILES" \
        --cluster_file "$CLUSTERF" \
        --cores 4 \
        > "$LOGDIR/03_pkn_eval.log" 2>&1
    echo "  exit=$? -> $LOGDIR/03_pkn_eval.log"
    # Propagate the .mvf interaction files PKN eval produced into ModuleViewer_files/
    # (heatmaps/overview read them from there).
    for src in "$RUN/results/ModuleViewer_files" "$RUN/ModuleViewer_files"; do
        [ -d "$src" ] && cp -f "$src"/*.mvf "$RUN/ModuleViewer_files/" 2>/dev/null || true
    done
fi

# =============================================================================
run_step 4 "MODULE_VIEWER_HEATMAPS"
if [ "$FIRST" -le 4 ] && [ "$LAST" -ge 4 ]; then
    REGFILES=""
    for pfx in Lovering Metabolites Lipids Phospho; do
        f=$(find "$RUN/ModuleViewer_files" -name "${pfx}.selected_regs_list.txt" | head -1)
        [ -n "$f" ] && REGFILES="${REGFILES:+$REGFILES,}${pfx}:${f}"
    done
    NMARG="$RUN/Preprocessing/name_mapping.tsv"
    sing env PYTHONNOUSERSITE=1 python3 "$SCR/module_viewer.py" \
        --input_dir "$RUN" --output_dir "$RUN/heatmaps" \
        --regulator_files "$REGFILES" \
        --regulator_types "$REG_TYPES" \
        --dpi 300 --show_regulator_scores \
        --annotation_types "$GROUP_COL" \
        --name_mapping "$NMARG" \
        > "$LOGDIR/04_heatmaps.log" 2>&1
    echo "  exit=$? -> $LOGDIR/04_heatmaps.log"
fi

# =============================================================================
run_step 5 "ENRICHMENT_ANALYSIS"
if [ "$FIRST" -le 5 ] && [ "$LAST" -ge 5 ]; then
    sing Rscript "$SCR/enrichment_analysis.R" \
        --input_dir "$RUN" --output_dir "$RUN" \
        --analysis_method "$ENRICH_METHOD" \
        --top_n_percent_regulators "$TOP_PCT" \
        --coherence_threshold "$COHERENCE" \
        --n_threads 4 \
        --regulator_types "$REG_TYPES" \
        --organism "$ORGANISM" \
        > "$LOGDIR/05_enrichment.log" 2>&1
    echo "  exit=$? -> $LOGDIR/05_enrichment.log"
fi

# =============================================================================
run_step 6 "MODULE_OVERVIEW (MegaGO reused from base run)"
if [ "$FIRST" -le 6 ] && [ "$LAST" -ge 6 ]; then
    REGFILES=""; REGSCORES=""
    for pfx in Lovering Metabolites Lipids Phospho; do
        f=$(find "$RUN/ModuleViewer_files" -name "${pfx}.selected_regs_list.txt" | head -1)
        s=$(find "$RUN/ModuleViewer_files" -name "${pfx}.selected_regulators_scores.txt" | head -1)
        [ -n "$f" ] && REGFILES="${REGFILES:+$REGFILES,}${pfx}:${f}"
        [ -n "$s" ] && REGSCORES="${REGSCORES:+$REGSCORES,}${pfx}:${s}"
    done
    EXPR="$RUN/Preprocessing/LemonPreprocessed_expression.txt"
    METAF="$RUN/Preprocessing/DESeq_groups.txt"
    # --n_clusters 1 disables MegaGO recomputation (functional clustering off);
    # base-run module-level MegaGO/rrvgo labels are merged in afterwards.
    sing python3 "$SCR/module_overview_interactive.py" \
        --input_dir "$RUN" --output_dir "$RUN" \
        --enrichment_method "$ENRICH_METHOD" \
        --coherence_threshold "$COHERENCE" \
        --group_column "$GROUP_COL" \
        --regulator_files "$REGFILES" \
        --regulator_score_files "$REGSCORES" \
        --n_clusters 1 \
        --organism "$ORGANISM" \
        --run_id "$RUN_ID" \
        --pkn_file "$PKN" \
        --metabolite_mapping "$RUN/Preprocessing/name_map.csv" \
        --name_mapping "$RUN/Preprocessing/name_mapping.tsv" \
        --prioritize_by_expression \
        --expression_file "$EXPR" \
        --metadata_file "$METAF" \
        > "$LOGDIR/06_overview.log" 2>&1
    echo "  overview exit=$? -> $LOGDIR/06_overview.log"

    # Inject base-run MegaGO/rrvgo module clustering (modules are identical).
    # Run inside the container (has pandas); bind the phospho script dir.
    INJECT=/home/borisvdm/repo/LemonIte/Wang_GBM/Lemonite_phospho/inject_base_megago.py
    OVERVIEW_CSV="$RUN/Module_Overview/Module_Overview.csv"
    singularity exec --cleanenv --bind "$RUN:$RUN" --bind "$BASE:$BASE" \
        --bind "$(dirname "$INJECT"):$(dirname "$INJECT")" "$SIF" \
        env PYTHONNOUSERSITE=1 python3 "$INJECT" \
            --overview "$OVERVIEW_CSV" \
            --base-overview "$BASE/Module_Overview_Complete.csv" \
            >> "$LOGDIR/06_overview.log" 2>&1
    echo "  megago-inject exit=$? (see 06_overview.log tail)"
fi

# =============================================================================
run_step 7 "SUMMARY_REPORT"
if [ "$FIRST" -le 7 ] && [ "$LAST" -ge 7 ]; then
    sing python3 "$SCR/generate_summary_report.py" \
        --input_dir "$RUN" --output_dir "$RUN" --run_id . \
        --regulator_types "$REG_TYPES" \
        --organism "$ORGANISM" \
        > "$LOGDIR/07_summary.log" 2>&1
    echo "  exit=$? -> $LOGDIR/07_summary.log"
fi

echo "DONE variant=$VARIANT steps $FIRST..$LAST"
