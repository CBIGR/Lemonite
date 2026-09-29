#!/usr/bin/env bash
# Resume the phospho-clustering runs. All expensive PREP is already on disk
# (Preprocessing/ + Activity/ for both variants), so this only (re)runs the LemonTree
# clustering + regulator assignment where missing. Safe to re-run: it skips finished steps.
set -uo pipefail

RES=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/phospho_clustering
HERE=$(dirname "$(readlink -f "$0")")

for V in top5000_hvg all_sites; do
    OUT="$RES/$V/Lemon_out"
    # count finished regulator layers
    ndone=$(ls "$OUT"/{TFs,Metabolites,Proteins,Lipids,KinaseActivity}.topreg.txt 2>/dev/null | wc -l)
    if [ "$ndone" -eq 5 ]; then
        echo "[$V] complete (5/5 regulator layers) — skipping"
        continue
    fi
    if [ ! -s "$OUT/tight_clusters.txt" ]; then
        echo "[$V] no tight_clusters yet -> full clustering + regulators"
        nohup bash "$HERE/run_lemontree.sh" "$RES/$V" 10 0.25 6 > "$RES/$V/lemontree_run.log" 2>&1 &
        echo "  launched (pid $!), log: $RES/$V/lemontree_run.log"
    else
        echo "[$V] clustering present, $ndone/5 regulator layers done -> (re)running regulators"
        # run_lemontree.sh re-does clustering; to only redo regulators, call the reg block:
        # simplest safe path is to re-run the whole script (clustering is idempotent-ish but
        # slow); instead we just re-launch regulators against existing tight_clusters.
        (
          SIF=/home/borisvdm/repo/LemonIte/nextflow/lemontree-pipeline_v1.0.0.sif
          JAR=/opt/lemontree/lemontree_v3.1.1.jar
          P="$RES/$V/Preprocessing"
          for pair in TFs:tfs Metabolites:metabolites Proteins:proteins Lipids:lipids KinaseActivity:kinaseactivity; do
            pfx="${pair%%:*}"; rl="$P/${pair##*:}.txt"
            [ -s "$OUT/$pfx.topreg.txt" ] && continue
            singularity exec --cleanenv --bind "$RES:$RES" --bind /home/borisvdm/repo/LemonIte:/home/borisvdm/repo/LemonIte "$SIF" \
              java -cp "$JAR" lemontree.modulenetwork.RunCli -task regulators \
              -data_file "$P/LemonPreprocessed_complete.txt" -reg_file "$rl" \
              -cluster_file "$OUT/tight_clusters.txt" -output_file "$OUT/$pfx" > "$OUT/reg_$pfx.log" 2>&1 &
          done
          wait
        ) &
        echo "  launched regulators (pid $!)"
    fi
done
echo "Resume launched. Check progress with: ls $RES/*/Lemon_out/*.topreg.txt"
