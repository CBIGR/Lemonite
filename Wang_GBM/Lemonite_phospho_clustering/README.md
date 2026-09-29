# Lemonite phospho-clustering variant (Wang GBM)

A second phospho extension of the Wang GBM Lemonite analysis. **Here phosphoproteomics is the
CLUSTERING INPUT** (LemonTree modules = co-regulated phosphosites), and TFs, metabolites,
proteins, lipids **and kinase activity** are assigned as regulators of those phospho-modules.
This is distinct from `../Lemonite_phospho/`, where phospho was a *regulator* of the existing
gene-expression modules.

Two clustering variants (like the earlier extension): **`all_sites`** (all 34,076 phosphosites
passing the ≥50%-valid filter) and **`top5000_hvg`** (top-5,000 most-variable sites).

All pipeline steps run with the pipeline's own JAR/scripts inside
`nextflow/lemontree-pipeline_v1.0.0.sif`. The nextflow `main.nf` is **not** used directly (it
can't ingest a phospho-as-primary matrix); see `PIPELINE_ISSUES.md`.

## Pipeline order (run these in sequence)

### A. Data manipulation (pre-Lemonite) — produces LemonTree inputs

| # | Script | What it does | Output |
|---|---|---|---|
| 1 | `build_phospho_clustering_input.py` | mmc3 `phosphoproteome_normalized` → subset samples, drop <50%-valid sites, NA→0, z-score per site, (opt) top-N variable | `Preprocessing/LemonPreprocessed_expression.txt` (the clustering primary), `phospho_sites.txt`, `samples.txt` |
| 2 | `kinase_and_tf_activity.R` | decoupleR `run_ulm`: **kinase activity** on the phospho matrix vs OmniPath KSN (rows re-keyed `GENE_residue`); **TF activity** vs CollecTRI | `Activity/kinase_activity.txt`, `TF_activity.txt`, `LemonPreprocessed_kinaseactivity.txt` (regulator layer) |
| 3 | `build_regulators.py` | build z-scored regulator layers (Proteins, Metabolites, Lipids, TFs=Lovering∩proteome, KinaseActivity) aligned to the phospho samples, then assemble `LemonPreprocessed_complete.txt` (phospho rows + all regulator rows) | `Preprocessing/LemonPreprocessed_<layer>.txt`, `<layer>.txt`, `LemonPreprocessed_complete.txt` |

### B. Lemonite call (clustering + regulator assignment)

| # | Script | What it does |
|---|---|---|
| 4 | `run_lemontree.sh <variant-dir> [n_clusters] [min_weight]` | inside the `.sif`, runs LemonTree `ganesh` (N Gibbs samplers) → `tight_clusters` → `regulators` per layer, exactly as the pipeline would after preprocessing. Outputs `Lemon_out/tight_clusters.txt` + `<layer>.topreg.txt` |

Run per variant, e.g.:
```bash
RES=/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/phospho_clustering
PY=/home/borisvdm/Software/miniconda3/envs/LemonIte/bin/python3
SIF=/home/borisvdm/repo/LemonIte/nextflow/lemontree-pipeline_v1.0.0.sif

for V in all_sites "top5000_hvg --top-var 5000"; do
  set -- $V; name=$1; shift
  # 1. clustering input
  $PY build_phospho_clustering_input.py --out-dir "$RES/$name/Preprocessing" "$@"
  # 2. kinase + TF activity (in the container)
  singularity exec --cleanenv --bind /home/borisvdm/repo/LemonIte:/home/borisvdm/repo/LemonIte \
    --bind "$RES:$RES" "$SIF" Rscript kinase_and_tf_activity.R \
    --expr "$RES/$name/Preprocessing/LemonPreprocessed_expression.txt" --out-dir "$RES/$name/Activity"
  # 3. regulator layers + complete matrix
  $PY build_regulators.py --variant-dir "$RES/$name"
  # 4. LemonTree clustering + regulator assignment
  bash run_lemontree.sh "$RES/$name" 10 0.25
done
```

## Outputs (results dir, per variant)

`…/results/LemonTree/phospho_clustering/<variant>/`
- `Preprocessing/` — `LemonPreprocessed_expression.txt` (clustering input),
  `LemonPreprocessed_complete.txt`, per-layer regulator matrices + lists.
- `Activity/` — `kinase_activity.txt`, `TF_activity.txt` (feature × sample), regulator files.
- `Lemon_out/` — `tight_clusters.txt` (the phospho modules) and
  `{TFs,Metabolites,Proteins,Lipids,KinaseActivity}.topreg.txt` (regulator→module scores).

## Files here

- `build_phospho_clustering_input.py`, `kinase_and_tf_activity.R`, `build_regulators.py` — data manipulation (A).
- `run_lemontree.sh` — the Lemonite call (B).
- `PLAN.md` — design/decisions/progress. `PIPELINE_ISSUES.md` — why `main.nf` isn't used.
- `FIRST_ANALYSIS.md` — first results + biological interpretation (written after the run).

## Notes / caveats
- Kinase activity uses decoupleR `run_ulm` (not the full `decouple()` consensus, which fails
  on colinear kinases sharing substrates). Kinases with ≥5 matched substrate sites are kept
  (89 in top5000, 259 in all_sites).
- Samples: phospho + proteome cover all 109; metabolome/lipids ~82–83 (LemonTree pads the rest).
- First-pass clustering uses 10 Gibbs samplers (pipeline default is 100) for tractable runtime;
  raise `n_clusters` for a production run.
- Module members are **phosphosites**, so the pipeline's gene-based network/enrichment/MegaGO
  steps are not run here (would need phosphosite→gene mapping; see PIPELINE_ISSUES.md #4).
