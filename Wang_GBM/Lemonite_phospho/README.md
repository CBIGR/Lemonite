# Adding phosphoproteomics as a regulator layer — Wang GBM Lemonite

This directory extends the existing **Wang GBM** Lemonite analysis by adding
**phosphoproteomics** as an additional regulator layer, **without touching the
existing gene co-expression modules or the other regulators** (TFs, metabolites,
lipids, proteins). Only LemonTree's *regulator-assignment* step is re-run, so the
new phospho regulators are scored against the exact same modules as before.

Everything here is kept separate from the original run:

| | Path |
|---|---|
| New scripts | `Wang_GBM/Lemonite_phospho/` (this folder) |
| New results | `…/results/LemonTree/noProteomics_percentile2_divide_by_sum_phospho/` |
| Base run (unmodified, read-only source of modules) | `…/results/LemonTree/noProteomics_percentile2_divide_by_sum/` |

## What "regulator assignment" is here

A LemonTree run has two conceptually separate phases:

1. **Clustering → tight clusters** — genes are grouped into co-expression *modules*.
   This produces `Lemon_out/tight_clusters.txt`. **We do NOT redo this.**
2. **Regulator assignment** (`-task regulators`) — for each candidate regulator,
   LemonTree scores how well its abundance profile explains the expression of each
   module and emits `<Prefix>.topreg.txt` (regulator → module → score). Each
   omics layer (TFs, metabolites, lipids, …) is one independent `-task regulators`
   call against the *same* `tight_clusters.txt`.

Adding phosphoproteomics is therefore just **one more `-task regulators` call**,
using the existing `tight_clusters.txt` and a phospho abundance matrix — no
re-clustering, existing regulators untouched.

## Input data

`phosphoproteome_normalized` sheet of the Wang 2021 CPTAC-GBM supplement:

```
/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/data/1-s2.0-S1535610821000507-mmc3.xlsx
```

- 70,330 phosphosites (peptide level) × 109 tumour samples.
- Row identity = `site_id`, e.g. `AAAS:NP_056480.1:S495`.
- Missing values encoded as the string `NA`.

## Normalization / scaling — what we did and why

**No re-normalization was needed.** Per the supplement's README (Tab 9), the
phosphoproteome is already:

> *log2 transformed, normalized by median polishing and batch corrected by ComBat;
> only entries with ≤ 90% NA across samples kept.*

The **only** transformation we apply is the same one the Lemonite pipeline applies
to every continuous feature so that all rows of `LemonPreprocessed_complete.txt`
share a common scale (verified: every gene/metabolite/lipid row in that file has
per-row mean ≈ 0, sd ≈ 1):

1. **Subset** to the 74 samples present in the existing LemonTree run
   (all 74 are present in the phospho data — full overlap, none padded).
2. **NA → 0** on the (already log2/median-polished) abundances.
3. **Z-score each phosphosite** (row-wise `scale`: subtract mean, divide by sd),
   i.e. `t(scale(t(x)))`.

This mirrors exactly the pipeline's PTM branch in
`nextflow/scripts/Preprocessing_TFA_Proteomics.R` (the `hPTM` block:
`hPTM[is.na(hPTM)] <- 0; hPTM <- t(scale(t(hPTM)))`), so the phospho layer is
treated identically to how the pipeline would treat a histone-PTM / PTM layer.

Sparsity filter: phosphosites measured in **< 50 %** of the 74 samples are dropped
(`--min-valid-frac 0.5`) — otherwise a mostly-`NA→0` row is dominated by imputed
zeros and adds noise. 34,255 of 70,330 sites pass this filter.

## Two runs

At the user's request, two independent regulator assignments were produced from
the same filtered, z-scored phospho matrix:

| Variant | Dir | # phospho regulators |
|---|---|---|
| **All phosphosites** | `all_phosphosites/Lemon_out/` | 34,255 (all sites with ≥ 50 % valid samples) |
| **Top-2000 most variable** | `top2000_variable/Lemon_out/` | 2,000 (highest pre-scaling variance) |

Variance for the top-2000 selection is computed on the NA→0 abundances **before**
z-scoring (z-scoring flattens every row's variance to ≈ 1, so it must be measured
first).

## How to reproduce

```bash
cd Wang_GBM/Lemonite_phospho
BASE=…/results/LemonTree/noProteomics_percentile2_divide_by_sum/Lemon_out
RD=…/results/LemonTree/noProteomics_percentile2_divide_by_sum_phospho
MMC3=…/Wang2021/data/1-s2.0-S1535610821000507-mmc3.xlsx
PY=…/miniconda3/envs/LemonIte/bin/python3   # any env with numpy + openpyxl

# 1. Build phospho abundance matrix + append to complete matrix (both variants)
for v in "all_phosphosites:0" "top2000_variable:2000"; do
  name=${v%%:*}; top=${v##*:}
  mkdir -p "$RD/$name/Lemon_out"
  cp "$BASE/tight_clusters.txt" "$RD/$name/Lemon_out/"      # existing modules, kept as-is
  $PY prepare_phospho_regulator.py \
      --mmc3 "$MMC3" --complete "$BASE/LemonPreprocessed_complete.txt" \
      --out-dir "$RD/$name/Lemon_out" --min-valid-frac 0.5 --top-var "$top"
done

# 2. Run ONLY the regulator-assignment step (LemonTree jar inside the pipeline .sif)
bash run_phospho_regulators.sh "$RD/all_phosphosites/Lemon_out"
bash run_phospho_regulators.sh "$RD/top2000_variable/Lemon_out"
```

LemonTree is executed inside the pipeline's Singularity image
(`nextflow/lemontree-pipeline_v1.0.0.sif`, jar
`/opt/lemontree/lemontree_v3.1.1.jar`) so the runtime is identical to the pipeline.

## Files

| File | Purpose |
|---|---|
| `prepare_phospho_regulator.py` | Parse mmc3 → subset samples → NA→0 → z-score → (optional) top-var → write reg list + phospho matrix + `LemonPreprocessed_complete_phospho.txt`. |
| `run_phospho_regulators.sh` | Run `lemontree -task regulators` in the Singularity image against the existing `tight_clusters.txt`. |
| `stage_downstream.sh` | Assemble a per-variant `run/` dir combining the base run's TF/Metabolite/Lipid regulator outputs with the new Phospho outputs. |
| `run_downstream.sh` | Run the nextflow post-regulator steps (network → subnetworks → PKN eval → heatmaps → enrichment → module overview → summary) over the phospho-augmented regulator set. |
| `inject_base_megago.py` | Merge the base run's per-module MegaGO/rrvgo functional-cluster labels into a variant's Module_Overview (modules are identical, so MegaGO is reused, not recomputed). |
| `README.md` | This file. |

## Results

Both regulator assignments completed against the same 63 existing modules (module
ids 0–62, 2,589 genes). `Phospho.topreg.txt` = `phosphosite <TAB> module_id <TAB> score`.

| Variant | phospho regulators fed in | assignments (topreg rows) | distinct phosphosites assigned | modules with ≥1 phospho regulator |
|---|---|---|---|---|
| Top-2000 variable | 2,000 | 275 | 216 | 62 / 63 |
| All phosphosites | 34,255 | 395 | 359 | 62 / 63 |

Sanity check — the assignments are biologically coherent. E.g. module 0 is a
neuronal/synaptic module (SYT4, SYNGR3, SYP, STX1B, GABRA5, KCNS1) and its top
phospho regulators are ion-channel / synaptic-protein phosphosites
(KCNQ2 S476, ANK2 S1858/S1891, GNG3 T5, HEPACAM S364;S368).

### Outputs (per variant, in `<variant>/Lemon_out/`)

- `phospho.txt` — regulator list (one phosphosite id per line).
- `LemonPreprocessed_phospho.txt` — z-scored phospho abundance matrix.
- `LemonPreprocessed_complete_phospho.txt` — base complete matrix **+** phospho rows
  (the `-data_file` LemonTree reads regulator profiles from).
- `tight_clusters.txt` — copy of the base run's modules (unchanged; input only).
- `Phospho`, `Phospho.topreg.txt`, `Phospho.allreg.txt`, `Phospho.xml.gz`, … —
  the phospho→module regulator assignments (LemonTree output). `Phospho.topreg.txt`
  is `phosphosite <TAB> module_id <TAB> score`, same format as
  `Metabolites.topreg.txt` etc. in the base run.

---

# Downstream pipeline (network → … → summary report)

After regulator assignment, the nextflow steps that run **after** POST_CLUSTERING
were reproduced over the phospho-augmented regulator set. Phospho is added as a
**4th regulator layer** alongside the base run's TFs (Lovering), Metabolites and
Lipids; the existing three layers are unchanged.

`run_downstream.sh <variant>` runs, inside the pipeline Singularity image, the same
scripts `main.nf` calls after POST_CLUSTERING:

| Step | nextflow process | script | notes |
|---|---|---|---|
| 1 | NETWORK_GENERATION | `lemontree_to_network.py` | selects regulators, builds `LemonNetwork_*`, `Phospho2targets_*`, coherence-filters modules |
| 2 | SUBNETWORK_GRAPHS | `create_subnetwork_graphs.py` | per-module graphviz subnetworks |
| 3 | PKN_EVALUATION | `evaluate_against_PKN.py` | in-silico validation vs `PKN/Lemonite_PKN.tsv`; PPI/metabolite enrichment; `.mvf` files |
| 4 | MODULE_VIEWER_HEATMAPS | `module_viewer.py` | per-module expression + regulator-score heatmaps |
| 5 | ENRICHMENT_ANALYSIS | `enrichment_analysis.R` | EnrichR per-module pathway enrichment |
| 6 | MODULE_OVERVIEW | `module_overview_interactive.py` | interactive overview; **MegaGO reused, not recomputed** (see below) |
| 7 | SUMMARY_REPORT | `generate_summary_report.py` | `Lemonite_Summary_Report.html` |

### Key parameters (match the original GBM run)

- `regulator_types = Lovering:lovering_TFs.txt, Metabolites:metabolites.txt, Lipids:lipids.txt, Phospho:phospho.txt`
- `regulator_selection_method = percentage`, `top_n_percent_regulators = 2.0`
- **`coherence_threshold = 0.5`** — the original GBM analysis kept **46 modules**,
  which corresponds to a 0.5 coherence cutoff (the pipeline *default* of 0.6 would
  keep only 36). The tight-cluster modules and their coherence scores are identical
  to the base run; only this threshold controls how many pass the filter.

### MegaGO functional clustering is reused, not recomputed

MegaGO clusters **modules** by GO-term semantic similarity. The modules are identical
to the base run, so MegaGO would give the same answer. Instead of recomputing it, the
overview is run with `--n_clusters 1` (which disables the MegaGO computation) and the
base run's per-module `MegaGO_cluster` / `MegaGO_label` / `MegaGO_representative_GO_ID`
are joined back in by module id (`inject_base_megago.py`). This is faster and
guarantees the functional clustering matches the initial analysis.

### Downstream results

Both variants reproduce the original **46-module** network
(`LemonNetwork_top2.0pct_46modules.txt`) with Phospho added as a regulator layer, and
all 7 steps complete successfully. Every module receives phospho regulators.

| | Top-2000 variable | All phosphosites |
|---|---|---|
| Modules in network | 46 | 46 |
| Selected Phospho→target edges (`Phospho2targets`) | 319 | 540 |
| Phospho-gene interactions (PKN eval) | 19,048 | 25,557 |
| Per-module heatmaps | 46 | 46 |
| Subnetwork graphs | 38 | 38 |
| MegaGO labels injected (from base run) | 44/46 | 44/46 |
| Summary report | ✓ | ✓ |

The per-module MegaGO labels (e.g. Cluster_5 = "synaptic signaling",
Cluster_2 = "immune system process") are carried over from the base run, and the
synaptic modules are regulated by synaptic/ion-channel phosphosites (SYN1, KCNQ2,
ANK2, GNG3 …), consistent with the regulator-assignment sanity check above.

### Outputs (per variant, in `<variant>/run/`)

- `Networks/` — `LemonNetwork_top2.0pct_46modules.txt`, `Phospho2targets_*`,
  `<Prefix>2targets_*`, `specific_modules.txt`, `Module_coherence_scores.txt`,
  ranked-regulator tables, `subnetworks/`.
- `ModuleViewer_files/` — `*.selected_regs_list.txt`, `evaluation_summary.txt`,
  `metabolite_LemoniteKG_interactions.mvf`, `PPI_interactions.mvf`, `HumanNet_interactions.mvf`.
- `heatmaps/` — per-module expression + regulator-score heatmaps.
- `Enrichment/` — EnrichR per-module pathway CSVs.
- `Module_Overview/Module_Overview.csv` — module overview incl. `Phospho_regulators`
  column and the reused `MegaGO_Cluster` / `MegaGO_Label`.
- `Lemonite_Summary_Report.html` — final HTML report.
- `downstream_logs/0N_*.log` — per-step logs.

### How to run the downstream

```bash
cd Wang_GBM/Lemonite_phospho
bash stage_downstream.sh top2000_variable       # stage run/ dir (base regs + Phospho)
bash stage_downstream.sh all_phosphosites
bash run_downstream.sh   top2000_variable 1 7   # steps 1..7
bash run_downstream.sh   all_phosphosites 1 7
```
