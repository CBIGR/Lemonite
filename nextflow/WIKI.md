# Lemonite Nextflow Pipeline

This file is the detailed user guide for the Nextflow pipeline under `nextflow/`. The shorter entrypoint lives in [../README.md](../README.md).

## Pipeline Overview

Lemonite is a DSL2 Nextflow workflow for transcriptomics-metabolomics centered multi-omics data integration. The current pipeline implementation is defined by:

- `main.nf` for workflow order and parameter validation
- `nextflow.config` for defaults and profiles
- `conf/base.config` for published output locations
- `modules/*.nf` and `scripts/*` for step-specific behavior

At a high level, the published workflow is:

1. **Parameter logging** — records the full effective parameter set to `pipeline_parameters_log.txt`.
2. **Preprocessing and TFA** — normalises the expression matrix (DESeq2 or pre-scaled), selects the top variable features, optionally infers transcription-factor activity (TFA) via decoupleR/CollecTRI, and writes LemonTree-ready input files.
3. **Parallel LemonTree clustering** — runs `n_clusters` independent Gibbs-sampling instances in parallel, each producing its own cluster assignment.
4. **Post-clustering** — consolidates the parallel cluster runs into tight (high-coherence) clusters and assigns regulators (TFs, metabolites, and any other configured regulator types) to each module.
5. **Network generation and subnetwork graphs** — filters modules by coherence threshold, selects top regulators per module, builds the regulator→target network, exports Cytoscape-compatible files, and renders PNG subnetwork graphs via pygraphviz.
6. **PKN-based evaluation** — compares the data-driven network against the Lemonite prior-knowledge network (PKN) to categorise edges as known or novel, and generates a metabolite-gene knowledge-graph MVF file.
7. **Module heatmaps** — generates per-module expression heatmaps for all omics layers, annotated with the selected metadata columns.
8. **Pathway enrichment analysis** — runs EnrichR and/or GSEA on module gene sets and regulator-target sets for all configured regulator types.
9. **Module overview** — builds the canonical MegaGO + rrvgo functional overview: BP enrichment terms are clustered into `overview_n_clusters` semantic groups and each module is labelled, producing a module-expression heatmap and CSV outputs.
10. **Summary report** — assembles all results into a self-contained `Lemonite_Summary_Report.html` with interactive network visualisations and key statistics.

## Requirements

### Software

- Nextflow `>= 23.04.0`
- Java 11+
- Singularity (required — the only supported execution backend)

### Hardware

- Minimum: 16 GB RAM, 4 CPU cores
- Recommended: 32 GB RAM, 8+ CPU cores

## Profiles

The pipeline defines these profiles in `nextflow.config` and `conf/*.config`:

| Profile | Purpose | Notes |
| --- | --- | --- |
| `singularity` | Required execution profile | Uses `lemontree-pipeline_v1.0.0.sif`; must be specified for every run |
| `local` | Lightweight local execution | Uses the local executor with reduced resource defaults |
| `hpc` | Resource overrides for clusters | Combine with `singularity`: `-profile singularity,hpc` |
| `test` | Bundled smoke-test dataset | Sets `input_dir` to `nextflow/test_dataset` and reduces cluster count; always combine with `singularity` |
| `dev` | Development bind mounts | Bind-mounts the host `scripts/` and `PKN/` directories so script changes are tested without rebuilding the image; combine with `singularity` |

## Installation

```bash
git clone https://github.com/CBIGR/Lemonite.git
cd Lemonite/nextflow

curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

Build the Singularity container image:

```bash
./build-singularity.sh
```

The build script creates `lemontree-pipeline_v1.0.0.sif` in the `nextflow/` directory. All pipeline runs require this image via the `singularity` profile.

### Building on HPC (sandbox containers)

On many HPC systems — including UGent's kyukon/VSC clusters — `apptainer build`
**cannot produce a `.sif` image**. These clusters run an *unprivileged* Apptainer
(no setuid `starter-suid`) and user accounts have no `/etc/subuid` mapping, so
`--fakeroot` falls back to `proot` for the final `squashfs` step. On the RHEL9
worker nodes `proot`'s `ptrace` is blocked by seccomp, and the build aborts with:

```
proot error: ptrace(TRACEME): Operation not permitted
FATAL:   While performing build: while creating squashfs: ... mksquashfs command failed
```

This affects every route through `apptainer build`'s SIF/squashfs step
(`PROOT_NO_SECCOMP=1`, `--userns`, sandbox→sif conversion, and pointing Apptainer
at the system `mksquashfs` all fail identically).

**Workaround — build and run a sandbox directory.** A `--sandbox` build produces a
plain root-filesystem *directory* instead of a squashfs `.sif`, so it never invokes
`proot`. Apptainer (and therefore Nextflow) can `exec` a sandbox directory exactly
like a `.sif`. Use the parallel build script, which does this automatically:

```bash
./build-singularity_parallel.sh
```

It:

1. Builds with `apptainer build --fakeroot --fix-perms --sandbox …` on the node's
   local `/tmp` (UGent recommends building on local disk, not scratch).
2. Relocates the finished sandbox to the `nextflow/` directory on scratch using
   `tar --xattrs` (a plain `mv`/`cp` fails on the fakeroot-owned `apt/.../partial`
   directories and the `user.rootlesscontainers` ownership-faking xattrs).
3. Produces `lemontree-pipeline_v1.0.0_parallel.sandbox/` (a ~7 GB directory) and
   smoke-tests it (`R --version`, `python3 --version`).

Point the container at the sandbox **directory path** (no `.sif` suffix):

```groovy
// conf/singularity.config or a dataset-specific config
process.container = "/abs/path/to/nextflow/lemontree-pipeline_v1.0.0_parallel.sandbox"
```

Notes:

- Run it like any other container: `nextflow run main.nf -profile singularity …`.
- A sandbox is a directory of many small files; keep it on scratch (not `$VSC_HOME`).
- The proper long-term fix is to ask HPC support to add `/etc/subuid` and
  `/etc/subgid` entries for your account, which makes native `--fakeroot` `.sif`
  builds work again.
- `git lfs` is required to fetch the `PKN/` and `scripts/` data the build copies in
  — on UGent HPC load it with `module load git-lfs/3.6.1` before `git lfs pull`.

## Running the Pipeline

### Smoke Test

```bash
nextflow run main.nf \
  -profile test,singularity
```

This uses the bundled dataset in `test_dataset/` and a reduced cluster count (`n_clusters = 5`, `top_n_genes = 1000`) from `conf/test.config`. The `singularity` profile is mandatory; the Singularity image (`lemontree-pipeline_v1.0.0.sif`) must be built before running.

### Standard Run

```bash
nextflow run main.nf \
  --input_dir /path/to/project \
  --organism human \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  --n_clusters 100 \
  --top_n_genes 5000 \
  --coherence_threshold 0.6 \
  --regulator_selection_method percentage \
  --top_n_percent_regulators 2.0 \
  -profile singularity \
  -resume
```

### Standard Run with Fold-Based Regulator Selection

```bash
nextflow run main.nf \
  --input_dir /path/to/project \
  --organism human \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:Proteomics.txt" \
  --n_clusters 5 \
  --top_n_genes 2000 \
  --coherence_threshold 0.6 \
  --regulator_selection_method fold_per_module \
  --regulator_fold_cutoff 2.0 \
  --metadata_columns "diagnosis,Survival,Sex" \
  -profile singularity \
  --max_cpus 6 \
  -resume
```

### HPC Run

```bash
nextflow run main.nf \
  --input_dir /path/to/project \
  --organism human \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  -profile singularity,hpc \
  --max_cpus 32 \
  --max_memory 128.GB \
  -resume
```

## Results Location and Run IDs

The published analysis results are written to:

```text
{output_dir or input_dir/results}/{computed_run_id}/
```

Two related defaults matter:

- `--output_dir` defaults to `{input_dir}/results`
- `--work_dir` is intended to default to `{input_dir}/work`, but **currently does not**:
  `workDir` is set at the top of `nextflow.config`, which Nextflow evaluates before
  `--input_dir` (or a profile's `input_dir`) has been merged into `params`, so it always
  falls back to `{launchDir}/work`. Pass `-w /path/to/work` on the command line if you
  need the work directory somewhere specific.

If you do not provide `--run_id`, the pipeline auto-generates one using:

```text
{top_n_genes}HVG_coherence{coherence_threshold}_{method_suffix}_clusters{n_clusters}
```

Examples:

- `5000HVG_coherence0.6_top2.0pct_clusters100`
- `2000HVG_coherence0.6_fold2.0x_clusters5`

## Input Directory Layout

The workflow expects an `input_dir` containing a `data/` folder.

```text
input_dir/
└── data/
    ├── Counts.tsv
    ├── Metadata.txt
    ├── Metabolomics.txt
    ├── Lipidomics.txt
    ├── Proteomics.txt
    ├── name_map.csv
    └── MyCustom_network.txt
```

Only the expression matrix and metadata are mandatory at pipeline start. Additional files are required only if you reference them through `--regulator_types` or want metabolite-aware validation outputs.

## Input File Rules

### Expression Matrix

Expression can be provided either by auto-detection in `data/` or explicitly with `--expression_file`.

Auto-detection looks for one of these patterns:

- `Counts.tsv`
- `*counts*.tsv`
- `*expression*.tsv`
- `*host_tx_counts.tsv`

Rules:

- One matching file must be found, otherwise preprocessing stops with an error.
- Rows are features and columns are samples.
- The pipeline is designed around transcriptomics-like inputs, but can also be used with other expression-like matrices when biologically appropriate.

### Metadata

Metadata can be auto-detected or explicitly provided with `--metadata_file`.

Auto-detection looks for:

- `Metadata.txt`
- `*metadata*.txt`

Expected columns:

- `Sample_ID` by default, configurable with `--sample_id_col`
- A primary grouping column such as `diagnosis`, configurable with `--deseq_contrast1`

### Regulator Files

Regulators are configured with a comma-separated `--regulator_types` string using this format:

```text
Prefix:DataFile[:DataType]
```

Where:

- `Prefix` is the regulator label used throughout the pipeline, for example `TFs`, `Metabolites`, `Lipids`, `Proteins`
- `DataFile` is the file name to look up
- `DataType` is optional and can be:
  - `c` for continuous data, the default
  - `d` for discrete or binary data

Examples:

```bash
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt"
--regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Proteins:Proteomics.txt"
--regulator_types "TFs:MyTFs.txt,ClinicalParams:clinical_binary.txt:d"
```

Important behavior:

- If a regulator entry looks like a TF list, Lemonite checks `data/` first and then falls back to the bundled copy in `nextflow/PKN/`.
- The built-in human default is `Lovering_TF_list.txt`.
- When `--organism mouse` is used, the preprocessing script automatically swaps the default TF file to `lovering_TF_list_mouse.txt`.
- Discrete regulators marked with `:d` are passed through without the continuous-data transformations used for abundance matrices.

### Choosing a Preprocessing Script

The pipeline ships two preprocessing scripts that are selected via `--preprocessing_type`:

| `preprocessing_type` | Script | Use when |
| --- | --- | --- |
| `rna` *(default)* | `Preprocessing_TFA_RNA.R` | Input is raw RNA-seq count data. DESeq2 is used for normalisation and variance-stabilising transformation. |
| `proteomics` | `Preprocessing_TFA_Proteomics.R` | Input is pre-scaled continuous data (proteomics, phospho-proteomics, multi-omics). No DESeq2 is applied; top variable features are selected by variance. |

Both scripts:

- Select the top `--top_n_genes` most variable features.
- Optionally run **TF Activity inference** via decoupleR/CollecTRI when `--perform_tfa true`.
- Produce **PCA plots** (per omics layer) in `LemonTree/Preprocessing/`.
- Write identical output files consumed by the rest of the pipeline:
  - `LemonPreprocessed_expression.txt` — primary expression matrix for clustering
  - `LemonPreprocessed_complete.txt` — all omics layers concatenated
  - `TFA_consensus.txt` + `TFA/*.pdf` — TFA results (when `perform_tfa = true`)
  - `DESeq_groups.txt`, `name_mapping.tsv`
  - `hptms.txt`, `metabolites.txt` — regulator name lists for post-clustering

**When to disable TFA** (`--perform_tfa false`): use this for proteomics or multi-omics datasets where no RNA-level TF expression is available (i.e., you cannot infer TF activity from the measured features). TFA will still run on proteomics data if transcription factors are present in the dataset.

**Example** — running the proteomics mode:

```bash
nextflow run main.nf \
  --preprocessing_type proteomics \
  --perform_tfa false \
  -c /path/to/your.config
```

### Custom Prior Network for TFA

There is no dedicated pipeline parameter for a prior network at the Nextflow level. The preprocessing module instead auto-detects a custom TF prior network by file name inside `data/`.

If you want to override the bundled CollecTRI network, place a file matching either of these patterns inside `data/`:

- `*CollecTRI*.txt`
- `*network*.txt`

If no custom file is found, Lemonite uses the bundled CollecTRI network from `nextflow/PKN/`, with human or mouse selection determined by `--organism`.

### Metabolite Name Mapping

`name_map.csv` is optional but recommended when metabolite regulators are used.

Expected columns:

- `Query`
- `HMDB`

It is used for:

- PKN-based metabolite validation outputs
- restoring original metabolite names in downstream visualizations

## Parameter Reference

### Core Parameters

| Parameter | Default | Notes |
| --- | --- | --- |
| `--input_dir` | required | Directory containing `data/` |
| `--output_dir` | `{input_dir}/results` | Parent directory for published runs |
| `--work_dir` | `{launchDir}/work` | Nextflow work directory. The `{input_dir}/work` default does not take effect — see the note above; use `-w` to override. |
| `--run_id` | auto-generated | If omitted, derived from analysis settings |
| `--organism` | `human` | Allowed values: `human`, `mouse` |
| `--expression_file` | auto-detect | Explicit expression matrix path |
| `--metadata_file` | auto-detect | Explicit metadata path |
| `--max_cpus` | `16` | Upper bound for process CPU allocation |
| `--max_memory` | `64.GB` | Upper bound for process memory allocation |
| `--max_time` | `24.h` | Upper bound for process runtime |

### Preprocessing and Metadata Parameters

| Parameter | Default | Notes |
| --- | --- | --- |
| `--preprocessing_type` | `rna` | Selects the preprocessing script: `rna` (DESeq2 + TFA, for RNA-seq counts) or `proteomics` (pre-scaled data, no DESeq2) |
| `--top_n_genes` | `5000` | Number of highly variable genes retained |
| `--perform_tfa` | `true` | Enables TF activity inference (works for both preprocessing types) |
| `--use_omics_specific_scaling` | `true` | Pareto scaling for metabolomics/lipidomics; z-score scaling for transcriptomics. Set `false` to force z-score for all layers (legacy behaviour). |
| `--gene_annotation_file` | `null` | Path to a pre-downloaded BioMart annotation TSV. Required on HPC nodes without internet access. |
| `--metabolomics_labels_file` | `null` | Optional path to a metabolite name→ID mapping TSV/CSV. Overrides `data/metabolomics_name_map.csv` when provided. |
| `--deseq_contrast1` | `diagnosis` | Main metadata grouping column, used for DESeq2 and PCA colouring |
| `--design_formula` | `~ diagnosis` | DESeq2 design formula; adjust for confounders, e.g. `"~ batch + diagnosis"` |
| `--metadata_columns` | `diagnosis` | Comma-separated metadata columns retained for downstream visualisations and heatmap annotations |
| `--heatmap_metadata_cols` | `null` | Optional override specifying which metadata columns appear on heatmaps (defaults to `metadata_columns`) |
| `--expression_col` | `count` | Column name holding expression values in the input matrix |
| `--sample_id_col` | `Sample_ID` | Sample identifier column in the metadata file |
| `--stop_after_network` | `false` | When `true`, the pipeline exits after network generation and subnetwork graphs (steps 1–5), skipping enrichment, heatmaps, overview, and the summary report. Useful for quick network inspection. |

### Clustering and Regulator Parameters

| Parameter | Default | Notes |
| --- | --- | --- |
| `--n_clusters` | `100` | Number of parallel LemonTree runs |
| `--random_seed` | `42` | Base seed used across clustering runs |
| `--coherence_threshold` | `0.6` | Module coherence cutoff used downstream |
| `--regulator_types` | `TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt` | Comma-separated regulator definitions |
| `--regulator_selection_method` | `percentage` | Allowed values: `percentage`, `fold_per_module` |
| `--top_n_percent_regulators` | `2.0` | Used by the `percentage` selector |
| `--regulator_fold_cutoff` | `2.0` | Used by the `fold_per_module` selector |

### Enrichment and Overview Parameters

| Parameter | Default | Notes |
| --- | --- | --- |
| `--enrichment_method` | `EnrichR` | Allowed values: `EnrichR`, `GSEA`, `both`. `auto` is treated as `both`. EnrichR requires internet access; use `GSEA` for offline runs. |
| `--enrichr_libraries` | `GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,Reactome_Pathways_2024` | Comma-separated EnrichR library names |
| `--prioritize_by_expression` | `true` | Ranks modules by differential expression magnitude in the overview |
| `--overview_n_clusters` | `5` | Number of canonical MegaGO functional clusters shown in the module overview |
| `--interactive_overview` | `false` | When `true`, additionally generates interactive HTML network visualisations in `Module_Overview/` |
| `--skip_megago` | `false` | Skip MegaGO semantic-similarity clustering in the module overview and use the pathway-similarity (Jaccard) fallback instead. **Defaults to `true` under `-profile test`**, since the smoke-test dataset doesn't need the (slower) canonical MegaGO workflow — pass `--skip_megago false` to exercise MegaGO on the test dataset. |
| `--pkn_network` | `PKN/Lemonite_PKN.tsv` | Prior-knowledge network TSV used for edge categorisation (known vs. novel) in the module overview. **Note:** `PKN_EVALUATION` and `SUBNETWORK_GRAPHS` do not read this parameter — they search `{projectDir}/PKN/Lemonite_PKN.tsv` first and fall back to the copy baked into the container at `/opt/PKN/Lemonite_PKN.tsv`. Pointing `--pkn_network` at a different file therefore only affects the overview stage. |

The overview stage uses a fixed canonical workflow: it selects `GO Biological Process` terms from a single enrichment source, retains the top 30 terms per module, clusters modules semantically with MegaGO (when available and `--skip_megago` is not set), and labels each cluster with rrvgo. When both EnrichR and GSEA outputs are present, the overview prefers EnrichR and falls back to GSEA only when EnrichR outputs are absent. Set `--interactive_overview true` to also produce interactive HTML network files. Set `--skip_megago true` to bypass the MegaGO binary entirely (e.g. when it isn't installed, or for a faster test run) — module clustering then falls back to Jaccard pathway similarity, the same fallback used automatically when the `megago` binary isn't found on `PATH`.

### Advanced Cluster and Network Knobs

These are present in `nextflow.config` and logged by the pipeline, but are not part of the minimal quick-start surface.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--use_deseq_priors` | `true` | Use DESeq2 priors during module discovery |
| `--min_cluster_size` | `10` | Minimum number of genes a cluster must have |
| `--tight_clusters_only` | `false` | Restrict all downstream analysis to tight (high co-occurrence) clusters only |
| `--lemontree_tight_min_weight` | `0.25` | Minimum co-occurrence weight for a module to be called a tight cluster (0.0–1.0) |
| `--max_n_iterations` | `1000` | Maximum Gibbs-sampling iterations per clustering run |
| `--random_seed` | `42` | Base random seed; each parallel clustering run adds its index to this value |
| `--min_regulator_size` | `3` | Minimum number of regulators a module must have to be included in the network |
| `--max_regulator_size` | `100` | Maximum regulators per module (excess regulators are ranked and trimmed) |
| `--min_module_size` | `10` | Minimum module size (in genes) for network generation |
| `--min_targets` | `3` | Minimum number of targets a regulator must have to be retained in the network |
| `--min_expression_fold_threshold` | `1.5` | Minimum expression fold-change used when filtering network edges |
| `--max_pvalue_threshold` | `0.05` | Maximum p-value threshold for network edge inclusion |

## Output Structure

The published run directory is assembled through `conf/base.config` publish rules.

```text
{output_parent}/{run_id}/
├── Lemonite_Summary_Report.html
├── pipeline_parameters_log.txt
└── LemonTree/
    ├── Preprocessing/
    │   ├── LemonPreprocessed_expression.txt
    │   ├── LemonPreprocessed_complete.txt
    │   ├── LemonPreprocessed_*.txt
    │   ├── DESeq_groups.txt
    │   ├── name_mapping.tsv
    │   ├── PCA_*.pdf
    │   └── PCA_*.png                   # embedded in the HTML summary report
    ├── Lemon_out/
    │   ├── Lemon_results/cluster_*/
    │   ├── *.allreg.txt
    │   ├── *.randomreg.txt
    │   ├── clusters_list.txt
    │   └── tight_clusters.txt
    ├── Networks/
    │   ├── LemonNetwork_*.txt
    │   ├── *2targets*.txt
    │   ├── Cytoscape_*.txt
    │   ├── Module_coherence_scores.txt
    │   ├── specific_modules.txt
    │   └── subnetworks_graphviz/graph_*_graphviz.{png,pdf}
    ├── ModuleViewer_files/
    │   ├── *.selected_regs_list.txt
    │   ├── *.selected_regulators_scores.txt
    │   ├── clusters_list.txt
    │   ├── sample_mapping.mvf
    │   ├── metabolite_LemoniteKG_interactions.mvf
    │   ├── PPI_interactions.mvf
    │   ├── HumanNet_interactions.mvf
    │   ├── PPI_enrichment_results.csv
    │   ├── Metabolite_Gene_enrichment_results.csv
    │   └── evaluation_summary.txt      # PKN evaluation metrics
    ├── module_heatmaps/
    │   ├── heatmaps/Module_*_heatmap.{png,pdf}
    │   └── module_viewer_summary.txt
    ├── Enrichment/
    │   ├── Modules_enrichr/            # or Modules_gsea/
    │   ├── TFs2targets/
    │   └── Metabolites2targets/
    └── Module_Overview/
        ├── Module_Overview.csv
        ├── Module_Overview.xlsx
        ├── module_expression_analysis.csv
        ├── interactive_module_network_movable.html
        ├── regulator_rankings.html
        ├── module_network_edges.txt
        ├── module_network_node_attributes.txt
        └── top_30/
            ├── bp_terms_top_30.csv
            ├── cluster_assignments_top_30.csv
            ├── rrvgo_cluster_labels_top_30.csv
            ├── rrvgo_module_labels_top_30.csv
            └── rrvgo_reduced_terms_top_30.csv

`PKN_Evaluation/` is created by the pipeline but `evaluate_against_PKN.py` writes its
outputs (including `evaluation_summary.txt`) into `ModuleViewer_files/` instead, so the
directory is normally empty and is not published.
```

Top-level notes:

- `Lemonite_Summary_Report.html` is the main report intended for end users.
- `pipeline_parameters_log.txt` records the effective parameter set and full command line.
- The `LemonTree/` directory contains the step-specific published artifacts.

### Nextflow Trace Artifacts

Execution timeline, trace, DAG, and report files are configured separately from the published analysis outputs. By default they are written under:

```text
./results/pipeline_info/
```

relative to the pipeline launch directory, unless the underlying Nextflow tracing parameters are overridden.

## Practical Notes

### EnrichR Requires Internet Access

`EnrichR` calls the remote service. If compute nodes do not have outbound internet access:

- use `--enrichment_method GSEA` for fully offline enrichment
- or use `--enrichment_method both` only when internet access is available

### Multiple Candidate Input Files

If `data/` contains more than one file matching the expression or metadata auto-detection rules, the pipeline stops and asks for explicit `--expression_file` or `--metadata_file` values.

### Human versus Mouse Analyses

`--organism` controls:

- TF list selection
- gene annotation lookup
- CollecTRI prior network fallback
- enrichment organism settings

If you run mouse data with `--organism human`, the default TF list and enrichment settings will not match your input.

### Development Mode

`-profile singularity,dev` bind-mounts the host `scripts/` and `PKN/` directories into the container so script changes can be tested without rebuilding the image. Each module's process block already checks for host-side scripts before falling back to the container path, so changes in `scripts/` take effect immediately.

### Singularity is the Only Supported Runtime

Docker support has been removed from this pipeline. All runs must use `-profile singularity` (optionally combined with `hpc`, `local`, or `dev`). Build the image once with `./build-singularity.sh` and reuse it across runs.

## Citation

Use the Lemonite preprint citation from [../README.md](../README.md#citation):

```bibtex
@article{vandemoortele2026lemonite,
  title={Lemonite: identification of regulatory metabolites through data-driven, interpretable integration of transcriptomics and metabolomics data},
  author={Vandemoortele, Boris and Devlies, Hilde and Michoel, Tom and Vanhaecke, Lynn and Vandenbroucke, Roosmarijn E. and Laukens, Debby and Vermeirssen, Vanessa},
  journal={bioRxiv},
  year={2026},
  doi={10.64898/2026.03.27.714373},
  url={https://doi.org/10.64898/2026.03.27.714373}
}
```


Original LemonTree citation:

```bibtex
@article{bonnet2015lemontree,
  title={Integrative Multi-omics Module Network Inference with Lemon-Tree},
  author={Bonnet, Eric and Calzone, Laurence and Michoel, Tom},
  journal={PLOS Computational Biology},
  volume={11},
  number={2},
  pages={e1003983},
  year={2015},
  doi={10.1371/journal.pcbi.1003983}
}
```
