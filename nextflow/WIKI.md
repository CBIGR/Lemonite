# Lemonite Nextflow Pipeline

This file is the detailed user guide for the Nextflow pipeline under `nextflow/`. The shorter entrypoint lives in [../README.md](../README.md).

## Pipeline Overview

Lemonite is a DSL2 Nextflow workflow for transcriptomics-metabolomics centered multi-omics data integration. The current pipeline implementation is defined by:

- `main.nf` for workflow order and parameter validation
- `nextflow.config` for defaults and profiles
- `conf/base.config` for published output locations
- `modules/*.nf` and `scripts/*` for step-specific behavior

At a high level, the published workflow is:

1. Preprocessing and optional TF activity inference.
2. Parallel LemonTree clustering runs.
3. Post-clustering consolidation and regulator assignment.
4. Network generation and regulator filtering.
5. Subnetwork graph generation and PKN-based evaluation.
6. Module heatmap generation.
7. Enrichment analysis.
8. Module overview assembly and top-level HTML summary report.

## Requirements

### Software

- Nextflow `>= 23.04.0`
- Java 11+
- Singularity runtime is the primary execution path in this repository

### Hardware

- Minimum: 16 GB RAM, 4 CPU cores
- Recommended: 32 GB RAM, 8+ CPU cores

## Profiles

The pipeline defines these profiles in `nextflow.config` and `conf/*.config`:

| Profile | Purpose | Notes |
| --- | --- | --- |
| `singularity` | Recommended default | Uses `lemontree-pipeline_v1.0.0.sif` and Singularity bind mounts |
| `docker` | Docker execution | Available, but repository scripts and examples focus on Singularity |
| `local` | Lightweight local execution | Uses the local executor with smaller defaults |
| `hpc` | Resource overrides for clusters | Typically combined with `singularity`, for example `-profile singularity,hpc` |
| `test` | Bundled smoke test dataset | Sets `input_dir` to `nextflow/test_dataset` and reduces runtime |
| `dev` | Development bind mounts | Intended for rapid iteration with host-mounted scripts and PKN files |

## Installation

```bash
git clone https://github.com/CBIGR/Lemonite.git
cd Lemonite/nextflow

curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

Build the container image used by the `singularity` profile:

```bash
./build-singularity.sh
```

The build script creates `lemontree-pipeline_v1.0.0.sif` in the `nextflow/` directory.

## Running the Pipeline

### Smoke Test

```bash
nextflow run main.nf \
  -profile test,singularity
```

This uses the bundled dataset in `test_dataset/` and a smaller cluster count from `conf/test.config`.

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
  --metadata_columns diagnosis \
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
  --metadata_columns diagnosis,Survival,Sex \
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
- `--work_dir` defaults to `{input_dir}/work`

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
| `--work_dir` | `{input_dir}/work` | Nextflow work directory |
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
| `--top_n_genes` | `5000` | Number of highly variable genes retained |
| `--perform_tfa` | `true` | Enables TF activity inference |
| `--use_omics_specific_scaling` | `true` | Uses transcriptomics-specific and metabolomics-style scaling |
| `--gene_annotation_file` | `null` | Optional offline annotation TSV for HPC/no-internet setups |
| `--deseq_contrast1` | `diagnosis` | Main metadata grouping column |
| `--design_formula` | `~ diagnosis` | DESeq2 design formula |
| `--metadata_columns` | `diagnosis` | Metadata columns retained for downstream use |
| `--heatmap_metadata_cols` | `null` | Optional override for heatmap annotations |
| `--expression_col` | `count` | Expression value column label used by scripts |
| `--sample_id_col` | `Sample_ID` | Sample identifier column in metadata |

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
| `--enrichment_method` | `EnrichR` | Allowed values: `EnrichR`, `GSEA`, `both`, `auto` |
| `--enrichr_libraries` | configured list | Default EnrichR libraries from `nextflow.config` |
| `--prioritize_by_expression` | `true` | Enables module prioritization from expression differences |
| `--overview_n_clusters` | `5` | Number of canonical MegaGO clusters in the overview |
| `--pkn_network` | `nextflow/PKN/Lemonite_PKN.tsv` | Prior knowledge network used for edge categorization and evaluation |

The overview stage now always uses one canonical workflow: select `BP` terms from a single enrichment source, keep the top `30` terms per module, cluster modules with MegaGO when available, and label each cluster with `rrvgo`. When both EnrichR and GSEA outputs are present, the overview prefers EnrichR and falls back to GSEA only when EnrichR outputs are absent.

### Advanced Cluster and Network Knobs

These are present in `nextflow.config` and logged by the pipeline, but are not part of the minimal quick-start surface.

| Parameter | Default | Purpose |
| --- | --- | --- |
| `--use_deseq_priors` | `true` | Use DESeq2 priors during module discovery |
| `--min_cluster_size` | `10` | Minimum cluster size |
| `--tight_clusters_only` | `false` | Restrict to tight clusters only |
| `--max_n_iterations` | `1000` | Maximum clustering iterations |
| `--min_regulator_size` | `3` | Minimum regulators per module |
| `--max_regulator_size` | `100` | Maximum regulators per module |
| `--min_module_size` | `10` | Minimum module size for network generation |
| `--min_targets` | `3` | Minimum targets per regulator |
| `--min_expression_fold_threshold` | `1.5` | Threshold used in network filtering |
| `--max_pvalue_threshold` | `0.05` | Threshold used in network filtering |

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
    │   └── PCA_*.pdf
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
    │   ├── subnetworks/graph_*.png
    │   └── cytoscape/*.tsv
    ├── ModuleViewer_files/
    │   ├── *.selected_regs_list.txt
    │   ├── *.selected_regulators_scores.txt
    │   ├── sample_mapping.mvf
    │   └── metabolite_LemoniteKG_interactions.mvf
    ├── PKN_Evaluation/
    ├── module_heatmaps/
    │   ├── heatmaps/*.png
    │   └── module_viewer_summary.txt
    ├── Enrichment/
    └── Module_Overview/
        ├── Module_Overview.csv
        ├── interactive_module_network.html
        ├── interactive_module_network_movable.html
        ├── module_network_edges.txt
        ├── module_network_node_attributes.txt
      ├── Module_Expression_Heatmap.png
      └── top_30/
        ├── bp_terms_top_30.csv
        ├── cluster_assignments_top_30.csv
        ├── rrvgo_cluster_labels_top_30.csv
        ├── rrvgo_module_labels_top_30.csv
        └── rrvgo_reduced_terms_top_30.csv
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

`-profile singularity,dev` bind-mounts the host `scripts/` and `PKN/` directories into the container so script changes can be tested without rebuilding the image.

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
