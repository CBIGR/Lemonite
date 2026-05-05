# Lemonite Nextflow Pipeline

This directory contains the executable Nextflow pipeline for Lemonite multi-omics integration.

Detailed documentation lives in [WIKI.md](WIKI.md). The repository-level overview lives in [../README.md](../README.md).

## Requirements

- Nextflow `>= 23.04.0`
- Java 11+
- Singularity runtime is the primary execution path for this repository

## Quick Start

### Smoke Test

```bash
./build-singularity.sh

nextflow run main.nf \
  -profile test,singularity
```

If you want to run the bundled test dataset without Singularity, use the existing `LemonIte` conda environment:

```bash
conda activate LemonIte
nextflow run main.nf \
  -profile conda,test
```

### Standard Run

```bash
nextflow run main.nf \
  --input_dir /path/to/project \
  --organism human \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  --n_clusters 100 \
  --top_n_genes 5000 \
  --coherence_threshold 0.6 \
  -profile singularity \
  -resume
```

## Profiles

| Profile | Purpose |
| --- | --- |
| `singularity` | Recommended default execution profile |
| `hpc` | Cluster resource overrides, usually combined with `singularity` |
| `local` | Local executor settings |
| `docker` | Docker execution profile |
| `test` | Bundled smoke-test dataset |
| `dev` | Development bind mounts for host scripts and PKN resources |

## Inputs and Outputs

Inputs are read from:

```text
{input_dir}/data/
```

Published analysis results are written to:

```text
{output_dir or input_dir/results}/{run_id}/
```

Top-level run outputs include:

- `Lemonite_Summary_Report.html`
- `pipeline_parameters_log.txt`
- `LemonTree/` with preprocessing, network, heatmap, enrichment, evaluation, and overview results

## More Documentation

- [WIKI.md](WIKI.md) for the detailed parameter and output reference
- [../README.md](../README.md) for the repository-level overview and citation
