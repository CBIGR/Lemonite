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

> **MegaGO clustering:** the module overview stage (step 9) normally clusters modules
> with MegaGO + rrvgo. This is skipped by default under `-profile test`
> (`--skip_megago true`) since it's slow and not needed to check pipeline function —
> pass `--skip_megago false` to exercise it on the test dataset. See
> [WIKI.md](WIKI.md#enrichment-and-overview-parameters) for details.

> **Building on HPC (e.g. UGent kyukon):** these nodes run an unprivileged
> Apptainer without a `subuid` mapping, so `apptainer build` cannot create a
> `.sif` (it falls back to `proot`, whose `ptrace` is blocked → `ptrace(TRACEME):
> Operation not permitted`). Use `./build-singularity_parallel.sh` instead, which
> builds a **sandbox directory** (`lemontree-pipeline_v1.0.0_parallel.sandbox`)
> and runs from it. Apptainer/Nextflow execute a sandbox directory exactly like a
> `.sif`; point `process.container` at the sandbox path. See
> [WIKI.md](WIKI.md#building-on-hpc-sandbox-containers) for details.

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
