# Lemonite

**Uncovering regulatory metabolites through interpretable, data-driven multi-omics integration**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Singularity](https://img.shields.io/badge/singularity-available-blue.svg)](https://sylabs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

## Overview

Lemonite is a framework for multi-omics data integration that identifies gene co-expression modules and their regulators, with a particular focus on regulatory metabolites. It builds on the [LemonTree algorithm](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983) and extends it with pipeline automation, support for additional omics and regulator types (proteomics, phosphoproteomics, epigenomics, microbiomics, copy number variation data...), extensive in silico validation using an in-house consctructed knowledge-graph, enrichment analyses, rich visualizations and interactive summaries.

Website: [www.lemonite.ugent.be](http://www.lemonite.ugent.be)

## Key Features

- Integrates transcriptomics with metabolomics and additional omics layers.
- Supports multiple regulator types, including continuous and discrete/binary regulators via `Prefix:File[:c|d]` entries in `--regulator_types`.
- Supports both human and mouse analyses through species-aware TF lists, annotation, and enrichment settings.
- Produces regulator-target networks, knowledge-graph-based in silico validation, module heatmaps, enrichment results, and an interactive module overview.
- Ships as a reproducible Nextflow pipeline with `singularity` profile.

## Repository Structure

```text
Lemonite/
├── nextflow/             # Nextflow pipeline and bundled resources
│   ├── main.nf           # Main pipeline workflow
│   ├── nextflow.config   # Default parameters and profiles
│   ├── conf/             # Base, singularity, local, hpc, and test configs
│   ├── modules/          # Workflow modules
│   ├── scripts/          # R, Python, and shell entrypoints used by modules
│   ├── PKN/              # Bundled prior-knowledge resources and default TF lists
│   ├── test_dataset/     # Minimal bundled dataset for tests
│   └── WIKI.md           # Detailed pipeline documentation
├── build_PKN/            # Prior-knowledge-network construction notebooks and pipeline
├── Lloyd-Price/          # Analysis scripts for the IBD use case
├── Wang_GBM/             # Analysis scripts for the GBM use case
└── LICENSE               # MIT License
```

## Nextflow Pipeline

### Requirements

- Nextflow `>= 23.04.0`
- Java 11+
- Singularity (required — the only supported execution backend; use the `singularity` command)
- Minimum 16 GB RAM and 4 CPU cores

Available profiles: `singularity` (required), `hpc` (cluster resource overrides), `local` (reduced defaults), `test` (bundled smoke-test dataset), `dev` (bind-mount host scripts for rapid iteration). All profiles must be combined with `singularity`.

### Quick Test

```bash
git clone https://github.com/CBIGR/Lemonite.git
cd Lemonite/nextflow

./build-singularity.sh

nextflow run main.nf \
  -profile test,singularity
```

The `test` profile sets `input_dir` to `nextflow/test_dataset/`, reduces the cluster count to 5 and the gene count to 1000, and publishes results under `nextflow/test_dataset/results/{run_id}/`.

### Run on Your Own Data

```bash
cd Lemonite/nextflow

nextflow run main.nf \
  --input_dir /path/to/project \
  --organism human \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  --regulator_selection_method fold_per_module \
  --regulator_fold_cutoff 2.0 \
  -profile singularity \
  -resume
```

`--input_dir` must contain a `data/` directory. By default the pipeline auto-detects the expression matrix and metadata file from that directory, then publishes results into:

```text
{output_dir or input_dir/results}/{run_id}/
```

If `--run_id` is omitted, the pipeline generates one from the main analysis settings.

### Input Directory

```text
input_dir/
└── data/
    ├── Counts.tsv or *counts*.tsv      # required, unless --expression_file is provided
    ├── Metadata.txt or *metadata*.txt  # required, unless --metadata_file is provided
    ├── Metabolomics.txt                # optional regulator abundance file
    ├── Lipidomics.txt / Proteomics.txt # optional additional regulator files
    ├── name_map.csv                    # maps metabolite names to HMDB IDs, required for validation with knowledge graph
    └── *network*.txt                   # optional custom TF prior network override
```

Notes:

- If you keep the default TF regulator entry `TFs:Lovering_TF_list.txt`, Lemonite checks `data/` first and then falls back to the bundled TF list in `nextflow/PKN/`.
- Custom regulator files are referenced through `--regulator_types`; discrete/binary regulators are marked with `:d`.
- If you place a custom `*network*.txt` or `*CollecTRI*.txt` file in `data/`, the preprocessing step uses it instead of the bundled CollecTRI prior network.

### Published Outputs

```text
{run_dir}/
├── Lemonite_Summary_Report.html
├── pipeline_parameters_log.txt
└── LemonTree/
    ├── Preprocessing/
    ├── Lemon_out/
    ├── Networks/
    ├── ModuleViewer_files/
    ├── PKN_Evaluation/
    ├── module_heatmaps/
    ├── Enrichment/
    └── Module_Overview/
```

`Lemonite_Summary_Report.html` is the top-level run summary. The detailed output layout, profiles, and parameter reference are documented in [nextflow/WIKI.md](nextflow/WIKI.md).

## Citation

If you use Lemonite in your research, cite the bioRxiv preprint:

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

Original LemonTree algorithm:

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

## Support

- Issues: [GitHub Issues](https://github.com/CBIGR/LemonIte/issues)
- Contact: boris.vandemoortele@ugent.be
- Lab: [Vermeirssen Lab](https://www.crig.ugent.be/en/prof-vanessa-vermeirssen-phd)

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
