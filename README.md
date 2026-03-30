# 🍋 Lemonite

**Lemonite: identification of regulatory metabolites through data-driven, interpretable integration of transcriptomics and metabolomics data**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A525.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-available-blue.svg)](https://www.docker.com/)
[![Singularity](https://img.shields.io/badge/singularity-available-blue.svg)](https://sylabs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

---

## Overview

**Lemonite** is a comprehensive framework for multi-omics data integration that identifies gene co-expression modules and their regulators, with a particular focus on discovering regulatory metabolites. Building on the [LemonTree algorithm](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003983), Lemonite extends it with support for multiple regulator types (transcription factors, metabolites, lipids...), enhanced regulator selection and network construction, prioritization of modules and regulators, validation against curated prior knowledge networks (PKNs), functional analyses, and rich visualizations.

🌐 **[Visit www.lemonite.ugent.be](http://www.lemonite.ugent.be)**

### Key Features

- **Multi-omics Integration**: Integrate transcriptomics, metabolomics, lipidomics, and other omics data types
- **Flexible Regulator Assignment**: Support for TFs, metabolites, lipids, and or other omics (binary data currently not supported)
- **Prior Knowledge Networks**: Lemonite knowledge graph for validation of metabolite-gene and protein-protein interactions
- **Rich Visualizations**: Interactive network graphs, module heatmaps, in silico validation with knowledge graph, and enrichment results

---

## 📦 Repository Structure

```
Lemonite/
├── nextflow/             # Automated Nextflow pipeline
│   ├── main.nf           # Main pipeline workflow
│   ├── modules/          # Modular pipeline components
│   ├── scripts/          # Pipeline scripts
│   ├── LemonTree/        # Source code for the LemonTree algorithm
│   ├── PKN/              # Lemonite prior knowledge graph
│   ├── conf/             # Configuration profiles (Docker, Singularity, HPC)
│
├── build_PKN/            # 🕸️ Prior Knowledge Network construction
│   ├── build_PKN/        # pipeline implementation to build the PKN - UNDER DEVELOPMENT!!!
│   ├── Collect_PKNdata_metabolites.ipynb
│   ├── Collect_PKNdata_proteins.ipynb
│   └── Build_final_PKN.ipynb
│
├── Wang_GBM/             # Scripts to recreate results on GBM dataset presented in manuscript
├── Lloyd-Price/          # Scripts to recreate results on IBD dataset presented in manuscript
│
├── web_app/              # Code for www.lemonite.ugent.be (keep on ugent github)
│   ├── Home.py           # Dashboard home page 
│   └── pages/            # Analysis visualization pages
│
└── LICENSE               # MIT License
```

---

## 🚀 Quick Start

### Prerequisites

- **For Nextflow Pipeline:**
  - Nextflow ≥ 25.04.0
  - Docker **OR** Singularity
  - Minimum 16 GB RAM, 4 CPUs

- **For R/Python Scripts:**
  - R ≥ 4.0 with Bioconductor packages
  - Python ≥ 3.8 with scanpy, pandas, numpy etc (see docker file for complete list of R packages and python modules)

### Nextflow Pipeline quick start

```bash
# Clone the repository
git clone https://github.com/CBIGR/Lemonite.git
cd Lemonite

# For Nextflow pipeline
cd nextflow
curl -s https://get.nextflow.io | bash # install nextflow if necessary
sudo mv nextflow /usr/local/bin/

# Build container (Docker or Singularity)
./build-docker.sh      # For Docker
# OR
./build-singularity.sh # For Singularity/HPC

# Run pipeline
nextflow run main.nf \
  --input_dir /path/to/data \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \ # You can add other omics to assign as regulators here
  --n_clusters 100 \
  -profile docker # or singularity
```

#### 🔄 Pipeline Workflow

```
📊 Input: Transcriptomics + Metabolomics/Lipidomics
         ↓
    1. PREPROCESSING & TFA
       • Normalization & log-transformation
       • Omics-specific scaling
       • TF activity inference (optional)
       • Highly variable gene selection
         ↓
    2. LEMONTREE CLUSTERING
       • Parallel Gibbs sampling (10-100 runs)
       • Module-based co-expression clustering
       • Regulator assignment (TFs, metabolites, etc.)
         ↓
    3. NETWORK GENERATION & FILTERING
       • Module coherence filtering
       • Regulator selection
       • Network construction
         ↓
    4. DOWNSTREAM ANALYSES
       ├─→ PKN validation
       ├─→ Subnetwork visualization
       ├─→ Module heatmaps
       └─→ Enrichment analysis
         ↓
    5. MODULE OVERVIEW
       • Interactive network visualization
       • Module prioritization
```


#### Input Data Requirements

Your input directory should contain:

```
input_directory/
└── data/
    ├── Counts.tsv               # Gene expression raw counts (REQUIRED)
    ├── Metabolomics.txt         # Metabolomics data (OPTIONAL but recommended)
    ├── Other_omics.txt          # Additional omics data for regulator assignment (OPTIONAL)
    ├── Metadata.txt             # Sample metadata (REQUIRED)
    ├── lovering_TF_list.txt     # TF list (REQUIRED, can be copied from this test data directory)
    └── name_map.csv             # Metabolite ID mapping (REQUIRED if using metabolites)
```

**File Format Requirements:**

| File | Format | Key Requirements |
|------|--------|------------------|
| **Counts.tsv** | Tab-separated | Rows: genes (HGNC), Columns: samples, Values: raw counts |
| **Metabolomics.txt** | Tab-separated | Rows: metabolites, Columns: samples (matching Counts) |
| **Metadata.txt** | Tab-separated | Must include `Sample_ID` and `diagnosis` columns |
| **lovering_TF_list.txt** | Plain text | One TF per line (HGNC symbols) |
| **name_map.csv** | CSV | Columns: `Query`, `HMDB` (metabolite name to HMDB ID) |

#### Output Structure

```
results/{run_id}/LemonTree/
├── Preprocessing/          # Normalized data & QC plots
├── Lemon_out/              # Raw LemonTree clustering results
├── Networks/               # Regulatory networks & subnetwork graphs
├── module_heatmaps/        # Expression heatmaps per module
├── Enrichment/             # Pathway enrichment results
├── PKN_evaluation/         # Prior knowledge validation
└── Module_Overview/        # Interactive overview & summary tables
```

**Key Output Files:**
- **`Networks/LemonNetwork_*.txt`**: Complete regulatory network
- **`Networks/subnetworks/graph_*.png`**: Per-module network visualizations with Lemonite knowledge graph
- **`Module_Overview/interactive_module_network.html`**: Interactive module overview figure
- **`Module_Overview/Module_Overview.csv`**: Master results table
- **`Enrichment/Modules_enrichr/`**: Functional enrichment per module

#### Advanced Usage

📖 **[Complete pipeline documentation & advanced usage →](https://github.com/CBIGR/Lemonite/wiki)**

---

## Reference

If you use Lemonite in your research, please cite:

```bibtex
@article{vandemoortele2025lemonite,
  title={Lemonite: identification of regulatory metabolites through data-driven, interpretable integration of transcriptomics and metabolomics data},
  author={Vandemoortele, Boris et al.},
  journal={[bioRxiv]},
  year={2026},
  note={https://doi.org/10.64898/2026.03.27.714373}
}
```



## Support & Contact

- **Report Issues**: [GitHub Issues](https://github.com/CBIGR/Lemonite/issues)

---

## Acknowledgments

- **Lemon-Tree**: Original algorithm by Eric Bonnet & Tom Michoel: doi:10.1371/journal.pcbi.1003983
- **Data Sources**: Lloyd-Price et al. (IBD), Wang et al. (GBM)
- **PKN Databases**: IntAct, STITCH, BioGRID, STRING, UniProtKB, ChEMBL, LINCS, MetalinksDB, Human1-GEM

---

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2026 CBIGR - Lab Vanessa Vermeirssen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction...
```

---




