# 🍋🌳 Lemonite Nextflow Pipeline

This directory contains the automated Nextflow pipeline for Lemonite multi-omics integration.

---

## 📖 Documentation

**Complete pipeline documentation, including:**
- Installation instructions
- Input data format specifications
- Configuration options and parameters
- Output file descriptions
- HPC/cluster execution guides
- Troubleshooting and FAQ

**→ [Visit the Lemonite Wiki for more extensive documentation](https://github.com/CBIGR/Lemonite/wiki)**

---

## 🚀 Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) ≥ 23.04.0
- [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (recommended container runtime)

### Build & Run

```bash
# Build Singularity container
./build-singularity.sh

# Run pipeline
nextflow run main.nf \
  --input_dir /path/to/data \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  --n_clusters 100 \
  -profile singularity
```

### Profiles

| Profile | Description |
|---------|-------------|
| `singularity` | Run with Singularity containers (recommended) |
| `hpc` | HPC-optimized resources (use with `singularity`, e.g. `-profile singularity,hpc`) |
| `dev` | Bind-mount host scripts/PKN for development (use with `singularity`) |
| `test` | Minimal test run with bundled test dataset |
| `local` | Local execution settings |

### Run Test Dataset

```bash
nextflow run main.nf -profile singularity,test
```

---

## 📧 Support

- 🐛 **Issues**: [GitHub Issues](https://github.com/CBIGR/Lemonite/issues)
- 📧 **Email**: boris.vandemoortele@ugent.be
- 📖 **Documentation**: [Wiki](https://github.com/CBIGR/Lemonite/wiki)

