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

```bash
# Build container
./build-docker.sh      # For Docker
# OR
./build-singularity.sh # For Singularity/HPC

# Run pipeline
nextflow run main.nf \
  --input_dir /path/to/data \
  --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
  --n_clusters 100 \
  -profile docker
```

---

## 📧 Support

- 🐛 **Issues**: [GitHub Issues](https://github.com/CBIGR/Lemonite/issues)
- 📧 **Email**: boris.vandemoortele@ugent.be
- 📖 **Documentation**: [Wiki](https://github.com/CBIGR/Lemonite/wiki)

