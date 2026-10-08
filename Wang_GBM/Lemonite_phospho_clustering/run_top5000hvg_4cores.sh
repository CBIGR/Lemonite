#!/bin/bash

cd /home/borisvdm/repo/LemonIte/nextflow

nextflow run main.nf \
  --input_dir /home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/phospho_clustering_nf/input \
  --output_dir /home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/phospho_clustering_nf/results \
  --preprocessing_type proteomics \
  --expression_file data/phospho_expression_full.tsv \
  --metadata_file data/metadata.tsv \
  --sample_id_col Sample_ID \
  --organism human \
  --perform_tfa false \
  --top_n_genes 5000 \
  --regulator_types "Metabolites:metabolome.csv,Proteins:proteins_raw.tsv,Lipids:lipids_raw.tsv,TFs:tf_activity_raw_top5000hvg.tsv,KinaseActivity:kinase_activity_raw_top5000hvg.tsv" \
  --n_clusters 100 \
  --stop_after_network \
  -profile singularity,dev,local \
  -c /tmp/claude-2544/-home-borisvdm-repo-LemonIte-Wang-GBM-Lemonite-phospho-clustering/09388d11-fc80-4a6e-80a9-e77b7e6740f5/scratchpad/nf_6cores_extended.config \
  -work-dir /home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/phospho_clustering_nf/work_top5000hvg \
  -resume 2>&1 | tee -a /home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/phospho_clustering_nf/rerun_top5000hvg_fix_v3.log
