#!/bin/bash -euo pipefail
cat > pipeline_parameters_log.txt << 'EOF'
================================================================================
                    LEMONITE PIPELINE PARAMETERS LOG
================================================================================
Run ID: minW0p15
Timestamp: 2026-05-29 11:42:39
Nextflow Version: 25.04.6
Profile: standard
================================================================================

INPUT/OUTPUT PARAMETERS:
------------------------
input_dir                     = /home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens
output_dir                    = /home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics_nextflow
workDir                       = /home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics_nextflow/work
run_id                        = minW0p15

PREPROCESSING PARAMETERS:
-------------------------
top_n_genes                   = 1000
perform_tfa                   = false
use_omics_specific_scaling    = true
expression_col                = count
sample_id_col                 = ID
metadata_columns              = Cell_line,FAB,MITO_subtype,Risk_classification
design_formula                = ~ Cell_line
deseq_contrast1               = Cell_line

CLUSTERING PARAMETERS:
----------------------
n_clusters                    = 5
coherence_threshold           = 0.5
use_deseq_priors              = true
min_cluster_size              = 10
tight_clusters_only           = false
max_n_iterations              = 1000

REGULATOR PARAMETERS:
---------------------
regulator_types               = hPTMs:Histone_proteomics.csv,Metabolites:metabolomics_combined.csv
regulator_selection_method    = fold_per_module
top_n_percent_regulators      = 2.0
regulator_fold_cutoff         = 2.0

NETWORK GENERATION PARAMETERS:
------------------------------
min_regulator_size            = 3
max_regulator_size            = 100
min_module_size               = 10
min_targets                   = 3
min_expression_fold_threshold = 1.5
max_pvalue_threshold          = 0.05

PRIOR KNOWLEDGE NETWORK (PKN):
------------------------------
pkn_network                   = /home/borisvdm/repo/LemonIte/nextflow/PKN/Lemonite_PKN.tsv

ENRICHMENT PARAMETERS:
----------------------
enrichment_method             = EnrichR
enrichr_libraries             = GO_Biological_Process_2025,GO_Molecular_Function_2025,GO_Cellular_Component_2025,KEGG_2021_Human,Reactome_Pathways_2024

OVERVIEW PARAMETERS:
--------------------
prioritize_by_expression      = true
overview_n_clusters           = 5
overview_labeling_workflow    = canonical MegaGO top_30 + rrvgo

RESOURCE LIMITS:
----------------
max_cpus                      = 5
max_memory                    = 32.GB
max_time                      = 12.h

================================================================================
COMMAND LINE:
nextflow run main.nf -c /home/borisvdm/repo/AML_integromics_BVM/Lemonite/Multiomics/aml_proteomics.config --output_dir /home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics_nextflow --outdir /home/borisvdm/Documents/PhD/24_AML_Dhaenens-Lammens/results/LemonTree/Proteomics_nextflow --run_id minW0p15 --lemontree_tight_min_weight 0.15 --max_cpus 5 -resume
================================================================================
EOF
