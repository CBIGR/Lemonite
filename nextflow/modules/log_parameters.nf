process LOG_PARAMETERS {
    tag "Logging pipeline parameters"
    publishDir "${(params.output_dir ?: (params.input_dir ? params.input_dir + '/results' : './results'))}/${(params.computed_run_id ?: (params.run_id ?: 'run_auto'))}", mode: 'copy'

    output:
    path "pipeline_parameters_log.txt", emit: parameters_log

    script:
    def timestamp = new Date().format('yyyy-MM-dd HH:mm:ss')
    """
    cat > pipeline_parameters_log.txt << 'EOF'
================================================================================
                    LEMONITE PIPELINE PARAMETERS LOG
================================================================================
Run ID: ${params.computed_run_id}
Timestamp: ${timestamp}
Nextflow Version: ${nextflow.version}
Profile: ${workflow.profile}
================================================================================

INPUT/OUTPUT PARAMETERS:
------------------------
input_dir                     = ${params.input_dir}
output_dir                    = ${params.output_dir}
workDir                       = ${params.output_dir}/work
run_id                        = ${params.run_id ?: '(auto-generated)'}

PREPROCESSING PARAMETERS:
-------------------------
top_n_genes                   = ${params.top_n_genes}
perform_tfa                   = ${params.perform_tfa}
use_omics_specific_scaling    = ${params.use_omics_specific_scaling}
expression_col                = ${params.expression_col}
sample_id_col                 = ${params.sample_id_col}
metadata_columns              = ${params.metadata_columns}
design_formula                = ${params.design_formula}
deseq_contrast1               = ${params.deseq_contrast1}

CLUSTERING PARAMETERS:
----------------------
n_clusters                    = ${params.n_clusters}
coherence_threshold           = ${params.coherence_threshold}
use_deseq_priors              = ${params.use_deseq_priors}
min_cluster_size              = ${params.min_cluster_size}
tight_clusters_only           = ${params.tight_clusters_only}
max_n_iterations              = ${params.max_n_iterations}

REGULATOR PARAMETERS:
---------------------
regulator_types               = ${params.regulator_types}
regulator_selection_method    = ${params.regulator_selection_method}
top_n_percent_regulators      = ${params.top_n_percent_regulators}
regulator_fold_cutoff         = ${params.regulator_fold_cutoff}

NETWORK GENERATION PARAMETERS:
------------------------------
min_regulator_size            = ${params.min_regulator_size}
max_regulator_size            = ${params.max_regulator_size}
min_module_size               = ${params.min_module_size}
min_targets                   = ${params.min_targets}
min_expression_fold_threshold = ${params.min_expression_fold_threshold}
max_pvalue_threshold          = ${params.max_pvalue_threshold}

PRIOR KNOWLEDGE NETWORK (PKN):
------------------------------
pkn_network                   = ${params.pkn_network}

ENRICHMENT PARAMETERS:
----------------------
use_megago                    = ${params.use_megago}
enrichment_method             = ${params.enrichment_method}
enrichr_libraries             = ${params.enrichr_libraries}

OVERVIEW PARAMETERS:
--------------------
prioritize_by_expression      = ${params.prioritize_by_expression}
overview_n_clusters           = ${params.overview_n_clusters}

RESOURCE LIMITS:
----------------
max_cpus                      = ${params.max_cpus}
max_memory                    = ${params.max_memory}
max_time                      = ${params.max_time}

================================================================================
COMMAND LINE:
${workflow.commandLine}
================================================================================
EOF
    """

    stub:
    """
    touch pipeline_parameters_log.txt
    """
}
