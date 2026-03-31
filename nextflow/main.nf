#!/usr/bin/env nextflow



nextflow.enable.dsl = 2

// Default parameters (minimal set, most defaults are in nextflow.config)
params.input_dir = ""
params.run_id = "" // Will be auto-generated based on parameters if not specified
params.perform_tfa = true  // Enable TFA analysis
params.use_megago = true  // Use megago for functional clustering
params.prioritize_by_expression = true  // Enable expression-based prioritization
params.overview_n_clusters = 5  // Number of functional clusters
params.expression_col = "count"  // Expression column name
params.sample_id_col = "Sample_ID"  // Sample ID column name
params.max_cpus = 16  // Maximum CPUs
params.max_memory = "64.GB"  // Maximum memory
params.max_time = "24.h"  // Maximum time
// Organism parameter controls species-specific annotation and enrichment (human|mouse)
params.organism = 'human'

// Compute the final run_id - either user provided or auto-generated
def method_suffix = params.regulator_selection_method == "percentage" ? "top${params.top_n_percent_regulators}pct" :
                    "fold${params.regulator_fold_cutoff}x"
def auto_run_id = "${params.top_n_genes}HVG_coherence${params.coherence_threshold}_${method_suffix}_clusters${params.n_clusters}"
params.computed_run_id = (params.run_id && params.run_id != "") ? params.run_id : auto_run_id

// Input validation
if (!params.input_dir) {
    error "Please provide an input directory with --input_dir"
}

// Validate key parameters
def valid_organisms = ['human', 'mouse']
if (!(params.organism in valid_organisms)) {
    error "Invalid --organism '${params.organism}'. Must be one of: ${valid_organisms.join(', ')}"
}

def valid_enrichment_methods = ['EnrichR', 'GSEA', 'both', 'auto']
if (!(params.enrichment_method in valid_enrichment_methods)) {
    error "Invalid --enrichment_method '${params.enrichment_method}'. Must be one of: ${valid_enrichment_methods.join(', ')}"
}

def valid_regulator_methods = ['percentage', 'fold_per_module']
if (!(params.regulator_selection_method in valid_regulator_methods)) {
    error "Invalid --regulator_selection_method '${params.regulator_selection_method}'. Must be one of: ${valid_regulator_methods.join(', ')}"
}

if (params.n_clusters < 1) {
    error "Invalid --n_clusters '${params.n_clusters}'. Must be a positive integer."
}

if (params.coherence_threshold < 0 || params.coherence_threshold > 1) {
    error "Invalid --coherence_threshold '${params.coherence_threshold}'. Must be between 0 and 1."
}

// Compute output/work directories (use def variables, not params mutation)
def finalOutputDir = params.output_dir ?: "${params.input_dir}/results"
def finalWorkDir = params.work_dir ?: "${params.input_dir}/work"

// Prepare display variables (computed after defaults are set)
def displayInputDir = params.input_dir.toString()
def displayOutputDir = finalOutputDir.toString()
def displayWorkDir = finalWorkDir.toString()

log.info """\

    ╔════════════════════════════════════════════════════════════════════════════════════════════════╗
    ║                                                                                                ║
    ║                L E M O N I T E      A N A L Y S I S       P I P E L I N E                      ║
    ║                                                                                                ║
    ║            ██╗     ███████╗███╗   ███╗ ██████╗ ███╗   ██╗██╗████████╗███████╗                  ║
    ║            ██║     ██╔════╝████╗ ████║██╔═══██╗████╗  ██║██║╚══██╔══╝██╔════╝                  ║
    ║            ██║     █████╗  ██╔████╔██║██║   ██║██╔██╗ ██║██║   ██║   █████╗                    ║
    ║            ██║     ██╔══╝  ██║╚██╔╝██║██║   ██║██║╚██╗██║██║   ██║   ██╔══╝                    ║
    ║            ███████╗███████╗██║ ╚═╝ ██║╚██████╔╝██║ ╚████║██║   ██║   ███████╗                  ║
    ║            ╚══════╝╚══════╝╚═╝     ╚═╝ ╚═════╝ ╚═╝  ╚═══╝╚═╝   ╚═╝   ╚══════╝                  ║
    ║                                                                                                ║
    ║            🧬 Advanced Multi-Omics Integration & Network Analysis Platform 🧬                  ║
    ║                                                                                                ║
    ╠════════════════════════════════════════════════════════════════════════════════════════════════╣
    ║                                                                                                ║
    ║  PIPELINE CONFIGURATION                                                                        ║
    ║  ┌──────────────────────────────────────────────────────────────────────────────────────────┐  ║
    ║  │  Input Directory        : ${displayInputDir.padRight(62)} │  ║
    ║  │  Output Directory       : ${displayOutputDir.padRight(62)} │  ║
    ║  │  Work Directory         : ${displayWorkDir.padRight(62)} │  ║
    ║  │  Run Identifier         : ${params.computed_run_id.padRight(62)} │  ║
    ║  │  Number of Clusters     : ${String.valueOf(params.n_clusters).padRight(62)} │  ║
    ║  │  Coherence Threshold    : ${String.valueOf(params.coherence_threshold).padRight(62)} │  ║
    ║  │  Top N% Regulators      : ${String.valueOf(params.top_n_percent_regulators).padRight(62)} │  ║
    ║  │  Top Variable Genes     : ${String.valueOf(params.top_n_genes).padRight(62)} │  ║
    ║  │  Enrichment Method      : ${String.valueOf(params.enrichment_method).padRight(62)} │  ║
    ║  │  Organism               : ${String.valueOf(params.organism).padRight(62)} │  ║
    ║  └──────────────────────────────────────────────────────────────────────────────────────────┘  ║
    ║                                                                                                ║
    ║  PIPELINE WORKFLOW STEPS                                                                       ║
    ║  ┌──────────────────────────────────────────────────────────────────────────────────────────┐  ║
    ║  │  1. Preprocessing & TFA Analysis                                                         │  ║
    ║  │  2. LemonTree Clustering - ${("${params.n_clusters} module-based Gibbs samplers").padRight(61)} │  ║
    ║  │  3. LemonTree regulator assignment                                                       │  ║
    ║  │  4. Module filtering, regulator selection and network generation                         │  ║
    ║  │  5. Comparison Lemonite data-driven network with Lemonite PKN                            │  ║
    ║  │  6. Create module heatmaps                                                               │  ║
    ║  │  7. Pathway Enrichment Analysis on module genes and regulator targets                    │  ║
    ║  │  8. Differential Expression Analysis & build module overview                             │  ║
    ║  └──────────────────────────────────────────────────────────────────────────────────────────┘  ║
    ║                                                                                                ║
    ║                                                                                                ║
    ║  💡 Pipeline developed by Vandemoortele Boris, CBIGR lab @ Ghent university                    ║
    ║                                                                                                ║
    ╚════════════════════════════════════════════════════════════════════════════════════════════════╝


    
    🌟 Initializing Lemonite Analysis Pipeline... 🌟
    ⏰ Start Time: ${new Date().format('yyyy-MM-dd HH:mm:ss')}

    ☕️ This will take a while, go grab a coffee ☕️

    """
    .stripIndent()

// Include process definitions
include { LOG_PARAMETERS } from './modules/log_parameters.nf'
include { PREPROCESSING_TFA } from './modules/preprocessing.nf'
include { CLUSTERING } from './modules/clustering.nf' 
include { POST_CLUSTERING; NETWORK_GENERATION; SUBNETWORK_GRAPHS } from './modules/network.nf'
include { PKN_EVALUATION } from './modules/evaluation.nf'
include { MODULE_VIEWER_HEATMAPS } from './modules/viewer_heatmaps.nf'
include { ENRICHMENT_ANALYSIS } from './modules/enrichment.nf'
include { MODULE_OVERVIEW_INTERACTIVE } from './modules/overview.nf'
include { SUMMARY_REPORT } from './modules/summary_report.nf'

workflow {
    // Create channels for input data
    input_ch = Channel.fromPath(params.input_dir, type: 'dir')
    
    // data_dir is always input_dir/data
    data_dir = "${params.input_dir}/data"
    data_ch = Channel.fromPath(data_dir, type: 'dir')
    name_map_ch = file("${data_dir}/name_map.csv").exists() ? Channel.fromPath("${data_dir}/name_map.csv") : Channel.empty()
    
    // PKN file for edge categorization in module overview
    pkn_ch = file(params.pkn_network).exists() ? Channel.fromPath(params.pkn_network) : Channel.empty()
    
    // Create channel for run_id
    run_id_ch = Channel.value(params.computed_run_id)
    
    // Log all parameters at the start of the pipeline
    LOG_PARAMETERS()
    
    // Step 1: Preprocessing and TFA
    PREPROCESSING_TFA(input_ch, data_ch, run_id_ch)
    
    // Step 2: Lemonite clustering (parallel jobs)
    cluster_ids = Channel.of(1..params.n_clusters)
    CLUSTERING(PREPROCESSING_TFA.out.preprocessed_data.combine(cluster_ids), run_id_ch)
    
    // Step 3a: Post-clustering analysis
    clustered_data = CLUSTERING.out.clusters.collect()
    POST_CLUSTERING(
        PREPROCESSING_TFA.out.lemontree_inputs,
        clustered_data,
        run_id_ch
    )
    
    // Step 3b: Network generation and filtering
    NETWORK_GENERATION(
        POST_CLUSTERING.out.lemontree_outputs,
        POST_CLUSTERING.out.preprocessing_files,
        run_id_ch
    )
    
    // Step 3c: Subnetwork visualization graphs
    SUBNETWORK_GRAPHS(
        NETWORK_GENERATION.out.viewer_files,
        POST_CLUSTERING.out.preprocessing_files,
        run_id_ch
    )
    
    // Step 4: PKN evaluation
    PKN_EVALUATION(
        NETWORK_GENERATION.out.all_regulator_targets.mix(NETWORK_GENERATION.out.main_network).collect(),
        NETWORK_GENERATION.out.viewer_files,
        data_ch,
        run_id_ch,
        params.regulator_types
    )
    
    // Step 5: Module viewer heatmaps (after PKN evaluation)
    MODULE_VIEWER_HEATMAPS(
        PKN_EVALUATION.out.viewer_files.ifEmpty(NETWORK_GENERATION.out.viewer_files),
        NETWORK_GENERATION.out.filtered_modules,
        input_ch,
        PREPROCESSING_TFA.out.preprocessed_data,   // expression file
        PREPROCESSING_TFA.out.complete_data,
        PREPROCESSING_TFA.out.omics_preprocessed.collect(),   // omics-specific preprocessed files
        run_id_ch,
        params.regulator_types
    )
    
    // Step 6: Enrichment analysis
    ENRICHMENT_ANALYSIS(
        NETWORK_GENERATION.out.filtered_modules,
        NETWORK_GENERATION.out.all_regulator_targets,
        NETWORK_GENERATION.out.network_files,
        PKN_EVALUATION.out.viewer_files.ifEmpty(NETWORK_GENERATION.out.viewer_files),
        PREPROCESSING_TFA.out.preprocessed_data,
        run_id_ch
    )
    
    // Step 7: Module overview generation (runs after enrichment analysis)
    MODULE_OVERVIEW_INTERACTIVE(
        PKN_EVALUATION.out.viewer_files.ifEmpty(NETWORK_GENERATION.out.viewer_files),
        NETWORK_GENERATION.out.filtered_modules,
        NETWORK_GENERATION.out.coherence_scores,
        ENRICHMENT_ANALYSIS.out.enrichment_results,
        PREPROCESSING_TFA.out.preprocessed_data,
        PREPROCESSING_TFA.out.metadata,
        pkn_ch.ifEmpty([]),
        name_map_ch.ifEmpty([]),
        run_id_ch,
        params.regulator_types
    )
    
    // Step 8: Generate comprehensive HTML summary report
    SUMMARY_REPORT(
        PKN_EVALUATION.out.viewer_files.ifEmpty(NETWORK_GENERATION.out.viewer_files),
        NETWORK_GENERATION.out.filtered_modules,
        NETWORK_GENERATION.out.network_files,
        NETWORK_GENERATION.out.coherence_scores,
        ENRICHMENT_ANALYSIS.out.enrichment_results,
        PKN_EVALUATION.out.pkn_results.ifEmpty([]),
        PREPROCESSING_TFA.out.preprocessing_results,
        CLUSTERING.out.clusters.collect(),
        MODULE_OVERVIEW_INTERACTIVE.out.overview_results,
        LOG_PARAMETERS.out.parameters_log,
        run_id_ch,
        params.regulator_types
    )
}

workflow.onComplete {
    if (workflow.success) {
        log.info """\
        
        ╔══════════════════════════════════════════════════════════════════════════════════╗
        ║                                                                                  ║
        ║         🎉✨ L E M O N I T E   P I P E L I N E   S U C C E S S ! ✨🎉            ║
        ║                                                                                  ║
        ║        ███████╗██╗   ██╗ ██████╗ ██████╗███████╗███████╗███████╗██╗              ║
        ║        ██╔════╝██║   ██║██╔════╝██╔════╝██╔════╝██╔════╝██╔════╝██║              ║
        ║        ███████╗██║   ██║██║     ██║     █████╗  ███████╗███████╗██║              ║
        ║        ╚════██║██║   ██║██║     ██║     ██╔══╝  ╚════██║╚════██║╚═╝              ║
        ║        ███████║╚██████╔╝╚██████╗╚██████╗███████╗███████║███████║██╗              ║
        ║        ╚══════╝ ╚═════╝  ╚═════╝ ╚═════╝╚══════╝╚══════╝╚══════╝╚═╝              ║
        ║                                                                                  ║
        ║          🌟 Your multi-omics analysis has completed successfully! 🌟             ║
        ║                                                                                  ║
        ╠══════════════════════════════════════════════════════════════════════════════════╣
        ║                                                                                  ║
        ║  🍋🌳 PIPELINE SUMMARY                                                           ║
        ║  ┌────────────────────────────────────────────────────────────────────────────┐  ║
        ║  │  ✅ All processes finished without errors                                  │  ║
        ║  │  ⏰ Start Time: ${workflow.start.toString().padRight(62)} │  ║
        ║  │  🏁 End Time: ${workflow.complete.toString().padRight(62)} │  ║
        ║  │  ⌛ Duration: ${workflow.duration.toString().padRight(62)} │  ║
        ║  │  📁 Results: ${(params.output_dir ?: params.outdir).toString().padRight(62)} │  ║
        ║  └────────────────────────────────────────────────────────────────────────────┘  ║
        ║                                                                                  ║
        ║  🎯 Your results include:                                                        ║
        ║  • 🧩 Gene modules with TF and metabolite regulators                             ║
        ║  • 🔗 Metabolite-gene regulatory networks with PKN evaluation                    ║
        ║  • 📊 Enrichment analysis and functional annotations                             ║
        ║  • 🎨 Interactive network visualizations and heatmaps                            ║
        ║  • 📋 Comprehensive module overview and statistics                               ║
        ║                                                                                  ║
        ║  🌟          Thank you for using Lemonite! Happy analyzing! 🌟                   ║
        ║  🧠  Remember: don't just run pipelines, understand what's happening! 🧠         ║
        ║                                                                                  ║
        ╚══════════════════════════════════════════════════════════════════════════════════╝
        
        """
    } else {
        log.info """\
        
        ╔══════════════════════════════════════════════════════════════════════════════════╗
        ║                                                                                  ║
        ║      ❌💥 L E M O N I T E   P I P E L I N E   F A I L E D ! 💥❌                 ║
        ║                                                                                  ║
        ║         ███████╗ █████╗ ██╗██╗     ███████╗██████╗ ██╗                           ║
        ║         ██╔════╝██╔══██╗██║██║     ██╔════╝██╔══██╗██║                           ║
        ║         █████╗  ███████║██║██║     █████╗  ██║  ██║██║                           ║
        ║         ██╔══╝  ██╔══██║██║██║     ██╔══╝  ██║  ██║╚═╝                           ║
        ║         ██║     ██║  ██║██║███████╗███████╗██████╔╝██╗                           ║
        ║         ╚═╝     ╚═╝  ╚═╝╚═╝╚══════╝╚══════╝╚═════╝ ╚═╝                           ║
        ║                                                                                  ║
        ║              😞 Something went wrong during the analysis...                      ║
        ║                                                                                  ║
        ╠══════════════════════════════════════════════════════════════════════════════════╣
        ║                                                                                  ║
        ║  🔍 TROUBLESHOOTING TIPS                                                         ║
        ║  ┌────────────────────────────────────────────────────────────────────────────┐  ║
        ║  │  1. 📋 Check the .nextflow.log file for detailed error messages            │  ║
        ║  │  2. 🔍 Look at work directories for intermediate results                   │  ║
        ║  │  3. 📁 Verify input data format and file permissions                       │  ║
        ║  │  4. 🐳 Ensure Docker/Singularity containers are working properly           │  ║
        ║  │  5. 💾 Check available disk space and memory                               │  ║
        ║  │  6. 📞 Contact support with error details if needed                        │  ║
        ║  └────────────────────────────────────────────────────────────────────────────┘  ║
        ║                                                                                  ║
        ║  💪 Don't give up! Science requires persistence! 💪                              ║
        ║                                                                                  ║
        ╚══════════════════════════════════════════════════════════════════════════════════╝
        
        """
    }
}
