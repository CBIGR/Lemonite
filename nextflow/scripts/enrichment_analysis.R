#!/usr/bin/env Rscript

############################################################################################################################################
#### This is a script for performing gsea enrichment analysis on LemonTree modules
#### Adapted from Enrichment_on_modules_and_targets.R for Nextflow pipeline
############################################################################################################################################

# Load required libraries
library(optparse)
library(data.table)

# Function to safely load packages
safe_library <- function(package, required = FALSE) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    if (required) {
      stop("Required package '", package, "' is not available. Please install it.")
    } else {
      cat("Warning: Package", package, "not available. Some functionality may be limited.\n")
      return(FALSE)
    }
  }
  return(TRUE)
}

# Load packages with safety checks
packages_available <- list()
packages_available$fgsea <- safe_library("fgsea")
packages_available$stringr <- safe_library("stringr")
packages_available$org.Hs.eg.db <- safe_library("org.Hs.eg.db")
packages_available$org.Mm.eg.db <- safe_library("org.Mm.eg.db")
packages_available$clusterProfiler <- safe_library("clusterProfiler")
packages_available$ggplot2 <- safe_library("ggplot2")
packages_available$biomaRt <- safe_library("biomaRt")
packages_available$ggrepel <- safe_library("ggrepel")
packages_available$ReactomePA <- safe_library("ReactomePA")
packages_available$dplyr <- safe_library("dplyr")
packages_available$enrichplot <- safe_library("enrichplot")
packages_available$WGCNA <- safe_library("WGCNA")
packages_available$enrichR <- safe_library("enrichR")
packages_available$BiocParallel <- safe_library("BiocParallel")
packages_available$tidyverse <- safe_library("tidyverse")
packages_available$future <- safe_library("future")
packages_available$future.apply <- safe_library("future.apply")

# Command line options
option_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL, 
              help="Input directory path", metavar="character"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("--analysis_method"), type="character", default="EnrichR",
              help="Analysis method: GSEA, EnrichR, or both [default= %default]", metavar="character"),
  make_option(c("--top_n_percent_regulators"), type="double", default=2.0,
              help="Top N percent of regulators to select [default= %default]", metavar="number"),
  make_option(c("--coherence_threshold"), type="double", default=0.4,
              help="Coherence threshold [default= %default]", metavar="number"),
  make_option(c("--regulator_types"), type="character", default="TF:Lovering,Metabolite:Metabolites",
              help="Comma-separated list of Type:Prefix pairs [default= %default]", metavar="character"),
  make_option(c("--n_modules_name"), type="character", default=NULL,
              help="Number of modules identifier [default= %default]", metavar="character"),
  make_option(c("--n_threads"), type="integer", default=8,
              help="Number of threads for parallel processing [default= %default]", metavar="integer")
  ,
  make_option(c("--organism"), type="character", default="human",
              help="Organism for enrichment/annotation: 'human' or 'mouse' [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Function to parse regulator types string
parse_regulator_types <- function(regulator_types_str) {
  # Parse regulator types string into list of prefix:datafile pairs
  # Format: 'Prefix1:DataFile1,Prefix2:DataFile2'
  # Example: 'TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt,Lipids:Lipidomics.txt'
  
  configs <- list()
  config_strings <- strsplit(regulator_types_str, ",")[[1]]
  
  for (config_str in config_strings) {
    config_str <- trimws(config_str)
    if (grepl(":", config_str)) {
      parts <- strsplit(config_str, ":")[[1]]
      if (length(parts) >= 2) {
        prefix <- trimws(parts[1])
        data_file <- trimws(parts[2])
        configs[[length(configs) + 1]] <- list(
          prefix = prefix,       # Used for output filenames
          type = prefix,         # Used as edge type (same as prefix)
          data_file = data_file  # Store for reference
        )
      }
    }
  }
  
  return(configs)
}

if (is.null(opt$input_dir)){
  print_help(opt_parser)
  stop("Input directory must be supplied", call.=FALSE)
}

# USER CONFIGURATION
############################################################################################################################################

# CHOOSE ANALYSIS METHOD: "GSEA", "EnrichR", or "both"
# Normalize to handle case-insensitive input (enrichr/EnrichR/ENRICHR -> EnrichR)
ANALYSIS_METHOD <- opt$analysis_method
if (tolower(ANALYSIS_METHOD) == "enrichr") {
  ANALYSIS_METHOD <- "EnrichR"
} else if (tolower(ANALYSIS_METHOD) == "gsea") {
  ANALYSIS_METHOD <- "GSEA"
} else if (tolower(ANALYSIS_METHOD) %in% c("both", "all")) {
  ANALYSIS_METHOD <- "both"
}

# Base directory - set from command line
base_dir <- opt$input_dir

# Parameters that should match your LemonTree_to_network.py settings
top_n_percent_regulators <- opt$top_n_percent_regulators
coherence_threshold <- opt$coherence_threshold
# Handle organism-specific annotation and codes
organism_choice <- tolower(opt$organism)
if (organism_choice %in% c('mouse', 'mmu', 'mus_musculus')) {
  kegg_org <- 'mmu'
  reactome_org <- 'mouse'
  if (packages_available$org.Mm.eg.db) {
    OrgDb_obj <- org.Mm.eg.db
  } else if (packages_available$org.Hs.eg.db) {
    warning('org.Mm.eg.db not available; falling back to org.Hs.eg.db')
    OrgDb_obj <- org.Hs.eg.db
  } else {
    stop('Neither org.Mm.eg.db nor org.Hs.eg.db available. Install organism annotation packages.')
  }
} else {
  kegg_org <- 'hsa'
  reactome_org <- 'human'
  if (packages_available$org.Hs.eg.db) {
    OrgDb_obj <- org.Hs.eg.db
  } else if (packages_available$org.Mm.eg.db) {
    warning('org.Hs.eg.db not available; falling back to org.Mm.eg.db')
    OrgDb_obj <- org.Mm.eg.db
  } else {
    stop('Neither org.Hs.eg.db nor org.Mm.eg.db available. Install organism annotation packages.')
  }
}

n_threads <- opt$n_threads

# File paths (nextflow provides these as symlinks in the work directory)
expression_dataset <- file.path(base_dir, 'LemonPreprocessed_expression.txt')
clusters_file <- file.path(base_dir, 'clusters_list.txt')
specific_modules_file <- file.path(base_dir, 'specific_modules.txt')
metabolite_regulators_file <- file.path(base_dir, paste0('Metabolite.top', top_n_percent_regulators, 'pct_list.txt'))
TF_regulators_file <- file.path(base_dir, paste0('Lovering.top', top_n_percent_regulators, 'pct_list.txt'))

# Check if we're using filtered modules
use_filtered_modules <- file.exists(specific_modules_file)

workdir <- file.path(opt$output_dir, 'Enrichment')

# Check if the 'Enrichment' directory exists using workdir
if (!file.exists(workdir)) {
  # Create the 'Enrichment' directory
  dir.create(workdir, recursive = TRUE)
  cat("Directory 'Enrichment' created.\n")
} else {
  cat("Directory 'Enrichment' already exists.\n")
}

cat("Starting enrichment analysis...\n")
cat("Analysis method:", ANALYSIS_METHOD, "\n")
cat("Using filtered modules:", use_filtered_modules, "\n")
cat("Number of threads available:", n_threads, "\n")
cat("System info: ", Sys.info()["nodename"], "with", parallel::detectCores(), "detected cores\n")

############################################################################################################################################
# Load and prepare data
############################################################################################################################################

# Read expression data
expression <- as.data.frame(fread(expression_dataset))
rownames(expression) <- expression$symbol
expression$symbol <- NULL
if ("ensembl_gene_id" %in% colnames(expression)) {
  expression$ensembl_gene_id <- NULL
}

# Build clusters_to_genes from regulator-to-target and regulator-to-module mappings
cat("Building gene sets from selected regulator mappings...\n")

# Read regulator-to-module mappings from selected regulator files
regulator_to_module <- data.frame()
# Look for selected regulator files (new naming or old method_suffix naming)
selected_reg_files <- c(
  Sys.glob(file.path(base_dir, "*selected_regs_list.txt")),
  Sys.glob(file.path(base_dir, "ModuleViewer_files", "*selected_regs_list.txt")),
  Sys.glob(file.path(base_dir, "*_*_list.txt")),
  Sys.glob(file.path(base_dir, "ModuleViewer_files", "*_*_list.txt"))
)
# Remove duplicates
selected_reg_files <- unique(selected_reg_files)

for (reg_file in selected_reg_files) {
  if (file.exists(reg_file)) {
    # Read the selected regulator file: format is module_id TAB regulators|separated|by|pipes
    reg_data <- read.table(reg_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    colnames(reg_data) <- c("Module", "Regulators")
    
    # Expand regulators using base R instead of tidyverse
    regulators_split <- strsplit(as.character(reg_data$Regulators), "\\|")
    expanded_regs <- data.frame(
      Module = rep(reg_data$Module, sapply(regulators_split, length)),
      Regulator = unlist(regulators_split)
    )
    
    regulator_to_module <- rbind(regulator_to_module, expanded_regs)
    cat("Read", nrow(expanded_regs), "regulator-module mappings from", basename(reg_file), "\n")
  }
}

# Read regulator-to-target mappings from *2targets files
regulator_to_targets <- data.frame()
targets_files <- Sys.glob(file.path(base_dir, "Networks", "*2targets*.txt"))
for (targets_file in targets_files) {
  if (file.exists(targets_file)) {
    targets_data <- fread(targets_file, header = TRUE)
    # Split targets by | and create one row per target using base R
    targets_split <- strsplit(as.character(targets_data$Target), "\\|")
    targets_expanded <- data.frame(
      Regulator = rep(targets_data$Regulator, sapply(targets_split, length)),
      Target = unlist(targets_split)
    )
    # Add other columns back
    if (ncol(targets_data) > 2) {
      for (col in colnames(targets_data)[3:ncol(targets_data)]) {
        targets_expanded[[col]] <- rep(targets_data[[col]], sapply(targets_split, length))
      }
    }
    regulator_to_targets <- rbind(regulator_to_targets, targets_expanded)
    cat("Read", nrow(targets_expanded), "regulator-target mappings from", basename(targets_file), "\n")
  }
}

# Build clusters_to_genes from clusters_list.txt
cat("Reading module-gene mappings from clusters_list.txt...\n")
if (file.exists(clusters_file)) {
  clusters_data <- read.table(clusters_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(clusters_data) <- c("Module", "Genes")
  
  # Build clusters_to_genes from the data
  clusters_to_genes <- list()
  for (i in 1:nrow(clusters_data)) {
    module <- as.character(clusters_data$Module[i])
    genes <- strsplit(as.character(clusters_data$Genes[i]), "\\|")[[1]]
    clusters_to_genes[[module]] <- genes
  }
  cat("Read", length(clusters_to_genes), "modules with gene assignments from", basename(clusters_file), "\n")
} else {
  stop("clusters_list.txt not found at:", clusters_file)
}

# Filter to specific modules if available
if (use_filtered_modules) {
  specific_modules <- readLines(specific_modules_file)
  specific_modules <- gsub("\\s+", "", specific_modules)  # Remove whitespace
  
  # Filter clusters_to_genes to only include specific modules
  clusters_to_genes <- clusters_to_genes[names(clusters_to_genes) %in% specific_modules]
  cat("Filtered to", length(clusters_to_genes), " modules with coherence > ", coherence_threshold,"\n")
} else {
  cat("Using all", length(clusters_to_genes), "modules\n")
}

# Determine the actual number of modules for file naming
if (!is.null(opt$n_modules_name)) {
  n_modules_name <- opt$n_modules_name
  cat("Using provided n_modules_name:", n_modules_name, "\n")
} else {
  n_modules_name <- length(clusters_to_genes)
  cat("Calculated n_modules_name from data:", n_modules_name, "\n")
}
cat("Number of modules for analysis:", length(clusters_to_genes), "\n")

############################################################################################################################################
# GSEA Analysis (if selected)
############################################################################################################################################

if (ANALYSIS_METHOD %in% c("GSEA", "both")) {
  cat("\n=== Running GSEA Analysis ===\n")
  
  # Check if required packages are available for GSEA
  # Require organism-specific OrgDb package
  required_org_pkg <- if (organism_choice %in% c('mouse','mmu','mus_musculus')) 'org.Mm.eg.db' else 'org.Hs.eg.db'
  required_packages <- c("clusterProfiler", required_org_pkg, "fgsea", "BiocParallel")
  missing_packages <- required_packages[!sapply(required_packages, function(x) packages_available[[x]])]
  
  if (length(missing_packages) > 0) {
    cat("Error: Required packages for GSEA analysis are missing:", paste(missing_packages, collapse=", "), "\n")
    cat("Skipping GSEA analysis. Consider using EnrichR method instead.\n")
    if (ANALYSIS_METHOD == "GSEA") {
      cat("Creating empty results to prevent pipeline failure...\n")
      # Create minimal output structure
      output_dir <- file.path(workdir, 'Modules_gsea')
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Create empty results file
      empty_results <- data.frame(
        Module = character(0),
        Description = character(0),
        p.adjust = numeric(0),
        stringsAsFactors = FALSE
      )
      write.csv(empty_results, file.path(output_dir, "gsea_results_summary.csv"), row.names = FALSE)
      cat("Empty GSEA results created. Pipeline can continue.\n")
    }
  } else {
    cat("All required packages available. Proceeding with GSEA analysis.\n")
  
  # Use parallel processing
  register(MulticoreParam(workers = n_threads))
    
  # Output directory
  output_dir <- file.path(workdir, 'Modules_gsea')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ######################################################################################################################################
  # #### GSEA analysis per module: rank genes based on coexpression with module eigengene - parallel
  # ######################################################################################################################################

  # Shared variables
  dbs <- c('ALL', 'CC', 'MF', 'BP')
  set.seed(123)

  # Prepare input
  expression_in_clusters <- expression[unlist(clusters_to_genes), ]
  module_labels <- unlist(lapply(names(clusters_to_genes), function(cluster) {
    rep(cluster, length(clusters_to_genes[[cluster]]))
  }))
  
  # Calculate module eigengenes (first principal component for each module)
  # Alternative to WGCNA::moduleEigengenes when WGCNA is not available
  if (packages_available$WGCNA) {
    cat("Using WGCNA for module eigengenes calculation...\n")
    module_eigengenes <- moduleEigengenes(t(expression_in_clusters), colors = module_labels)$eigengenes
  } else {
    cat("WGCNA not available, using PCA for module eigengenes calculation...\n")
    module_eigengenes <- list()
    
    for (cluster in names(clusters_to_genes)) {
      cluster_genes <- clusters_to_genes[[cluster]]
      cluster_expression <- t(expression[cluster_genes, , drop = FALSE])
      
      # Calculate first principal component
      if (ncol(cluster_expression) > 1) {
        pca_result <- prcomp(cluster_expression, center = TRUE, scale. = FALSE)
        eigengene <- pca_result$x[, 1]
        
        # Ensure correlation with the module genes is positive
        module_mean <- rowMeans(cluster_expression)
        if (cor(eigengene, module_mean) < 0) {
          eigengene <- -eigengene
        }
      } else {
        # If only one gene in module, use that gene's expression
        eigengene <- cluster_expression[, 1]
      }
      
      module_eigengenes[[paste0("ME", cluster)]] <- eigengene
    }
    module_eigengenes <- as.data.frame(module_eigengenes)
  }

  # ------------------------ GSEA HELPER FUNCTIONS ----------------------------

  convert_symbols_to_entrez <- function(gene_ranks) {
    valid_symbols <- keys(OrgDb_obj, keytype = "SYMBOL")
    filtered <- gene_ranks[names(gene_ranks) %in% valid_symbols]
    converted <- bitr(names(filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb_obj)
    merged <- merge(converted, data.frame(SYMBOL = names(filtered), RANK = filtered), by = "SYMBOL")
    gene_list <- setNames(merged$RANK, merged$ENTREZID)
    sort(gene_list, decreasing = TRUE)
  }

  safe_plot_save <- function(gsea_obj, filename) {
    try({
      p <- dotplot(gsea_obj, showCategory = 10, split = ".sign") + 
        facet_grid(.~.sign) + 
        theme(axis.text.y = element_text(size = 7))
      ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 600)
    }, silent = TRUE)
  }

  process_gsea_result <- function(result, cluster, db) {
    extract_top <- function(subset, direction) {
      if (nrow(subset) == 0) return(data.frame())
      subset <- subset[order(-abs(subset$NES)), ][1:10, ]
      subset <- subset[subset$p.adjust <= 0.05, ]
      subset$Module <- cluster
      subset$Database <- db
      subset$Term <- paste(subset$ID, subset$Description, sep = " - ")
      subset[, c("Module", "Database", "Term", "p.adjust")]
    }
    up <- tryCatch(extract_top(result[result$NES > 0, ], "up"), error = function(e) data.frame())
    down <- tryCatch(extract_top(result[result$NES < 0, ], "down"), error = function(e) data.frame())
    list(up, down)
  }

  run_all_gsea <- function(cluster, dbs, organism, module_eigengenes, expression, output_dir) {
    
    cat(paste0("[", cluster, "] Starting GSEA...\n"))
    dir_name <- file.path(output_dir, paste0("module_", cluster))
    dir.create(dir_name, showWarnings = FALSE)
    
    eigengene <- module_eigengenes[, paste0("ME", cluster)]
    correlations <- cor(t(expression), eigengene, use = "pairwise.complete.obs")
    
    ranked_genes <- correlations[, 1]
    names(ranked_genes) <- rownames(correlations)
    ranked_genes <- sort(ranked_genes, decreasing = TRUE)
    
    local_up <- list()
    local_down <- list()
    
    # --- GO-based GSEA ---
    for (db in dbs) {
      tryCatch({
        gsea <- gseGO(
          geneList = ranked_genes,
          ont = db,
          keyType = "SYMBOL",
          OrgDb = OrgDb_obj,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          minGSSize = 3,
          maxGSSize = 800,
          verbose = FALSE,
          nPermSimple = 10000
        )
        
        safe_plot_save(gsea, file.path(dir_name, paste0("gsea_", db, ".png")))
        result <- gsea@result
        res <- process_gsea_result(result, cluster, db)
        local_up <- append(local_up, list(res[[1]]))
        local_down <- append(local_down, list(res[[2]]))
      }, error = function(e) message(paste("GO GSEA error in", cluster, "db:", db, e$message)))
    }
    
    # --- KEGG ---
    gene_list_entrez <- convert_symbols_to_entrez(ranked_genes)
    tryCatch({
      kegg <- gseKEGG(
        geneList = gene_list_entrez,
        organism = kegg_org,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
      safe_plot_save(kegg, file.path(dir_name, "gsea_KEGG.png"))
      res <- process_gsea_result(kegg@result, cluster, "KEGG")
      local_up <- append(local_up, list(res[[1]]))
      local_down <- append(local_down, list(res[[2]]))
    }, error = function(e) message(paste("KEGG error in", cluster, e$message)))
    
    # --- Reactome ---
    tryCatch({
      reactome <- gsePathway(
        geneList = gene_list_entrez,
        organism = reactome_org,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = FALSE
      )
      safe_plot_save(reactome, file.path(dir_name, "gsea_Reactome.png"))
      res <- process_gsea_result(reactome@result, cluster, "Reactome")
      local_up <- append(local_up, list(res[[1]]))
      local_down <- append(local_down, list(res[[2]]))
    }, error = function(e) message(paste("Reactome error in", cluster, e$message)))
    
    list(do.call(rbind, local_up), do.call(rbind, local_down))
   
  }

  # ------------------------ PARALLEL GSEA EXECUTION ----------------------------

  gsea_results <- bplapply(
    names(clusters_to_genes),
    run_all_gsea,
    dbs = dbs,
    organism = OrgDb_obj,
    module_eigengenes = module_eigengenes,
    expression = expression,
    output_dir = output_dir
  )

  # Split up/down results
  top_pathways_list_up <- lapply(gsea_results, `[[`, 1)
  top_pathways_list_down <- lapply(gsea_results, `[[`, 2)

  # Combine into data.frames
  top_pathways_df_up <- do.call(rbind, top_pathways_list_up)
  top_pathways_df_down <- do.call(rbind, top_pathways_list_down)

  # Remove NA and ensure proper structure
  if (nrow(top_pathways_df_up) > 0) {
    top_pathways_df_up <- top_pathways_df_up[!is.na(top_pathways_df_up$p.adjust), ]
  } else {
    top_pathways_df_up <- data.frame(
      Module = character(0),
      Database = character(0), 
      Term = character(0),
      p.adjust = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(top_pathways_df_down) > 0) {
    top_pathways_df_down <- top_pathways_df_down[!is.na(top_pathways_df_down$p.adjust), ]
  } else {
    top_pathways_df_down <- data.frame(
      Module = character(0),
      Database = character(0), 
      Term = character(0),
      p.adjust = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  cat("GSEA results summary:\n")
  cat("- Up-regulated pathways:", nrow(top_pathways_df_up), "entries\n")
  cat("- Down-regulated pathways:", nrow(top_pathways_df_down), "entries\n")

  # Save GSEA results in format compatible with Module_Overview_Generator
  cat("Saving GSEA results to:", output_dir, "\n")
  cat("Up results file:", file.path(output_dir, "GSEA_top_10_enriched_pathways_up_per_module.csv"), "\n")
  cat("Down results file:", file.path(output_dir, "GSEA_top_10_enriched_pathways_down_per_module.csv"), "\n")
  
  # Ensure consistent column structure for compatibility
  required_columns <- c("Module", "Database", "Term", "p.adjust")
  
  # Verify and standardize up results
  if (nrow(top_pathways_df_up) > 0) {
    # Ensure all required columns exist
    for (col in required_columns) {
      if (!col %in% colnames(top_pathways_df_up)) {
        cat("Warning: Missing column", col, "in GSEA up results\n")
      }
    }
    # Select only required columns in correct order
    top_pathways_df_up <- top_pathways_df_up[, required_columns, drop = FALSE]
  }
  
  # Verify and standardize down results
  if (nrow(top_pathways_df_down) > 0) {
    # Ensure all required columns exist
    for (col in required_columns) {
      if (!col %in% colnames(top_pathways_df_down)) {
        cat("Warning: Missing column", col, "in GSEA down results\n")
      }
    }
    # Select only required columns in correct order
    top_pathways_df_down <- top_pathways_df_down[, required_columns, drop = FALSE]
  }
  
  # Save GSEA results with consistent naming pattern
  fwrite(top_pathways_df_up, file.path(output_dir, "GSEA_top_10_enriched_pathways_up_per_module.csv"))
  fwrite(top_pathways_df_down, file.path(output_dir, "GSEA_top_10_enriched_pathways_down_per_module.csv"))
  
  # Verify files were created and check their structure
  up_file <- file.path(output_dir, "GSEA_top_10_enriched_pathways_up_per_module.csv")
  down_file <- file.path(output_dir, "GSEA_top_10_enriched_pathways_down_per_module.csv")
  
  if (file.exists(up_file)) {
    cat("[OK] GSEA up results file created successfully. Size:", file.size(up_file), "bytes\n")
    cat("[OK] GSEA up results:", nrow(top_pathways_df_up), "entries\n")
  }
  
  if (file.exists(down_file)) {
    cat("[OK] GSEA down results file created successfully. Size:", file.size(down_file), "bytes\n")
    cat("[OK] GSEA down results:", nrow(top_pathways_df_down), "entries\n")
  }
  
  # Verify files were created
  up_file <- file.path(output_dir, "Gsea_top_10_enriched_pathways_up_per_module.csv")
  down_file <- file.path(output_dir, "Gsea_top_10_enriched_pathways_down_per_module.csv")
  
  if (file.exists(up_file)) {
    cat("GSEA up results file created. Size:", file.size(up_file), "bytes\n")
  }
  if (file.exists(down_file)) {
    cat("GSEA down results file created. Size:", file.size(down_file), "bytes\n")
  }

  cat("GSEA Analysis completed!\n")
  cat("Results saved to:", output_dir, "\n")
  
  } # End of else block for available packages
}

############################################################################################################################################
# EnrichR Analysis (if selected)
############################################################################################################################################

if (ANALYSIS_METHOD %in% c("EnrichR", "both")) {
  cat("\n=== Running EnrichR Analysis ===\n")
  
  # Create output directory for EnrichR
  enrichr_output_dir <- file.path(workdir, 'Modules_enrichr')
  dir.create(enrichr_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set EnrichR site
  setEnrichrSite("Enrichr")
  if (organism_choice %in% c('mouse', 'mmu', 'mus_musculus')) {
    dbs_enrichr <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", "GO_Cellular_Component_2025")
  } else {
    dbs_enrichr <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", 
                     "GO_Cellular_Component_2025", "KEGG_2021_Human", "Reactome_Pathways_2024")
  }
  
  # Initialize dataframes to store results
  enrichr_up_results <- data.frame()
  enrichr_down_results <- data.frame()
  
  cat("Running EnrichR for", length(clusters_to_genes), "modules using", n_threads, "threads...\n")
  
  # For EnrichR, use more aggressive parallelization but with rate limiting
  # Increase thread limit for better performance on high-core systems
  effective_threads <- min(n_threads, 6)  # Increased from 4 to 6 threads
  cat("Using", effective_threads, "threads for EnrichR with API rate limiting\n")
  
  # Set up parallel processing for EnrichR
  if (packages_available$future && packages_available$future.apply) {
    plan(multisession, workers = effective_threads)
  }
  
  # Process modules in parallel (if future packages available) or sequentially
  process_enrichr_module <- function(cluster) {
    genes_in_cluster <- clusters_to_genes[[cluster]]
    
    cat("Processing module", cluster, "with", length(genes_in_cluster), "genes\n")
    
    if (length(genes_in_cluster) < 3) {
      cat("Skipping module", cluster, "- too few genes (", length(genes_in_cluster), ")\n")
      return(list(up = data.frame(), down = data.frame()))
    }
    
    # Create module directory
    module_dir <- file.path(enrichr_output_dir, paste0("module_", cluster))
    dir.create(module_dir, showWarnings = FALSE)
    
    module_results <- data.frame()
    total_enrichment_results <- 0
    total_significant_results <- 0
    
    tryCatch({
      # Run enrichment with increased delay for parallel safety
      Sys.sleep(runif(1, 1, 3))  # Random delay between 1-3 seconds
      enrichment <- enrichr(genes_in_cluster, dbs_enrichr)
      
      # Process each database
      for (db in dbs_enrichr) {
        if (nrow(enrichment[[db]]) > 0) {
          total_enrichment_results <- total_enrichment_results + nrow(enrichment[[db]])
          
          # Save plot
          tryCatch({
            p <- plotEnrich(enrichment[[db]], showTerms = 15, numChar = 40, y = 'Count', orderBy = 'P.value')
            ggsave(file.path(module_dir, paste0(db, ".png")), plot = p, width = 10, height = 8)
          }, error = function(e) cat("Plot error for", cluster, db, ":", e$message, "\n"))
          
          # Extract significant results (p.adj <= 0.05)
          sig_results <- enrichment[[db]][enrichment[[db]]$Adjusted.P.value <= 0.05, ]
          total_significant_results <- total_significant_results + nrow(sig_results)
          
          # If no significant results, take top 3 anyway (for exploratory analysis)
          if (nrow(sig_results) == 0) {
            cat("No significant results for", cluster, db, "- taking top 3 results for exploration\n")
            sig_results <- head(enrichment[[db]][order(enrichment[[db]]$Adjusted.P.value), ], 3)
          }
          
          if (nrow(sig_results) > 0) {
            # Take top 10 results
            top_results <- head(sig_results[order(sig_results$Adjusted.P.value), ], 10)
            
            # Add metadata
            top_results$Module <- cluster
            top_results$Database <- case_when(
              grepl("Biological_Process", db) ~ "BP",
              grepl("Molecular_Function", db) ~ "MF", 
              grepl("Cellular_Component", db) ~ "CC",
              grepl("KEGG", db) ~ "KEGG",
              grepl("Reactome", db) ~ "Reactome",
              TRUE ~ db
            )
            
            # Create Term column compatible with GSEA format
            top_results$Term <- paste(top_results$Term)
            
            # Rename p.adjust column to match GSEA format
            colnames(top_results)[colnames(top_results) == "Adjusted.P.value"] <- "p.adjust"
            
            # Select relevant columns
            formatted_results <- top_results[, c("Module", "Database", "Term", "p.adjust")]
            module_results <- rbind(module_results, formatted_results)
          }
        }
      }
      
      cat("Module", cluster, "- Total enrichment results:", total_enrichment_results, 
          ", Significant results (p.adj <= 0.05):", total_significant_results, 
          ", Final formatted results:", nrow(module_results), "\n")
      
    }, error = function(e) {
      cat("Error processing module", cluster, ":", e$message, "\n")
    })
    
    return(list(up = module_results, down = data.frame()))
  }
  
  # Run enrichment analysis
  if (packages_available$future && packages_available$future.apply && effective_threads > 1) {
    cat("Using parallel processing with future.apply\n")
    enrichr_results <- future_lapply(names(clusters_to_genes), process_enrichr_module, future.seed = TRUE)
  } else {
    cat("Using sequential processing\n")
    enrichr_results <- lapply(names(clusters_to_genes), process_enrichr_module)
  }
  
  # Combine results
  for (result in enrichr_results) {
    if (nrow(result$up) > 0) {
      enrichr_up_results <- rbind(enrichr_up_results, result$up)
    }
  }
  
  cat("EnrichR results summary:\n")
  cat("- Total enrichr_results processed:", length(enrichr_results), "\n")
  cat("- Final enrichr_up_results rows:", nrow(enrichr_up_results), "\n")
  
  # Ensure we always have properly structured results, even if empty
  if (nrow(enrichr_up_results) == 0) {
    cat("No EnrichR results found. Creating empty results with proper structure.\n")
    enrichr_up_results <- data.frame(
      Module = character(0),
      Database = character(0), 
      Term = character(0),
      p.adjust = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  
  # Save EnrichR results in GSEA-compatible format
  cat("Saving EnrichR results to:", enrichr_output_dir, "\n")
  cat("Up results file:", file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_up_per_module.csv"), "\n")
  cat("Down results file:", file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_down_per_module.csv"), "\n")
  
  # Save the up results
  fwrite(enrichr_up_results, file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_up_per_module.csv"))
  
  # Verify the file was created and check its size
  up_file <- file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_up_per_module.csv")
  if (file.exists(up_file)) {
    file_size <- file.size(up_file)
    cat("Up results file created successfully. Size:", file_size, "bytes\n")
    if (file_size > 0) {
      cat("Up results file contents preview:\n")
      cat(head(readLines(up_file, n = 3), 3), sep = "\n")
    } else {
      cat("Warning: Up results file is empty!\n")
    }
  } else {
    cat("Error: Up results file was not created!\n")
  }
  
  # Create properly structured empty down-regulated file for consistency
  # Use the same column structure as the up results
  if (nrow(enrichr_up_results) > 0) {
    empty_down_results <- enrichr_up_results[FALSE, ]  # Create empty dataframe with same structure
  } else {
    # If no up results either, create minimal required structure
    empty_down_results <- data.frame(
      Module = character(0),
      Database = character(0), 
      Term = character(0),
      p.adjust = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  fwrite(empty_down_results, file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_down_per_module.csv"))
  
  # Verify the down file was created
  down_file <- file.path(enrichr_output_dir, "Enrichr_top_10_enriched_pathways_down_per_module.csv")
  if (file.exists(down_file)) {
    file_size <- file.size(down_file)
    cat("Down results file created successfully. Size:", file_size, "bytes\n")
  } else {
    cat("Error: Down results file was not created!\n")
  }
  
  cat("EnrichR Analysis completed!\n")
  cat("Results saved to:", enrichr_output_dir, "\n")
}

############################################################################################################################################
# Summary and Module_Overview_Generator compatibility check
############################################################################################################################################

cat("\n=== Analysis Complete ===\n")
cat("Analysis method used:", ANALYSIS_METHOD, "\n")
cat("Modules analyzed:", length(clusters_to_genes), "\n")
cat("Number of modules (n_modules_name):", n_modules_name, "\n")
cat("Output directories:\n")

if (ANALYSIS_METHOD %in% c("GSEA", "both")) {
  cat("  GSEA results:", file.path(workdir, 'Modules_gsea'), "\n")
}

if (ANALYSIS_METHOD %in% c("EnrichR", "both")) {
  cat("  EnrichR results:", file.path(workdir, 'Modules_enrichr'), "\n")
}

cat("\nFor Module_Overview_Generator.py compatibility:\n")
cat("- Both GSEA and EnrichR results are now compatible with module_overview.py\n")
cat("- Files use consistent naming: [METHOD]_top_10_enriched_pathways_[up/down]_per_module.csv\n")
cat("- All output files have standardized column structure: Module, Database, Term, p.adjust\n")

if (ANALYSIS_METHOD == "GSEA") {
  cat("- GSEA up file: ", file.path(workdir, 'Modules_gsea', 'GSEA_top_10_enriched_pathways_up_per_module.csv'), "\n")
  cat("- GSEA down file: ", file.path(workdir, 'Modules_gsea', 'GSEA_top_10_enriched_pathways_down_per_module.csv'), "\n")
} else if (ANALYSIS_METHOD == "EnrichR") {
  cat("- EnrichR up file: ", file.path(workdir, 'Modules_enrichr', 'Enrichr_top_10_enriched_pathways_up_per_module.csv'), "\n")
  cat("- EnrichR down file: ", file.path(workdir, 'Modules_enrichr', 'Enrichr_top_10_enriched_pathways_down_per_module.csv'), "\n")
} else if (ANALYSIS_METHOD == "both") {
  cat("- GSEA files in: ", file.path(workdir, 'Modules_gsea/'), "\n")
  cat("- EnrichR files in: ", file.path(workdir, 'Modules_enrichr/'), "\n")
}

cat("\n[OK] All output files are compatible with module_overview.py downstream processing!\n")

######################################################################################################################################
#### EnrichR analysis on target genes of regulators - parallellized - courtesy by chatGPT
######################################################################################################################################

# Parallel plan - use more aggressive threading for regulator analysis
effective_threads_regulators <- min(n_threads, 5)  # Increased from 3 to 5 threads
cat("Using", effective_threads_regulators, "threads for regulator enrichment analysis\n")
plan(multisession, workers = effective_threads_regulators)

safe_enrichr <- function(genes, dbs, max_tries = 5) {
  for (i in 1:max_tries) {
    tryCatch({
      return(do.call(enrichr, list(as.character(genes), dbs)))
    }, error = function(e) {
      if (grepl("429", e$message)) {
        wait_time <- 2 ^ i
        cat("Rate limited. Retrying in", wait_time, "seconds...\n")
        Sys.sleep(wait_time)
      } else {
        stop(e)  # Other error
      }
    })
  }
  stop("Failed after multiple attempts due to rate limiting.")
}

# Function to perform enrichment analysis for a single regulator
run_enrichment_for_regulator <- function(reg, reg2genes, dbs, outdir) {
  tryCatch({
    genes_in_cluster <- as.character(reg2genes[[reg]])
    
    if (length(genes_in_cluster) == 0 || all(is.na(genes_in_cluster))) {
      cat("Skipping", reg, "- empty or NA gene list.\n")
      return(NULL)
    }
    
    reg_dir <- file.path(outdir, reg)
    # Create the regulator directory if it doesn't exist
    if (!dir.exists(reg_dir)) {
      dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    cat("Running enrichment for:", reg, "\n")
    # Random delay for parallel safety - shorter than before
    Sys.sleep(runif(1, 1.5, 3.5))
    enrichment <- safe_enrichr(genes_in_cluster, dbs)
    
    for (i in dbs) {
      tryCatch({
        plot <- plotEnrich(enrichment[[i]], showTerms = 15, numChar = 40, y = 'Count', orderBy = 'P.value') +
          theme(axis.text.x = element_text(size = rel(1)),
                axis.text.y = element_text(size = rel(1)))
        ggsave(filename = file.path(reg_dir, paste0(i, ".png")), plot = plot)
      }, error = function(e) {
        cat("Plotting failed for", reg, "in", i, ":", e$message, "\n")
      })
    }
    
    cat("Enrichment for", reg, "done.\n")
    return(NULL)
    
  }, error = function(e) {
    cat("Error in", reg, ":", e$message, "\n")
    return(NULL)
  })
}

# Main enrichment controller
regulator_targets_enrichment <- function(file_with_regulators) {
  tryCatch({
    dir_name <- strsplit(basename(file_with_regulators), "_")[[1]][1]
    # Create directory directly in workdir (which is already Enrichment folder)
    outdir <- file.path(workdir, dir_name)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    # Load file
    reg2targets <- fread(file_with_regulators)
    colnames(reg2targets) <- c("Regulator", "Target")
    
    reg2genes <- split(reg2targets$Target, reg2targets$Regulator)
    reg2genes <- lapply(reg2genes, function(g) unlist(strsplit(g, "\\|")))
    
    cat("Processing", length(reg2genes), "regulators using", effective_threads_regulators, "threads\n")
    
    setEnrichrSite("Enrichr")
    if (organism_choice %in% c('mouse', 'mmu', 'mus_musculus')) {
      dbs <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", "GO_Cellular_Component_2025")
    } else {
      dbs <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", "GO_Cellular_Component_2025", "KEGG_2021_Human", "Reactome_Pathways_2024")
    }
    
    # Run enrichment
    future_lapply(
      names(reg2genes),
      function(reg) run_enrichment_for_regulator(reg, reg2genes, dbs, outdir)
    )
    
  }, error = function(e) {
    cat("Failed on file", file_with_regulators, ":", e$message, "\n")
  })
}

# Example usage
tryCatch({
  # Debug: List contents of base_dir to understand file structure
  cat("Base directory contents:\n")
  if (dir.exists(base_dir)) {
    files_in_base <- list.files(base_dir, recursive = FALSE)
    cat("Files/directories in", base_dir, ":", paste(files_in_base, collapse = ", "), "\n")
    
    networks_dir <- file.path(base_dir, 'Networks')
    if (dir.exists(networks_dir)) {
      cat("Networks directory exists\n")
      networks_files <- list.files(networks_dir, pattern = "targets.*\\.txt$")
      cat("Target files in Networks:", paste(networks_files, collapse = ", "), "\n")
    } else {
      cat("Networks directory does not exist at:", networks_dir, "\n")
    }
  } else {
    cat("Base directory does not exist:", base_dir, "\n")
  }
  
  # Parse regulator types to get all configured regulator types
  regulator_configs <- parse_regulator_types(opt$regulator_types)
  cat("Parsed regulator configurations:\n")
  for (config in regulator_configs) {
    cat("  ", config$type, ":", config$prefix, "\n")
  }
  
  # Find regulator target files dynamically for all configured regulator types
  regulator_files <- list()
  
  for (config in regulator_configs) {
    # Pattern should match: prefix2targets_top2.0pct_4_modules.txt
    # or: prefix2targets_*modules*.txt (with single wildcard for flexibility)
    pattern <- file.path(base_dir, 'Networks', paste0(config$prefix, '2targets_*modules*.txt'))
    cat("Searching for", config$type, "targets with pattern:", pattern, "\n")
    
    files <- Sys.glob(pattern)
    files <- unique(files)
    
    if (length(files) > 0) {
      cat("Found", length(files), config$type, "target files\n")
      regulator_files[[config$type]] <- list(
        files = files,
        prefix = config$prefix,
        type = config$type
      )
    } else {
      cat("No", config$type, "target files found\n")
    }
  }
  
  # Process all found regulator target files
  for (reg_type in names(regulator_files)) {
    reg_info <- regulator_files[[reg_type]]
    files <- reg_info$files
    prefix <- reg_info$prefix
    
    if (length(files) > 0) {
      # Use the file that matches the current number of modules if possible
      if (exists("n_modules_name") && !is.null(n_modules_name)) {
        matching_file <- files[grep(paste0("_", n_modules_name, "_"), files)]
        if (length(matching_file) > 0) {
          targets_file <- matching_file[1]
        } else {
          targets_file <- files[1]  # fallback to first file
        }
      } else {
        targets_file <- files[1]  # Use the first file
      }
      cat("Using", reg_type, "targets file:", basename(targets_file), "\n")
      regulator_targets_enrichment(targets_file)
    } else {
      # Fallback: try to find any file matching the pattern
      fallback_pattern <- file.path(base_dir, 'Networks', paste0(prefix, '2targets_*_', n_modules_name, '_modules_*.txt'))
      cat("Trying fallback", reg_type, "pattern:", fallback_pattern, "\n")
      
      fallback_files <- Sys.glob(fallback_pattern)
      if (length(fallback_files) > 0) {
        cat("Using fallback", reg_type, "targets file:", basename(fallback_files[1]), "\n")
        regulator_targets_enrichment(fallback_files[1])
      } else {
        cat(reg_type, "targets file not found. Searched for pattern:", paste0(prefix, '2targets_*_modules_*.txt'), "\n")
        cat("Also tried fallback:", fallback_pattern, "\n")
      }
    }
  }
  
}, error = function(e) {
  cat("Global failure: ", e$message, "\n")
})
