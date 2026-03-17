#!/usr/bin/Rscript

############################################################################################################################################
#### This is a script for performing gsea enrichment analysis on LemonTree modules
############################################################################################################################################

library(data.table)
library(fgsea)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(biomaRt)
library(ggrepel)
library(ReactomePA)
library(dplyr)
library(enrichplot)
library(WGCNA)
library(enrichR)
library(BiocParallel)
library(tidyverse)

# ------------------------ Plotting options ----------------------------
# Use custom EnrichR plotting function if available (plots top pathways, no background colour)
USE_CUSTOM_ENRICHR_PLOTTING <- TRUE

# Source plotting utilities from the canonical absolute path; fall back gracefully
plot_utils_path <- "/home/borisvdm/repo/gsea_plotting_utils.R"
found_plot_utils <- FALSE
if (file.exists(plot_utils_path)) {
  try({ source(plot_utils_path); found_plot_utils <- TRUE })
} else {
  message("plot_enrichr not found at ", plot_utils_path, "; EnrichR plotting will fall back to default plotEnrich if available.")
} 

# USER CONFIGURATION
############################################################################################################################################

# CHOOSE ANALYSIS METHOD: "GSEA", "EnrichR", or "both"
ANALYSIS_METHOD <- "EnrichR"  # Change this to choose your analysis method

# Base directory - UPDATE THIS PATH for IBD Lloyd-Price data
base_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum/'

# Parameters that should match your LemonTree_to_network.ipynb settings
fold <- 2
n_modules_name <- '46'  # This should match the output from LemonTree_to_network.ipynb
coherence_threshold <- 0.6
organism <- org.Hs.eg.db
n_threads <- 8
set.seed(1234)

# File paths (these should match LemonTree_to_network.ipynb outputs)
expression_dataset <- paste0(base_dir, 'Preprocessing/LemonPreprocessed_complete.txt')
clusters_file <- paste0(base_dir, 'ModuleViewer_files/clusters_list.txt')
specific_modules_file <- paste0(base_dir, 'Networks/specific_modules.txt')
metabolite_regulators_file <- paste0(base_dir, 'ModuleViewer_files/Metabolite.fold', fold, '_list.txt')
TF_regulators_file <- paste0(base_dir, 'ModuleViewer_files/Lovering.fold', fold, '_list.txt')

# Check if we're using filtered modules
use_filtered_modules <- file.exists(specific_modules_file)

workdir <- paste0(base_dir, 'Enrichment')

# Check if the 'Enrichment' directory exists using workdir
if (!file.exists(workdir)) {
  # Create the 'Enrichment' directory
  dir.create(workdir, recursive = TRUE)
  cat("Directory 'Enrichment' created.\n")
} else {
  cat("Directory 'Enrichment' already exists.\n")
}

setwd(workdir)

cat("Starting enrichment analysis...\n")
cat("Analysis method:", ANALYSIS_METHOD, "\n")
cat("Using filtered modules:", use_filtered_modules, "\n")

############################################################################################################################################
# Load and prepare data
############################################################################################################################################

# Read expression data
expression <- as.data.frame(fread(expression_dataset))
# remove rows with duplicate gene symbols
expression <- expression[!duplicated(expression$symbol), ]
# remove rows with missing value in symbol
expression <- expression[!is.na(expression$symbol), ]
rownames(expression) <- expression$symbol
expression$symbol <- NULL
if ("ensembl_gene_id" %in% colnames(expression)) {
  expression$ensembl_gene_id <- NULL
}

# Read clusters file
clusters <- fread(clusters_file)
clusters_to_genes <- vector('list', nrow(clusters))
names(clusters_to_genes) <- as.character(clusters$V1)

# Add genes to the named list
for (i in (1:nrow(clusters))){
  genes <- as.character(clusters[i, 2])
  genes <- str_split(genes, "\\|")
  clusters_nr <- as.character(clusters[i,1])
  clusters_to_genes[clusters_nr] <- genes
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


############################################################################################################################################
# GSEA Analysis (if selected)
############################################################################################################################################

if (ANALYSIS_METHOD %in% c("GSEA", "both")) {
  cat("\n=== Running GSEA Analysis ===\n")
  # Use 6 cores
  #plan(multisession, workers = 1)
  register(MulticoreParam(workers = n_threads))
    
  # Output directory
  output_dir <- paste0(base_dir, 'Enrichment/Modules_gsea')
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
  module_eigengenes <- moduleEigengenes(t(expression_in_clusters), colors = module_labels)$eigengenes

  # ------------------------ GSEA HELPER FUNCTIONS ----------------------------

  convert_symbols_to_entrez <- function(gene_ranks) {
    valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    filtered <- gene_ranks[names(gene_ranks) %in% valid_symbols]
    converted <- bitr(names(filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    merged <- merge(converted, data.frame(SYMBOL = names(filtered), RANK = filtered), by = "SYMBOL")
    gene_list <- setNames(merged$RANK, merged$ENTREZID)
    sort(gene_list, decreasing = TRUE)
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

    #cluster <- '20'

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
          OrgDb = organism,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          minGSSize = 3,
          maxGSSize = 800,
          verbose = FALSE,
          nPermSimple = 10000,
          seed = TRUE
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
        organism = "hsa",
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = TRUE
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
        organism = "human",
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        verbose = FALSE,
        seed = TRUE
      )
      safe_plot_save(reactome, file.path(dir_name, "gsea_Reactome.png"))
      res <- process_gsea_result(reactome@result, cluster, "Reactome")
      local_up <- append(local_up, list(res[[1]]))
      local_down <- append(local_down, list(res[[2]]))
    }, error = function(e) message(paste("Reactome error in", cluster, e$message)))

    list(do.call(rbind, local_up), do.call(rbind, local_down))

  }



  # ------------------------ PARALLEL GSEA EXECUTION ----------------------------
  #gsea_results <- bplapply(names(clusters_to_genes), run_all_gsea)

  gsea_results <- bplapply(
    names(clusters_to_genes),
    run_all_gsea,
    dbs = dbs,
    organism = organism,
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

  # Remove NA
  top_pathways_df_up <- top_pathways_df_up[!is.na(top_pathways_df_up$p.adjust), ]
  top_pathways_df_down <- top_pathways_df_down[!is.na(top_pathways_df_down$p.adjust), ]

  # Save GSEA results in format compatible with Module_Overview_Generator
  fwrite(top_pathways_df_up, file.path(output_dir, "Gsea_top_10_enriched_pathways_up_per_module.csv"))
  fwrite(top_pathways_df_down, file.path(output_dir, "Gsea_top_10_enriched_pathways_down_per_module.csv"))

  cat("GSEA Analysis completed!\n")
  cat("Results saved to:", output_dir, "\n")
} else {
  cat("Skipping GSEA (ANALYSIS_METHOD =", ANALYSIS_METHOD, ")\n")
}



############################################################################################################################################
# EnrichR Analysis (if selected)
############################################################################################################################################

if (ANALYSIS_METHOD %in% c("EnrichR", "both")) {
  cat("\n=== Running EnrichR Analysis ===\n")
  
  # Create output directory for EnrichR
  enrichr_output_dir <- paste0(base_dir, 'Enrichment/Modules_enrichr')
  dir.create(enrichr_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Set EnrichR site
  setEnrichrSite("Enrichr")
  dbs_enrichr <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", 
                   "GO_Cellular_Component_2025", "KEGG_2021_Human", "Reactome_Pathways_2024")
  
  # Initialize dataframes to store results
  enrichr_up_results <- data.frame()
  enrichr_down_results <- data.frame()
  
  cat("Running EnrichR for", length(clusters_to_genes), "modules...\n")
  
  for (cluster in names(clusters_to_genes)) {
    cat("Processing module:", cluster, "\n")
    
    genes_in_cluster <- clusters_to_genes[[cluster]]
    
    if (length(genes_in_cluster) < 3) {
      cat("Skipping module", cluster, "- too few genes\n")
      next
    }
    
    # Create module directory
    module_dir <- file.path(enrichr_output_dir, paste0("module_", cluster))
    dir.create(module_dir, showWarnings = FALSE)
    
    tryCatch({
      # Run enrichment
      Sys.sleep(1)  # Rate limiting
      enrichment <- enrichr(genes_in_cluster, dbs_enrichr)
      
      # Process each database
      for (db in dbs_enrichr) {
        if (nrow(enrichment[[db]]) > 0) {
          # Save plot
          tryCatch({
            if (exists("plot_enrichr") && USE_CUSTOM_ENRICHR_PLOTTING) {
              p <- plot_enrichr(enrichment[[db]], top_n = 10, numChar = 60, title = db) +
            theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 12), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
            } else {
              p <- plotEnrich(enrichment[[db]], showTerms = 10, numChar = 80, y = 'Count', orderBy = 'P.value') +
                theme(axis.text.y = element_text(size = 14), axis.text.x = element_text(size = 12), panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
            }
            ggsave(file.path(module_dir, paste0(db, ".png")), plot = p, width = 12, height = 5)
          }, error = function(e) cat("Plot error for", cluster, db, ":", e$message, "\n"))
          
          # Extract significant results
          sig_results <- enrichment[[db]][enrichment[[db]]$Adjusted.P.value <= 0.05, ]
          
          if (nrow(sig_results) > 0) {
            # Take top 10 results
            # top_results <- head(sig_results[order(sig_results$Adjusted.P.value), ], 10)
            
            # Take all results
            top_results <- sig_results[order(sig_results$Adjusted.P.value), ]
            
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
            
            # For EnrichR, we'll treat all results as "up" since it doesn't distinguish direction
            enrichr_up_results <- rbind(enrichr_up_results, formatted_results)
          }
        }
      }
      
    }, error = function(e) {
      cat("Error processing module", cluster, ":", e$message, "\n")
    })
  }
  
  # Save EnrichR results in GSEA-compatible format
  fwrite(enrichr_up_results, file.path(enrichr_output_dir, "Enrichr_all_enriched_pathways_up_per_module.csv"))
  # Create empty down-regulated file for consistency
  fwrite(data.frame(), file.path(enrichr_output_dir, "Enrichr_all_enriched_pathways_down_per_module.csv"))
  
  cat("EnrichR Analysis completed!\n")
  cat("Results saved to:", enrichr_output_dir, "\n")
}

############################################################################################################################################
# Summary and Module_Overview_Generator compatibility check
############################################################################################################################################

cat("\n=== Analysis Complete ===\n")
cat("Analysis method used:", ANALYSIS_METHOD, "\n")
cat("Modules analyzed:", length(clusters_to_genes), "\n")
cat("Output directories:\n")

if (ANALYSIS_METHOD %in% c("GSEA", "both")) {
  cat("  GSEA results:", paste0(base_dir, 'Enrichment/Modules_gsea'), "\n")
}

if (ANALYSIS_METHOD %in% c("EnrichR", "both")) {
  cat("  EnrichR results:", paste0(base_dir, 'Enrichment/Modules_enrichr'), "\n")
}

cat("\nFor Module_Overview_Generator.ipynb compatibility:\n")
cat("- Update enrichment_up_file path to point to your chosen analysis results\n")
cat("- Update enrichment_down_file path accordingly\n")

if (ANALYSIS_METHOD == "GSEA") {
  cat("- Use: enrichment_up_file = base_dir + '/Enrichment/Modules_gsea/Gsea_top_10_enriched_pathways_up_per_module.csv'\n")
  cat("- Use: enrichment_down_file = base_dir + '/Enrichment/Modules_gsea/Gsea_top_10_enriched_pathways_down_per_module.csv'\n")
} else if (ANALYSIS_METHOD == "EnrichR") {
  cat("- Use: enrichment_up_file = base_dir + '/Enrichment/Modules_enrichr/Enrichr_top_10_enriched_pathways_up_per_module.csv'\n")
  cat("- Use: enrichment_down_file = base_dir + '/Enrichment/Modules_enrichr/Enrichr_top_10_enriched_pathways_down_per_module.csv'\n")
}

cat("\nAll output files are compatible with Module_Overview_Generator.ipynb!\n")


######################################################################################################################################
#### EnrichR analysis on target genes of regulators - parallellized - courtesy by chatGPT
######################################################################################################################################


library(future)
library(future.apply)

# Parallel plan
plan(multisession, workers = 6)

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
    reg_dir <- file.path(outdir, chartr("/", "_", reg))
    #reg_dir <- file.path(outdir, reg)
    
    if (!dir.exists(reg_dir)) {
      dir.create(reg_dir, recursive = TRUE)
    } else {
      cat("Enrichment for", reg, "already exists. Skipping.\n")
      return(NULL)
    }

    cat("Running enrichment for:", reg, "\n")
    Sys.sleep(2.5)
    enrichment <- safe_enrichr(genes_in_cluster, dbs)
    
    for (i in dbs) {
      tryCatch({
        if (exists("plot_enrichr") && USE_CUSTOM_ENRICHR_PLOTTING) {
          plot <- plot_enrichr(enrichment[[i]], top_n = 10, numChar = 80, title = i) +
            theme(axis.text.x = element_text(size = rel(1)),
                  axis.text.y = element_text(size = rel(1.2)),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
        } else {
          plot <- plotEnrich(enrichment[[i]], showTerms = 10, numChar = 80, y = 'Count', orderBy = 'P.value') +
            theme(axis.text.x = element_text(size = rel(1)),
                  axis.text.y = element_text(size = rel(1.2)),
                  panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
        }
        ggsave(filename = file.path(reg_dir, paste0(i, ".png")), plot = plot, width = 14, height = 5) 
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
    outdir <- file.path(base_dir, "Enrichment", dir_name)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    
    # Load file
    reg2targets <- fread(file_with_regulators)
    colnames(reg2targets) <- c("Regulator", "Target")
    
    reg2genes <- split(reg2targets$Target, reg2targets$Regulator)
    reg2genes <- lapply(reg2genes, function(g) unlist(str_split(g, "\\|")))
    
    setEnrichrSite("Enrichr")
    dbs <- c("GO_Biological_Process_2025", "GO_Molecular_Function_2025", 
                     "GO_Cellular_Component_2025", "KEGG_2021_Human", "Reactome_Pathways_2024")
    
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
  #regulator_targets_enrichment(file.path(base_dir, 'Networks', paste0('hPTMs2targets_fold', fold, '_selected_', n_modules, '_modules.txt')))
  regulator_targets_enrichment(file.path(base_dir, 'Networks', paste0('Metabolites2targets_percentile', fold, '_', n_modules_name, '_modules.txt')))
  regulator_targets_enrichment(file.path(base_dir, 'Networks', paste0('TFs2targets_percentile', fold, '_', n_modules_name, '_modules.txt')))
  regulator_targets_enrichment(file.path(base_dir, 'Networks', paste0('Lipids2targets_percentile', fold, '_', n_modules_name, '_modules.txt')))
  
}, error = function(e) {
  cat("Global failure: ", e$message, "\n")
})


