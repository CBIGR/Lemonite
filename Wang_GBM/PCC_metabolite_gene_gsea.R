#!/usr/bin/Rscript


############################################################################################################################################
#### This is a script for performing gsea enrichment analysis on correlation between hPTM abundance and protein abundance
############################################################################################################################################

library(data.table)
library(fgsea)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
# library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
#library(enrichR)
library(biomaRt)
library(parallel)
library(ggrepel)
library(ReactomePA)
library(dplyr)

# Load necessary libraries explicitly for each core
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)


base_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree_bulk/no_controls/Variability_0.7_5817genes_old/'
setwd(paste0(base_dir, 'Enrichment/'))

organism <- org.Hs.eg.db

# Load Wang GBM data (preprocessed without TFA)
cat("Loading data...\n")
metabolomics <- fread(paste0(base_dir, 'Preprocessing/LemonPreprocessed_metabolomics.txt'), data.table = FALSE)
lipidomics <- fread(paste0(base_dir, 'Preprocessing/LemonPreprocessed_lipidomics.txt'), data.table = FALSE)
gene_expression <- fread(paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression.txt'), data.table = FALSE)

# Clean lipid names BEFORE setting as rownames: replace special characters with underscores
# This handles characters like (), /, :, -, spaces, etc.
cat("Cleaning lipid names...\n")
lipidomics$symbol_cleaned <- gsub("[^A-Za-z0-9_]", "_", lipidomics$symbol)
# Remove multiple consecutive underscores
lipidomics$symbol_cleaned <- gsub("_{2,}", "_", lipidomics$symbol_cleaned)
# Remove leading/trailing underscores
lipidomics$symbol_cleaned <- gsub("^_|_$", "", lipidomics$symbol_cleaned)

# Check for and remove duplicates in lipid names
cat("Checking for duplicate lipid names...\n")
n_lipids_before <- nrow(lipidomics)
duplicated_lipids <- lipidomics$symbol_cleaned[duplicated(lipidomics$symbol_cleaned)]
if (length(duplicated_lipids) > 0) {
  cat(paste0("Warning: Found ", length(duplicated_lipids), " duplicate lipid names after cleaning\n"))
  cat(paste0("Example duplicates: ", paste(head(duplicated_lipids, 3), collapse=", "), "\n"))
  # Keep only first occurrence of each lipid name
  lipidomics <- lipidomics[!duplicated(lipidomics$symbol_cleaned), ]
  cat(paste0("Removed ", n_lipids_before - nrow(lipidomics), " duplicate entries\n"))
} else {
  cat("No duplicate lipid names found\n")
}

# set column 'symbol' as rownames, then drop identifier columns
rownames(metabolomics) <- metabolomics$symbol
rownames(lipidomics) <- lipidomics$symbol_cleaned
rownames(gene_expression) <- gene_expression$symbol
metabolomics$symbol <- NULL; metabolomics$ensembl_gene_id <- NULL
lipidomics$symbol <- NULL; lipidomics$ensembl_gene_id <- NULL; lipidomics$symbol_cleaned <- NULL
gene_expression$symbol <- NULL; gene_expression$ensembl_gene_id <- NULL

cat("Data loaded successfully!\n")
cat(paste0("- Metabolites: ", nrow(metabolomics), " features x ", ncol(metabolomics), " samples\n"))
cat(paste0("- Lipids: ", nrow(lipidomics), " features x ", ncol(lipidomics), " samples\n"))
cat(paste0("- Genes: ", nrow(gene_expression), " features x ", ncol(gene_expression), " samples\n"))

# List of GO databases to run GSEA on
dbs <- c("BP", "MF", "CC", "ALL", "KEGG", "Reactome")  # Biological Process, Molecular Function, Cellular Component, All GO terms, KEGG pathways, Reactome pathways


##############################################################################################################
#### FUNCTION TO RUN CORRELATION + GSEA ANALYSIS
##############################################################################################################

run_correlation_gsea <- function(omics_data, gene_expression, omics_type, output_prefix, output_dir) {
  cat("\n========================================\n")
  cat(paste0("Starting analysis for: ", omics_type, "\n"))
  cat("========================================\n")
  
  # Make sure column names match and have the same order
  common_samples <- intersect(colnames(omics_data), colnames(gene_expression))
  omics_subset <- omics_data[, common_samples]
  genes_subset <- gene_expression[, common_samples]
  
  cat(paste0("Common samples: ", length(common_samples), "\n"))
  
  # create omics-gene correlation table
  cat(paste0("Creating correlation table between ", omics_type, " and genes...\n"))
  cor_table <- cor(t(omics_subset), t(genes_subset), method = "pearson", use = "pairwise.complete.obs")
  cat(paste0("Correlation table created with dimensions: ", nrow(cor_table), " x ", ncol(cor_table), "\n"))
  
  # Set the data for GSEA analysis (transpose so genes are rows)
  dat <- t(cor_table)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # write to file in the specified output directory
  output_file <- file.path(output_dir, paste0('PCC_', output_prefix, '_transcripts.csv'))
  write.table(dat, file = output_file, sep = '\t', quote = FALSE, row.names = TRUE, col.names = NA)
  cat(paste0("Correlation table saved to: ", output_file, "\n"))
  
  # Print data summary
  cat(paste0("\nData summary for GSEA:\n"))
  cat(paste0("- Number of genes (rows): ", nrow(dat), "\n"))
  cat(paste0("- Number of ", omics_type, " features (columns): ", ncol(dat), "\n"))
  cat(paste0("- Sample gene names: ", paste(head(rownames(dat), 3), collapse=", "), "...\n"))
  cat(paste0("- Sample ", omics_type, " feature names: ", paste(head(colnames(dat), 3), collapse=", "), "...\n\n"))
  
  return(dat)
}


# Safe plot saving function (from working script)
safe_plot_save <- function(gsea_obj, filename) {
  try({
    p <- dotplot(gsea_obj, showCategory = 10, split = ".sign") + 
      facet_grid(.~.sign) + 
      theme(axis.text.y = element_text(size = 7))
    ggsave(filename = filename, plot = p, width = 8, height = 6, dpi = 600)
  }, silent = TRUE)
}

# Function to convert gene symbols to Entrez IDs (from working script)
convert_symbols_to_entrez <- function(gene_ranks) {
  valid_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
  filtered <- gene_ranks[names(gene_ranks) %in% valid_symbols]
  converted <- bitr(names(filtered), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  merged <- merge(converted, data.frame(SYMBOL = names(filtered), RANK = filtered), by = "SYMBOL")
  gene_list <- setNames(merged$RANK, merged$ENTREZID)
  sort(gene_list, decreasing = TRUE)
}

# Function to run GSEA for a single column
run_gsea_for_column <- function(i, dat, dbs, organism) {
  # Load necessary libraries for each parallel worker
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(data.table)
  #i <- 3
  # Get and sanitize column name
  col_name <- colnames(dat)[i]
  col_name_sanitized <- gsub(" ", "_", col_name)
  col_name_sanitized <- gsub("/", "_", col_name)
  
  # Progress message
  cat(paste0("Starting GSEA for column ", i, "/", ncol(dat), ": ", col_name, "\n"))
  
  # Create output directory
  output_dir <- file.path('.', col_name_sanitized)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
    cat(paste0("Directory '", col_name_sanitized, "' created.\n"))
    # stop function
    #return(paste0("GSEA for column: ", col_name), " already completed")
  } else {
    cat(paste0("Directory '", col_name_sanitized, "' already exists.\n"))
  }
  
  # Extract gene vector and remove NAs
  vec <- dat[, i]
  names(vec) <- rownames(dat)
  vec <- vec[!is.na(vec)]
  vec <- sort(vec, decreasing = TRUE)
  
  # Check if we have enough valid genes
  if (length(vec) < 10) {
    cat(paste0("Warning: Only ", length(vec), " valid genes for column ", col_name, ". Skipping.\n"))
    return(paste0("Skipped column ", col_name, " - insufficient genes"))
  }
  
  cat(paste0("Processing ", length(vec), " genes for column ", col_name, "\n"))
  
  # Function to run GSEA for a single database
  run_gsea_single_db <- function(db) {
    cat(paste0("  Running GSEA for ", db, " ontology in column ", col_name, "...\n"))
    tryCatch({
      if (db %in% c("KEGG", "Reactome")) {
        # For KEGG and Reactome, convert symbols to Entrez IDs first
        gene_list_entrez <- convert_symbols_to_entrez(vec)
        
        if (db == "KEGG") {
          gse_all <- gseKEGG(geneList = gene_list_entrez, organism = "hsa", minGSSize = 3, maxGSSize = 800,
                            pvalueCutoff = 0.05, pAdjustMethod = "BH", verbose = FALSE)
        } else if (db == "Reactome") {
          gse_all <- gsePathway(geneList = gene_list_entrez, organism = "human", pvalueCutoff = 0.05,
                               pAdjustMethod = "BH", verbose = FALSE)
        }
      } else {
        # For GO databases, use the original approach
        gse_all <- gseGO(geneList = vec,
                         ont = db,
                         keyType = "SYMBOL",
                         nPermSimple = 10000,
                         minGSSize = 3,
                         maxGSSize = 800,
                         pvalueCutoff = 0.05,
                         verbose = FALSE,  # Set to FALSE to reduce output in parallel
                         OrgDb = organism,
                         pAdjustMethod = "BH")
      }
      
      # Use safe plot saving (same as working script)
      plot_file <- file.path(output_dir, paste0('gsea_dotplot_', db, '.png'))
      safe_plot_save(gse_all, plot_file)
      
      if (nrow(gse_all@result) == 0) {
        cat(paste0('  No enriched terms for ', db, ' in column ', col_name, '\n'))
      } else {
        cat(paste0('  Found ', nrow(gse_all@result), ' enriched terms for ', db, ' in column ', col_name, '\n'))
        cat(paste0('  Saved plot: ', plot_file, '\n'))
      }
      
      return(paste0("Completed GSEA for ", db, " in column ", col_name))
      
    }, error = function(e) {
      cat(paste0('  Error in GSEA for ', db, ' in column ', col_name, ': ', e$message, '\n'))
      return(NULL)
    })
  }
  
  # Run GSEA for all databases in this column 
  # Note: Could use mclapply here too for nested parallelization if memory allows
  db_results <- lapply(dbs, run_gsea_single_db)
  
  return(paste0("Completed GSEA for column: ", col_name))
}


##############################################################################################################
#### MAIN EXECUTION: Run analysis for both metabolomics and lipidomics
##############################################################################################################

# Determine number of cores to use (use all available cores minus 6, or at least 1)
n_cores <- max(1, detectCores() - 10)
cat(paste0("\nUsing ", n_cores, " cores for parallel processing\n\n"))

# List to store results from both analyses
all_results <- list()

##############################################################################################################
# ANALYSIS 1: METABOLOMICS
##############################################################################################################
cat("\n\n#######################################################\n")
cat("### STARTING METABOLOMICS ANALYSIS ###\n")
cat("#######################################################\n\n")

# Define output directory for metabolomics
metabolomics_dir <- file.path(paste0(base_dir, 'Enrichment/'), "PCC_metabolites_transcripts")

# Run correlation analysis for metabolomics
dat_metabolomics <- run_correlation_gsea(metabolomics, gene_expression, "metabolomics", "metabolites", metabolomics_dir)

# Change to metabolomics directory for GSEA results
setwd(metabolomics_dir)

cat(paste0("Total metabolite columns to process: ", ncol(dat_metabolomics), "\n"))
cat("Starting parallel GSEA analysis for metabolomics...\n")

# Run GSEA analysis in parallel across metabolite columns
results_metabolomics <- mclapply(seq_len(ncol(dat_metabolomics)), function(i) {
  run_gsea_for_column(i, dat_metabolomics, dbs, organism)
}, mc.cores = n_cores)

cat("\nMetabolomics GSEA analysis completed!\n")
cat(paste0("Results saved in: ", getwd(), "\n"))
all_results[["metabolomics"]] <- results_metabolomics

# Return to base enrichment directory
setwd(paste0(base_dir, 'Enrichment/'))


##############################################################################################################
# ANALYSIS 2: LIPIDOMICS
##############################################################################################################
cat("\n\n#######################################################\n")
cat("### STARTING LIPIDOMICS ANALYSIS ###\n")
cat("#######################################################\n\n")

# Define output directory for lipidomics
lipidomics_dir <- file.path(paste0(base_dir, 'Enrichment/'), "PCC_lipids_transcripts")

# Run correlation analysis for lipidomics
dat_lipidomics <- run_correlation_gsea(lipidomics, gene_expression, "lipidomics", "lipids", lipidomics_dir)

# Change to lipidomics directory for GSEA results
setwd(lipidomics_dir)

cat(paste0("Total lipid columns to process: ", ncol(dat_lipidomics), "\n"))
cat("Starting parallel GSEA analysis for lipidomics...\n")

# Run GSEA analysis in parallel across lipid columns
results_lipidomics <- mclapply(seq_len(ncol(dat_lipidomics)), function(i) {
  run_gsea_for_column(i, dat_lipidomics, dbs, organism)
}, mc.cores = n_cores)

cat("\nLipidomics GSEA analysis completed!\n")
cat(paste0("Results saved in: ", getwd(), "\n"))
all_results[["lipidomics"]] <- results_lipidomics


##############################################################################################################
# FINAL SUMMARY
##############################################################################################################
cat("\n\n#######################################################\n")
cat("### ALL ANALYSES COMPLETED! ###\n")
cat("#######################################################\n\n")
cat("Summary:\n")
cat(paste0("- Metabolomics: ", ncol(dat_metabolomics), " features analyzed\n"))
cat(paste0("  Results in: ", metabolomics_dir, "\n"))
cat(paste0("- Lipidomics: ", ncol(dat_lipidomics), " features analyzed\n"))
cat(paste0("  Results in: ", lipidomics_dir, "\n"))
cat("\nCorrelation tables saved:\n")
cat(paste0("- ", metabolomics_dir, "/PCC_metabolites_transcripts.csv\n"))
cat(paste0("- ", lipidomics_dir, "/PCC_lipids_transcripts.csv\n"))




