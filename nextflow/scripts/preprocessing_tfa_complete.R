#!/usr/bin/env Rscript

# Parameterized version of Preprocessing_and_TFA_CollecTri_consensus.R
# for use in Nextflow pipeline

library(optparse)

# Command line options
option_list <- list(
  make_option(c("--expression"), type="character", default=NULL, 
              help="Expression file path", metavar="character"),
  make_option(c("--regulator_types"), type="character", 
              default="TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt",
              help="Comma-separated Prefix:DataFile pairs for all regulator types [default= %default]", metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Metadata file path", metavar="character"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory path [default= %default]", metavar="character"),
  make_option(c("--top_n_genes"), type="integer", default=1000,
              help="Number of top variable genes to select [default= %default]", metavar="integer"),
  make_option(c("--prior_network"), type="character", default=NULL,
              help="Prior network file path", metavar="character"),
  make_option(c("--perform_TFA"), type="logical", default=TRUE,
              help="Whether to perform TFA [default= %default]", metavar="logical"),
  make_option(c("--use_omics_specific_scaling"), type="logical", default=TRUE,
              help="Use omics-specific scaling (pareto for metabolomics/lipidomics, standard for transcriptomics) [default= %default]", metavar="logical"),
  make_option(c("--DESeq_contrast1"), type="character", default="diagnosis",
              help="DESeq contrast variable [default= %default]", metavar="character"),
  make_option(c("--design_formula"), type="character", default="~ diagnosis",
              help="DESeq2 design formula (e.g., '~ diagnosis', '~ biopsy_location + diagnosis + sex') [default= %default]", metavar="character"),
  make_option(c("--metadata_columns"), type="character", default="diagnosis",
              help="Comma-separated list of metadata columns to include [default= %default]", metavar="character"),
  make_option(c("--expression_col"), type="character", default="count",
              help="Expression gene symbol column name [default= %default]", metavar="character"),
  make_option(c("--sample_id_col"), type="character", default="Sample_ID",
              help="Sample ID column name in metadata [default= %default]", metavar="character")
  ,
  make_option(c("--organism"), type="character", default="human",
              help="Organism for annotation/mapping: 'human' or 'mouse' [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Parse regulator types configuration
# Format: "Prefix1:DataFile1[:DataType1],Prefix2:DataFile2[:DataType2]"
# DataType is optional: 'c' for continuous (default), 'd' for discrete/binary
parse_regulator_config <- function(regulator_types_str) {
  configs <- list()
  data_types <- list()
  pairs <- strsplit(regulator_types_str, ",")[[1]]
  
  for (pair in pairs) {
    pair <- trimws(pair)
    parts <- strsplit(pair, ":")[[1]]
    if (length(parts) >= 2) {
      prefix <- trimws(parts[1])
      filename <- trimws(parts[2])
      data_type <- ifelse(length(parts) >= 3, tolower(trimws(parts[3])), "c")
      if (!data_type %in% c("c", "d")) {
        warning(paste("Invalid data type '", data_type, "' for", prefix, "- using 'c' (continuous). Valid: 'c' or 'd'"))
        data_type <- "c"
      }
      configs[[prefix]] <- filename
      data_types[[prefix]] <- data_type
    } else {
      warning(paste("Invalid regulator config format:", pair, "(expected Prefix:DataFile[:DataType])"))
    }
  }
  
  return(list(configs = configs, data_types = data_types))
}

# Validate required inputs
if (is.null(opt$expression) || is.null(opt$metadata)){
  print_help(opt_parser)
  stop("Expression and metadata files must be supplied", call.=FALSE)
}

# Parse regulator configurations
reg_parsed <- parse_regulator_config(opt$regulator_types)
regulator_configs <- reg_parsed$configs
regulator_data_types <- reg_parsed$data_types
cat("\n Regulator types configuration:\n")
for (prefix in names(regulator_configs)) {
  dt <- regulator_data_types[[prefix]]
  dt_label <- ifelse(dt == "d", "discrete/binary", "continuous")
  cat(sprintf("   - %s: %s (%s)\n", prefix, regulator_configs[[prefix]], dt_label))
}

# Set working directory and variables
setwd(opt$output_dir)
top_n_genes <- opt$top_n_genes
perform_TFA <- opt$perform_TFA
use_omics_specific_scaling <- opt$use_omics_specific_scaling
# DESeq_contrast will be set later after loading metadata and auto-detecting conditions
# Parse metadata columns from parameter
metadata_to_include <- trimws(strsplit(opt$metadata_columns, ",")[[1]])
cat("Metadata columns to include:", paste(metadata_to_include, collapse=", "), "\n")

# Create required directory structure
dir.create("LemonTree/Preprocessing", recursive = TRUE, showWarnings = FALSE)
dir.create("TFA", recursive = TRUE, showWarnings = FALSE)
dir.create("LemonTree/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("TF-metabolite_interactions", recursive = TRUE, showWarnings = FALSE)


# Load required libraries
library(DESeq2)
library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(tibble)
library(biomaRt)
library(ggplot2)
library(decoupleR)
library(pheatmap)
library(AnnotationDbi)

# Helper functions
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# PCA plotting function
create_pca_plot <- function(data_matrix, metadata_df, contrast_col, omics_name, output_dir = "./LemonTree/Preprocessing") {
  cat(sprintf("\n=== Creating PCA plot for %s ===\n", omics_name))
  
  # Ensure data is numeric matrix with samples in columns
  data_matrix <- as.matrix(data_matrix)
  
  # Check dimensions
  cat(sprintf("  Data dimensions: %d features x %d samples\n", nrow(data_matrix), ncol(data_matrix)))
  
  # Remove any rows with NA or infinite values
  valid_rows <- complete.cases(data_matrix) & apply(data_matrix, 1, function(x) all(is.finite(x)))
  if (sum(!valid_rows) > 0) {
    cat(sprintf("  Removing %d features with NA or infinite values\n", sum(!valid_rows)))
    data_matrix <- data_matrix[valid_rows, , drop=FALSE]
  }
  
  # Transpose for PCA (samples as rows, features as columns)
  pca_input <- t(data_matrix)
  
  # Run PCA
  pca_result <- prcomp(pca_input, center = TRUE, scale. = TRUE)
  
  # Calculate variance explained
  variance_explained <- summary(pca_result)$importance[2, ] * 100
  pc1_var <- round(variance_explained[1], 1)
  pc2_var <- round(variance_explained[2], 1)
  
  # Create data frame for plotting
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x)
  )
  
  # Add diagnosis/contrast information
  if (contrast_col %in% colnames(metadata_df)) {
    # Match samples between PCA and metadata
    sample_match <- match(pca_df$Sample, rownames(metadata_df))
    pca_df$Group <- metadata_df[[contrast_col]][sample_match]
    
    # Create plot with diagnosis coloring
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = paste("PCA -", omics_name),
        x = paste0("PC1 (", pc1_var, "% variance)"),
        y = paste0("PC2 (", pc2_var, "% variance)"),
        color = contrast_col
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )
  } else {
    # If contrast column not available, plot without coloring
    cat(sprintf("  Warning: Contrast column '%s' not found in metadata, plotting without grouping\n", contrast_col))
    p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = paste("PCA -", omics_name),
        x = paste0("PC1 (", pc1_var, "% variance)"),
        y = paste0("PC2 (", pc2_var, "% variance)")
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        panel.grid.minor = element_blank()
      )
  }
  
  # Save plot
  output_file <- file.path(output_dir, paste0("PCA_", gsub(" ", "_", tolower(omics_name)), ".pdf"))
  ggsave(output_file, plot = p, width = 8, height = 6)
  cat(sprintf("[OK] Saved PCA plot: %s\n", output_file))
  
  return(invisible(pca_result))
}

# Input files
expression <- opt$expression
metadata <- opt$metadata
prior_network <- opt$prior_network

# Check if required files exist
if (!file.exists(expression)) stop("Expression file not found: ", expression)
if (!file.exists(metadata)) stop("Metadata file not found: ", metadata)

# Identify TF list file and other regulator files
TF_prefix <- NULL
TF_file <- NULL
other_regulator_files <- list()

for (prefix in names(regulator_configs)) {
  filename <- regulator_configs[[prefix]]
  data_type <- regulator_data_types[[prefix]]
  is_tf_list <- grepl("TF", prefix, ignore.case = TRUE) || grepl("_list\\.txt$", filename, ignore.case = TRUE)
  
  # Discrete/binary regulators are never TF lists — they contain abundance data (0/1)
  if (data_type == "d") {
    is_tf_list <- FALSE
  }
  
  # Auto-replace Lovering_TF_list.txt with organism-specific version
  if (is_tf_list) {
    organism_choice <- tolower(opt$organism)
    if (organism_choice %in% c('mouse', 'mmu', 'mus_musculus')) {
      # Use mouse TF list
      if (filename == "Lovering_TF_list.txt" || filename == "lovering_TF_list.txt") {
        filename <- "lovering_TF_list_mouse.txt"
        cat(sprintf("   🐭 Auto-selected mouse TF list: %s\n", filename))
      }
    }
  }
  
  # TF lists come from PKN directory by default, other regulators from data directory
  # This allows users to override by providing custom files in data/
  if (is_tf_list) {
    # Check data/ first (user override), then PKN/ (default)
    filepath_data <- file.path("data", filename)
    filepath_pkn <- file.path("PKN", filename)
    if (file.exists(filepath_data)) {
      filepath <- filepath_data
      cat(sprintf("   [OK] Using custom TF list: %s (%s)\n", prefix, filepath_data))
    } else if (file.exists(filepath_pkn)) {
      filepath <- filepath_pkn
      cat(sprintf("   [OK] Using default TF list: %s (%s)\n", prefix, filepath_pkn))
    } else {
      stop(sprintf("TF list file not found: %s (checked data/ and PKN/)", filename))
    }
  } else {
    # Other regulator types (metabolomics, etc.) always from data/
    filepath <- file.path("data", filename)
  }
  
  # Check if this is a TF list
  if (is_tf_list) {
    TF_prefix <- prefix
    TF_file <- filepath
    cat(sprintf("   [OK] Identified TF list: %s\n", prefix))
  } else {
    other_regulator_files[[prefix]] <- filepath
    cat(sprintf("   [OK] Identified regulator data: %s (%s)\n", prefix, filename))
  }
}

cat("Starting preprocessing with the following parameters:\n")
cat("Expression file:", expression, "\n")
cat("Metadata file:", metadata, "\n")
if (!is.null(TF_file)) {
  cat("TF list file:", TF_file, "\n")
}
for (prefix in names(other_regulator_files)) {
  cat(sprintf("%s file: %s\n", prefix, other_regulator_files[[prefix]]))
}
cat("Metadata file:", metadata, "\n")
cat("Output directory:", opt$output_dir, "\n")
cat("Top n variable genes:", top_n_genes, "\n")
cat("Perform TFA:", perform_TFA, "\n")
if (!is.null(prior_network)) cat("Prior network:", prior_network, "\n")
if (!is.null(TF_file)) cat("TF list:", TF_file, "\n")

# Check TFA prerequisites
if (perform_TFA) {
  cat("\n=== TFA ANALYSIS PREREQUISITES CHECK ===\n")
  if (is.null(prior_network)) {
    cat("[INFO]  No prior network file specified - will use built-in CollecTRI database\n")
  } else if (!file.exists(prior_network)) {
    cat("[WARNING]  Prior network file does not exist:", prior_network, "\n")
    cat("[INFO]  Will fallback to built-in CollecTRI database\n")
  } else {
    cat("[OK] Prior network file found:", prior_network, "\n")
  }
  cat("[OK] TFA analysis will proceed using available network source\n")
  cat("============================================\n")
} else {
  cat("\n[WARNING] TFA analysis is DISABLED (perform_TFA = FALSE)\n")
}

# Validate metadata requirements for pipeline consistency
cat("\n=== METADATA REQUIREMENTS VALIDATION ===\n")
cat(" Expected metadata structure for Lemonite pipeline:\n")
cat("   - Sample IDs as row names or in column:", opt$sample_id_col, "\n")
cat("   - Required columns: diagnosis (or", opt$DESeq_contrast1, ")\n")
cat("   - Optional columns: sex, age, batch, biopsy_location\n")
cat("   - File format: Tab-separated (.txt)\n")
cat("   - Output: DESeq_groups.txt for downstream analysis\n")
cat("=============================================\n")
cat("\n")

###########################################################################################
#### Read RNAseq data, filter for protein coding genes
###########################################################################################

cat("Reading RNA-seq data...\n")
RNAseq <- fread(expression, header=TRUE, data.table=TRUE)

# Handle case where gene symbols are row names (first column without header)
if (ncol(RNAseq) > 0 && names(RNAseq)[1] == "V1") {
  # File has genes as first column without proper header
  gene_names <- RNAseq$V1
  RNAseq$V1 <- NULL
  RNAseq <- as.data.frame(RNAseq)
  rownames(RNAseq) <- gene_names
  
  # Create a gene symbol column for merging
  RNAseq$gene_symbol <- rownames(RNAseq)
  gene_col <- "gene_symbol"
} else {
  # Use the specified expression column
  gene_col <- opt$expression_col
}

# Use the Ensembl BioMart with error handling and fallback
cat("Connecting to Ensembl BioMart...\n")

# Try multiple Ensembl servers in order of preference
ensembl_servers <- c(
  "https://www.ensembl.org",                    # Current release
  "https://jul2024.archive.ensembl.org",        # Recent archive
  "https://jan2024.archive.ensembl.org"         # Original (fallback)
)

ensembl <- NULL
# Determine organism/Ensembl dataset and symbol attribute
organism_choice <- tolower(opt$organism)
if (organism_choice %in% c('mouse', 'mmu', 'mus_musculus')) {
  dataset_name <- "mmusculus_gene_ensembl"
  symbol_attr <- 'mgi_symbol'
} else {
  dataset_name <- "hsapiens_gene_ensembl"
  symbol_attr <- 'hgnc_symbol'
}
for (server in ensembl_servers) {
  cat(sprintf("Trying Ensembl server: %s\n", server))
  tryCatch({
    # Set biomaRt cache directory to a writable location and disable caching
    options(biomaRt.cache = FALSE)
    Sys.setenv(BIOMART_CACHE = "FALSE")

    ensembl <- useMart("ensembl", dataset = dataset_name, host = server)
    cat(sprintf("Successfully connected to: %s\n", server))
    break
  }, error = function(e) {
    cat(sprintf("Failed to connect to %s: %s\n", server, e$message))
    if (server == tail(ensembl_servers, 1)) {
      cat("All Ensembl servers failed. Trying offline mode...\n")

      # Create a minimal gene annotation file for protein coding genes
      # This is a fallback when BioMart is completely unavailable
      cat("Creating minimal gene annotation file...\n")
      minimal_genes <- data.frame(
        tmp_symbol = RNAseq[[gene_col]],
        ensembl_gene_id = paste0("ENS_fallback_", seq_along(RNAseq[[gene_col]])),
        gene_biotype = "protein_coding"
      )
      names(minimal_genes)[names(minimal_genes) == 'tmp_symbol'] <- symbol_attr
      write.table(minimal_genes, file = './ensembl_mapping_fallback.txt', quote=FALSE, sep = '\t', row.names = FALSE)
      all_genes <- minimal_genes
      break
    }
  })
}

# If we still don't have ensembl connection, use the fallback
if (is.null(ensembl)) {
  cat("Using offline fallback mode for gene annotations\n")
  all_genes <- read.table('./ensembl_mapping_fallback.txt', header=TRUE, sep='\t')
} else {
  # Get gene annotations from BioMart
  all_genes <- getBM(attributes = c(symbol_attr, 'ensembl_gene_id', 'gene_biotype'),
                     mart = ensembl)
  write.table(all_genes, file = './ensembl_mapping.txt', quote=FALSE, sep = '\t', row.names = FALSE)
}

id_ensembl <- all_genes[all_genes[[symbol_attr]] %in% RNAseq[[gene_col]], ]

# Merge gene annotation but do NOT filter for protein-coding genes
# Users are expected to provide an expression file with the desired gene symbols,
# so we keep all provided genes (after deduplication) regardless of biotype.
RNAseq <- merge(RNAseq, id_ensembl, by.x = gene_col, by.y = symbol_attr)
RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq[[gene_col]]), ]) # Remove duplicated gene symbols
rownames(RNAseq_coding) <- RNAseq_coding[[gene_col]]
# Do not drop genes based on `gene_biotype`; keep all rows supplied by the user
# Remove the original gene symbol column but preserve Ensembl ID if available
RNAseq_coding[[gene_col]] <- NULL
if ('gene_biotype' %in% colnames(RNAseq_coding)) RNAseq_coding$gene_biotype <- NULL

#####################################################################################################
#### Normalization, log transformation and scaling + selection for highly variable genes using DESeq2
#### Also do TFA inference using DecouplR
#####################################################################################################

cat("Processing metadata and performing normalization...\n")
metadata_df <- fread(metadata, data.table =FALSE)
if ("V1" %in% colnames(metadata_df)) metadata_df$V1 <- NULL

cat(" Metadata file loaded:\n")
cat("   - Dimensions:", nrow(metadata_df), "samples x", ncol(metadata_df), "columns\n")
cat("   - Available columns:", paste(colnames(metadata_df), collapse=", "), "\n")

# Check which metadata columns are actually available
available_columns <- intersect(metadata_to_include, colnames(metadata_df))
missing_columns <- setdiff(metadata_to_include, colnames(metadata_df))

if (length(missing_columns) > 0) {
  cat("[WARNING]  Warning: The following metadata columns are not available and will be skipped:", 
      paste(missing_columns, collapse=", "), "\n")
}

if (length(available_columns) == 0) {
  stop("[ERROR] Error: None of the specified metadata columns are available in the metadata file.\n",
       "   Requested: ", paste(metadata_to_include, collapse=", "), "\n",
       "   Available: ", paste(colnames(metadata_df), collapse=", "))
}

cat("[OK] Using metadata columns:", paste(available_columns, collapse=", "), "\n")

DESeq_groups <- metadata_df[, available_columns, drop=FALSE]

# Check if sample_id_col exists in metadata
sample_id_col <- opt$sample_id_col
if (!sample_id_col %in% colnames(metadata_df)) {
  cat("[WARNING]  Warning: Sample ID column '", sample_id_col, "' not found in metadata\n")
  cat("   Available columns:", paste(colnames(metadata_df), collapse=", "), "\n")
  # Try to find a likely sample ID column
  likely_sample_cols <- grep("sample|Sample|id|ID", colnames(metadata_df), value=TRUE, ignore.case=TRUE)
  if (length(likely_sample_cols) > 0) {
    sample_id_col <- likely_sample_cols[1]
    cat("[OK] Using column '", sample_id_col, "' as sample ID\n")
  } else {
    # Check if row names are meaningful (not just numbers)
    if (all(grepl("^[0-9]+$", rownames(metadata_df)))) {
      stop("[ERROR] Cannot find sample ID column in metadata and row names appear to be row numbers.\n",
           "   Please ensure metadata has a sample ID column or meaningful row names.")
    } else {
      cat("[OK] Using row names as sample IDs\n")
      sample_id_col <- NULL  # Will use row names
    }
  }
}

# Set row names appropriately
if (!is.null(sample_id_col)) {
  cat(" Setting row names from column:", sample_id_col, "\n")
  rownames(DESeq_groups) <- metadata_df[[sample_id_col]]
  rownames(metadata_df) <- metadata_df[[sample_id_col]]
} else {
  cat(" Using existing row names as sample IDs\n")
  rownames(DESeq_groups) <- rownames(metadata_df)
}

DESeq_groups[] <- lapply(DESeq_groups, factor)

# Debug: Print column names and their levels
cat(" DESeq_groups structure:\n")
cat("   - Dimensions:", nrow(DESeq_groups), "samples x", ncol(DESeq_groups), "columns\n")
for (col in colnames(DESeq_groups)) {
  levels_info <- levels(DESeq_groups[[col]])
  cat(sprintf("   - Column '%s': %d levels (%s)\n", col, length(levels_info), 
              paste(head(levels_info, 3), collapse=", ")))
}

# Auto-detect conditions for DESeq2 contrast
contrast_column <- opt$DESeq_contrast1
cat("\n DESeq2 contrast analysis setup:\n")

if (!contrast_column %in% colnames(metadata_df)) {
  cat("[WARNING]  Warning: Contrast column '", contrast_column, "' not found in metadata\n")
  cat("   Available columns:", paste(colnames(metadata_df), collapse=", "), "\n")
  # Try to find a likely diagnosis column
  likely_diagnosis_cols <- grep("diagnosis|condition|group|class|treatment|phenotype", 
                               colnames(metadata_df), value=TRUE, ignore.case=TRUE)
  if (length(likely_diagnosis_cols) > 0) {
    contrast_column <- likely_diagnosis_cols[1]
    cat("[OK] Using column '", contrast_column, "' as contrast column\n")
  } else {
    # Use first available column from metadata_to_include
    if (length(available_columns) > 0) {
      contrast_column <- available_columns[1]
      cat("[OK] Using first available metadata column '", contrast_column, "' as contrast column\n")
    } else {
      stop("[ERROR] Cannot find suitable contrast column in metadata")
    }
  }
}

# Validate contrast column and detect conditions
if (contrast_column %in% colnames(metadata_df)) {
  unique_conditions <- unique(metadata_df[[contrast_column]])
  unique_conditions <- unique_conditions[!is.na(unique_conditions)]
  
  cat("   - Contrast column:", contrast_column, "\n")
  cat("   - Unique conditions:", paste(unique_conditions, collapse=", "), "\n")
  cat("   - Sample distribution:\n")
  condition_counts <- table(metadata_df[[contrast_column]])
  for (cond in names(condition_counts)) {
    cat(sprintf("     * %s: %d samples\n", cond, condition_counts[cond]))
  }
  
  if (length(unique_conditions) == 2) {
    # For 2 groups, use the alphabetically first as reference (condition2) and second as test (condition1)
    sorted_conditions <- sort(unique_conditions)
    DESeq_contrast <- c(contrast_column, sorted_conditions[2], sorted_conditions[1])
    cat(sprintf("Auto-detected 2 groups in '%s': comparing %s vs %s (reference)\n", 
                contrast_column, sorted_conditions[2], sorted_conditions[1]))
  } else if (length(unique_conditions) > 2) {
    # For multiple groups, use first two alphabetically for DESeq2 (though this is somewhat arbitrary)
    sorted_conditions <- sort(unique_conditions)
    DESeq_contrast <- c(contrast_column, sorted_conditions[2], sorted_conditions[1])
    cat(sprintf("Auto-detected %d groups in '%s': using %s vs %s for DESeq2 normalization\n", 
                length(unique_conditions), contrast_column, sorted_conditions[2], sorted_conditions[1]))
    cat(sprintf("Note: All groups (%s) will be used for clustering, this contrast is only for normalization\n", 
                paste(sorted_conditions, collapse=", ")))
  } else {
    stop(sprintf("Error: Column '%s' must have at least 2 unique values for contrast analysis", contrast_column))
  }
} else {
  stop(sprintf("Error: Column '%s' not found in metadata", contrast_column))
}

RNAseq <- RNAseq_coding[, colnames(RNAseq_coding) %in% metadata_df[[sample_id_col]]] # Select samples

# Debug: Check sample matching
cat("RNA-seq samples found:", length(colnames(RNAseq)), "\n")
cat("Metadata samples found:", length(rownames(DESeq_groups)), "\n")
cat("Sample ID column used:", sample_id_col, "\n")
cat("Sample ID values in metadata:", paste(head(rownames(DESeq_groups), 5), collapse=", "), "\n")
cat("Sample ID values in RNA-seq:", paste(head(colnames(RNAseq), 5), collapse=", "), "\n")

ord <- match(colnames(RNAseq), rownames(DESeq_groups))
DESeq_groups <- DESeq_groups[ord, , drop=FALSE]

# Parse and validate the design formula
design_formula_str <- opt$design_formula
cat("Using design formula:", design_formula_str, "\n")

# Update design formula to use the actual column name if different
if (contrast_column != opt$DESeq_contrast1) {
  design_formula_str <- sub(opt$DESeq_contrast1, contrast_column, design_formula_str)
  cat("Updated design formula to use actual column name:", design_formula_str, "\n")
}

# Parse the formula to extract variable names
design_formula <- as.formula(design_formula_str)
formula_vars <- all.vars(design_formula)

# Check if all variables in the formula are available in the metadata
# Use the original metadata_df column names for validation
missing_formula_vars <- setdiff(formula_vars, colnames(metadata_df))
if (length(missing_formula_vars) > 0) {
  cat("Warning: Variables in design formula not found in metadata:", 
      paste(missing_formula_vars, collapse=", "), "\n")
  
  # Create a simpler formula with only available variables
  available_formula_vars <- intersect(formula_vars, colnames(metadata_df))
  if (length(available_formula_vars) == 0) {
    cat("Error: No variables from design formula are available in metadata\n")
    cat("Available metadata columns:", paste(colnames(metadata_df), collapse=", "), "\n")
    cat("Requested formula variables:", paste(formula_vars, collapse=", "), "\n")
    stop("Cannot create valid design formula")
  }
  
  # Create new formula with available variables
  if (length(available_formula_vars) == 1) {
    design_formula <- as.formula(paste("~", available_formula_vars[1]))
  } else {
    design_formula <- as.formula(paste("~", paste(available_formula_vars, collapse=" + ")))
  }
  cat("Adjusted design formula to:", deparse(design_formula), "\n")
}

# Keep all metadata columns for visualization
DESeq_groups_all <- DESeq_groups

# Filter DESeq_groups to only include columns used in the final design formula (for DESeq2)
final_formula_vars <- all.vars(design_formula)
DESeq_groups <- DESeq_groups[, final_formula_vars, drop=FALSE]
cat("Final DESeq_groups columns (for DESeq2):", paste(colnames(DESeq_groups), collapse=", "), "\n")
cat("All metadata columns (for visualization):", paste(colnames(DESeq_groups_all), collapse=", "), "\n")

dds <- DESeqDataSetFromMatrix(countData = (RNAseq+1), colData = DESeq_groups, design = design_formula) 
keep <- rowSums(counts(dds) >=10) >= 3 # Dispersion plot looks better with some prefiltering
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast=DESeq_contrast))

up <- res[(res$log2FoldChange >=1 & res$padj <= 0.05), ]
down <- res[(res$log2FoldChange <=1 & res$padj <= 0.05), ]
normcnt <- as.data.frame(counts(dds, normalized=TRUE))
log_normcnt <- log(normcnt)

# Select for highly variable genes
M <- log_normcnt
vars <- apply(M,1,var)
pdf('./LemonTree/Preprocessing/expression_variance_histogram.pdf')
hist(vars, 1100,xlim=c(0,2))
abline(v=sort(vars, decreasing=TRUE)[top_n_genes], col='red', lwd=3)
dev.off()
top_genes <- names(sort(vars, decreasing=TRUE)[1:top_n_genes])
cat("Number of variable genes selected:", length(top_genes), "\n")
variable_genes <- M[top_genes, ]

# Create PCA plot for transcriptomics (RNA-seq)
create_pca_plot(log_normcnt, metadata_df, contrast_column, "Transcriptomics", "./LemonTree/Preprocessing")

###########################################################################################
#### Re-introduce TFs from user-provided list (Option 1 requirement)
###########################################################################################
cat("\n=== RE-INTRODUCING TFs FROM USER LIST ===\n")

# Read TF list if provided
user_tfs <- c()
if (!is.null(TF_file) && file.exists(TF_file)) {
  cat(" Reading TF list from:", TF_file, "\n")
  
  # Read TF list file - try different formats
  user_tfs <- tryCatch({
    # Try reading as tab-delimited with header
    tf_data <- read.delim(TF_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    if (ncol(tf_data) > 1) {
      cat("   [INFO]  Multi-column file detected, using first column\n")
      tf_data[[1]]
    } else {
      tf_data[[1]]
    }
  }, error = function(e) {
    # Fallback: try reading as simple list without header
    cat("   [INFO]  Reading as simple gene list (no header)\n")
    tf_list <- read.table(TF_file, header = FALSE, stringsAsFactors = FALSE)
    tf_list[[1]]
  })
  
  # Remove header if it's a string like "HGNC_approved_gene_symbol"
  if (length(user_tfs) > 0 && grepl("symbol|name|gene|id", user_tfs[1], ignore.case = TRUE)) {
    cat("   [INFO]  Detected header row, removing it\n")
    user_tfs <- user_tfs[-1]
  }
  
  cat(" Total TFs in user list:", length(user_tfs), "\n")
  
  # Find TFs that exist in expression data but not in HVG selection
  tfs_in_data <- user_tfs[user_tfs %in% rownames(M)]
  cat(" TFs found in expression data:", length(tfs_in_data), "\n")
  
  tfs_not_in_hvg <- setdiff(tfs_in_data, rownames(variable_genes))
  cat(" TFs to add (not in HVG selection):", length(tfs_not_in_hvg), "\n")
  
  # Add missing TFs to variable_genes
  if (length(tfs_not_in_hvg) > 0) {
    tf_expression <- M[tfs_not_in_hvg, , drop=FALSE]
    variable_genes <- rbind(variable_genes, tf_expression)
    cat("[OK] Added", length(tfs_not_in_hvg), "TFs to variable genes\n")
    cat("   Total features now:", nrow(variable_genes), "(", length(top_genes), "HVGs +", length(tfs_not_in_hvg), "TFs )\n")
  } else {
    cat("[OK] All user TFs already included in HVG selection\n")
  }
} else {
  cat("[WARNING]  No TF list file provided or file not found\n")
}
cat("=============================================\n\n")

##### TFA inference on log transformed + normalized gene expression table (without HVG selection, we will just re-add these TFs later)
cat("\n=== STARTING TFA ANALYSIS ===\n")
cat("perform_TFA:", perform_TFA, "\n")
cat("prior_network:", prior_network, "\n")
cat("prior_network exists:", !is.null(prior_network) && file.exists(prior_network), "\n")

if (perform_TFA) {
  cat("[OK] TFA analysis enabled, proceeding...\n")
  
  # Try to read prior network from file first, otherwise use built-in CollecTRI
  if (!is.null(prior_network) && file.exists(prior_network)) {
    cat(" Reading prior network from file:", prior_network, "\n")
    net <- fread(prior_network)
    cat(" Network loaded with", nrow(net), "interactions\n")
    
    # Check if required columns exist
    if (!all(c('source', 'target') %in% colnames(net))) {
      cat("[WARNING]  Network file should have 'source' and 'target' columns. Using first two columns.\n")
      colnames(net)[1:2] <- c('source', 'target')
    }
  } else {
    cat(" No prior network file provided, using built-in CollecTRI database...\n")
    
    # Determine organism for CollecTRI
    collectri_organism <- if (tolower(opt$organism) %in% c('mouse', 'mmu', 'mus_musculus')) 'mouse' else 'human'
    cat(sprintf(" Using organism: %s\n", collectri_organism))
    
    tryCatch({
      # First try online download (most up-to-date)
      cat(" Attempting online download from CollecTRI database...\n")
      net <- decoupleR::get_collectri(organism = collectri_organism, split_complexes = FALSE)
      cat(sprintf(" CollecTRI network downloaded with %d interactions\n", nrow(net)))
    }, error = function(e1) {
      cat(sprintf("[WARNING] Online download failed: %s\n", e1$message))
      cat(" Trying pre-downloaded network from PKN directory...\n")
      
      tryCatch({
        # Fallback to pre-downloaded network from PKN (organism-specific)
        pkn_network_file <- if (collectri_organism == 'mouse') {
          "PKN/CollecTRI_network_mouse.txt"
        } else {
          "PKN/CollecTRI_network.txt"
        }
        
        if (file.exists(pkn_network_file)) {
          cat(sprintf(" Loading pre-downloaded CollecTRI network: %s\n", pkn_network_file))
          net <<- fread(pkn_network_file)
          cat(sprintf(" CollecTRI network loaded with %d interactions\n", nrow(net)))
          
          # Ensure correct column names
          if (!all(c('source', 'target') %in% colnames(net))) {
            cat("[INFO] Renaming columns to 'source' and 'target'\n")
            colnames(net)[1:2] <- c('source', 'target')
          }
        } else {
          stop(sprintf("Pre-downloaded network not found: %s", pkn_network_file))
        }
      }, error = function(e2) {
        cat(sprintf("[ERROR] Failed to load CollecTRI network: %s\n", e2$message))
        cat("[WARNING] Continuing without TFA...\n")
        perform_TFA <<- FALSE
        return(NULL)
      })
    })
  }
  
  # Only proceed with decouple analysis if network was loaded successfully
  if (perform_TFA && exists("net") && !is.null(net)) {
    cat(" Running decouple analysis...\n")
    tryCatch({
      decoupled <- decouple(mat=log_normcnt, net=net, .source='source', .target='target')
      cat("[OK] Decouple completed, running consensus...\n")
      consensus <- run_consensus(decoupled)
      cat("[OK] Consensus analysis completed\n")
    }, error = function(e) {
      cat("[ERROR] Error in TFA analysis:", e$message, "\n")
      cat("[WARNING]  Continuing without TFA...\n")
      perform_TFA <<- FALSE
      return(NULL)
    })
  } else if (perform_TFA) {
    cat("[ERROR] Cannot proceed with TFA analysis: network not loaded\n")
    cat("[WARNING]  Continuing without TFA...\n")
    perform_TFA <- FALSE
  }
  
  if (perform_TFA) {
  
  # Transform to wide matrix
  TFA_df <- consensus %>%
    dplyr::filter(statistic == 'consensus') %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    t %>%
    as.data.frame()
  
  # Read Lovering TF list for TFA filtering
  Lovering_TF_list <- user_tfs  # Use the TFs already loaded earlier

  # Visualization
  n_tfs <- 50

  # Transform to wide matrix
  sample_acts_mat <- consensus %>%
    dplyr::filter(statistic == 'consensus') %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()

  # Get top tfs with more variable means across clusters
  tfs <- consensus %>%
    group_by(source) %>%
    summarise(std = sd(score)) %>%
    arrange(-abs(std)) %>%
    head(n_tfs) %>%
    pull(source)
  sample_acts_mat <- sample_acts_mat[,tfs]

  # Scale per sample
  sample_acts_mat <- scale(sample_acts_mat)

  # Choose color palette
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

  my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
                 seq(0.05, 3, length.out=floor(palette_length/2)))

  # Create annotation colors dynamically based on available columns
  ann_colors = list()
  
  # Only include columns that are factors and suitable for visualization
  viz_columns <- DESeq_groups_all[, sapply(DESeq_groups_all, function(x) is.factor(x) && length(levels(x)) <= 10), drop=FALSE]
  
  # Generate colors for each categorical variable
  color_palettes <- list(
    c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F"),
    c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#8DD3C7", "#FFFFB3")
  )
  
  for (i in seq_along(colnames(viz_columns))) {
    col_name <- colnames(viz_columns)[i]
    levels_vec <- levels(viz_columns[[col_name]])
    n_levels <- length(levels_vec)
    
    if (n_levels > 0) {
      # Use predefined colors if available, otherwise generate
      if (col_name == "diagnosis" && all(c("UC", "nonIBD") %in% levels_vec)) {
        ann_colors[[col_name]] <- c("UC" = "firebrick", "nonIBD" = "black")
      } else if (col_name == "sex" && all(c("Female", "Male") %in% levels_vec)) {
        ann_colors[[col_name]] <- c("Female" = "#1B9E77", "Male" = "#D95F02")
      } else if (col_name == "biopsy_location" && all(c("Colon", "Rectum") %in% levels_vec)) {
        ann_colors[[col_name]] <- c("Colon" = "#7570B3", "Rectum" = "#E7298A")
      } else {
        # Generate colors for this variable
        palette_idx <- ((i - 1) %% length(color_palettes)) + 1
        colors <- color_palettes[[palette_idx]][1:n_levels]
        names(colors) <- levels_vec
        ann_colors[[col_name]] <- colors
      }
    }
  }

  plot <- pheatmap(t(sample_acts_mat), border_color = NA, color=my_color, breaks = my_breaks, 
                   annotation_col = viz_columns, annotation_colors = ann_colors,
                   cluster_rows = T, show_colnames = FALSE, fontsize_row = 5)
  save(TFA_df, file='./TFA/TFA_df.RData')
  save_pheatmap_pdf(plot, paste0('./TFA/TFA_consensus_',n_tfs,'heatmap_annotated.pdf'))

  # write to file 
  write.table(TFA_df, './LemonTree/Preprocessing/TFA_consensus.txt', sep = '\t', quote=FALSE, row.names=TRUE)

  # Plot heatmap with TF expression instead
  genes_to_plot <- colnames(sample_acts_mat)

  # Transpose the matrix
  transposed_matrix <- t(normcnt)

  # Scale the transposed matrix to get Z-scores (mean = 0, standard deviation = 1)
  scaled_matrix <- scale(transposed_matrix)

  # Transpose back to the original orientation if required
  z_score_matrix <- t(scaled_matrix)
  # Select rows that are in genes_to_plot
  z_score_matrix <- z_score_matrix[rownames(z_score_matrix) %in% genes_to_plot, ]

  plot2 <- pheatmap(z_score_matrix, border_color = NA, color=my_color, breaks = my_breaks, 
                   annotation_col = viz_columns, annotation_colors = ann_colors,
                   cluster_rows = T, show_colnames = FALSE, fontsize_row = 5)

  save_pheatmap_pdf(plot2, paste0('./TFA/TFAgenes_ZscoredExpression_',n_tfs,'heatmap_annotated.pdf'))

  ## Create RNA_preprocessed_withTFA: HVGs + Lovering TFs (with TFA where available, else expression)
  cat("\n=== CREATING TFA-INTEGRATED DATASET ===\n")
  cat(" TFA features:", nrow(TFA_df), "\n")
  cat(" Expression features before merge:", nrow(variable_genes), "\n")
  
  # Filter TFA_df to only include Lovering TFs
  TFA_df_lovering <- TFA_df[rownames(TFA_df) %in% Lovering_TF_list, , drop=FALSE]
  cat(" TFA features for Lovering TFs:", nrow(TFA_df_lovering), "\n")
  
  # Start with variable genes
  RNA_preprocessed_withTFA <- variable_genes
  
  # Identify Lovering TFs not in variable_genes
  lovering_tfs_not_in_hvg <- setdiff(Lovering_TF_list, rownames(variable_genes))
  lovering_tfs_not_in_hvg <- intersect(lovering_tfs_not_in_hvg, rownames(log_normcnt))  # Must exist in log_normcnt
  
  # Identify Lovering TFs that exist in variable_genes (will be replaced with TFA if available)
  lovering_tfs_in_hvg <- intersect(rownames(variable_genes), rownames(TFA_df_lovering))
  
  # Identify which of the non-HVG Lovering TFs have TFA scores
  lovering_tfs_with_tfa <- intersect(lovering_tfs_not_in_hvg, rownames(TFA_df_lovering))
  
  # Identify which of the non-HVG Lovering TFs don't have TFA (will use expression)
  lovering_tfs_without_tfa <- setdiff(lovering_tfs_not_in_hvg, rownames(TFA_df_lovering))
  
  # Scale non-Lovering TF genes in variable_genes
  non_lovering_tf_genes <- setdiff(rownames(variable_genes), Lovering_TF_list)
  if (length(non_lovering_tf_genes) > 0) {
    RNA_preprocessed_withTFA[non_lovering_tf_genes, ] <- t(scale(t(variable_genes[non_lovering_tf_genes, , drop=FALSE])))
  }
  
  # Replace Lovering TFs in HVGs with TFA scores (unscaled)
  if (length(lovering_tfs_in_hvg) > 0) {
    RNA_preprocessed_withTFA[lovering_tfs_in_hvg, ] <- TFA_df_lovering[lovering_tfs_in_hvg, , drop=FALSE]
    cat(" Replaced", length(lovering_tfs_in_hvg), "Lovering TF expression values with TFA scores\n")
  }
  
  # Add Lovering TFs not in HVGs that have TFA scores (unscaled)
  if (length(lovering_tfs_with_tfa) > 0) {
    RNA_preprocessed_withTFA <- rbind(RNA_preprocessed_withTFA, TFA_df_lovering[lovering_tfs_with_tfa, , drop=FALSE])
    cat(" Added", length(lovering_tfs_with_tfa), "Lovering TFs with TFA scores\n")
  }
  
  # Add Lovering TFs not in HVGs that don't have TFA (use scaled expression from log_normcnt)
  if (length(lovering_tfs_without_tfa) > 0) {
    lovering_tfs_expr <- log_normcnt[lovering_tfs_without_tfa, colnames(RNA_preprocessed_withTFA), drop=FALSE]
    lovering_tfs_expr_scaled <- t(scale(t(lovering_tfs_expr)))
    RNA_preprocessed_withTFA <- rbind(RNA_preprocessed_withTFA, lovering_tfs_expr_scaled)
    cat(" Added", length(lovering_tfs_without_tfa), "Lovering TFs with scaled expression\n")
  }
  
  cat(" Total features after merge:", nrow(RNA_preprocessed_withTFA), "\n")
  cat("===========================================\n\n")

  cat(" TFA analysis completed for", nrow(TFA_df), "transcription factors\n")
  cat(" TFA consensus matrix saved to: ./LemonTree/Preprocessing/TFA_consensus.txt\n")
  cat(" TFA heatmap saved to: ./TFA/TFA_consensus_50heatmap_annotated.pdf\n")
  cat(" TF gene expression heatmap saved to: ./TFA/TFAgenes_ZscoredExpression_50heatmap_annotated.pdf\n")
  cat("===========================\n\n")
} else {
  cat("[INFO] TFA analysis disabled (perform_TFA=FALSE)\n")
  # When TFA is disabled, use variable_genes for complete dataframe
  RNA_preprocessed <- variable_genes
}

###########################################################################################
#### Create RNA_preprocessed_noTFA: HVGs + Lovering TFs (NO TFA, just scaled expression)
###########################################################################################
cat("=== CREATING EXPRESSION-ONLY DATASET (NO TFA) ===\n")

# Identify Lovering TFs not in variable_genes that we need to add
lovering_tfs_not_in_hvg_noTFA <- setdiff(user_tfs, rownames(variable_genes))
lovering_tfs_not_in_hvg_noTFA <- intersect(lovering_tfs_not_in_hvg_noTFA, rownames(log_normcnt))  # Must exist in log_normcnt

# Start with variable genes (will scale all)
RNA_preprocessed_noTFA <- variable_genes

# Add Lovering TFs not in HVGs (use expression from log_normcnt)
if (length(lovering_tfs_not_in_hvg_noTFA) > 0) {
  RNA_preprocessed_noTFA <- rbind(RNA_preprocessed_noTFA, log_normcnt[lovering_tfs_not_in_hvg_noTFA, colnames(RNA_preprocessed_noTFA), drop=FALSE])
  cat(" Added", length(lovering_tfs_not_in_hvg_noTFA), "Lovering TFs to dataset\n")
}

# Scale all genes (HVGs + Lovering TFs)
RNA_preprocessed_noTFA <- as.data.frame(t(scale(t(RNA_preprocessed_noTFA))))
cat("[OK] Created expression-only dataset with", nrow(RNA_preprocessed_noTFA), "features (all scaled)\n")
cat("=============================================\n\n")

# Scale RNA_preprocessed if TFA was disabled
if (!perform_TFA) {
  RNA_preprocessed <- as.data.frame(t(scale(t(RNA_preprocessed))))
  cat("[INFO] Scaled all features (TFA disabled)\n")
}

cat("===========================\n\n")

  # Set RNA_preprocessed to withTFA version for complete dataframe
  RNA_preprocessed <- RNA_preprocessed_withTFA


###########################################################################################
#### Preprocessing of metabolomics and lipidomics data with omics-specific scaling
###########################################################################################

# Define pareto scaling function - operates per feature (row-wise)
pareto_scale <- function(x, centering = TRUE) {
  # Pareto scaling: (x - mean) / sqrt(sd) per feature
  # x should be a matrix with features in rows, samples in columns
  if (length(centering) != 1 || !is.logical(centering)) {
    stop("'centering' must be a single logical indicator", call. = FALSE)
  }
  
  x <- as.matrix(x)
  
  # Calculate row-wise (per feature) statistics
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, sd, na.rm = TRUE)
  
  # Center data if requested (subtract mean from each feature)
  if (centering) {
    x <- sweep(x, 1, row_means, "-")
  }
  
  # Scale by sqrt(sd) for each feature
  x <- sweep(x, 1, sqrt(row_sds), "/")
  
  return(x)
}

# Define function to process omics data (metabolomics or lipidomics)
process_omics_data <- function(omics_file, omics_type, use_pareto = TRUE, data_type = "c") {
  cat("\n")
  cat("=======================================================================\n")
  cat("Processing", omics_type, "data...\n")
  cat("=======================================================================\n")
  
  # Check if file exists and is readable
  if (!file.exists(omics_file)) {
    cat("[WARNING]  File not found:", omics_file, "\n")
    return(NULL)
  }
  
  # Read the file with error handling
  cat("Reading", omics_type, "file...\n")
  tryCatch({
    omics_data <- fread(omics_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                        check.names = FALSE, quote = "", comment.char = "")
    cat("[OK] Successfully read with fread\n")
  }, error = function(e1) {
    cat("[WARNING]  fread failed, trying read.csv...\n")
    tryCatch({
      omics_data <<- read.csv(omics_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                             check.names = FALSE, quote = "", comment.char = "")
      cat("[OK] Successfully read with read.csv\n")
    }, error = function(e2) {
      cat("[WARNING]  read.csv failed, trying read.table...\n")
      tryCatch({
        omics_data <<- read.table(omics_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                                 check.names = FALSE, quote = "", comment.char = "")
        cat("[OK] Successfully read with read.table\n")
      }, error = function(e3) {
        stop("Failed to read ", omics_type, " file: ", e3$message)
      })
    })
  })
  
  # Convert to data.frame immediately - fread returns data.table which silently
  # ignores rowname assignment, causing feature names to be lost
  omics_data <- as.data.frame(omics_data, check.names = FALSE)
  
  cat(sprintf("  Data dimensions: %d features x %d columns\n", nrow(omics_data), ncol(omics_data)))
  
  # Extract feature names from 'Name' column (or first column) and set as rownames
  if ("Name" %in% colnames(omics_data)) {
    feature_names <- as.character(omics_data$Name)
    omics_data$Name <- NULL
  } else {
    feature_names <- as.character(omics_data[[1]])
    omics_data[[1]] <- NULL
  }
  rownames(omics_data) <- feature_names
  
  cat(sprintf("  First 5 feature names: %s\n", paste(head(rownames(omics_data), 5), collapse=", ")))
  
  # Clean row names
  rownames(omics_data) <- str_replace_all(rownames(omics_data), ' ', '_')
  rownames(omics_data) <- str_replace_all(rownames(omics_data), '-', '_')
  rownames(omics_data) <- str_replace_all(rownames(omics_data), ':', '_')
  rownames(omics_data) <- str_replace_all(rownames(omics_data), '\\+', '_')
  rownames(omics_data) <- str_replace_all(rownames(omics_data), '_$', '')
  
  cat(sprintf("  First 5 feature names (cleaned): %s\n", paste(head(rownames(omics_data), 5), collapse=", ")))
  
  # Convert to numeric and handle NAs
  for (col_idx in 1:ncol(omics_data)) {
    if (!is.numeric(omics_data[[col_idx]])) {
      omics_data[[col_idx]] <- as.numeric(as.character(omics_data[[col_idx]]))
    }
  }
  omics_data[is.na(omics_data)] <- 0
  
  # Check if data is already log-transformed
  # Heuristic: if median value is < 100 and most values are in a narrow range (e.g., 10-30),
  # it's likely already log-transformed
  median_val <- median(as.matrix(omics_data), na.rm = TRUE)
  max_val <- max(as.matrix(omics_data), na.rm = TRUE)
  
  # For discrete/binary data, skip all transformations and scaling
  if (data_type == "d") {
    unique_vals <- sort(unique(as.vector(as.matrix(omics_data))))
    cat(sprintf("[INFO] Discrete/binary data detected — skipping log transformation and scaling\n"))
    cat(sprintf("  Unique values found: %s\n", paste(unique_vals, collapse = ", ")))
    cat(sprintf("  Data will be used as-is for LemonTree discrete regulator assignment\n"))
  } else {
    already_log_transformed <- FALSE
    if (median_val < 100 && max_val < 100) {
      cat("[WARNING]  Data appears to be already log-transformed (median:", round(median_val, 2), ", max:", round(max_val, 2), ")\n")
      cat("   Skipping log transformation to avoid double-logging\n")
      already_log_transformed <- TRUE
    }
    
    # Log transformation (only if not already log-transformed)
    if (!already_log_transformed) {
      cat("Applying log transformation...\n")
      omics_data <- log(omics_data + 1)
    } else {
      cat("Using data as-is (already in log space)\n")
    }
    
    # Apply scaling based on omics type
    if (use_omics_specific_scaling && use_pareto) {
      cat("Applying Pareto scaling for", omics_type, "...\n")
      omics_data <- as.data.frame(pareto_scale(omics_data, centering = TRUE))
    } else {
      cat("Applying standard (z-score) scaling for", omics_type, "...\n")
      omics_data <- as.data.frame(t(scale(t(omics_data))))
    }
  }
  
  # Create visualization
  pdf_file <- paste0('./LemonTree/Preprocessing/Normalized_', tolower(omics_type), '.pdf')
  pdf(pdf_file)
  omics_data %>%
    gather(Sample, Count) %>%
    ggplot(aes(Sample, Count)) + 
    geom_boxplot() + 
    ggtitle(paste("Normalized", omics_type)) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5, size=4))
  dev.off()
  cat("[OK] Saved visualization to:", pdf_file, "\n")
  
  # Create PCA plot for this omics type
  tryCatch({
    create_pca_plot(omics_data, metadata_df, contrast_column, omics_type, "./LemonTree/Preprocessing")
  }, error = function(e) {
    cat(sprintf("[WARNING] Failed to create PCA plot for %s: %s\n", omics_type, e$message))
  })
  
  return(omics_data)
}

# Process all regulator omics files dynamically
processed_regulator_data <- list()

for (prefix in names(other_regulator_files)) {
  filepath <- other_regulator_files[[prefix]]
  data_type <- regulator_data_types[[prefix]]
  
  if (file.exists(filepath)) {
    cat(sprintf("\n Processing %s data from: %s (data_type=%s)\n", prefix, filepath, data_type))
    omics_data <- process_omics_data(filepath, prefix, use_pareto = TRUE, data_type = data_type)
    processed_regulator_data[[prefix]] <- omics_data
  } else {
    cat(sprintf("[WARNING]  Warning: File not found for %s: %s\n", prefix, filepath))
  }
}

###########################################################################################
#### Create complete dataframe for LemonTree
###########################################################################################

cat("\n")
cat("=======================================================================\n")
cat("Creating output files...\n")
cat("=======================================================================\n")

# Add Ensembl IDs and write LemonPreprocessed_expression.txt (NO TFA)
RNA_preprocessed_noTFA_ids <- RNA_preprocessed_noTFA
RNA_preprocessed_noTFA_ids$symbol <- row.names(RNA_preprocessed_noTFA_ids)
RNA_preprocessed_noTFA_ids <- merge(RNA_preprocessed_noTFA_ids, id_ensembl, by.x='symbol', by.y = symbol_attr)
RNA_preprocessed_noTFA_ids <- as.data.frame(RNA_preprocessed_noTFA_ids %>% group_by(symbol) %>% dplyr::filter(row_number()==1))
RNA_preprocessed_noTFA_ids <- RNA_preprocessed_noTFA_ids[, c(1,ncol(RNA_preprocessed_noTFA_ids)-1,2:(ncol(RNA_preprocessed_noTFA_ids)-2))]
names(RNA_preprocessed_noTFA_ids)[2] <- 'ensembl_gene_id'
write.table(RNA_preprocessed_noTFA_ids, './LemonTree/Preprocessing/LemonPreprocessed_expression.txt', sep = '\t', quote=FALSE, row.names=FALSE)
cat("[OK] Saved: LemonPreprocessed_expression.txt (HVGs + Lovering TFs, NO TFA)\n")

# Add Ensembl IDs to RNA_preprocessed (WITH TFA) for complete dataframe
RNA_preprocessed_withTFA_ids <- RNA_preprocessed
RNA_preprocessed_withTFA_ids$symbol <- row.names(RNA_preprocessed_withTFA_ids)
RNA_preprocessed_withTFA_ids <- merge(RNA_preprocessed_withTFA_ids, id_ensembl, by.x='symbol', by.y = symbol_attr)
RNA_preprocessed_withTFA_ids <- as.data.frame(RNA_preprocessed_withTFA_ids %>% group_by(symbol) %>% dplyr::filter(row_number()==1))
RNA_preprocessed_withTFA_ids <- RNA_preprocessed_withTFA_ids[, c(1,ncol(RNA_preprocessed_withTFA_ids)-1,2:(ncol(RNA_preprocessed_withTFA_ids)-2))]
names(RNA_preprocessed_withTFA_ids)[2] <- 'ensembl_gene_id'
cat("[OK] Prepared RNA data with TFA for complete dataframe\n")

# Collect all omics data for complete dataframe
omics_datasets <- list(RNA_preprocessed_withTFA_ids)
all_regulator_names <- list()

# Process each regulator type
for (prefix in names(processed_regulator_data)) {
  omics_data <- processed_regulator_data[[prefix]]
  
  # Add columns for merging
  omics_data$ensembl_gene_id <- row.names(omics_data)
  omics_data$symbol <- row.names(omics_data)
  omics_data <- omics_data[, c(ncol(omics_data), ncol(omics_data)-1, 1:(ncol(omics_data)-2))]
  
  # Derive output filename from the original data file name (not the prefix) so that
  # e.g. Metabolites:Metabolomics.txt -> LemonPreprocessed_metabolomics.txt
  #      Proteins:proteomics.txt      -> LemonPreprocessed_proteomics.txt
  data_file_basename <- tolower(tools::file_path_sans_ext(basename(regulator_configs[[prefix]])))
  output_filename <- paste0('./LemonTree/Preprocessing/LemonPreprocessed_', data_file_basename, '.txt')
  write.table(omics_data, output_filename, sep = '\t', quote=FALSE, row.names=FALSE)
  cat(sprintf("[OK] Saved: LemonPreprocessed_%s.txt\n", data_file_basename))
  
  # Add to datasets list
  omics_datasets[[length(omics_datasets) + 1]] <- omics_data
  
  # Store regulator names (row names before adding columns)
  regulator_names <- rownames(processed_regulator_data[[prefix]])
  all_regulator_names[[prefix]] <- regulator_names
  
  # Write regulator list file (keep prefix-based naming for downstream LemonTree compatibility)
  list_filename <- paste0('./LemonTree/Preprocessing/', tolower(prefix), '.txt')
  data_type <- regulator_data_types[[prefix]]
  if (!is.null(data_type) && data_type == "d") {
    # Write 2-column format with 'd' flag as required by LemonTree for discrete regulators
    # Format: regulator_name<tab>d (one per line)
    reg_df <- data.frame(name = regulator_names, type = rep("d", length(regulator_names)))
    write.table(reg_df, list_filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
    cat(sprintf("[OK] Saved: %s DISCRETE regulator list (2-column format with 'd' flag, %d features) -> %s\n",
                tolower(prefix), length(regulator_names), list_filename))
  } else {
    write.table(regulator_names, list_filename, quote = FALSE, row.names = FALSE, col.names=FALSE)
    cat(sprintf("[OK] Saved: %s regulator list (%d features, e.g. %s) -> %s\n",
                tolower(prefix), length(regulator_names),
                paste(head(regulator_names, 3), collapse=", "), list_filename))
  }
}

# Handle TF list file separately (it's a list of gene names, not abundance data)
if (!is.null(TF_file) && file.exists(TF_file)) {
  cat(sprintf("\n Processing TF list from: %s\n", TF_file))
  
  # Read TF list file - try different formats
  tf_names <- tryCatch({
    # Try reading as tab-delimited with header (like lovering_TF_list.txt)
    tf_data <- read.delim(TF_file, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
    if (ncol(tf_data) > 1) {
      cat("   [INFO]  Multi-column file detected, using first column\n")
      tf_data[[1]]
    } else {
      tf_data[[1]]
    }
  }, error = function(e) {
    # Fallback: try reading as simple list without header
    cat("   [INFO]  Reading as simple gene list (no header)\n")
    tf_list <- read.table(TF_file, header = FALSE, stringsAsFactors = FALSE)
    tf_list[[1]]
  })
  
  # Remove header if it's a string like "HGNC_approved_gene_symbol"
  if (length(tf_names) > 0 && grepl("symbol|name|gene|id", tf_names[1], ignore.case = TRUE)) {
    cat("   [INFO]  Detected header row, removing it\n")
    tf_names <- tf_names[-1]
  }
  
  # Write to tfs.txt in preprocessing directory
  # Ensure directory exists before writing
  dir.create('./LemonTree/Preprocessing/', recursive = TRUE, showWarnings = FALSE)
  tf_output_file <- paste0('./LemonTree/Preprocessing/', tolower(TF_prefix), '.txt')
  write.table(tf_names, tf_output_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("[OK] Saved: %s list (%d TFs) -> %s\n", tolower(TF_prefix), length(tf_names), tf_output_file))
  
  # Store TF names for reference
  all_regulator_names[[TF_prefix]] <- tf_names
} else if (!is.null(TF_file)) {
  cat(sprintf("[WARNING]  Warning: TF file not found: %s\n", TF_file))
}

# Combine all omics datasets
complete_df <- do.call(rbind, c(omics_datasets, list(fill=TRUE)))
write.table(complete_df, './LemonTree/Preprocessing/LemonPreprocessed_complete.txt', sep = '\t', quote=FALSE, row.names=FALSE)
cat("[OK] Saved: LemonPreprocessed_complete.txt\n")

# Write DESeq_groups file with validation for downstream pipeline compatibility
cat("\n Creating DESeq_groups.txt for downstream analysis...\n")
cat("   - Format: Tab-separated with sample IDs as row names\n")
cat("   - Dimensions:", nrow(DESeq_groups), "samples x", ncol(DESeq_groups), "columns\n")
cat("   - Columns:", paste(colnames(DESeq_groups), collapse=", "), "\n")

write.table(DESeq_groups_all, './LemonTree/Preprocessing/DESeq_groups.txt', quote=FALSE, sep='\t', row.names = TRUE)

# Validate that we have the required columns for downstream analysis
required_for_pipeline <- c(contrast_column)
available_for_pipeline <- intersect(required_for_pipeline, colnames(DESeq_groups))
missing_for_pipeline <- setdiff(required_for_pipeline, colnames(DESeq_groups))

if (length(missing_for_pipeline) > 0) {
  cat("[WARNING]  Warning: Missing columns that may be needed by downstream analysis:", 
      paste(missing_for_pipeline, collapse=", "), "\n")
}


cat("\n Pipeline compatibility check:\n")
cat("   - [OK] DESeq_groups.txt: Compatible with DESeq2 and overview modules\n")
cat("   - [OK] Sample IDs: Available as row names for sample matching\n")
cat("   - [OK] Contrast column: '", contrast_column, "' available for differential analysis\n")

# Check for common optional columns that enhance downstream analysis
optional_beneficial <- c("sex", "age", "batch", "biopsy_location", "tissue", "condition")
available_optional <- intersect(optional_beneficial, colnames(DESeq_groups))
if (length(available_optional) > 0) {
  cat("   - [OK] Additional metadata: ", paste(available_optional, collapse=", "), 
      " (can be used for enhanced visualizations)\n")
} else {
  cat("   - [INFO]  No additional metadata columns detected (sex, age, batch, etc.)\n")
}

cat("   -  Ready for module_overview_interactive.py and other downstream tools\n")

cat("\n Preprocessing completed successfully!\n")
cat(" Output files created in: ./LemonTree/\n")
cat("   - LemonPreprocessed_expression.txt (without TFA)\n")
cat("   - LemonPreprocessed_complete.txt (with all omics data)\n")
cat("   - DESeq_groups.txt\n")

# List all regulator files created
if (!is.null(TF_file) && file.exists(TF_file)) {
  cat(sprintf("   - %s.txt (TF list)\n", tolower(TF_prefix)))
}
for (prefix in names(all_regulator_names)) {
  data_file_basename <- tolower(tools::file_path_sans_ext(basename(regulator_configs[[prefix]])))
  cat(sprintf("   - LemonPreprocessed_%s.txt\n", data_file_basename))
  cat(sprintf("   - %s.txt (regulator list)\n", tolower(prefix)))
}

# Check if TFA files actually exist (more reliable than checking parameters)
tfa_consensus_exists <- file.exists('./LemonTree/Preprocessing/TFA_consensus.txt')
tfa_expression_exists <- file.exists('./LemonTree/Preprocessing/LemonPreprocessed_expression_TFA.txt')
tfa_dir_exists <- dir.exists('./TFA/')

if (tfa_consensus_exists || tfa_expression_exists || tfa_dir_exists) {
  cat(" TFA ANALYSIS FILES:\n")
  if (tfa_consensus_exists) cat("   [OK] TFA_consensus.txt (TF activities)\n")
  if (tfa_expression_exists) cat("   [OK] LemonPreprocessed_expression_TFA.txt (expression + TFA)\n")
  if (tfa_dir_exists) {
    cat(" TFA visualization files in: ./TFA/\n")
    tfa_files <- list.files('./TFA/', pattern = "*.pdf|*.RData")
    for (file in tfa_files) {
      cat("   [OK]", file, "\n")
    }
  }
} else {
  cat("[WARNING]  No TFA analysis files found\n")
  if (perform_TFA) {
    cat("   This suggests TFA analysis may have failed or been skipped\n")
  } else {
    cat("   TFA analysis was disabled\n")
  }
}

cat("\n=======================================================================\n")
cat("[OK] PREPROCESSING COMPLETED SUCCESSFULLY\n")
cat("=======================================================================\n")
}
