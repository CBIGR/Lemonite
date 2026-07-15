#!/usr/bin/env Rscript

# ============================================================================
# Preprocessing script for pre-scaled proteomics/metabolomics data
# ============================================================================
# Adapted from Preprocessing_and_TFA_clean.R for use in Nextflow pipeline
# Handles pre-scaled continuous data without DESeq2
# Supports optional hPTM and metabolomics preprocessing

library(optparse)
library(data.table)
library(tidyverse)
library(ggplot2)

# Command line options
option_list <- list(
  make_option(c("--expression"), type="character", default=NULL,
              help="Expression file path (required)", metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Metadata file path (required)", metavar="character"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory [default= %default]", metavar="character"),
  make_option(c("--top_n_genes"), type="integer", default=5000,
              help="Number of top variable genes to select [default= %default]", metavar="integer"),
  make_option(c("--id_mapping"), type="character", default=NULL,
              help="Gene ID mapping file (format: old_id\\tnew_id)", metavar="character"),
  make_option(c("--hptm_file"), type="character", default=NULL,
              help="hPTM file path (optional)", metavar="character"),
  make_option(c("--metabolomics_file"), type="character", default=NULL,
              help="Metabolomics file path (optional)", metavar="character"),
  make_option(c("--metabolomics_labels"), type="character", default=NULL,
              help="Metabolomics name mapping file (optional)", metavar="character"),
  make_option(c("--sample_id_col"), type="character", default="ID",
              help="Sample ID column in metadata [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$expression)) {
  print_help(opt_parser)
  stop("--expression is required")
}
if (is.null(opt$metadata)) {
  print_help(opt_parser)
  stop("--metadata is required")
}

cat("=== Preprocessing Pre-scaled Omics Data ===\n")
cat("Expression file:", opt$expression, "\n")
cat("Metadata file:", opt$metadata, "\n")
cat("Output dir:", opt$output_dir, "\n")

# Create output directory
output_dir <- opt$output_dir
preprocessing_dir <- file.path(output_dir, "LemonTree/Preprocessing")
dir.create(preprocessing_dir, recursive = TRUE, showWarnings = FALSE)

# Helper function for Pareto scaling (for metabolomics)
pareto_scale <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, sd, na.rm = TRUE)
  return((x - row_means) / sqrt(row_sds))
}

# ============================================================================
# Read and preprocess PROTEOMICS data
# ============================================================================

cat("\nReading proteomics data...\n")
proteomics <- fread(opt$expression, header=TRUE, data.table=FALSE)

# Set first column as row names if it contains gene symbols
rownames(proteomics) <- proteomics[, 1]
proteomics[, 1] <- NULL
cat("Proteomics dimensions:", nrow(proteomics), "proteins x", ncol(proteomics), "samples\n")

# ============================================================================
# Read metadata and format sample IDs
# ============================================================================

cat("\nReading metadata...\n")
metadata <- fread(opt$metadata, data.table=FALSE)
if ("V1" %in% colnames(metadata)) metadata$V1 <- NULL

# Set sample ID column as row names
sample_col <- opt$sample_id_col
if (sample_col %in% colnames(metadata)) {
  rownames(metadata) <- metadata[[sample_col]]
  metadata[[sample_col]] <- NULL
} else {
  rownames(metadata) <- metadata[[1]]
  metadata[[1]] <- NULL
}

# Remove date prefix from sample IDs if present (e.g., "240221_FKH1_1" -> "FKH1_1")
if (any(grepl("^[0-9]+_", rownames(metadata)))) {
  rownames(metadata) <- sub("^[0-9]+_", "", rownames(metadata))
}
cat("Metadata samples:", nrow(metadata), "\n")

# Match samples between expression and metadata
proteomics <- proteomics[, colnames(proteomics) %in% rownames(metadata), drop=FALSE]
ord <- match(colnames(proteomics), rownames(metadata))
metadata <- metadata[ord, , drop=FALSE]
cat("Matched samples:", ncol(proteomics), "\n")

# Convert metadata columns to factors
metadata[] <- lapply(metadata, factor)

# ============================================================================
# Select top variable genes
# ============================================================================

cat("\nSelecting top", opt$top_n_genes, "variable genes...\n")
vars <- apply(proteomics, 1, var, na.rm=TRUE)
n_genes <- min(opt$top_n_genes, nrow(proteomics))
top_genes <- names(sort(vars, decreasing=TRUE)[1:n_genes])
proteomics <- proteomics[top_genes, ]
cat("Selected", nrow(proteomics), "genes\n")

# ============================================================================
# Add gene IDs for LemonTree (Protein_id column)
# ============================================================================

cat("\nAdding gene ID mapping...\n")
proteomics$Gene_symbol <- rownames(proteomics)

if (!is.null(opt$id_mapping) && file.exists(opt$id_mapping)) {
  cat("Reading ID mapping from:", opt$id_mapping, "\n")
  id_mapping <- fread(opt$id_mapping, header=TRUE, data.table=FALSE)
  id_mapping <- id_mapping[!duplicated(id_mapping[[1]]), ]
  proteomics <- merge(proteomics, id_mapping, by.x="Gene_symbol", by.y=colnames(id_mapping)[1], all.x=TRUE)
  proteomics$Protein_id <- if ("Protein_id" %in% colnames(proteomics)) proteomics$Protein_id else proteomics$Gene_symbol
} else {
  # Fallback: use gene symbol as ID
  proteomics$Protein_id <- proteomics$Gene_symbol
}

# Reorder columns: Gene_symbol, Protein_id, then samples
all_cols <- colnames(proteomics)
sample_cols <- all_cols[!(all_cols %in% c("Gene_symbol", "Protein_id"))]
proteomics <- proteomics[, c("Gene_symbol", "Protein_id", sample_cols)]

# Save preprocessed proteomics
output_file <- file.path(preprocessing_dir, "LemonPreprocessed_proteomics.txt")
write.table(proteomics, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("[OK] Saved:", output_file, "\n")

# Keep a version without gene ID columns for downstream use
RNA_preprocessed <- proteomics[, -(1:2)]
rownames(RNA_preprocessed) <- proteomics$Gene_symbol

# ============================================================================
# Process hPTM data (if provided)
# ============================================================================

if (!is.null(opt$hptm_file) && file.exists(opt$hptm_file)) {
  cat("\n=== Processing hPTM data ===\n")
  
  hPTM <- fread(opt$hptm_file, data.table=FALSE, header=TRUE)
  rownames(hPTM) <- hPTM[, 1]
  hPTM[, 1] <- NULL
  cat("hPTM dimensions:", nrow(hPTM), "x", ncol(hPTM), "\n")
  
  # Clean row names
  rownames(hPTM) <- str_replace_all(rownames(hPTM), ' ', '_')
  rownames(hPTM) <- str_replace_all(rownames(hPTM), '-', '_')
  
  # Fix sample ID formatting in column names — only strip date prefix if present
  if (any(grepl("^[0-9]+_", colnames(hPTM)))) {
    colnames(hPTM) <- sub("^[0-9]+_", "", colnames(hPTM))
  }
  
  # Match to proteomics samples (use intersection to handle missing samples)
  common_samples <- intersect(colnames(hPTM), colnames(RNA_preprocessed))
  cat("hPTM samples matched to proteomics:", length(common_samples), "\n")
  hPTM <- hPTM[, common_samples, drop=FALSE]
  RNA_hptm <- RNA_preprocessed[, common_samples, drop=FALSE]
  
  # Handle missing values and normalize
  hPTM[is.na(hPTM)] <- 0
  hPTM <- as.data.frame(t(scale(t(hPTM))))
  
  # Format for LemonTree
  hPTM_formatted <- hPTM
  hPTM_formatted$Gene_symbol <- rownames(hPTM_formatted)
  hPTM_formatted$Protein_id <- rownames(hPTM_formatted)
  hPTM_formatted <- hPTM_formatted[, c("Gene_symbol", "Protein_id", common_samples)]
  
  # Save
  output_file <- file.path(preprocessing_dir, "LemonPreprocessed_hPTM.txt")
  write.table(hPTM_formatted, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", output_file, "\n")

  # Save regulator name list (used as -reg_file by LemonTree regulators task)
  hptm_list_file <- file.path(preprocessing_dir, "hptms.txt")
  writeLines(rownames(hPTM), hptm_list_file)
  cat("[OK] Saved hPTM regulator list:", hptm_list_file, "\n")
}

# ============================================================================
# Process metabolomics data (if provided)
# ============================================================================

if (!is.null(opt$metabolomics_file) && file.exists(opt$metabolomics_file)) {
  cat("\n=== Processing metabolomics data ===\n")
  
  metabolomics <- fread(opt$metabolomics_file, data.table=FALSE, header=TRUE)
  rownames(metabolomics) <- metabolomics[, 1]
  metabolomics[, 1] <- NULL
  cat("Metabolomics dimensions:", nrow(metabolomics), "x", ncol(metabolomics), "\n")
  
  # Fix sample ID formatting — only strip date prefix if present
  if (any(grepl("^[0-9]+_", colnames(metabolomics)))) {
    colnames(metabolomics) <- sub("^[0-9]+_", "", colnames(metabolomics))
  }
  
  # Match to proteomics samples (use intersection)
  common_metabo_samples <- intersect(colnames(metabolomics), colnames(RNA_preprocessed))
  cat("Metabolomics samples matched to proteomics:", length(common_metabo_samples), "\n")
  metabolomics <- metabolomics[, common_metabo_samples, drop=FALSE]
  
  # Pareto scaling for metabolomics
  metabolomics <- as.data.frame(t(pareto_scale(t(metabolomics))))
  
  # Apply name mapping if provided
  if (!is.null(opt$metabolomics_labels) && file.exists(opt$metabolomics_labels)) {
    metabolomics_labels <- fread(opt$metabolomics_labels, data.table=FALSE, header=TRUE)
    metabolomics <- merge(metabolomics, metabolomics_labels, by.x="row.names", by.y=1, all.x=TRUE)
    metabolomics <- metabolomics[!duplicated(metabolomics$Row.names), ]
    rownames(metabolomics) <- metabolomics$Row.names
    metabolomics$Row.names <- NULL
    # Move ID column to second position
    id_col <- colnames(metabolomics)[ncol(metabolomics)]
    metabolomics <- metabolomics[, c(id_col, setdiff(colnames(metabolomics), id_col))]
  }
  
  # Format for LemonTree
  metabolomics_formatted <- metabolomics
  metabolomics_formatted$Gene_symbol <- rownames(metabolomics_formatted)
  metabolomics_formatted$Protein_id <- if (ncol(metabolomics_formatted) > 1) metabolomics_formatted[[1]] else rownames(metabolomics_formatted)
  
  # Ensure proper column order
  other_cols <- setdiff(colnames(metabolomics_formatted), c("Gene_symbol", "Protein_id"))
  metabolomics_formatted <- metabolomics_formatted[, c("Gene_symbol", "Protein_id", other_cols)]
  
  # Save
  output_file <- file.path(preprocessing_dir, "LemonPreprocessed_metabolomics.txt")
  write.table(metabolomics_formatted, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
  cat("[OK] Saved:", output_file, "\n")

  # Save regulator name list (used as -reg_file by LemonTree regulators task)
  metabolites_list_file <- file.path(preprocessing_dir, "metabolites.txt")
  writeLines(rownames(metabolomics), metabolites_list_file)
  cat("[OK] Saved metabolomics regulator list:", metabolites_list_file, "\n")
}

# ============================================================================
# Create combined LemonTree input file (all omics together)
# ============================================================================

cat("\n=== Creating combined LemonTree input ===\n")

# All sample columns from proteomics (the primary dataset)
all_sample_cols <- colnames(RNA_preprocessed)

# Helper: pad a formatted data frame to include all proteomics sample columns (NAs for missing)
pad_to_samples <- function(df_formatted, all_cols) {
  sample_cols_present <- intersect(all_cols, colnames(df_formatted))
  sample_cols_missing <- setdiff(all_cols, colnames(df_formatted))
  if (length(sample_cols_missing) > 0) {
    df_formatted[, sample_cols_missing] <- NA
  }
  df_formatted[, c("Gene_symbol", "Protein_id", all_cols)]
}

combined <- pad_to_samples(proteomics, all_sample_cols)

if (exists("hPTM_formatted")) {
  combined <- rbind(combined, pad_to_samples(hPTM_formatted, all_sample_cols))
}

if (exists("metabolomics_formatted")) {
  combined <- rbind(combined, pad_to_samples(metabolomics_formatted, all_sample_cols))
}

output_file <- file.path(preprocessing_dir, "LemonPreprocessed_complete.txt")
write.table(combined, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("[OK] Saved combined file:", output_file, "\n")

# Also save proteomics as LemonPreprocessed_expression.txt (expected by the pipeline)
file.copy(file.path(preprocessing_dir, "LemonPreprocessed_proteomics.txt"),
          file.path(preprocessing_dir, "LemonPreprocessed_expression.txt"), overwrite=TRUE)
cat("[OK] Saved:", file.path(preprocessing_dir, "LemonPreprocessed_expression.txt"), "\n")

# ============================================================================
# Save metadata for downstream use
# ============================================================================

output_file <- file.path(preprocessing_dir, "DESeq_groups.txt")
write.table(metadata, file=output_file, sep="\t", quote=FALSE)
cat("[OK] Saved metadata:", output_file, "\n")

cat("\n=== Preprocessing complete ===\n")
cat("Output directory:", preprocessing_dir, "\n")

# Create minimal name_mapping.tsv (expected by the pipeline)
name_map_file <- file.path(preprocessing_dir, "name_mapping.tsv")
write.table(data.frame(original=character(0), mapped=character(0)),
            file=name_map_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("[OK] Saved:", name_map_file, "\n")
