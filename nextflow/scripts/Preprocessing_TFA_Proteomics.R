#!/usr/bin/env Rscript

# ============================================================================
# Preprocessing_TFA_Proteomics.R
# ============================================================================
# Preprocessing script for pre-scaled proteomics/multi-omics data.
# Handles pre-scaled continuous data (no DESeq2 normalisation).
# Supports optional hPTM and metabolomics co-preprocessing.
# Optionally runs TF Activity inference (decoupleR / CollecTRI).
# Produces the same outputs as Preprocessing_TFA_RNA.R.

library(optparse)
library(data.table)
library(tidyverse)
library(ggplot2)
library(decoupleR)

# Command line options
option_list <- list(
  make_option(c("--expression"), type="character", default=NULL,
              help="Expression / proteomics file path (required)", metavar="character"),
  make_option(c("--metadata"), type="character", default=NULL,
              help="Metadata file path (required)", metavar="character"),
  make_option(c("--output_dir"), type="character", default=".",
              help="Output directory [default= %default]", metavar="character"),
  make_option(c("--top_n_genes"), type="integer", default=5000,
              help="Number of top variable features to select [default= %default]", metavar="integer"),
  make_option(c("--id_mapping"), type="character", default=NULL,
              help="Gene/protein ID mapping file (tab-separated: symbol <TAB> id)", metavar="character"),
  make_option(c("--hptm_file"), type="character", default=NULL,
              help="hPTM data file path (optional)", metavar="character"),
  make_option(c("--metabolomics_file"), type="character", default=NULL,
              help="Metabolomics data file path (optional)", metavar="character"),
  make_option(c("--metabolomics_labels"), type="character", default=NULL,
              help="Metabolomics name mapping file (optional)", metavar="character"),
  make_option(c("--sample_id_col"), type="character", default="ID",
              help="Sample ID column in metadata [default= %default]", metavar="character"),
  make_option(c("--perform_TFA"), type="logical", default=TRUE,
              help="Run TF Activity inference with decoupleR/CollecTRI [default= %default]", metavar="logical"),
  make_option(c("--prior_network"), type="character", default=NULL,
              help="Prior network file for TFA (source/target columns). Falls back to CollecTRI if NULL.", metavar="character"),
  make_option(c("--organism"), type="character", default="human",
              help="Organism for CollecTRI TFA: 'human' or 'mouse' [default= %default]", metavar="character"),
  make_option(c("--deseq_contrast1"), type="character", default=NULL,
              help="Metadata column to use for PCA colouring [default: first factor column]", metavar="character"),
  make_option(c("--regulator_types"), type="character", default=NULL,
              help="Comma-separated Prefix:DataFile[:DataType] pairs for additional regulator layers (continuous 'c', default, or discrete 'd'). Files are read from ./data/. Prefixes 'hPTMs'/'Metabolites' are skipped here (handled by --hptm_file/--metabolomics_file above).", metavar="character")
)

# Parse --regulator_types into { prefix -> (filename, data_type) }, mirroring
# Preprocessing_TFA_RNA.R's parse_regulator_config(). Unlike the RNA script, every prefix here
# is treated as a value matrix (feature x sample), never a bare gene-name list: this script's
# regulator layers (e.g. protein/kinase/TF activity, lipids) are abundance/activity scores, not
# names to look up in the primary matrix.
parse_regulator_config <- function(regulator_types_str) {
  configs <- list(); data_types <- list()
  if (is.null(regulator_types_str) || !nzchar(regulator_types_str)) return(list(configs=configs, data_types=data_types))
  for (pair in strsplit(regulator_types_str, ",")[[1]]) {
    pair <- trimws(pair)
    if (!nzchar(pair)) next
    parts <- strsplit(pair, ":")[[1]]
    if (length(parts) < 2) { warning(paste("Invalid regulator config format:", pair)); next }
    prefix <- trimws(parts[1]); filename <- trimws(parts[2])
    data_type <- if (length(parts) >= 3) tolower(trimws(parts[3])) else "c"
    if (!data_type %in% c("c", "d")) data_type <- "c"
    configs[[prefix]] <- filename; data_types[[prefix]] <- data_type
  }
  list(configs = configs, data_types = data_types)
}

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$expression)) { print_help(opt_parser); stop("--expression is required") }
if (is.null(opt$metadata))   { print_help(opt_parser); stop("--metadata is required") }

cat("=== Preprocessing Pre-scaled Proteomics/Multi-Omics Data ===\n")
cat("Expression file:", opt$expression, "\n")
cat("Metadata file:", opt$metadata, "\n")
cat("Output dir:", opt$output_dir, "\n")
cat("TFA enabled:", opt$perform_TFA, "\n")

perform_TFA <- opt$perform_TFA

# Create output directories
output_dir        <- opt$output_dir
preprocessing_dir <- file.path(output_dir, "LemonTree/Preprocessing")
tfa_dir           <- file.path(output_dir, "TFA")
dir.create(preprocessing_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tfa_dir,           recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# Helper functions
# ============================================================================

pareto_scale <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds   <- apply(x, 1, sd, na.rm = TRUE)
  (x - row_means) / sqrt(row_sds)
}

save_pheatmap_pdf <- function(x, filename, width = 10, height = 8) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

create_pca_plot <- function(data_matrix, metadata_df, contrast_col, omics_name,
                            out_dir = preprocessing_dir) {
  cat(sprintf("\n=== Creating PCA plot for %s ===\n", omics_name))
  data_matrix <- as.matrix(data_matrix)
  valid_rows <- complete.cases(data_matrix) & apply(data_matrix, 1, function(x) all(is.finite(x)))
  if (sum(!valid_rows) > 0) {
    cat(sprintf("  Removing %d features with NA/Inf values\n", sum(!valid_rows)))
    data_matrix <- data_matrix[valid_rows, , drop = FALSE]
  }
  if (nrow(data_matrix) < 2 || ncol(data_matrix) < 2) {
    cat("  [SKIP] Not enough features/samples for PCA\n"); return(invisible(NULL))
  }
  pca_result <- prcomp(t(data_matrix), center = TRUE, scale. = TRUE)
  var_exp    <- summary(pca_result)$importance[2, ] * 100
  pca_df     <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2],
                            Sample = rownames(pca_result$x))
  if (!is.null(contrast_col) && contrast_col %in% colnames(metadata_df)) {
    idx          <- match(pca_df$Sample, rownames(metadata_df))
    pca_df$Group <- metadata_df[[contrast_col]][idx]
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = paste("PCA -", omics_name),
           x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
           y = sprintf("PC2 (%.1f%% variance)", var_exp[2]),
           color = contrast_col) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  } else {
    p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(title = paste("PCA -", omics_name),
           x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
           y = sprintf("PC2 (%.1f%% variance)", var_exp[2])) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  }
  out_file <- file.path(out_dir, paste0("PCA_", gsub(" ", "_", tolower(omics_name)), ".pdf"))
  ggsave(out_file, plot = p, width = 8, height = 6)
  cat(sprintf("[OK] Saved PCA plot: %s\n", out_file))
  return(invisible(pca_result))
}

# ============================================================================
# Read metadata
# ============================================================================

cat("\nReading metadata...\n")
metadata <- fread(opt$metadata, data.table = FALSE)
if ("V1" %in% colnames(metadata)) metadata$V1 <- NULL

sample_col <- opt$sample_id_col
if (sample_col %in% colnames(metadata)) {
  rownames(metadata) <- metadata[[sample_col]]; metadata[[sample_col]] <- NULL
} else {
  rownames(metadata) <- metadata[[1]]; metadata[[1]] <- NULL
}
if (any(grepl("^[0-9]+_", rownames(metadata)))) {
  rownames(metadata) <- sub("^[0-9]+_", "", rownames(metadata))
}
metadata[] <- lapply(metadata, factor)
cat("Metadata samples:", nrow(metadata), "\n")

# Determine PCA contrast column
contrast_col <- opt$deseq_contrast1
if (is.null(contrast_col) || !contrast_col %in% colnames(metadata)) {
  factor_cols  <- names(which(sapply(metadata, is.factor)))
  contrast_col <- if (length(factor_cols) > 0) factor_cols[1] else NULL
  if (!is.null(contrast_col)) cat("PCA contrast column:", contrast_col, "\n")
}

# ============================================================================
# Read and preprocess proteomics data
# ============================================================================

cat("\nReading proteomics data...\n")
proteomics <- fread(opt$expression, header = TRUE, data.table = FALSE)
rownames(proteomics) <- proteomics[, 1]; proteomics[, 1] <- NULL
cat("Proteomics dimensions:", nrow(proteomics), "x", ncol(proteomics), "\n")

proteomics <- proteomics[, colnames(proteomics) %in% rownames(metadata), drop = FALSE]
ord        <- match(colnames(proteomics), rownames(metadata))
metadata   <- metadata[ord, , drop = FALSE]
cat("Matched samples:", ncol(proteomics), "\n")

# ============================================================================
# Select top N variable features
# ============================================================================

cat("\nSelecting top", opt$top_n_genes, "variable features...\n")
vars       <- apply(proteomics, 1, var, na.rm = TRUE)
n_genes    <- min(opt$top_n_genes, nrow(proteomics))
proteomics <- proteomics[names(sort(vars, decreasing = TRUE)[1:n_genes]), ]
cat("Selected", nrow(proteomics), "features\n")

# Numeric matrix for TFA / PCA (features × samples)
RNA_preprocessed <- proteomics
create_pca_plot(RNA_preprocessed, metadata, contrast_col, "Proteomics")

# ============================================================================
# Add protein IDs for LemonTree
# ============================================================================

cat("\nAdding ID mapping...\n")
prot_with_ids <- as.data.frame(proteomics)
prot_with_ids$Gene_symbol <- rownames(prot_with_ids)

if (!is.null(opt$id_mapping) && file.exists(opt$id_mapping)) {
  id_map <- fread(opt$id_mapping, header = TRUE, data.table = FALSE)
  id_map <- id_map[!duplicated(id_map[[1]]), ]
  prot_with_ids <- merge(prot_with_ids, id_map, by.x = "Gene_symbol", by.y = colnames(id_map)[1], all.x = TRUE)
  if (!"Protein_id" %in% colnames(prot_with_ids)) prot_with_ids$Protein_id <- prot_with_ids$Gene_symbol
} else {
  prot_with_ids$Protein_id <- prot_with_ids$Gene_symbol
}

sample_cols   <- setdiff(colnames(prot_with_ids), c("Gene_symbol", "Protein_id"))
prot_with_ids <- prot_with_ids[, c("Gene_symbol", "Protein_id", sample_cols)]

write.table(prot_with_ids, file.path(preprocessing_dir, "LemonPreprocessed_proteomics.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[OK] Saved: LemonPreprocessed_proteomics.txt\n")

# ============================================================================
# TF Activity Inference (decoupleR / CollecTRI)
# ============================================================================

if (perform_TFA) {
  cat("\n=== Running TF Activity Inference ===\n")

  net <- NULL
  if (!is.null(opt$prior_network) && file.exists(opt$prior_network)) {
    cat("Reading prior network from:", opt$prior_network, "\n")
    net <- fread(opt$prior_network)
    if (!all(c("source", "target") %in% colnames(net))) colnames(net)[1:2] <- c("source", "target")
    cat("Network loaded:", nrow(net), "interactions\n")
  } else {
    collectri_org <- if (tolower(opt$organism) %in% c("mouse", "mmu", "mus_musculus")) "mouse" else "human"
    cat("Loading CollecTRI (organism:", collectri_org, ")...\n")
    net <- tryCatch({
      r <- decoupleR::get_collectri(organism = collectri_org, split_complexes = FALSE)
      cat("CollecTRI downloaded:", nrow(r), "interactions\n"); r
    }, error = function(e1) {
      cat("[WARNING] Online download failed:", e1$message, "\n")
      pkn_file <- if (collectri_org == "mouse") "PKN/CollecTRI_network_mouse.txt" else "PKN/CollecTRI_network.txt"
      tryCatch({
        r <- fread(pkn_file)
        if (!all(c("source", "target") %in% colnames(r))) colnames(r)[1:2] <- c("source", "target")
        cat("Loaded PKN fallback:", nrow(r), "interactions\n"); r
      }, error = function(e2) {
        cat("[WARNING] Fallback also failed:", e2$message, "\n"); NULL
      })
    })
  }

  if (is.null(net)) {
    cat("[WARNING] No prior network available — skipping TFA\n")
    perform_TFA <- FALSE
  } else {
    mat <- as.matrix(RNA_preprocessed)
    storage.mode(mat) <- "double"

    # decouple()'s default statistics = c('mlm','ulm','wsum'). MLM fits ALL sources
    # simultaneously as covariates of one linear model per sample, so regulators whose
    # target sets are (near-)identical after intersecting with mat's features make that
    # design matrix rank-deficient ("N sources ... colinear ... cannot fit a linear
    # model"), aborting the whole decouple() call incl. ulm/wsum. decoupleR ships
    # check_corr() (see ?decoupleR::check_corr, "Check correlation (colinearity)")
    # specifically to find these pairs beforehand -- restrict to the features mlm will
    # actually see (rows of mat) and drop the redundant half of each near-perfectly
    # correlated pair so mlm can fit.
    net_for_corr <- net[net$target %in% rownames(mat), ]
    corr_threshold <- 0.99
    collinear_pairs <- tryCatch({
      decoupleR::check_corr(net_for_corr, .source = "source", .target = "target") %>%
        dplyr::filter(abs(correlation) >= corr_threshold)
    }, error = function(e) {
      cat("[WARNING] check_corr() failed, skipping collinearity pre-filter:", e$message, "\n")
      NULL
    })
    if (!is.null(collinear_pairs) && nrow(collinear_pairs) > 0) {
      sources_to_drop <- unique(collinear_pairs$source.2)
      cat(sprintf("[INFO] %d source(s) are collinear (|r|>=%.2f) with another regulator's target set -- dropping the redundant ones so MLM can fit: %s%s\n",
                  length(sources_to_drop), corr_threshold,
                  paste(head(sources_to_drop, 20), collapse = ", "),
                  if (length(sources_to_drop) > 20) ", ..." else ""))
      net <- net[!(net$source %in% sources_to_drop), ]
    }

    tfa_failure_reason <- NULL
    decouple_ok <- tryCatch({
      cat("Running decouple analysis...\n")
      decoupled <- decouple(mat = mat, net = net, .source = "source", .target = "target")
      consensus <- run_consensus(decoupled)
      cat("[OK] Consensus TFA completed\n")

      TFA_df <- consensus %>%
        dplyr::filter(statistic == "consensus") %>%
        pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        column_to_rownames("condition") %>%
        t() %>% as.data.frame()

      save(TFA_df, file = file.path(tfa_dir, "TFA_df.RData"))
      write.table(TFA_df, file.path(preprocessing_dir, "TFA_consensus.txt"),
                  sep = "\t", quote = FALSE, row.names = TRUE)
      cat("[OK] Saved: TFA_consensus.txt\n")

      # Heatmap: top 50 most variable TFs
      n_tfs   <- 50
      top_tfs <- consensus %>% group_by(source) %>%
        summarise(std = sd(score)) %>% arrange(-abs(std)) %>% head(n_tfs) %>% pull(source)
      act_mat <- consensus %>%
        dplyr::filter(statistic == "consensus") %>%
        pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
        column_to_rownames("condition") %>% as.matrix()
      act_mat <- scale(act_mat[, top_tfs, drop = FALSE])

      palette_length <- 100
      my_color  <- colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
      my_breaks <- c(seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
                     seq(0.05, 3, length.out = floor(palette_length / 2)))

      viz_cols <- metadata[rownames(metadata) %in% rownames(act_mat), , drop = FALSE]
      viz_cols <- viz_cols[, sapply(viz_cols, function(x) is.factor(x) && nlevels(x) <= 10), drop = FALSE]

      suppressPackageStartupMessages(library(pheatmap))
      plot_tfa <- pheatmap(t(act_mat), border_color = NA, color = my_color, breaks = my_breaks,
                           annotation_col = if (ncol(viz_cols) > 0) viz_cols else NULL,
                           cluster_rows = TRUE, show_colnames = FALSE, fontsize_row = 5,
                           main = paste("TF Activity (top", n_tfs, "variable TFs)"))
      save_pheatmap_pdf(plot_tfa, file.path(tfa_dir,
                        paste0("TFA_consensus_", n_tfs, "heatmap_annotated.pdf")))
      cat("[OK] Saved TFA heatmap\n")
      TRUE
    }, error = function(e) {
      # Collinearity screening above catches near-duplicate regulon pairs, but MLM's rank
      # deficiency can also come from >2-way (multivariate) collinearity that a pairwise
      # correlation check won't catch. Rather than lose TFA entirely, fall back to ULM
      # (per-source univariate regression -- immune to collinearity between sources since
      # each source is fit independently) + WSUM, dropping MLM from the run.
      cat("[WARNING] decouple() with MLM failed:", e$message, "\n")
      cat("Retrying without MLM (ULM + WSUM only, both robust to collinear regulons)...\n")
      tryCatch({
        decoupled <<- decouple(mat = mat, net = net, .source = "source", .target = "target",
                                statistics = c("ulm", "wsum"), consensus_stats = c("ulm", "norm_wsum"))
        consensus <<- run_consensus(decoupled)
        cat("[OK] Consensus TFA completed without MLM\n")

        TFA_df <- consensus %>%
          dplyr::filter(statistic == "consensus") %>%
          pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
          column_to_rownames("condition") %>%
          t() %>% as.data.frame()

        save(TFA_df, file = file.path(tfa_dir, "TFA_df.RData"))
        write.table(TFA_df, file.path(preprocessing_dir, "TFA_consensus.txt"),
                    sep = "\t", quote = FALSE, row.names = TRUE)
        cat("[OK] Saved: TFA_consensus.txt\n")

        n_tfs   <- 50
        top_tfs <- consensus %>% group_by(source) %>%
          summarise(std = sd(score)) %>% arrange(-abs(std)) %>% head(n_tfs) %>% pull(source)
        act_mat <- consensus %>%
          dplyr::filter(statistic == "consensus") %>%
          pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
          column_to_rownames("condition") %>% as.matrix()
        act_mat <- scale(act_mat[, top_tfs, drop = FALSE])

        palette_length <- 100
        my_color  <- colorRampPalette(c("Darkblue", "white", "red"))(palette_length)
        my_breaks <- c(seq(-3, 0, length.out = ceiling(palette_length / 2) + 1),
                       seq(0.05, 3, length.out = floor(palette_length / 2)))

        viz_cols <- metadata[rownames(metadata) %in% rownames(act_mat), , drop = FALSE]
        viz_cols <- viz_cols[, sapply(viz_cols, function(x) is.factor(x) && nlevels(x) <= 10), drop = FALSE]

        suppressPackageStartupMessages(library(pheatmap))
        plot_tfa <- pheatmap(t(act_mat), border_color = NA, color = my_color, breaks = my_breaks,
                             annotation_col = if (ncol(viz_cols) > 0) viz_cols else NULL,
                             cluster_rows = TRUE, show_colnames = FALSE, fontsize_row = 5,
                             main = paste("TF Activity (top", n_tfs, "variable TFs, MLM skipped)"))
        save_pheatmap_pdf(plot_tfa, file.path(tfa_dir,
                          paste0("TFA_consensus_", n_tfs, "heatmap_annotated.pdf")))
        cat("[OK] Saved TFA heatmap\n")
        TRUE
      }, error = function(e2) {
        cat("[ERROR] TFA failed even without MLM:", e2$message, "\n[WARNING] Continuing without TFA\n")
        tfa_failure_reason <<- e2$message
        FALSE
      })
    })
    if (!decouple_ok) {
      perform_TFA <- FALSE
      if (is.null(tfa_failure_reason)) tfa_failure_reason <- "TFA analysis failed (see log above for details)"
    }
  }
}

tfa_status <- if (!opt$perform_TFA) {
  "DISABLED: perform_TFA=FALSE was requested"
} else if (perform_TFA) {
  "SUCCESS"
} else {
  paste0("FAILED: ", if (exists("tfa_failure_reason") && !is.null(tfa_failure_reason)) tfa_failure_reason else "no prior network available")
}
writeLines(tfa_status, file.path(tfa_dir, "TFA_status.txt"))
cat("[INFO] TFA status:", tfa_status, "\n")

# ============================================================================
# Process hPTM data (if provided)
# ============================================================================

if (!is.null(opt$hptm_file) && file.exists(opt$hptm_file)) {
  cat("\n=== Processing hPTM data ===\n")
  hPTM <- fread(opt$hptm_file, data.table = FALSE, header = TRUE)
  rownames(hPTM) <- str_replace_all(str_replace_all(hPTM[, 1], " ", "_"), "-", "_")
  hPTM[, 1] <- NULL
  cat("hPTM dimensions:", nrow(hPTM), "x", ncol(hPTM), "\n")

  if (any(grepl("^[0-9]+_", colnames(hPTM)))) colnames(hPTM) <- sub("^[0-9]+_", "", colnames(hPTM))
  common_samples <- intersect(colnames(hPTM), colnames(RNA_preprocessed))
  cat("hPTM samples matched:", length(common_samples), "\n")
  hPTM <- hPTM[, common_samples, drop = FALSE]
  hPTM[is.na(hPTM)] <- 0
  hPTM <- as.data.frame(t(scale(t(hPTM))))

  hPTM_formatted <- hPTM
  hPTM_formatted$Gene_symbol <- rownames(hPTM_formatted)
  hPTM_formatted$Protein_id  <- rownames(hPTM_formatted)
  hPTM_formatted <- hPTM_formatted[, c("Gene_symbol", "Protein_id", common_samples)]

  write.table(hPTM_formatted, file.path(preprocessing_dir, "LemonPreprocessed_hPTM.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(rownames(hPTM), file.path(preprocessing_dir, "hptms.txt"))
  cat("[OK] Saved LemonPreprocessed_hPTM.txt and hptms.txt\n")
  create_pca_plot(hPTM, metadata, contrast_col, "hPTM")
}

# ============================================================================
# Process metabolomics data (if provided)
# ============================================================================

if (!is.null(opt$metabolomics_file) && file.exists(opt$metabolomics_file)) {
  cat("\n=== Processing metabolomics data ===\n")
  metabolomics <- fread(opt$metabolomics_file, data.table = FALSE, header = TRUE)
  # LemonTree splits data_file rows on whitespace, so a name like "2-hydroxybutyric acid" would
  # shift columns and crash it (ArrayIndexOutOfBoundsException) -- sanitize like hPTM does above.
  rownames(metabolomics) <- str_replace_all(as.character(metabolomics[, 1]), " ", "_")
  metabolomics[, 1] <- NULL
  cat("Metabolomics dimensions:", nrow(metabolomics), "x", ncol(metabolomics), "\n")

  if (any(grepl("^[0-9]+_", colnames(metabolomics)))) {
    colnames(metabolomics) <- sub("^[0-9]+_", "", colnames(metabolomics))
  }
  common_metabo_samples <- intersect(colnames(metabolomics), colnames(RNA_preprocessed))
  cat("Metabolomics samples matched:", length(common_metabo_samples), "\n")
  metabolomics <- metabolomics[, common_metabo_samples, drop = FALSE]
  # pareto_scale() is row-wise (per feature), same as Preprocessing_TFA_RNA.R's version --
  # call it directly, do NOT wrap in t(...(t(x))) (that idiom is only correct for column-wise
  # functions like scale(), and was previously copied from the neighboring hPTM line, which
  # inverted the axis: it centered/scaled per SAMPLE across metabolites instead of per
  # METABOLITE across samples, flattening real cross-sample variation in heatmaps).
  metabolomics <- as.data.frame(pareto_scale(metabolomics))

  if (!is.null(opt$metabolomics_labels) && file.exists(opt$metabolomics_labels)) {
    met_labels   <- fread(opt$metabolomics_labels, data.table = FALSE, header = TRUE)
    metabolomics <- merge(metabolomics, met_labels, by.x = "row.names", by.y = 1, all.x = TRUE)
    metabolomics <- metabolomics[!duplicated(metabolomics$Row.names), ]
    rownames(metabolomics) <- metabolomics$Row.names; metabolomics$Row.names <- NULL
    id_col       <- colnames(metabolomics)[ncol(metabolomics)]
    metabolomics <- metabolomics[, c(id_col, setdiff(colnames(metabolomics), id_col))]
  }

  metabolomics_formatted <- metabolomics
  metabolomics_formatted$Gene_symbol <- rownames(metabolomics_formatted)
  metabolomics_formatted$Protein_id  <- if (ncol(metabolomics_formatted) > 1) metabolomics_formatted[[1]] else rownames(metabolomics_formatted)
  other_cols <- setdiff(colnames(metabolomics_formatted), c("Gene_symbol", "Protein_id"))
  metabolomics_formatted <- metabolomics_formatted[, c("Gene_symbol", "Protein_id", other_cols)]

  write.table(metabolomics_formatted, file.path(preprocessing_dir, "LemonPreprocessed_metabolomics.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(rownames(metabolomics), file.path(preprocessing_dir, "metabolites.txt"))
  cat("[OK] Saved LemonPreprocessed_metabolomics.txt and metabolites.txt\n")
  create_pca_plot(metabolomics[, common_metabo_samples, drop = FALSE], metadata, contrast_col, "Metabolomics")
}

# ============================================================================
# Generic additional regulator layers (--regulator_types), beyond hPTM/Metabolites
# ============================================================================
# Preprocessing_TFA_RNA.R has a fully generic Prefix:File[:type] regulator loader; this
# proteomics script only had the two fixed hooks above (hptm_file/metabolomics_file). This
# block adds the same generic capability so any number of extra continuous/discrete regulator
# layers (e.g. protein/kinase/TF activity, lipids) can be supplied via --regulator_types with
# no further code changes, matching the RNA path. Every prefix here is a value matrix (feature
# x sample) -- unlike the RNA script's TF-list special case, there is no bare-name-list mode,
# since these layers are activity/abundance scores rather than lookups into the primary matrix.

reg_parsed           <- parse_regulator_config(opt$regulator_types)
regulator_configs    <- reg_parsed$configs
regulator_data_types <- reg_parsed$data_types
generic_layers  <- list()
seen_names <- rownames(RNA_preprocessed)
if (exists("hPTM_formatted"))         seen_names <- c(seen_names, rownames(hPTM))
if (exists("metabolomics_formatted")) seen_names <- c(seen_names, rownames(metabolomics))

for (prefix in names(regulator_configs)) {
  if (tolower(prefix) %in% c("hptms", "hptm", "metabolites")) next  # handled by the hooks above
  filename  <- regulator_configs[[prefix]]
  data_type <- regulator_data_types[[prefix]]
  filepath  <- file.path("data", filename)
  if (!file.exists(filepath)) {
    cat(sprintf("[WARNING] Regulator file not found for %s: %s -- skipping\n", prefix, filepath))
    next
  }
  cat(sprintf("\n=== Processing %s regulator layer from: %s (data_type=%s) ===\n", prefix, filepath, data_type))
  layer <- fread(filepath, data.table = FALSE, header = TRUE)
  ids <- str_replace_all(str_replace_all(as.character(layer[, 1]), " ", "_"), "-", "_")
  layer[, 1] <- NULL
  n_dup_within <- sum(duplicated(ids))
  if (n_dup_within > 0) {
    cat(sprintf("  [INFO] %d duplicate feature names within %s -- keeping first occurrence of each\n",
                n_dup_within, prefix))
    keep <- !duplicated(ids)
    layer <- layer[keep, , drop = FALSE]
    ids   <- ids[keep]
  }
  rownames(layer) <- ids
  cat(sprintf("  %s dimensions: %d x %d\n", prefix, nrow(layer), ncol(layer)))

  if (any(grepl("^[0-9]+_", colnames(layer)))) colnames(layer) <- sub("^[0-9]+_", "", colnames(layer))
  common <- intersect(colnames(layer), colnames(RNA_preprocessed))
  cat(sprintf("  %s samples matched: %d\n", prefix, length(common)))
  layer <- layer[, common, drop = FALSE]
  layer[] <- lapply(layer, function(col) suppressWarnings(as.numeric(as.character(col))))
  layer[is.na(layer)] <- 0

  # LemonTree indexes rows by name; a name colliding with an earlier layer/primary corrupts its
  # array indexing (ArrayIndexOutOfBoundsException) -- disambiguate rather than allow that.
  dup <- rownames(layer) %in% seen_names
  if (any(dup)) {
    suffix <- paste0("_", tolower(gsub("[^A-Za-z0-9]", "", prefix)))
    rownames(layer)[dup] <- paste0(rownames(layer)[dup], suffix)
    cat(sprintf("  [INFO] %d %s feature names collided with an earlier layer -- suffixed with '%s'\n",
                sum(dup), prefix, suffix))
  }
  seen_names <- c(seen_names, rownames(layer))

  if (data_type == "d") {
    cat(sprintf("  %s: discrete/binary -- used as-is (no scaling)\n", prefix))
  } else {
    layer <- as.data.frame(t(scale(t(layer))))
    cat(sprintf("  %s: continuous -- z-scored per feature\n", prefix))
  }

  layer_formatted <- layer
  layer_formatted$Gene_symbol <- rownames(layer_formatted)
  layer_formatted$Protein_id  <- rownames(layer_formatted)
  layer_formatted <- layer_formatted[, c("Gene_symbol", "Protein_id", common)]

  data_file_basename <- tolower(tools::file_path_sans_ext(basename(filename)))
  write.table(layer_formatted, file.path(preprocessing_dir, paste0("LemonPreprocessed_", data_file_basename, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("[OK] Saved: LemonPreprocessed_%s.txt\n", data_file_basename))

  list_filename <- file.path(preprocessing_dir, paste0(tolower(prefix), ".txt"))
  if (data_type == "d") {
    reg_df <- data.frame(name = rownames(layer), type = "d")
    write.table(reg_df, list_filename, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
  } else {
    writeLines(rownames(layer), list_filename)
  }
  cat(sprintf("[OK] Saved: %s regulator list (%d features) -> %s\n", tolower(prefix), nrow(layer), list_filename))

  generic_layers[[prefix]] <- layer_formatted
}

# ============================================================================
# Create combined LemonTree input file
# ============================================================================

cat("\n=== Creating combined LemonTree input ===\n")

all_sample_cols <- colnames(RNA_preprocessed)

pad_to_samples <- function(df_formatted, all_cols) {
  missing <- setdiff(all_cols, colnames(df_formatted))
  if (length(missing) > 0) df_formatted[, missing] <- NA
  df_formatted[, c("Gene_symbol", "Protein_id", all_cols)]
}

combined <- pad_to_samples(prot_with_ids, all_sample_cols)
if (exists("hPTM_formatted"))         combined <- rbind(combined, pad_to_samples(hPTM_formatted, all_sample_cols))
if (exists("metabolomics_formatted")) combined <- rbind(combined, pad_to_samples(metabolomics_formatted, all_sample_cols))
for (prefix in names(generic_layers)) combined <- rbind(combined, pad_to_samples(generic_layers[[prefix]], all_sample_cols))

write.table(combined, file.path(preprocessing_dir, "LemonPreprocessed_complete.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("[OK] Saved: LemonPreprocessed_complete.txt\n")

# LemonPreprocessed_expression.txt = proteomics (primary data), used by CLUSTERING
file.copy(file.path(preprocessing_dir, "LemonPreprocessed_proteomics.txt"),
          file.path(preprocessing_dir, "LemonPreprocessed_expression.txt"), overwrite = TRUE)
cat("[OK] Saved: LemonPreprocessed_expression.txt\n")

# ============================================================================
# Metadata and name mapping
# ============================================================================

write.table(metadata, file.path(preprocessing_dir, "DESeq_groups.txt"), sep = "\t", quote = FALSE)
cat("[OK] Saved: DESeq_groups.txt\n")

write.table(data.frame(original = character(0), mapped = character(0)),
            file.path(preprocessing_dir, "name_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("[OK] Saved: name_mapping.tsv\n")

# ============================================================================
# Summary
# ============================================================================

cat("\n=== Preprocessing complete ===\n")
cat("Output directory:", preprocessing_dir, "\n")
cat("Files produced:\n")
for (f in list.files(preprocessing_dir)) cat(sprintf("  - LemonTree/Preprocessing/%s\n", f))
if (perform_TFA) {
  cat("TFA files:\n")
  for (f in list.files(tfa_dir)) cat(sprintf("  - TFA/%s\n", f))
} else {
  cat("TFA: skipped (perform_TFA=FALSE or network unavailable)\n")
}
