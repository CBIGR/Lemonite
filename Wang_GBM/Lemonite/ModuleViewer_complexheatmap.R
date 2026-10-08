#!/usr/bin/env Rscript
# ModuleViewer (ComplexHeatmap version) for the Wang GBM case study.
#
# Rebuild of nextflow/scripts/module_viewer.R (the generic, CLI-driven pipeline
# version) as a standalone script for this case study: same ComplexHeatmap
# rendering approach (independent color scale per omics block, packed legends,
# capped per-block row heights) but with the analysis's own hardcoded file
# layout instead of --regulator_files/--regulator_types CLI parsing, matching
# the conventions of ModuleViewer.ipynb and ChEA_TF_target_enrichment.R in
# this folder.
#
# Adds three right-annotation KG blocks versus the pipeline version, drawn
# alongside the existing metabolite-gene PKN interaction panel:
#   - TF-target ORA          (TF_target_ORA_interactions.mvf, ChEA_TF_target_enrichment.R)
#   - TF motif, predicted    (HOMER_motif_interactions.mvf, HOMER_motif_enrichment.R)
#   - TF motif, unpredicted  (HOMER_motif_discovery_interactions.mvf, same script)
# The two HOMER panels are kept separate on purpose: the first says a predicted
# regulator's motif is enriched in the module's promoters, the second says some
# other TF's motif is - evidence about the module, not about the prediction.
#
# One more block marks genes linked to an ENHANCER that carries an IRF-TF motif
# (HOMER_enhancer_IRF_interactions.mvf, made by make_IRF_enhancer_mvf.R from the extended-library
# HOMER_enhancer_enrichment.R results). Links are FANTOM5 expression-correlation enhancer-TSS
# associations and a "hit" is a motif match, so this is suggestive, not evidence of binding.
#
# Usage: Rscript ModuleViewer_complexheatmap.R [module,module,...]   (default: all modules)

suppressMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
ht_opt$message <- FALSE
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ------------------------------------------------------------------------------------------------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------------------------------------------------------------------------------------------------

base_dir   <- '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/'
viewer_dir <- file.path(base_dir, 'ModuleViewer_files')
output_dir <- file.path(base_dir, 'Figures_complexheatmap')
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expression_file <- file.path(base_dir, 'Preprocessing', 'LemonPreprocessed_expression.txt')
complete_file   <- file.path(base_dir, 'Preprocessing', 'LemonPreprocessed_complete.txt')
name_mapping_file <- file.path(base_dir, 'Preprocessing', 'name_mapping.tsv')
specific_modules_file <- file.path(base_dir, 'Networks', 'specific_modules.txt')

percentile <- 2      # matches the *.percentileN_list.txt / *.percentileN.txt naming used throughout this analysis
res <- 300           # PNG/PDF resolution (dpi)
SHOW_REGULATOR_SCORES <- TRUE

# Regulator blocks stacked above the Expression heatmap, in this order (mirrors
# ModuleViewer.ipynb's REGULATORS list; hPTM/Proteins are commented out there too).
# is_tf = TRUE draws straight from the primary expression matrix (TF symbols are
# gene symbols, already present there); is_tf = FALSE looks up a dedicated
# LemonPreprocessed_<omics_basename>.txt file for that layer's own values.
regulator_specs <- list(
  list(reg_type = "TFs",         list_file = "Lovering.percentile2_list.txt",  prefix = "Lovering",   is_tf = TRUE),
  list(reg_type = "Metabolites", list_file = "Metabolite.percentile2_list.txt", prefix = "Metabolite", is_tf = FALSE, omics_basename = "metabolomics"),
  list(reg_type = "Lipids",      list_file = "Lipid.percentile2_list.txt",      prefix = "Lipid",      is_tf = FALSE, omics_basename = "lipidomics")
)

# Sample annotation track to draw beneath the Expression heatmap (::TYPE= in sample_mapping.mvf)
annotation_type <- "Diagnosis"

CORUM_MAX_COLS <- 10

## ============================== Parsing helpers ============================

parse_sample_mapping <- function(path) {
  content <- readLines(path)
  txt <- paste(content, collapse = "\n")
  sections <- strsplit(txt, "\n---\n")[[1]]
  out <- list()
  for (sec in sections) {
    lines <- strsplit(sec, "\n")[[1]]
    type_line <- grep("^::TYPE=", lines, value = TRUE)
    legend_line <- grep("^::LEGEND=", lines, value = TRUE)
    map_line <- grep("^\\|", lines, value = TRUE)
    if (length(type_line) == 0 || length(legend_line) == 0 || length(map_line) == 0) next
    type_name <- sub("^::TYPE=", "", type_line)
    legend_part <- sub("^::LEGEND=", "", legend_line)
    parts <- strsplit(legend_part, "\t")[[1]]
    legend_items <- strsplit(parts[1], "\\|")[[1]]
    label <- if (length(parts) > 1) parts[2] else type_name
    legend <- list()
    for (item in legend_items) {
      kv <- strsplit(item, ":")[[1]]
      if (length(kv) == 2) legend[[trimws(kv[1])]] <- tolower(trimws(kv[2]))
    }
    sample_map <- list()
    map_items <- strsplit(gsub("^\\||\\|$", "", map_line[1]), "\\|")[[1]]
    for (item in map_items) {
      kv <- strsplit(item, ":")[[1]]
      if (length(kv) == 2) sample_map[[trimws(kv[1])]] <- tolower(trimws(kv[2]))
    }
    out[[type_name]] <- list(label = label, legend = legend, samples = sample_map)
  }
  out
}

parse_list_file <- function(path) {
  lines <- readLines(path)
  lines <- lines[nzchar(lines)]
  mod <- sub("\t.*$", "", lines)
  vals <- sub("^[^\t]*\t", "", lines)
  setNames(strsplit(vals, "\\|"), mod)
}

# module<TAB>gene1|gene2|...<TAB>value  (value optional -- PPI/HumanNet only have 2 cols)
parse_kg_mvf <- function(path) {
  empty <- data.frame(Module = character(), Genes = character(), Value = character())
  if (!file.exists(path)) return(empty)
  lines <- readLines(path)
  lines <- lines[!grepl("^::", lines)]
  lines <- lines[nzchar(lines)]
  if (length(lines) == 0) return(empty)
  rows <- strsplit(lines, "\t")
  data.frame(
    Module = sapply(rows, `[`, 1),
    Genes  = sapply(rows, `[`, 2),
    Value  = sapply(rows, function(r) if (length(r) >= 3) r[3] else NA),
    stringsAsFactors = FALSE
  )
}

parse_mvf_metadata <- function(path) {
  meta <- list()
  if (!file.exists(path)) return(meta)
  for (line in readLines(path)) {
    if (startsWith(line, "::")) {
      kv <- strsplit(sub("^::", "", line), "=|:", perl = TRUE)[[1]]
      if (length(kv) >= 2) meta[[trimws(kv[1])]] <- trimws(kv[2])
    }
  }
  meta
}

web_to_r_color <- function(x) {
  map <- c(teal = "darkcyan", lime = "green2", olive = "darkolivegreen",
           navy = "navy", maroon = "maroon", coral = "coral")
  out <- ifelse(x %in% names(map), map[x], x)
  names(out) <- names(x)
  out
}

truncate_label <- function(x, max_chars = 40) {
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

rank_by_shared_genes <- function(module_data, max_cols = NULL) {
  values <- unique(module_data$Value)
  if (is.null(max_cols) || length(values) <= max_cols) return(values)
  n_shared <- sapply(values, function(v) {
    gl <- module_data$Genes[module_data$Value == v]
    max(sapply(gl, function(g) length(strsplit(g, "\\|")[[1]])), 0)
  })
  values[order(n_shared, decreasing = TRUE)][seq_len(max_cols)]
}

read_mat <- function(path) {
  d <- read.delim(path, check.names = FALSE)
  # ID column name varies by preprocessing path: RNA-seq writes "symbol",
  # proteomics/phospho/kinase/lipid variants write "Gene_symbol".
  id_col <- intersect(c("symbol", "Gene_symbol", "gene_symbol"), colnames(d))
  if (length(id_col) > 0) {
    ids <- d[[id_col[1]]]
    dup <- duplicated(ids)
    if (any(dup)) {
      cat(sprintf("Warning: %d duplicate row identifier(s) in %s, keeping first occurrence (e.g. %s)\n",
                  sum(dup), path, paste(unique(ids[dup]), collapse = ", ")))
      d <- d[!dup, , drop = FALSE]
      ids <- ids[!dup]
    }
    rownames(d) <- ids
    d[[id_col[1]]] <- NULL
  }
  # Drop any stray non-numeric/all-NA columns (e.g. gene_id/ensembl_gene_id)
  junk <- vapply(d, function(col) !is.numeric(col) || all(is.na(col)), logical(1))
  as.matrix(d[, !junk, drop = FALSE])
}

## ============================== Load shared inputs ==========================

name_lookup <- c()
if (file.exists(name_mapping_file)) {
  nm <- read.delim(name_mapping_file)
  if (all(c("cleaned", "original") %in% colnames(nm)) && nrow(nm) > 0) {
    name_lookup <- setNames(nm$original, nm$cleaned)
    cat(sprintf("Loaded %d name mappings from %s\n", length(name_lookup), name_mapping_file))
  }
}
restore_names <- function(x) {
  if (length(name_lookup) == 0) return(x)
  ifelse(x %in% names(name_lookup), unname(name_lookup[x]), x)
}

cat(sprintf("Loading expression data from: %s\n", expression_file))
expr_all <- read_mat(expression_file)

clusters <- parse_list_file(file.path(viewer_dir, "clusters_list.txt"))
sample_mapping <- parse_sample_mapping(file.path(viewer_dir, "sample_mapping.mvf"))

## specific_modules.txt (written by LemonTree_to_network after post-clustering
## consolidation) restricts which modules actually get rendered.
if (file.exists(specific_modules_file)) {
  specific_modules <- trimws(readLines(specific_modules_file))
  specific_modules <- specific_modules[nzchar(specific_modules)]
  before <- length(clusters)
  clusters <- clusters[names(clusters) %in% specific_modules]
  cat(sprintf("Modules before filtering: %d\nModules after filtering: %d\n", before, length(clusters)))
} else {
  cat("Warning: specific_modules.txt not found - processing all modules\n")
}

## ---- regulator blocks (config-driven from regulator_specs above) ----
regulators <- list()
for (spec in regulator_specs) {
  regs_path <- file.path(viewer_dir, spec$list_file)
  if (!file.exists(regs_path)) {
    cat(sprintf("Warning: regulator file not found: %s\n", regs_path))
    next
  }
  regs_list <- parse_list_file(regs_path)

  score_path <- file.path(base_dir, "Lemon_out", sprintf("%s.percentile%d.txt", spec$prefix, percentile))
  score_df <- if (file.exists(score_path)) read.delim(score_path) else NULL
  if (SHOW_REGULATOR_SCORES) {
    if (!is.null(score_df)) cat(sprintf("Loaded scores for %s from %s\n", spec$reg_type, score_path))
    else cat(sprintf("Warning: score file not found for %s: %s\n", spec$reg_type, score_path))
  }

  if (spec$is_tf) {
    omics_data <- expr_all
    cat(sprintf("Using expression data for %s regulators\n", spec$reg_type))
  } else {
    omics_path <- file.path(base_dir, "Preprocessing", sprintf("LemonPreprocessed_%s.txt", spec$omics_basename))
    if (file.exists(omics_path)) {
      omics_data <- read_mat(omics_path)
      cat(sprintf("Loaded omics-specific data for %s from %s\n", spec$reg_type, omics_path))
    } else {
      cat(sprintf("Warning: omics file not found for %s (%s), falling back to complete data\n", spec$reg_type, omics_path))
      omics_data <- read_mat(complete_file)
    }
  }

  regulators[[length(regulators) + 1]] <- list(
    name = spec$reg_type, is_tf = spec$is_tf, regs_list = regs_list,
    omics_data = omics_data, score_df = score_df
  )
}
cat(sprintf("Configured %d regulator types: %s\n", length(regulators), paste(sapply(regulators, `[[`, "name"), collapse = ", ")))

## ---- KG / interaction sources (all optional, guarded by file existence) ----
metabo_kg   <- parse_kg_mvf(file.path(viewer_dir, "metabolite_LemonIteKG_interactions.mvf"))
metabo_meta <- parse_mvf_metadata(file.path(viewer_dir, "metabolite_LemonIteKG_interactions.mvf"))
tfora_kg    <- parse_kg_mvf(file.path(viewer_dir, "TF_target_ORA_interactions.mvf"))
tfora_meta  <- parse_mvf_metadata(file.path(viewer_dir, "TF_target_ORA_interactions.mvf"))
homer_kg    <- parse_kg_mvf(file.path(viewer_dir, "HOMER_motif_interactions.mvf"))
homer_meta  <- parse_mvf_metadata(file.path(viewer_dir, "HOMER_motif_interactions.mvf"))
# HOMER_motif_discovery (promoter, unpredicted motifs) is intentionally not read here: it's
# real data (kept in the CSV/mvf for transparency) but not TF-validation evidence, so it does
# not belong on the same per-gene panel as the ORA/HOMER-predicted validation columns.
enhancer_kg   <- parse_kg_mvf(file.path(viewer_dir, "HOMER_enhancer_interactions.mvf"))
enhancer_meta <- parse_mvf_metadata(file.path(viewer_dir, "HOMER_enhancer_interactions.mvf"))
irf_enh_kg   <- parse_kg_mvf(file.path(viewer_dir, "HOMER_enhancer_IRF_interactions.mvf"))
irf_enh_meta <- parse_mvf_metadata(file.path(viewer_dir, "HOMER_enhancer_IRF_interactions.mvf"))
ppi_ints    <- parse_kg_mvf(file.path(viewer_dir, "PPI_interactions.mvf"))
ppi_meta    <- parse_mvf_metadata(file.path(viewer_dir, "PPI_interactions.mvf"))
hn_ints     <- parse_kg_mvf(file.path(viewer_dir, "HumanNet_interactions.mvf"))
hn_meta     <- parse_mvf_metadata(file.path(viewer_dir, "HumanNet_interactions.mvf"))

## ---- sample annotation (single track: Diagnosis) ----
if (!annotation_type %in% names(sample_mapping)) {
  cat(sprintf("Warning: annotation type '%s' not found in sample_mapping.mvf. Available: %s\n",
              annotation_type, paste(names(sample_mapping), collapse = ", ")))
  annotation_type <- names(sample_mapping)[1]
}
selected_metadata <- sample_mapping[annotation_type]
cat(sprintf("Displaying annotation track: %s\n", annotation_type))

cat(sprintf("%d modules to process, PPI rows=%d, HumanNet rows=%d, Metabolite-KG rows=%d, TF-target ORA rows=%d, HOMER promoter-motif rows=%d, HOMER enhancer-motif rows=%d\n",
            length(clusters), nrow(ppi_ints), nrow(hn_ints), nrow(metabo_kg), nrow(tfora_kg),
            nrow(homer_kg), nrow(enhancer_kg)))

## ============================== Sizing helpers ==============================

show_names <- function(n) n <= 150
fontsize_for <- function(n) max(7, min(9, 200 / n))
row_h_cm <- function(n) {
  max_total_cm <- 65
  if (show_names(n)) {
    fs <- fontsize_for(n)
    target <- (fs * 1.5) / 28.35
  } else {
    target <- 0.05
  }
  if (n * target > max_total_cm) target <- max_total_cm / n
  max(target, 0.01)
}

robust_col_fun <- function(mat, fixed_max = NULL) {
  if (is.null(fixed_max)) fixed_max <- max(0.5, quantile(abs(mat), 0.99, na.rm = TRUE))
  colorRamp2(c(-fixed_max, 0, fixed_max), c("blue", "black", "yellow"))
}

# raw_names/display_names must be the SAME length and in the same order: raw_names is
# matched against the score file's Regulator column (which uses the internal cleaned
# name, e.g. "L_lysine", not the restored original "L-lysine"), and the resulting score
# suffix is appended to the corresponding (possibly restored) display_names entry.
score_suffix_labels <- function(raw_names, display_names, score_df, module_id) {
  if (is.null(score_df)) return(display_names)
  mod_scores <- score_df[as.character(score_df$Target) == as.character(module_id), ]
  score_map <- setNames(mod_scores$Score, mod_scores$Regulator)
  mapply(function(raw, disp) {
    if (raw %in% names(score_map)) sprintf("%s (%d)", disp, round(score_map[[raw]])) else disp
  }, raw_names, display_names, USE.NAMES = FALSE)
}

## ============================ Per-module render =============================

process_module <- function(module_id) {
  genes <- clusters[[module_id]]
  if (is.null(genes)) { cat(sprintf("Module %s: not found in clusters_list, skipping\n", module_id)); return(invisible(NULL)) }

  expr_mat <- expr_all[rownames(expr_all) %in% genes, , drop = FALSE]
  if (nrow(expr_mat) == 0) { cat(sprintf("Module %s: no genes matched expression data, skipping\n", module_id)); return(invisible(NULL)) }
  if (nrow(expr_mat) < length(genes)) {
    missing <- setdiff(genes, rownames(expr_mat))
    cat(sprintf("Warning: Module %s requested %d genes, found %d in expression data. Missing examples: %s\n",
                module_id, length(genes), nrow(expr_mat), paste(head(missing, 10), collapse = ", ")))
  }
  cat(sprintf("Processing module %s: %d/%d genes matched\n", module_id, nrow(expr_mat), length(genes)))

  if (nrow(expr_mat) < 2) {
    eigengene <- colMeans(expr_mat, na.rm = TRUE)
  } else {
    pca_input <- expr_mat
    if (anyNA(pca_input)) {
      row_means <- rowMeans(pca_input, na.rm = TRUE)
      for (i in seq_len(nrow(pca_input))) {
        na_j <- is.na(pca_input[i, ])
        if (any(na_j)) pca_input[i, na_j] <- if (is.nan(row_means[i])) 0 else row_means[i]
      }
    }
    pca <- prcomp(t(pca_input), scale. = FALSE, center = TRUE)
    eigengene <- pca$x[, 1]
  }
  sorted_samples <- names(sort(eigengene))
  expr_mat <- expr_mat[, sorted_samples, drop = FALSE]

  ## ---- regulator blocks, in the order given by regulator_specs ----
  reg_blocks <- list()
  for (reg in regulators) {
    feats <- reg$regs_list[[module_id]]
    if (is.null(feats)) next
    present <- feats[feats %in% rownames(reg$omics_data)]
    if (length(present) == 0) next
    avail_samples <- intersect(sorted_samples, colnames(reg$omics_data))
    mat <- matrix(NA_real_, nrow = length(present), ncol = length(sorted_samples),
                  dimnames = list(present, sorted_samples))
    mat[, avail_samples] <- reg$omics_data[present, avail_samples, drop = FALSE]
    if (reg$is_tf) {
      mat <- t(scale(t(mat)))
      mat[is.na(mat)] <- 0
    }
    display <- restore_names(rownames(mat))
    if (SHOW_REGULATOR_SCORES) display <- score_suffix_labels(rownames(mat), display, reg$score_df, module_id)
    reg_blocks[[length(reg_blocks) + 1]] <- list(mat = mat, title = reg$name, row_labels = display)
  }

  expr_col_fun <- robust_col_fun(expr_mat, fixed_max = 2.0)
  labels_expr <- restore_names(rownames(expr_mat))

  ## ---- column annotation (sample metadata) ----
  anno_df <- data.frame(row.names = sorted_samples)
  anno_colors <- list()
  for (type_name in names(selected_metadata)) {
    entry <- selected_metadata[[type_name]]
    label <- entry$label
    color_to_label <- setNames(names(entry$legend), unlist(entry$legend))
    labs <- sapply(sorted_samples, function(s) {
      col <- entry$samples[[s]]
      if (is.null(col)) return(NA_character_)
      lab <- color_to_label[[col]]
      if (is.null(lab)) col else lab
    })
    anno_df[[label]] <- labs
    anno_colors[[label]] <- web_to_r_color(unlist(entry$legend))
  }
  col_anno <- HeatmapAnnotation(df = anno_df, col = anno_colors,
                                 annotation_name_side = "left",
                                 simple_anno_size = unit(0.35, "cm"))

  ## ---- combined KG match panel (metabolite + TF-target ORA) as
  ## ---- right_annotation on the Expression heatmap.
  build_kg_cols <- function(kg_df, max_cols, color, legend_label) {
    mod_data <- kg_df[kg_df$Module == as.character(module_id), ]
    if (nrow(mod_data) == 0) return(NULL)
    values <- rank_by_shared_genes(mod_data, max_cols)
    mat <- matrix(0L, nrow = length(labels_expr), ncol = length(values),
                  dimnames = list(labels_expr, truncate_label(values)))
    for (v in values) {
      gl <- mod_data$Genes[mod_data$Value == v]
      hit <- intersect(unlist(strsplit(gl, "\\|")), labels_expr)
      if (length(hit) > 0) mat[hit, truncate_label(v)] <- 1L
    }
    list(mat = mat, color = color, legend_label = legend_label)
  }
  kg_blocks <- list(
    build_kg_cols(metabo_kg, NULL, web_to_r_color(tolower(metabo_meta$COLOR %||% "yellow")), "Metabolite interaction (PKN)"),
    build_kg_cols(tfora_kg, NULL, web_to_r_color(tolower(tfora_meta$COLOR %||% "purple")), "TF-target ORA (validated)"),
    build_kg_cols(homer_kg, NULL, web_to_r_color(tolower(homer_meta$COLOR %||% "teal")), "TF motif, promoter (HOMER)"),
    build_kg_cols(enhancer_kg, NULL, web_to_r_color(tolower(enhancer_meta$COLOR %||% "darkcyan")), "TF motif, enhancer (HOMER)"),
    build_kg_cols(irf_enh_kg, NULL, web_to_r_color(tolower(irf_enh_meta$COLOR %||% "magenta")), "Enhancer with IRF motif (FANTOM5 link)")
  )
  kg_blocks <- kg_blocks[!sapply(kg_blocks, is.null)]

  kg_row_anno <- NULL
  kg_legend <- NULL
  kg_col_list <- list()
  if (length(kg_blocks) > 0) {
    kg_mat_all <- do.call(cbind, lapply(kg_blocks, `[[`, "mat"))
    kg_df <- as.data.frame(kg_mat_all)
    for (cn in colnames(kg_df)) kg_df[[cn]] <- as.character(kg_df[[cn]])
    col_idx <- 1
    for (blk in kg_blocks) {
      for (j in seq_len(ncol(blk$mat))) {
        cn <- colnames(kg_df)[col_idx]
        kg_col_list[[cn]] <- c("0" = "white", "1" = blk$color)
        col_idx <- col_idx + 1
      }
    }
    kg_row_anno <- rowAnnotation(df = kg_df, col = kg_col_list,
                                  show_legend = FALSE,
                                  simple_anno_size = unit(0.35, "cm"),
                                  annotation_name_rot = 90,
                                  annotation_name_gp = gpar(fontsize = 6),
                                  gp = gpar(col = "grey70", lwd = 0.5))
    kg_legend <- Legend(labels = sapply(kg_blocks, `[[`, "legend_label"),
                         legend_gp = gpar(fill = sapply(kg_blocks, `[[`, "color")),
                         title = "Gene-KG interactions")
  }

  ## ---- PPI/HumanNet arc connectors: left_annotation on the Expression
  ## ---- heatmap via anno_empty() + a post-draw decorate_annotation callback.
  resolve_pairs <- function(kg_df) {
    mod_data <- kg_df[kg_df$Module == as.character(module_id), ]
    if (nrow(mod_data) == 0) return(list())
    pairs <- list()
    for (i in seq_len(nrow(mod_data))) {
      gp <- strsplit(mod_data$Genes[i], "\\|")[[1]]
      if (length(gp) != 2) next
      i1 <- match(gp[1], labels_expr); i2 <- match(gp[2], labels_expr)
      if (!is.na(i1) && !is.na(i2)) pairs[[length(pairs) + 1]] <- c(i1, i2)
    }
    pairs
  }
  ppi_pairs <- resolve_pairs(ppi_ints)
  hn_pairs  <- resolve_pairs(hn_ints)

  fix_color <- function(col, bad, fallback) if (tolower(col) %in% bad) fallback else col
  ppi_color <- fix_color(web_to_r_color(tolower(ppi_meta$COLOR %||% "darkgreen")), c("blue", "darkblue"), "darkgreen")
  hn_color  <- fix_color(web_to_r_color(tolower(hn_meta$COLOR %||% "saddlebrown")), c("orange", "darkorange"), "saddlebrown")

  arc_anno <- NULL
  arc_width_cm <- 0
  n_rows <- length(labels_expr)
  if (length(ppi_pairs) > 0 || length(hn_pairs) > 0) {
    all_pairs <- c(hn_pairs, ppi_pairs)  # HumanNet drawn first (behind), PPI on top
    max_dist <- max(sapply(all_pairs, function(p) abs(p[2] - p[1])))
    arc_width_cm <- max(1.0, min(4.0, 1.0 + max_dist / n_rows * 3.0))
    arc_anno <- HeatmapAnnotation(arcs = anno_empty(border = FALSE, width = unit(arc_width_cm, "cm")),
                                   which = "row")
  }

  draw_arc_set <- function(pairs, color) {
    if (length(pairs) == 0) return(invisible(NULL))
    for (p in pairs) {
      i1 <- p[1]; i2 <- p[2]
      y1 <- 1 - (i1 - 0.5) / n_rows
      y2 <- 1 - (i2 - 0.5) / n_rows
      dist <- abs(i2 - i1)
      ctrl_x <- max(0.05, 1 - 0.9 * (dist / n_rows))
      t <- seq(0, 1, length.out = 15)
      xs <- (1 - t)^2 * 1 + 2 * (1 - t) * t * ctrl_x + t^2 * 1
      ys <- (1 - t)^2 * y1 + 2 * (1 - t) * t * ((y1 + y2) / 2) + t^2 * y2
      grid.lines(x = xs, y = ys, default.units = "npc",
                 gp = gpar(col = color, lwd = 1, alpha = 0.5))
    }
  }

  ## ---- build the stacked heatmaps: regulator blocks (in configured order), then Expression ----
  build_ht <- function(mat, col_fun, title, row_labels_ = NULL, col_anno_ = NULL, right_anno_ = NULL, left_anno_ = NULL) {
    n <- nrow(mat)
    Heatmap(mat, name = title, col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = show_names(n), row_names_side = "right",
            row_labels = row_labels_ %||% rownames(mat),
            row_names_gp = gpar(fontsize = fontsize_for(n)),
            show_column_names = FALSE,
            height = unit(row_h_cm(n) * n, "cm"),
            use_raster = n > 100,
            border = TRUE,
            column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
            bottom_annotation = col_anno_,
            right_annotation = right_anno_,
            left_annotation = left_anno_,
            heatmap_legend_param = list(title = paste0(title, " value")))
  }

  ht_list <- NULL
  for (blk in reg_blocks) {
    col_fun <- robust_col_fun(blk$mat)
    h <- build_ht(blk$mat, col_fun, blk$title, row_labels_ = blk$row_labels)
    ht_list <- if (is.null(ht_list)) h else ht_list %v% h
  }
  expr_ht <- build_ht(expr_mat, expr_col_fun, "Expression", row_labels_ = labels_expr,
                       col_anno_ = col_anno, right_anno_ = kg_row_anno, left_anno_ = arc_anno)
  ht_list <- if (is.null(ht_list)) expr_ht else ht_list %v% expr_ht

  ## ---- canvas sizing: sum the ACTUAL (capped) per-block heights ----
  cm_to_px <- function(cm) round(cm * res / 2.54)
  blocks_n <- c(sapply(reg_blocks, function(b) nrow(b$mat)), nrow(expr_mat))
  content_height_cm <- sum(sapply(blocks_n, function(n) row_h_cm(n) * n)) +
    length(blocks_n) * 1.0 + (length(blocks_n) - 1) * 0.4 +
    length(selected_metadata) * 0.45 + 1.4
  n_discrete_legend_rows <- sum(sapply(selected_metadata, function(e) length(e$legend) + 1)) +
    (if (!is.null(kg_legend)) length(kg_blocks) + 1 else 0)
  n_continuous_legends <- length(blocks_n)
  legend_height_cm <- n_discrete_legend_rows * 0.5 + n_continuous_legends * 3.5
  height_cm <- max(content_height_cm, legend_height_cm) + 2
  height_px <- min(cm_to_px(height_cm), 9000)
  width_px <- 2400 + length(kg_col_list) * 22 + cm_to_px(arc_width_cm)

  out_png <- file.path(output_dir, sprintf("Module_%s_heatmap.png", module_id))
  out_pdf <- file.path(output_dir, sprintf("Module_%s_heatmap.pdf", module_id))

  render <- function(dev_open, dev_close) {
    dev_open()
    draw(ht_list, column_title = paste("Module", module_id),
         column_title_gp = gpar(fontsize = 16, fontface = "bold"),
         heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = FALSE, ht_gap = unit(4, "mm"),
         annotation_legend_list = if (!is.null(kg_legend)) list(kg_legend) else list())
    if (!is.null(arc_anno)) {
      decorate_annotation("arcs", {
        draw_arc_set(hn_pairs, hn_color)
        draw_arc_set(ppi_pairs, ppi_color)
      })
    }
    dev_close()
  }
  render(function() png(out_png, width = width_px, height = height_px, res = res), dev.off)
  render(function() pdf(out_pdf, width = width_px / res, height = height_px / res), dev.off)

  cat(sprintf("  Saved module %s (%d x %d px)\n", module_id, width_px, height_px))
}

## ============================== Main loop ===================================

cli_modules <- commandArgs(trailingOnly = TRUE)
only_modules <- if (length(cli_modules) >= 1) strsplit(cli_modules[1], ",")[[1]] else NULL

modules_processed <- 0
for (mid in names(clusters)) {
  if (!is.null(only_modules) && !(mid %in% only_modules)) next
  result <- try(process_module(mid), silent = FALSE)
  if (!inherits(result, "try-error")) modules_processed <- modules_processed + 1
}

cat(sprintf("\nProcessing complete! Generated heatmaps for %d modules\n", modules_processed))
cat(sprintf("Output directory: %s\n", output_dir))
