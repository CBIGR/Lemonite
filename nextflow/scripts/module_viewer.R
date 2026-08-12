#!/usr/bin/env Rscript
# ModuleViewer: ComplexHeatmap-based module heatmap visualizations for LemonTree modules.
#
# Replaces the matplotlib-based module_viewer.py. Same CLI contract (see
# viewer_heatmaps.nf), same input files (ModuleViewer_files/*.mvf + selected-regulator
# lists, Preprocessing/LemonPreprocessed_*.txt), same output naming
# (heatmaps/Module_<id>_heatmap.{png,pdf}).
#
# Why the rewrite: the matplotlib version applied one hardcoded +/-2 color scale to
# every omics block (crushing metabolite/hPTM/etc. blocks that don't live on that
# scale), sized the whole figure from figsize=(15, total_rows/2+8) with no cap (a
# ~250-feature module produced a 4545x30276px image), and positioned legends with
# manual bbox_to_anchor fractions that collide once a track has many categories.
# ComplexHeatmap gives each omics block its own color scale, packs legends without
# manual coordinate math, and lets per-block row height be capped independently of
# canvas size.
#
# Usage: Rscript module_viewer.R --input_dir . --output_dir heatmaps \
#          --regulator_files "TFs:TFs.selected_regs_list.txt,Metabolites:Metabolites.selected_regs_list.txt" \
#          --regulator_types "TFs:Lovering_TF_list.txt,Metabolites:Metabolomics.txt" \
#          --dpi 300 --show_regulator_scores --annotation_types diagnosis \
#          --name_mapping Preprocessing/name_mapping.tsv [--modules 0,1,2]

suppressMessages({
  library(optparse)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
ht_opt$message <- FALSE
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

CORUM_MAX_COLS <- 10

## ============================== CLI args ===================================

opt_list <- list(
  make_option("--input_dir", type = "character"),
  make_option("--output_dir", type = "character"),
  make_option("--regulator_files", type = "character"),
  make_option("--expression_file", type = "character", default = "LemonPreprocessed_expression.txt"),
  make_option("--complete_file", type = "character", default = "LemonPreprocessed_complete.txt"),
  make_option("--regulator_types", type = "character", default = NULL),
  make_option("--show_regulator_scores", action = "store_true", default = FALSE),
  make_option("--dpi", type = "integer", default = 300),
  make_option("--modules", type = "character", default = NULL),
  make_option("--annotation_types", type = "character", default = "diagnosis"),
  make_option("--annotation_labels", type = "character", default = NULL),
  make_option("--name_mapping", type = "character", default = NULL)
)
args <- parse_args(OptionParser(option_list = opt_list))
for (req in c("input_dir", "output_dir", "regulator_files")) {
  if (is.null(args[[req]])) stop(sprintf("--%s is required", req))
}

input_dir  <- args$input_dir
output_dir <- args$output_dir
viewer_dir <- file.path(input_dir, "ModuleViewer_files")
if (!dir.exists(viewer_dir)) {
  stop(sprintf("ModuleViewer files directory not found: %s\nEnsure the pipeline produced ModuleViewer_files with sample_mapping.mvf and metabolite_LemoniteKG_interactions.mvf", viewer_dir))
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
res <- args$dpi

cat(sprintf("Processing ModuleViewer data from: %s\n", viewer_dir))
cat(sprintf("Output directory: %s\n", output_dir))
cat(sprintf("Regulator files: %s\n", args$regulator_files))

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

# Mirrors the R preprocessing script's basename derivation
# (tolower(file_path_sans_ext(basename(...)))) so regulator types map onto their
# LemonPreprocessed_{basename}.txt file, e.g. "Proteins:proteomics.txt" ->
# LemonPreprocessed_proteomics.txt. A trailing ":c"/":d" continuous/discrete
# suffix (see --regulator_types docs) rides along inside what looks like the
# extension and drops out along with it, same as the old Python version's
# os.path.splitext behavior.
basename_no_ext <- function(path) {
  b <- basename(path)
  sub("\\.[^.]*$", "", b)
}

## Search several candidate locations for a data/score file, mirroring
## module_viewer.py's find_data_file()/score_file_locations fallbacks -- the
## exact preprocessing output layout varies (bind-mounted dev run vs. staged
## Nextflow work dir vs. a results/ tree).
find_file <- function(candidates) {
  for (loc in candidates) if (file.exists(loc)) return(loc)
  hits <- Sys.glob(file.path(input_dir, "results", "*", "LemonTree", "Preprocessing", basename(candidates[1])))
  if (length(hits) > 0) return(hits[1])
  NULL
}

## ============================== Load shared inputs ==========================

name_lookup <- c()
if (!is.null(args$name_mapping) && file.exists(args$name_mapping)) {
  nm <- read.delim(args$name_mapping)
  if (all(c("original", "mapped") %in% colnames(nm)) && nrow(nm) > 0) {
    name_lookup <- setNames(nm$original, nm$mapped)
    cat(sprintf("Loaded %d name mappings from %s\n", length(name_lookup), args$name_mapping))
  } else if (all(c("cleaned", "original") %in% colnames(nm)) && nrow(nm) > 0) {
    name_lookup <- setNames(nm$original, nm$cleaned)
    cat(sprintf("Loaded %d name mappings from %s\n", length(name_lookup), args$name_mapping))
  }
}
restore_names <- function(x) {
  if (length(name_lookup) == 0) return(x)
  ifelse(x %in% names(name_lookup), unname(name_lookup[x]), x)
}

read_mat <- function(path) {
  d <- read.delim(path, check.names = FALSE)
  if ("symbol" %in% colnames(d)) {
    rownames(d) <- d$symbol
    d$symbol <- NULL
  }
  # Drop any stray non-numeric/all-NA columns (mirrors create_subset()'s junk-column guard)
  junk <- vapply(d, function(col) !is.numeric(col) || all(is.na(col)), logical(1))
  as.matrix(d[, !junk, drop = FALSE])
}

expr_path <- find_file(c(
  file.path(viewer_dir, "..", "Preprocessing", args$expression_file),
  file.path(input_dir, args$expression_file),
  file.path(input_dir, "Preprocessing", args$expression_file)
))
if (is.null(expr_path)) stop(sprintf("Expression data file not found: %s", args$expression_file))
cat(sprintf("Loading expression data from: %s\n", expr_path))
expr_all <- read_mat(expr_path)

clusters <- parse_list_file(file.path(viewer_dir, "clusters_list.txt"))
sample_mapping <- parse_sample_mapping(file.path(viewer_dir, "sample_mapping.mvf"))

## specific_modules.txt (written by LemonTree_to_network after post-clustering
## consolidation) restricts which modules actually get rendered -- ports the
## same filter module_viewer.py applied so standalone runs of this script
## behave the same way even without the Nextflow --modules passthrough.
specific_modules_file <- file.path(input_dir, "Networks", "specific_modules.txt")
if (file.exists(specific_modules_file)) {
  specific_modules <- trimws(readLines(specific_modules_file))
  specific_modules <- specific_modules[nzchar(specific_modules)]
  before <- length(clusters)
  clusters <- clusters[names(clusters) %in% specific_modules]
  cat(sprintf("Modules before filtering: %d\nModules after filtering: %d\n", before, length(clusters)))
} else {
  cat("Warning: specific_modules.txt not found - processing all modules\n")
}
if (!is.null(args$modules)) {
  requested <- trimws(strsplit(args$modules, ",")[[1]])
  requested <- requested[nzchar(requested)]
  clusters <- clusters[names(clusters) %in% requested]
  cat(sprintf("Processing only requested modules: %s\n", paste(requested, collapse = ", ")))
}

## ---- regulator configuration (config-driven: parses --regulator_files /
## ---- --regulator_types the same way module_viewer.py does, so any
## ---- case-study's regulator_types set -- TFs, Metabolites, Proteins,
## ---- hPTMs, whatever --regulator_types "Prefix:File[:c|d]" pairs were
## ---- configured for this run -- works without code changes here.
regtype_to_basename <- list()
if (!is.null(args$regulator_types)) {
  for (entry in strsplit(args$regulator_types, ",")[[1]]) {
    entry <- trimws(entry)
    if (grepl(":", entry)) {
      parts <- strsplit(entry, ":")[[1]]
      rtype <- trimws(parts[1])
      rfile <- trimws(paste(parts[-1], collapse = ":"))
      regtype_to_basename[[rtype]] <- tolower(basename_no_ext(rfile))
    }
  }
  cat("Regulator type -> data file basename mapping:\n")
  for (n in names(regtype_to_basename)) cat(sprintf("  %s -> %s\n", n, regtype_to_basename[[n]]))
}

regulators <- list()
for (entry in strsplit(args$regulator_files, ",")[[1]]) {
  entry <- trimws(entry)
  if (!grepl(":", entry) || !nzchar(entry)) next
  parts <- strsplit(entry, ":")[[1]]
  reg_type <- trimws(parts[1])
  reg_path <- trimws(paste(parts[-1], collapse = ":"))
  reg_filename <- basename(reg_path)
  prefix <- strsplit(reg_filename, "\\.")[[1]][1]

  reg_type_lower <- tolower(reg_type)
  is_tf <- grepl("tf", reg_type_lower) || reg_type_lower %in% c("tfs", "transcription_factors")

  regs_path <- file.path(viewer_dir, reg_filename)
  if (!file.exists(regs_path)) {
    cat(sprintf("Warning: regulator file not found: %s\n", regs_path))
    next
  }
  regs_list <- parse_list_file(regs_path)

  score_path <- find_file(c(
    file.path(viewer_dir, sprintf("%s.selected_regulators_scores.txt", prefix)),
    file.path(input_dir, sprintf("%s.selected_regulators_scores.txt", prefix)),
    file.path(input_dir, "Lemon_out", sprintf("%s.selected_regulators_scores.txt", prefix)),
    file.path(viewer_dir, "..", "Lemon_out", sprintf("%s.selected_regulators_scores.txt", prefix))
  ))
  score_df <- if (!is.null(score_path)) read.delim(score_path) else NULL
  if (args$show_regulator_scores) {
    if (!is.null(score_path)) cat(sprintf("Loaded scores for %s from %s\n", reg_type, score_path))
    else cat(sprintf("Warning: score file not found for %s\n", reg_type))
  }

  if (is_tf) {
    omics_data <- expr_all
    cat(sprintf("Using expression data for %s regulators\n", reg_type))
  } else {
    data_basename <- regtype_to_basename[[reg_type]] %||% reg_type_lower
    omics_filename <- sprintf("LemonPreprocessed_%s.txt", data_basename)
    omics_path <- find_file(c(
      file.path(viewer_dir, "..", "Preprocessing", omics_filename),
      file.path(input_dir, omics_filename),
      file.path(input_dir, "Preprocessing", omics_filename)
    ))
    if (!is.null(omics_path)) {
      omics_data <- read_mat(omics_path)
      cat(sprintf("Loaded omics-specific data for %s from %s\n", reg_type, omics_path))
    } else {
      cat(sprintf("Falling back to complete data for %s\n", reg_type))
      complete_path <- find_file(c(
        file.path(viewer_dir, "..", "Preprocessing", args$complete_file),
        file.path(input_dir, args$complete_file),
        file.path(input_dir, "Preprocessing", args$complete_file)
      ))
      omics_data <- if (!is.null(complete_path)) read_mat(complete_path) else expr_all
    }
  }

  regulators[[length(regulators) + 1]] <- list(
    name = reg_type, is_tf = is_tf, regs_list = regs_list,
    omics_data = omics_data, score_df = score_df
  )
}
cat(sprintf("Configured %d regulator types: %s\n", length(regulators), paste(sapply(regulators, `[[`, "name"), collapse = ", ")))

## ---- KG / interaction sources (all optional, guarded by file existence) ----
metabo_kg  <- parse_kg_mvf(file.path(viewer_dir, "metabolite_LemoniteKG_interactions.mvf"))
metabo_meta <- parse_mvf_metadata(file.path(viewer_dir, "metabolite_LemoniteKG_interactions.mvf"))
phospho_kg <- parse_kg_mvf(file.path(viewer_dir, "phospho_LemoniteKG_interactions.mvf"))
phospho_meta <- parse_mvf_metadata(file.path(viewer_dir, "phospho_LemoniteKG_interactions.mvf"))
corum_kg   <- parse_kg_mvf(file.path(viewer_dir, "protein_complex_CORUM_interactions.mvf"))
corum_meta <- parse_mvf_metadata(file.path(viewer_dir, "protein_complex_CORUM_interactions.mvf"))
ppi_ints  <- parse_kg_mvf(file.path(viewer_dir, "PPI_interactions.mvf"))
ppi_meta  <- parse_mvf_metadata(file.path(viewer_dir, "PPI_interactions.mvf"))
hn_ints   <- parse_kg_mvf(file.path(viewer_dir, "HumanNet_interactions.mvf"))
hn_meta   <- parse_mvf_metadata(file.path(viewer_dir, "HumanNet_interactions.mvf"))

## ---- sample annotation types: honor --annotation_types / --annotation_labels,
## ---- falling back the same way module_viewer.py did if a requested type is
## ---- missing from sample_mapping.mvf.
requested_annotations <- trimws(strsplit(args$annotation_types, ",")[[1]])
selected_metadata <- sample_mapping[names(sample_mapping) %in% requested_annotations]
missing <- setdiff(requested_annotations, names(sample_mapping))
if (length(missing) > 0) {
  cat(sprintf("Warning: requested annotation type(s) not found in sample_mapping.mvf: %s\n", paste(missing, collapse = ", ")))
  cat(sprintf("Available types: %s\n", paste(names(sample_mapping), collapse = ", ")))
}
if (length(selected_metadata) == 0) {
  cat("Error: no valid annotation types found. Falling back to 'diagnosis' or the first available type.\n")
  fallback <- if ("diagnosis" %in% names(sample_mapping)) "diagnosis" else names(sample_mapping)[1]
  selected_metadata <- sample_mapping[fallback]
}
custom_labels <- list()
if (!is.null(args$annotation_labels)) {
  labels_list <- trimws(strsplit(args$annotation_labels, ",")[[1]])
  if (length(labels_list) == length(selected_metadata)) {
    custom_labels <- setNames(as.list(labels_list), names(selected_metadata))
  } else {
    cat(sprintf("Warning: number of custom labels (%d) doesn't match annotations (%d), using defaults\n",
                length(labels_list), length(selected_metadata)))
  }
}
cat(sprintf("Displaying %d annotation track(s): %s\n", length(selected_metadata), paste(names(selected_metadata), collapse = ", ")))

cat(sprintf("%d modules to process, PPI rows=%d, HumanNet rows=%d, CORUM rows=%d, Metabolite-KG rows=%d, Phospho-KG rows=%d\n",
            length(clusters), nrow(ppi_ints), nrow(hn_ints), nrow(corum_kg), nrow(metabo_kg), nrow(phospho_kg)))

## ============================== Sizing helpers ==============================

show_names <- function(n) n <= 150
fontsize_for <- function(n) max(3, min(9, 200 / n))
row_h_cm <- function(n) {
  max_total_cm <- 34
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

score_suffix_labels <- function(feature_names, score_df, module_id) {
  if (is.null(score_df)) return(feature_names)
  mod_scores <- score_df[as.character(score_df$Target) == as.character(module_id), ]
  score_map <- setNames(mod_scores$Score, mod_scores$Regulator)
  sapply(feature_names, function(f) {
    if (f %in% names(score_map)) sprintf("%s (%d)", f, round(score_map[[f]])) else f
  })
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
    eigengene <- colMeans(expr_mat)
  } else {
    pca <- prcomp(t(expr_mat), scale. = FALSE, center = TRUE)
    eigengene <- pca$x[, 1]
  }
  sorted_samples <- names(sort(eigengene))
  expr_mat <- expr_mat[, sorted_samples, drop = FALSE]

  ## ---- regulator blocks, in the order given by --regulator_files ----
  reg_blocks <- list()
  for (reg in regulators) {
    feats <- reg$regs_list[[module_id]]
    if (is.null(feats)) next
    present <- feats[feats %in% rownames(reg$omics_data)]
    if (length(present) == 0) next
    mat <- reg$omics_data[present, sorted_samples, drop = FALSE]
    # TFA scores aren't pre-scaled in LemonPreprocessed_complete.txt -- row-wise
    # z-score them (same special case module_viewer.py applied to TF blocks).
    if (reg$is_tf) {
      mat <- t(scale(t(mat)))
      mat[is.na(mat)] <- 0
    }
    display <- restore_names(rownames(mat))
    if (args$show_regulator_scores) display <- score_suffix_labels(display, reg$score_df, module_id)
    reg_blocks[[length(reg_blocks) + 1]] <- list(mat = mat, title = reg$name, row_labels = display)
  }

  expr_col_fun <- robust_col_fun(expr_mat, fixed_max = 2.0)
  labels_expr <- restore_names(rownames(expr_mat))

  ## ---- column annotation (sample metadata) ----
  anno_df <- data.frame(row.names = sorted_samples)
  anno_colors <- list()
  for (type_name in names(selected_metadata)) {
    entry <- selected_metadata[[type_name]]
    label <- custom_labels[[type_name]] %||% entry$label
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

  ## ---- combined KG match panel (metabolite + phospho + CORUM) as
  ## ---- right_annotation on the Expression heatmap. Matching uses the
  ## ---- restored (original) gene names, since the KG mvf files reference
  ## ---- genes by their real symbols, not the preprocessing-internal cleaned
  ## ---- ones -- same order of operations module_viewer.py used.
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
    build_kg_cols(metabo_kg, NULL, web_to_r_color(tolower(metabo_meta$COLOR %||% "yellow")), "Metabolite interaction"),
    build_kg_cols(phospho_kg, NULL, web_to_r_color(tolower(phospho_meta$COLOR %||% "orange")), "Phosphosite interaction"),
    build_kg_cols(corum_kg, CORUM_MAX_COLS, web_to_r_color(tolower(corum_meta$COLOR %||% "orange")), "CORUM complex")
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
  ## ---- Row i maps to npc y = 1 - (i-0.5)/n (verified against a standalone
  ## ---- grid test before first use).
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

  ## ---- canvas sizing: sum the ACTUAL (capped) per-block heights, not a
  ## ---- linear-in-rows guess -- the guess overshoots once row_h_cm's cap
  ## ---- kicks in for a large block, leaving ComplexHeatmap to center a much
  ## ---- shorter drawing inside an oversized device (a mostly-blank PNG).
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

modules_processed <- 0
for (mid in names(clusters)) {
  result <- try(process_module(mid), silent = FALSE)
  if (!inherits(result, "try-error")) modules_processed <- modules_processed + 1
}

cat(sprintf("\nProcessing complete! Generated heatmaps for %d modules\n", modules_processed))
cat(sprintf("Output directory: %s\n", output_dir))
