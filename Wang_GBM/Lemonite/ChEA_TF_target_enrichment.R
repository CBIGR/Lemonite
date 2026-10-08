#!/usr/bin/Rscript

############################################################################################################################################
#### Overrepresentation analysis (ORA) validating predicted TF regulators:
#### for each module, test whether its genes are enriched for the known target genes of
#### the TF(s) LemonTree assigned to that module, using several complementary EnrichR
#### TF-target libraries (ChEA/ENCODE ChIP-seq, TRRUST curated targets, TF perturbation
#### signatures, and JASPAR binding-motif matches).
####
#### Many of the TFs LemonTree assigns here are tissue-specific homeobox/developmental
#### factors that ChEA alone rarely has ChIP-seq data for, so relying on ChEA only leaves
#### most predictions untestable rather than "not significant" - pooling several
#### experimental libraries closes most of that gap. Motif (JASPAR) evidence is weaker
#### (a predicted binding site, not an observed regulatory event) and is kept in a
#### separate column rather than merged into the experimental verdict.
####
#### Output is written as a ModuleViewer-compatible .mvf annotation file so the result can
#### be drawn as an extra right-side panel next to the existing metabolite-gene interaction
#### panel in ModuleViewer.ipynb.
############################################################################################################################################

library(data.table)
library(stringr)
library(enrichR)

set.seed(1234)

# ------------------------------------------------------------------------------------------------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------------------------------------------------------------------------------------------------

base_dir <- '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/'

percentile <- 2               # matches the *.percentileN_list.txt naming used by ModuleViewer_files
n_modules_name <- '46'        # for naming the output CSV only

# EnrichR TF-target libraries to pool as "experimental" evidence (ChIP-seq/ChIP-chip,
# curated literature targets, and TF-perturbation-followed-by-expression signatures), plus
# one motif-based library kept separate because it reflects predicted binding sites rather
# than an observed regulatory event.
EXPERIMENTAL_DBS <- c(
  "ChEA_2022",
  "ENCODE_TF_ChIP-seq_2015",
  "TRRUST_Transcription_Factors_2019",
  "TF_Perturbations_Followed_by_Expression"
)
MOTIF_DBS <- c("JASPAR_PWM_Human_2025")
ALL_DBS <- c(EXPERIMENTAL_DBS, MOTIF_DBS)

padj_threshold <- 0.05
sleep_between_calls <- 1      # seconds; EnrichR rate limiting

clusters_file         <- paste0(base_dir, 'ModuleViewer_files/clusters_list.txt')
tf_regulators_file    <- paste0(base_dir, 'ModuleViewer_files/Lovering.percentile', percentile, '_list.txt')
specific_modules_file <- paste0(base_dir, 'Networks/specific_modules.txt')

output_dir      <- paste0(base_dir, 'Enrichment/ChEA_TF_validation')
mvf_output_file <- paste0(base_dir, 'ModuleViewer_files/TF_target_ORA_interactions.mvf')

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------------------------------------------------------------------
# LOAD MODULES AND THEIR PREDICTED TF REGULATORS
# ------------------------------------------------------------------------------------------------------------------------------------------

read_module_list <- function(file, value_col) {
  dt <- fread(file, header = FALSE, col.names = c("Module", value_col))
  dt[[value_col]] <- str_split(dt[[value_col]], "\\|")
  dt
}

clusters <- read_module_list(clusters_file, "Genes")
tf_regs  <- read_module_list(tf_regulators_file, "TFs")

if (file.exists(specific_modules_file)) {
  specific_modules <- gsub("\\s+", "", readLines(specific_modules_file))
  before <- nrow(clusters)
  clusters <- clusters[clusters$Module %in% specific_modules, ]
  cat("Filtered to", nrow(clusters), "of", before, "modules using specific_modules.txt\n")
} else {
  cat("Warning: specific_modules.txt not found - testing all modules\n")
}

module2genes <- setNames(clusters$Genes, as.character(clusters$Module))
module2tfs   <- setNames(tf_regs$TFs, as.character(tf_regs$Module))

# ------------------------------------------------------------------------------------------------------------------------------------------
# ENRICHR HELPERS
# ------------------------------------------------------------------------------------------------------------------------------------------

setEnrichrSite("Enrichr")

safe_enrichr <- function(genes, dbs, max_tries = 5) {
  for (i in 1:max_tries) {
    result <- tryCatch(enrichr(genes, dbs), error = function(e) e)
    if (!inherits(result, "error")) return(result)
    if (grepl("429", result$message)) {
      wait_time <- 2 ^ i
      cat("Rate limited. Retrying in", wait_time, "seconds...\n")
      Sys.sleep(wait_time)
    } else {
      stop(result)
    }
  }
  stop("Failed after multiple attempts due to rate limiting.")
}

# TF-target library terms are named "<TF_SYMBOL> <source id/description...>", so a TF's
# terms are found by matching the leading token of Term (case-insensitive). Returns the
# best (lowest padj) row among a term subset, or NULL if empty.
best_hit <- function(terms) {
  if (nrow(terms) == 0) return(NULL)
  terms[order(terms$Adjusted.P.value), ][1, ]
}

# For one module: test each of its predicted TFs against the module's own ORA results
# across all configured databases (already fetched for this module in `res_by_db`).
#
# Returns one row per module-TF pair with:
#   Status                   - "no_data" / "tested_not_significant" / "tested_significant",
#                               derived from EXPERIMENTAL_DBS only (ChIP-seq/curated/perturbation)
#   Significant_Databases    - which experimental database(s) drove a "tested_significant" call
#   Target_genes_in_module   - union of module genes overlapping the significant experimental term(s)
#   Motif_Support            - TRUE/FALSE/NA(no data): independent JASPAR motif-match verdict,
#                               kept separate since a motif match is weaker evidence than an
#                               observed ChIP-seq/perturbation/curated target relationship
# Collapsing "no data" into "not significant" would misrepresent untestable TFs as failed
# validations, so the distinction is kept explicit end-to-end.
test_module_tf_enrichment <- function(module, genes, tfs, res_by_db) {
  rows <- lapply(tfs, function(tf) {
    tf_upper <- toupper(tf)

    sig_dbs <- character(0)
    overlap_genes <- character(0)
    any_tested <- FALSE
    best_row <- NULL

    for (db in EXPERIMENTAL_DBS) {
      res <- res_by_db[[db]]
      if (is.null(res) || nrow(res) == 0) next
      tf_terms <- res[res$TF_symbol == tf_upper, ]
      if (nrow(tf_terms) == 0) next
      any_tested <- TRUE

      hits <- tf_terms[tf_terms$Adjusted.P.value <= padj_threshold, ]
      if (nrow(hits) > 0) {
        sig_dbs <- c(sig_dbs, db)
        db_overlap <- unique(unlist(str_split(hits$Genes, ";")))
        overlap_genes <- union(overlap_genes, intersect(db_overlap, genes))
        top <- best_hit(hits)
        if (is.null(best_row) || top$Adjusted.P.value < best_row$Adjusted.P.value) best_row <- top
      }
    }

    status <- if (length(sig_dbs) > 0) "tested_significant" else if (any_tested) "tested_not_significant" else "no_data"

    # Independent motif-based verdict (not folded into `status`)
    motif_support <- NA
    motif_db <- NA
    motif_padj <- NA
    for (db in MOTIF_DBS) {
      res <- res_by_db[[db]]
      if (is.null(res) || nrow(res) == 0) next
      tf_terms <- res[res$TF_symbol == tf_upper, ]
      if (nrow(tf_terms) == 0) next
      if (is.na(motif_support)) motif_support <- FALSE
      hits <- tf_terms[tf_terms$Adjusted.P.value <= padj_threshold, ]
      if (nrow(hits) > 0) {
        top <- best_hit(hits)
        if (is.na(motif_padj) || top$Adjusted.P.value < motif_padj) {
          motif_support <- TRUE
          motif_db <- db
          motif_padj <- top$Adjusted.P.value
        }
      }
    }

    data.frame(
      Module = module,
      TF = tf,
      Status = status,
      Significant_Databases = paste(sig_dbs, collapse = ";"),
      Best_Term = if (is.null(best_row)) NA else best_row$Term,
      P.value = if (is.null(best_row)) NA else best_row$P.value,
      Adjusted.P.value = if (is.null(best_row)) NA else best_row$Adjusted.P.value,
      N_target_genes_in_module = length(overlap_genes),
      Target_genes_in_module = paste(overlap_genes, collapse = "|"),
      Motif_Support = motif_support,
      Motif_Database = motif_db,
      Motif_Adjusted.P.value = motif_padj,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# ------------------------------------------------------------------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------------------------------------------------------------------

cat("Testing TF-target enrichment for", length(module2genes), "modules against:\n")
cat("  experimental:", paste(EXPERIMENTAL_DBS, collapse = ", "), "\n")
cat("  motif:", paste(MOTIF_DBS, collapse = ", "), "\n\n")

all_results <- list()
for (module in names(module2genes)) {
  tfs <- module2tfs[[module]]
  genes <- unique(module2genes[[module]])
  if (is.null(tfs) || length(tfs) == 0 || length(genes) < 3) next

  cat("Module", module, "- predicted TFs:", paste(tfs, collapse = ", "), "\n")

  Sys.sleep(sleep_between_calls)
  enrichment <- tryCatch(safe_enrichr(genes, ALL_DBS), error = function(e) {
    cat("EnrichR failed for module", module, ":", e$message, "\n")
    NULL
  })
  if (is.null(enrichment)) next

  res_by_db <- lapply(enrichment, function(res) {
    if (is.null(res) || nrow(res) == 0) return(res)
    res$TF_symbol <- toupper(sub("^([^ ]+).*", "\\1", res$Term))
    res
  })

  res <- test_module_tf_enrichment(module, genes, tfs, res_by_db)
  if (!is.null(res)) all_results[[module]] <- res
}

results_df <- do.call(rbind, all_results)
if (is.null(results_df)) {
  results_df <- data.frame(Module = character(), TF = character(), Status = character(),
                            Significant_Databases = character(), Best_Term = character(),
                            P.value = numeric(), Adjusted.P.value = numeric(),
                            N_target_genes_in_module = integer(), Target_genes_in_module = character(),
                            Motif_Support = logical(), Motif_Database = character(),
                            Motif_Adjusted.P.value = numeric())
}

results_csv <- file.path(output_dir, paste0("TF_target_enrichment_", n_modules_name, "_modules.csv"))
fwrite(results_df, results_csv)
cat("Saved enrichment table to:", results_csv, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# WRITE MODULEVIEWER-COMPATIBLE .mvf ANNOTATION FILE
# Only "tested_significant" rows (>=1 experimental database) are drawn as annotation hits
# (same 3-column Module / Genes / <label> format as metabolite_LemonIteKG_interactions.mvf);
# "no_data", "tested_not_significant", and motif-only support stay in the CSV for
# transparency but are not plotted, since they are not validated regulatory evidence.
# ------------------------------------------------------------------------------------------------------------------------------------------

sig_df <- results_df[results_df$Status == "tested_significant", ]

mvf_lines <- c(
  "::TYPE=TF_target_ORA",
  "::TITLE:TF_target_ORA",
  "::OBJECT=GENES",
  "::COLOR=PURPLE"
)
if (nrow(sig_df) > 0) {
  mvf_lines <- c(mvf_lines, apply(sig_df, 1, function(r) {
    paste(r[["Module"]], r[["Target_genes_in_module"]], r[["TF"]], sep = "\t")
  }))
}
writeLines(mvf_lines, mvf_output_file)
cat("Saved ModuleViewer annotation file to:", mvf_output_file, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# COVERAGE / RESULT SUMMARY
# ------------------------------------------------------------------------------------------------------------------------------------------

n_pairs <- nrow(results_df)
n_no_data <- sum(results_df$Status == "no_data")
n_tested <- sum(results_df$Status != "no_data")
n_sig <- sum(results_df$Status == "tested_significant")
n_motif_only <- sum(results_df$Status != "tested_significant" & !is.na(results_df$Motif_Support) & results_df$Motif_Support, na.rm = TRUE)

cat("\n=== TF-target ORA complete ===\n")
cat("Module-TF pairs evaluated:", n_pairs, "\n")
cat("  - no data in any experimental database:", n_no_data, "\n")
cat("  - tested (>=1 experimental database has the TF):", n_tested, "\n")
cat("      - significant (padj <=", padj_threshold, "):", n_sig, "\n")
cat("      - not significant:", n_tested - n_sig, "\n")
cat("  - motif-only support (JASPAR) without experimental significance:", n_motif_only, "\n")
cat("Modules with >=1 experimentally-validated predicted TF:",
    length(unique(sig_df$Module)), "out of", length(module2genes), "\n")
