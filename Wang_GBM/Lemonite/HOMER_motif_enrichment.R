#!/usr/bin/Rscript

############################################################################################################################################
#### HOMER known-motif enrichment validating predicted TF regulators.
####
#### Companion to ChEA_TF_target_enrichment.R. Same two-tier logic, different evidence type:
####   tier 1 (per module) - are the promoters (TSS +/-1kb) of a module's genes enriched for the
####                         binding motif of the TF LemonTree assigned to that module, relative to
####                         the promoters of all other expressed genes?
####   tier 2 (per gene)   - for module-TF pairs that pass tier 1, which individual promoters actually
####                         carry the motif? Those genes populate the .mvf annotation panel.
####
#### Tier 1 is what makes the call; tier 2 only decides which cells get coloured. A motif match on its
#### own is weak evidence (a short PWM hits a large fraction of any promoter set by chance), which is
#### why genes are never plotted unless their module-TF pair cleared the enrichment test first - the
#### same reason ChEA_TF_target_enrichment.R keeps JASPAR motif hits in a separate column instead of
#### merging them into the experimental verdict.
####
#### Only known motifs are used (-nomotif disables de novo discovery). Two motif libraries:
####   homer    (default) HOMER's known.motifs only. ~75% of module-TF pairs have no motif there.
####   extended HOMER's known.motifs PLUS JASPAR 2024 / HOCOMOCO v12 matrices for the predicted TFs HOMER
####            lacks, plus paralog-proxy motifs for four TFs with no motif anywhere (BATF2, NEUROD6,
####            TFAP2D, NFE2L3). Built by build_extended_motif_library.py. A paralog motif is used ONLY
####            for a TF with no motif of its own, and every such call is flagged: Relation =
####            "paralog_proxy" plus Motif_TF (which paralog) in the CSV, "TF (via PARALOG)" in the .mvf.
####            q-values are Benjamini over the whole extended library, so they are slightly more
####            conservative than in `homer` mode; the p-values of HOMER's own motifs are identical.
####   Usage: Rscript HOMER_motif_enrichment.R [homer|extended] [module,module,...]
####
#### Background: promoters of every gene in the expression matrix, minus the module's own genes. Using
#### HOMER's default random-genomic background instead would make GC-rich motifs (SP/KLF/E2F) look
#### significant in nearly every module, since promoters are GC-rich relative to the genome at large.
############################################################################################################################################

suppressMessages({
  library(data.table)
  library(stringr)
  library(org.Hs.eg.db)
})

# ------------------------------------------------------------------------------------------------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------------------------------------------------------------------------------------------------

base_dir <- '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/'

cli_args      <- commandArgs(trailingOnly = TRUE)
motif_library <- if (length(cli_args) >= 1) cli_args[1] else "homer"
only_modules  <- if (length(cli_args) >= 2) strsplit(cli_args[2], ",")[[1]] else NULL   # smoke tests
stopifnot(motif_library %in% c("homer", "extended"))
extended <- motif_library == "extended"
lib_tag  <- paste0(if (extended) "_extended" else "", if (!is.null(only_modules)) "_subset" else "")

extended_lib <- paste0(base_dir, 'Enrichment/motif_library/extended_known.motifs')   # HOMER known + JASPAR/HOCOMOCO
custom_map   <- paste0(base_dir, 'Enrichment/motif_library/custom_motif_map.tsv')    # motif <-> predicted TF (+ Relation)

homer_root    <- '/home/borisvdm/Bioinformatics_software/homer'
homer_bin     <- file.path(homer_root, 'bin')
known_motifs  <- file.path(homer_root, 'data/knownTFs/vertebrates/known.motifs')   # HOMER's own library (alias-matched to TFs)
promoter_set  <- 'human'      # HOMER promoter set (human-p); sequence-based, so genome build is irrelevant here

tss_start <- -1000            # bp relative to TSS
tss_end   <-  1000
qval_threshold <- 0.05        # HOMER "q-value (Benjamini)" column
n_cpus <- 8
max_discovery_motifs <- 5     # per module, cap on the discovery panel's width

percentile <- 2               # matches the *.percentileN_list.txt naming used by ModuleViewer_files

clusters_file         <- paste0(base_dir, 'ModuleViewer_files/clusters_list.txt')
tf_regulators_file    <- paste0(base_dir, 'ModuleViewer_files/Lovering.percentile', percentile, '_list.txt')
specific_modules_file <- paste0(base_dir, 'Networks/specific_modules.txt')
expression_file       <- paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression.txt')

output_dir      <- paste0(base_dir, 'Enrichment/HOMER_motif', lib_tag)
work_dir        <- file.path(output_dir, 'runs')
input_dir       <- file.path(output_dir, 'inputs')
mvf_output_file <- paste0(base_dir, 'ModuleViewer_files/HOMER_motif', lib_tag, '_interactions.mvf')
discovery_mvf_file <- paste0(base_dir, 'ModuleViewer_files/HOMER_motif', lib_tag, '_discovery_interactions.mvf')

for (d in c(output_dir, work_dir, input_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(PATH = paste(homer_bin, Sys.getenv("PATH"), sep = ":"))

# ------------------------------------------------------------------------------------------------------------------------------------------
# LOAD MODULES, PREDICTED TFs, AND THE BACKGROUND GENE UNIVERSE
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

if (!is.null(only_modules)) clusters <- clusters[as.character(clusters$Module) %in% only_modules, ]

module2genes <- setNames(clusters$Genes, as.character(clusters$Module))
module2tfs   <- setNames(tf_regs$TFs, as.character(tf_regs$Module))

universe <- unique(fread(expression_file, select = 1)[[1]])
universe <- universe[nzchar(universe)]
cat("Background universe:", length(universe), "expressed genes\n")

bg_all_file <- file.path(input_dir, "background_all_expressed.txt")
writeLines(sort(universe), bg_all_file)

# ------------------------------------------------------------------------------------------------------------------------------------------
# MAP HOMER KNOWN MOTIFS TO THE PREDICTED TF SYMBOLS
#
# HOMER labels motifs with historical/common protein names rather than gene symbols: TFAP2A's
# motif is "AP-2alpha(AP2)/Hela-AP2alpha-ChIP-Seq", TP63's is "p63(p53)/...", SNAI1's is
# "Snail1(Zf)/...". Matching on the motif name alone drops most predicted TFs and would mislabel
# them "no_motif_in_library" - i.e. report them as untestable when they were in fact tested.
#
# Resolution therefore runs outwards from the predicted TFs: each predicted symbol's own alias set
# is pulled from org.Hs.eg.db and HOMER's name tokens are matched against it. The reverse direction
# (alias -> symbol, e.g. limma::alias2SymbolTable) is unsafe here - "p63" resolves to RPE65, whose
# historical name it also is, silently attributing TP63's motif to an unrelated gene. An alias
# claimed by more than one predicted TF is dropped rather than guessed at.
# ------------------------------------------------------------------------------------------------------------------------------------------

norm_sym <- function(x) gsub("[^A-Z0-9]", "", toupper(x))

predicted_tfs <- unique(unlist(module2tfs[names(module2genes)]))
predicted_tfs <- predicted_tfs[nzchar(predicted_tfs)]

alias_tbl <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = predicted_tfs,
                                                    keytype = "SYMBOL", columns = "ALIAS"))
alias_dt <- unique(rbind(
  data.table(SYMBOL = predicted_tfs, ALIAS = predicted_tfs),
  as.data.table(alias_tbl)[!is.na(ALIAS)]
))
alias_dt[, KEY := norm_sym(ALIAS)]
alias_dt <- alias_dt[nzchar(KEY)]
ambig <- alias_dt[, .(n = uniqueN(SYMBOL)), by = KEY][n > 1, KEY]
if (length(ambig) > 0) cat("Dropping", length(ambig), "alias(es) claimed by multiple predicted TFs\n")
alias_dt <- alias_dt[!KEY %in% ambig]
alias2tf <- setNames(alias_dt$SYMBOL, alias_dt$KEY)

# Candidate symbols for one motif: the motif's own name plus the halves of dimer names
# ("FOXA1:AR"), and - from the source experiment field - ONLY the token immediately preceding
# "ChIP", which is where HOMER names the profiled factor ("MCF7-TFAP2C-ChIP-Seq" -> TFAP2C).
#
# Harvesting every source token instead is actively unsafe: it matched the predicted TF IRF6 to
# "NFkB-p65-Rel(RHD)/ThioMac-LPS-Expression" because IRF6's alias list contains "LPS" (Lip Pit
# Syndrome) and that source field's "LPS" means lipopolysaccharide, the treatment. HOMER has no
# IRF6 motif at all, so the pair must read no_motif_in_library rather than a confirmed prediction.
# Sources with no "ChIP" token (expression/treatment experiments) therefore contribute nothing.
motif_tokens <- function(name) {
  parts <- strsplit(name, "/")[[1]]
  common <- trimws(sub("\\(.*$", "", parts[1]))
  toks <- c(common, trimws(strsplit(common, "[:+]")[[1]]))
  if (length(parts) > 1) {
    src_toks <- trimws(strsplit(gsub("\\(.*?\\)", "", parts[2]), "[-_.]")[[1]])
    chip_at <- which(toupper(src_toks) == "CHIP")
    if (length(chip_at) > 0 && chip_at[1] > 1) {
      factor_tok <- src_toks[chip_at[1] - 1]
      # strip epitope-tag suffixes HOMER appends to the factor ("Sp5.Flag", "KLF14.GFP")
      toks <- c(toks, sub("\\.(FLAG|GFP|HA|V5|BIOTIN|MYC)$", "", factor_tok, ignore.case = TRUE))
    }
  }
  toks <- norm_sym(toks)
  setdiff(toks[nzchar(toks)],
          c("CHIP", "SEQ", "CHIPSEQ", "GSE", "HOMER", "EXPRESSION", "PROMOTERS",
            "GRO", "RNA", "DNASE", "FLAG", "GFP", "HA", "V5", "ENCODE"))
}

# "ISRE(IRF)/ThioMac-LPS-Expression(GSE23622)/Homer" -> "ISRE" (panel labels must stay short)
motif_label <- function(name) trimws(sub("\\(.*$", "", strsplit(name, "/")[[1]][1]))

motif_headers <- grep("^>", readLines(known_motifs), value = TRUE)
motif_index <- rbindlist(lapply(motif_headers, function(h) {
  f <- strsplit(sub("^>", "", h), "\t")[[1]]
  name <- f[2]
  hits <- unique(unname(alias2tf[intersect(motif_tokens(name), names(alias2tf))]))
  if (length(hits) == 0) return(NULL)
  data.table(Motif_Name = name, Symbol = hits)
}))

# rbindlist() of all-NULL is a column-less table: happens when none of the predicted TFs (e.g. a
# single-module test) has a HOMER-library motif
if (nrow(motif_index) == 0) motif_index <- data.table(Motif_Name = character(), Symbol = character())
motif_index[, `:=`(Relation = "direct", Motif_TF = Symbol, Motif_DB = "HOMER", Matrix_ID = NA_character_)]
cat("HOMER known motifs:", length(motif_headers), "| matched to",
    length(unique(motif_index$Symbol)), "of", length(predicted_tfs), "predicted TFs\n")

# Extended library: JASPAR/HOCOMOCO motifs are linked to predicted TFs through the mapping table
# written by build_extended_motif_library.py (never by alias guessing). Each row says whether the motif
# belongs to the TF itself ("direct") or to a same-family paralog standing in for it ("paralog_proxy").
if (extended) {
  cmap <- fread(custom_map)[Predicted_TF %in% predicted_tfs]
  motif_index <- rbind(motif_index,
                       cmap[, .(Motif_Name, Symbol = Predicted_TF, Relation, Motif_TF, Motif_DB, Matrix_ID)])
  cat("Extended library:", uniqueN(cmap$Motif_Name), "JASPAR/HOCOMOCO motifs |",
      uniqueN(cmap[Relation == "direct", Predicted_TF]), "TFs with a direct custom motif |",
      uniqueN(cmap[Relation == "paralog_proxy", Predicted_TF]), "TFs via paralog proxy\n")
}
tfs_with_motif <- unique(motif_index$Symbol)
cat("Predicted TFs testable with this library:", length(tfs_with_motif), "of", length(predicted_tfs), "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# HOMER WRAPPERS
# ------------------------------------------------------------------------------------------------------------------------------------------

run_findmotifs <- function(gene_file, bg_file, out_dir) {
  args <- c(gene_file, promoter_set, out_dir,
            "-start", tss_start, "-end", tss_end,
            "-nomotif", "-bg", bg_file, "-p", n_cpus)
  if (extended) args <- c(args, "-mknown", extended_lib)
  log <- file.path(out_dir, "homer.log")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  status <- system2(file.path(homer_bin, "findMotifs.pl"), args = as.character(args),
                    stdout = log, stderr = log)
  status == 0
}

# Per-gene attribution: rerun the same promoter set scanning ONLY the motif that passed tier 1.
# -find writes one row per motif occurrence; column 1 is the sequence (gene) identifier.
find_motif_genes <- function(gene_file, motif_file, out_dir) {
  hits_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(motif_file)), ".hits.txt"))
  args <- c(gene_file, promoter_set, file.path(out_dir, "find_tmp"),
            "-start", tss_start, "-end", tss_end,
            "-find", motif_file, "-p", n_cpus)
  status <- system2(file.path(homer_bin, "findMotifs.pl"), args = as.character(args),
                    stdout = hits_file, stderr = file.path(out_dir, "find.log"))
  if (status != 0 || !file.exists(hits_file)) return(character(0))
  hits <- tryCatch(fread(hits_file), error = function(e) NULL)
  if (is.null(hits) || nrow(hits) == 0) return(character(0))
  # Column 1 of -find output is the Entrez GeneID, NOT the symbol; the symbol is in "Name".
  # Intersecting column 1 with gene symbols would silently return nothing every time.
  sym_col <- intersect(c("Name", "name"), colnames(hits))
  if (length(sym_col) == 0) {
    cat("    Warning: no 'Name' column in", hits_file, "- cannot map motif hits to symbols\n")
    return(character(0))
  }
  unique(as.character(hits[[sym_col[1]]]))
}

# ------------------------------------------------------------------------------------------------------------------------------------------
# RUN: one HOMER job per module, then read off the rows for that module's predicted TFs
# ------------------------------------------------------------------------------------------------------------------------------------------

cat("\nTesting", length(module2genes), "modules | promoters TSS", tss_start, "to", tss_end,
    "| known motifs only | q <=", qval_threshold, "\n\n")

all_results <- list()
all_discovery <- list()
for (module in names(module2genes)) {
  tfs <- module2tfs[[module]]
  genes <- unique(module2genes[[module]])
  if (is.null(tfs) || length(tfs) == 0 || length(genes) < 3) next

  cat("Module", module, "-", length(genes), "genes | predicted TFs:", paste(tfs, collapse = ", "), "\n")

  gene_file <- file.path(input_dir, sprintf("module_%s.txt", module))
  bg_file   <- file.path(input_dir, sprintf("bg_minus_%s.txt", module))
  writeLines(sort(genes), gene_file)
  writeLines(sort(setdiff(universe, genes)), bg_file)

  out_dir <- file.path(work_dir, sprintf("module_%s", module))
  known_res_file <- file.path(out_dir, "knownResults.txt")
  if (!file.exists(known_res_file)) {
    ok <- run_findmotifs(gene_file, bg_file, out_dir)
    if (!ok || !file.exists(known_res_file)) {
      cat("  HOMER failed for module", module, "- see", file.path(out_dir, "homer.log"), "\n")
      next
    }
  } else {
    cat("  (reusing existing HOMER output)\n")
  }

  kr <- fread(known_res_file)
  setnames(kr, 1:9, c("Motif_Name", "Consensus", "P.value", "Log.P.value", "q.value",
                      "N_target_with_motif", "Pct_target_with_motif",
                      "N_bg_with_motif", "Pct_bg_with_motif"))
  kr[, Rank := .I]   # knownResults/known<Rank>.motif follows knownResults.txt row order

  rows <- lapply(tfs, function(tf) {
    # A paralog proxy stands in only for a TF with no motif of its own; it never competes with a direct one.
    tf_map <- motif_index[Symbol == tf]
    if (any(tf_map$Relation == "direct")) tf_map <- tf_map[Relation == "direct"]
    tf_motifs <- unique(tf_map$Motif_Name)
    hits <- kr[Motif_Name %in% tf_motifs]

    if (nrow(tf_map) == 0 || nrow(hits) == 0) {
      return(data.frame(Module = module, TF = tf, Status = "no_motif_in_library",
                        Relation = NA_character_, Motif_TF = NA_character_,
                        Motif_DB = NA_character_, Matrix_ID = NA_character_, N_motifs_for_TF = 0L,
                        Best_Motif = NA_character_, Consensus = NA_character_,
                        P.value = NA_real_, q.value = NA_real_,
                        Pct_target_with_motif = NA_character_, Pct_bg_with_motif = NA_character_,
                        N_target_genes_with_motif = 0L, Target_genes_with_motif = "",
                        stringsAsFactors = FALSE))
    }

    best <- hits[order(q.value, P.value)][1]
    info <- tf_map[Motif_Name == best$Motif_Name][1]
    significant <- !is.na(best$q.value) && best$q.value <= qval_threshold

    genes_with_motif <- character(0)
    if (significant) {
      motif_file <- file.path(out_dir, "knownResults", sprintf("known%d.motif", best$Rank))
      if (file.exists(motif_file)) {
        ids <- find_motif_genes(gene_file, motif_file, out_dir)
        genes_with_motif <- intersect(toupper(ids), toupper(genes))
        genes_with_motif <- genes[toupper(genes) %in% genes_with_motif]
      } else {
        cat("    Warning: motif file not found for", tf, ":", motif_file, "\n")
      }
    }

    data.frame(Module = module, TF = tf,
               Status = if (significant) "tested_significant" else "tested_not_significant",
               Relation = info$Relation, Motif_TF = info$Motif_TF, Motif_DB = info$Motif_DB,
               Matrix_ID = info$Matrix_ID, N_motifs_for_TF = length(tf_motifs),
               Best_Motif = best$Motif_Name, Consensus = best$Consensus,
               P.value = best$P.value, q.value = best$q.value,
               Pct_target_with_motif = as.character(best$Pct_target_with_motif),
               Pct_bg_with_motif = as.character(best$Pct_bg_with_motif),
               N_target_genes_with_motif = length(genes_with_motif),
               Target_genes_with_motif = paste(genes_with_motif, collapse = "|"),
               stringsAsFactors = FALSE)
  })

  res <- do.call(rbind, rows)
  all_results[[module]] <- res

  ## ---- discovery: significant motifs that are NOT attributable to a predicted TF ----
  ## Motif enrichment frequently identifies the actual DNA-binding driver of a module while the
  ## expression-based prediction names a partner/co-regulated factor (module 33: ISRE/IRF family
  ## at 65% vs 4% background, where the predictions were BATF2 and VDR). Those hits are recorded
  ## and plotted separately so they never get mistaken for a confirmation of the prediction.
  predicted_motifs <- motif_index[Symbol %in% tfs, unique(Motif_Name)]
  disc <- kr[!is.na(q.value) & q.value <= qval_threshold & !Motif_Name %in% predicted_motifs]
  disc <- disc[order(q.value, P.value)]
  if (nrow(disc) > max_discovery_motifs) disc <- disc[seq_len(max_discovery_motifs)]

  if (nrow(disc) > 0) {
    disc_rows <- lapply(seq_len(nrow(disc)), function(i) {
      d <- disc[i]
      motif_file <- file.path(out_dir, "knownResults", sprintf("known%d.motif", d$Rank))
      genes_with_motif <- character(0)
      if (file.exists(motif_file)) {
        ids <- find_motif_genes(gene_file, motif_file, out_dir)
        genes_with_motif <- genes[toupper(genes) %in% intersect(toupper(ids), toupper(genes))]
      }
      data.frame(Module = module, Motif_Label = motif_label(d$Motif_Name),
                 Motif_Name = d$Motif_Name, Consensus = d$Consensus,
                 P.value = d$P.value, q.value = d$q.value,
                 Pct_target_with_motif = as.character(d$Pct_target_with_motif),
                 Pct_bg_with_motif = as.character(d$Pct_bg_with_motif),
                 N_genes_with_motif = length(genes_with_motif),
                 Genes_with_motif = paste(genes_with_motif, collapse = "|"),
                 stringsAsFactors = FALSE)
    })
    all_discovery[[module]] <- do.call(rbind, disc_rows)
    cat("  -> discovery:", nrow(disc), "significant unpredicted motif(s):",
        paste(sapply(disc$Motif_Name, motif_label), collapse = ", "), "\n")
  }
  n_sig <- sum(res$Status == "tested_significant")
  if (n_sig > 0) cat("  ->", n_sig, "significant:", paste(res$TF[res$Status == "tested_significant"], collapse = ", "), "\n")
}

results_df <- do.call(rbind, all_results)
if (is.null(results_df)) {
  results_df <- data.frame(Module = character(), TF = character(), Status = character(),
                           Relation = character(), Motif_TF = character(), Motif_DB = character(),
                           Matrix_ID = character(), N_motifs_for_TF = integer(),
                           Best_Motif = character(), Consensus = character(),
                           P.value = numeric(), q.value = numeric(),
                           Pct_target_with_motif = character(), Pct_bg_with_motif = character(),
                           N_target_genes_with_motif = integer(), Target_genes_with_motif = character())
}

results_csv <- file.path(output_dir, paste0("HOMER_motif_enrichment_46_modules", lib_tag, ".csv"))
fwrite(results_df, results_csv)
cat("\nSaved enrichment table to:", results_csv, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# WRITE MODULEVIEWER-COMPATIBLE .mvf ANNOTATION FILE
# Only tier-1-significant module-TF pairs are drawn, and within those only the genes whose own
# promoter carries the motif (tier 2). "no_motif_in_library" and "tested_not_significant" stay in
# the CSV for transparency but are not plotted.
# ------------------------------------------------------------------------------------------------------------------------------------------

sig_df <- results_df[results_df$Status == "tested_significant" & nzchar(results_df$Target_genes_with_motif), ]

mvf_type <- paste0("HOMER_motif", lib_tag)
mvf_lines <- c(paste0("::TYPE=", mvf_type), paste0("::TITLE:", mvf_type), "::OBJECT=GENES", "::COLOR=TEAL")
if (nrow(sig_df) > 0) {
  # paralog-derived calls carry the paralog in the label so they stay distinguishable in the viewer
  sig_df$Label <- ifelse(sig_df$Relation == "paralog_proxy" & !is.na(sig_df$Relation),
                         sprintf("%s (via %s)", sig_df$TF, sig_df$Motif_TF), sig_df$TF)
  mvf_lines <- c(mvf_lines, apply(sig_df, 1, function(r) {
    paste(r[["Module"]], r[["Target_genes_with_motif"]], r[["Label"]], sep = "\t")
  }))
}
writeLines(mvf_lines, mvf_output_file)
cat("Saved ModuleViewer annotation file to:", mvf_output_file, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# DISCOVERY OUTPUTS (significant motifs with no predicted-TF match)
# ------------------------------------------------------------------------------------------------------------------------------------------

discovery_df <- do.call(rbind, all_discovery)
if (is.null(discovery_df)) {
  discovery_df <- data.frame(Module = character(), Motif_Label = character(), Motif_Name = character(),
                             Consensus = character(), P.value = numeric(), q.value = numeric(),
                             Pct_target_with_motif = character(), Pct_bg_with_motif = character(),
                             N_genes_with_motif = integer(), Genes_with_motif = character())
}
discovery_csv <- file.path(output_dir, paste0("HOMER_motif_discovery_46_modules", lib_tag, ".csv"))
fwrite(discovery_df, discovery_csv)
cat("Saved discovery table to:", discovery_csv, "\n")

disc_plot <- discovery_df[nzchar(discovery_df$Genes_with_motif), ]
disc_lines <- c(paste0("::TYPE=HOMER_motif_discovery", lib_tag), paste0("::TITLE:HOMER_motif_discovery", lib_tag),
                "::OBJECT=GENES", "::COLOR=MAROON")
if (nrow(disc_plot) > 0) {
  disc_lines <- c(disc_lines, apply(disc_plot, 1, function(r) {
    paste(r[["Module"]], r[["Genes_with_motif"]], r[["Motif_Label"]], sep = "\t")
  }))
}
writeLines(disc_lines, discovery_mvf_file)
cat("Saved discovery annotation file to:", discovery_mvf_file, "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# COVERAGE / RESULT SUMMARY
# ------------------------------------------------------------------------------------------------------------------------------------------

n_pairs <- nrow(results_df)
n_no_motif <- sum(results_df$Status == "no_motif_in_library")
n_tested <- n_pairs - n_no_motif
n_sig <- sum(results_df$Status == "tested_significant")

cat("\n=== HOMER known-motif enrichment complete ===\n")
cat("Module-TF pairs evaluated:", n_pairs, "\n")
cat("  - no motif in HOMER known library:", n_no_motif, "\n")
cat("  - tested:", n_tested, "\n")
cat("      - significant (q <=", qval_threshold, "):", n_sig, "\n")
cat("      - not significant:", n_tested - n_sig, "\n")
if (extended) {
  by_rel <- results_df[results_df$Status != "no_motif_in_library", ]
  cat("  tested by relation:", paste(names(table(by_rel$Relation)), table(by_rel$Relation), collapse = ", "), "\n")
  cat("  significant via paralog proxy:", sum(sig_df$Relation == "paralog_proxy", na.rm = TRUE), "\n")
}
cat("Modules with >=1 motif-validated predicted TF:", length(unique(sig_df$Module)),
    "out of", length(module2genes), "\n")
cat("Discovery: ", nrow(discovery_df), " significant unpredicted motif(s) across ",
    length(unique(discovery_df$Module)), " module(s)\n", sep = "")
