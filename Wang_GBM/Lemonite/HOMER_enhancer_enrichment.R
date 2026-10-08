#!/usr/bin/Rscript

############################################################################################################################################
#### HOMER known-motif enrichment on ENHANCERS (distal regulatory regions), per gene module.
####
#### Companion to HOMER_motif_enrichment.R, which runs the same two-tier logic on promoters
#### (TSS +/-1kb). This script asks the same questions of distal elements instead:
####   tier 1 (per module) - are the enhancers linked to a module's genes enriched for the binding
####                         motif of the TF LemonTree assigned to that module, relative to the
####                         enhancers of all other expressed genes?
####   tier 2 (per enhancer) - for pairs that pass tier 1, which enhancers actually carry the motif?
####                         Their linked genes populate the .mvf annotation panel.
####
#### WHY AN EXTERNAL ENHANCER CATALOGUE IS USED
#### The acetylome in this dataset (Acetylome_normalized.csv) cannot define enhancers. It is
#### mass-spectrometry acetylproteomics: 12,455 lysine sites keyed by protein residue
#### ("H1F0:NP_005309.1:K12"), with no genomic coordinates. Even its 165 histone sites are bulk
#### nuclear averages - MS reports how much H3/H2A acetylation exists, never at which locus - so it
#### cannot distinguish an acetylated enhancer from an acetylated gene body. Locus-resolved data
#### (H3K27ac ChIP-seq/CUT&RUN, ATAC-seq) would be required, and this dataset has none.
####
#### Enhancers and their gene links come from FANTOM5's own enhancer-TSS association file
#### (Andersson et al. 2014 Nature): enhancer/promoter CAGE expression is correlated across all
#### FANTOM5 libraries and every pair with FDR < 1e-5 is kept (66,942 associations, 11,607 genes).
#### This is the expression-validated link the analysis originally called for - NOT a proximity
#### heuristic. The file itself (enhancer.binf.ku.dk, where FANTOM5 originally hosted it) is no
#### longer served and no longer mirrored by RIKEN or Zenodo; it was recovered from a 2022-10-06
#### Wayback Machine snapshot. It is hg19; enhancer coordinates are lifted to hg38 with pyliftover
#### + UCSC's hg19ToHg38.over.chain (target-gene TSS coordinates are not needed downstream, so only
#### the enhancer interval is lifted). 66,936/66,942 associations (99.99%) lift cleanly.
####
#### Background: enhancers of every expressed gene, minus the module's own. Enhancer-vs-enhancer
#### comparison controls sequence composition far better than a random-genomic background would.
####
#### MOTIF LIBRARIES (same switch as HOMER_motif_enrichment.R)
####   homer    (default) HOMER's known.motifs only - ~75% of module-TF pairs have no motif there.
####   extended HOMER's known.motifs PLUS JASPAR 2024 / HOCOMOCO v12 matrices for predicted TFs HOMER
####            lacks, plus paralog-proxy motifs for four TFs with no motif anywhere (BATF2, NEUROD6,
####            TFAP2D, NFE2L3), from build_extended_motif_library.py. A paralog motif is used ONLY for a
####            TF with no motif of its own and is always flagged: Relation = "paralog_proxy" + Motif_TF
####            in the CSV, "TF (via PARALOG)" in the .mvf. q-values are Benjamini over the whole
####            extended library (slightly more conservative); p-values of HOMER's motifs are identical.
####   Usage: Rscript HOMER_enhancer_enrichment.R [homer|extended] [module,module,...]
############################################################################################################################################

suppressMessages({
  library(data.table)
  library(stringr)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------------------------------------------------------------------------------------------------------------------
# USER CONFIGURATION
# ------------------------------------------------------------------------------------------------------------------------------------------

base_dir   <- '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/'

cli_args      <- commandArgs(trailingOnly = TRUE)
motif_library <- if (length(cli_args) >= 1) cli_args[1] else "homer"
only_modules  <- if (length(cli_args) >= 2) strsplit(cli_args[2], ",")[[1]] else NULL   # smoke tests
stopifnot(motif_library %in% c("homer", "extended"))
extended <- motif_library == "extended"
lib_tag  <- paste0(if (extended) "_extended" else "", if (!is.null(only_modules)) "_subset" else "")

extended_lib <- paste0(base_dir, 'Enrichment/motif_library/extended_known.motifs')   # HOMER known + JASPAR/HOCOMOCO
custom_map   <- paste0(base_dir, 'Enrichment/motif_library/custom_motif_map.tsv')    # motif <-> predicted TF (+ Relation)
homer_bin  <- '/home/borisvdm/Bioinformatics_software/homer/bin'
genome     <- 'hg38'
known_motifs <- '/home/borisvdm/Bioinformatics_software/homer/data/knownTFs/vertebrates/known.motifs'   # HOMER's own library (alias-matched to TFs)

percentile     <- 2
MIN_ENHANCERS  <- 10        # skip modules with fewer linked enhancers (no power)
q_threshold    <- 0.05
n_cpus         <- 4
N_DISCOVERY    <- 5         # top unpredicted motifs reported per module

clusters_file         <- paste0(base_dir, 'ModuleViewer_files/clusters_list.txt')
tf_regulators_file    <- paste0(base_dir, 'ModuleViewer_files/Lovering.percentile', percentile, '_list.txt')
specific_modules_file <- paste0(base_dir, 'Networks/specific_modules.txt')
expression_file       <- paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression.txt')

output_dir  <- paste0(base_dir, 'Enrichment/HOMER_enhancer', lib_tag)
runs_dir    <- file.path(output_dir, 'runs')
# built once by scripts/build_fantom5_enhancer_links.py from the recovered FANTOM5
# enhancer-TSS association file (see header) - columns PeakID/Chr/Start/End/Gene/R/FDR,
# already restricted to FDR < 1e-5 and lifted to hg38.
links_file  <- paste0(base_dir, 'Enrichment/HOMER_enhancer/fantom5_enhancer_gene_links_hg38.tsv')   # shared by both libraries
mvf_valid   <- paste0(base_dir, 'ModuleViewer_files/HOMER_enhancer', lib_tag, '_interactions.mvf')
mvf_disc    <- paste0(base_dir, 'ModuleViewer_files/HOMER_enhancer', lib_tag, '_discovery_interactions.mvf')

dir.create(runs_dir, recursive = TRUE, showWarnings = FALSE)

# findMotifsGenome.pl invokes its helpers (homer2, findKnownMotifs.pl, ...) by bare name,
# so HOMER's bin must be on PATH for the child processes, not just addressed absolutely here.
Sys.setenv(PATH = paste(homer_bin, Sys.getenv("PATH"), sep = ":"))

if (!file.exists(links_file)) stop("Missing ", links_file, " - see header for how it is built")

# ------------------------------------------------------------------------------------------------------------------------------------------
# MODULES, PREDICTED TFs, EXPRESSED UNIVERSE
# ------------------------------------------------------------------------------------------------------------------------------------------

read_module_list <- function(file, value_col) {
  dt <- fread(file, header = FALSE, col.names = c("Module", value_col))
  dt[[value_col]] <- str_split(dt[[value_col]], "\\|")
  dt
}

clusters <- read_module_list(clusters_file, "Genes")
tf_regs  <- read_module_list(tf_regulators_file, "TFs")

if (file.exists(specific_modules_file)) {
  specific <- gsub("\\s+", "", readLines(specific_modules_file))
  clusters <- clusters[clusters$Module %in% specific, ]
  cat("Filtered to", nrow(clusters), "modules using specific_modules.txt\n")
}
if (!is.null(only_modules)) clusters <- clusters[as.character(clusters$Module) %in% only_modules, ]
module2genes <- setNames(clusters$Genes, as.character(clusters$Module))
module2tfs   <- setNames(tf_regs$TFs, as.character(tf_regs$Module))

expressed <- unique(fread(expression_file, select = 1)[[1]])
cat("Expressed-gene universe:", length(expressed), "\n")

# ------------------------------------------------------------------------------------------------------------------------------------------
# ENHANCER -> GENE ASSIGNMENT (FANTOM5 expression-correlation links, FDR < 1e-5)
# ------------------------------------------------------------------------------------------------------------------------------------------

enh <- fread(links_file, sep = "\t", showProgress = FALSE)
enh <- enh[nzchar(Gene) & Gene %in% expressed]
cat("FANTOM5 enhancer-gene associations landing on an expressed gene:", nrow(enh), "\n")

peak2gene <- setNames(enh$Gene, enh$PeakID)

write_peaks <- function(dt, path) {
  fwrite(dt[, .(PeakID, Chr, Start, End, Strand = "+")], path, sep = "\t", col.names = FALSE)
}

# ------------------------------------------------------------------------------------------------------------------------------------------
# MOTIF <-> PREDICTED TF RESOLUTION
#
# Identical to HOMER_motif_enrichment.R: HOMER names motifs by historical protein name, so each
# predicted symbol's alias set is pulled from org.Hs.eg.db and HOMER's tokens matched against it.
# Only the source token immediately preceding "ChIP" is used, because harvesting every token
# matched IRF6 to an NF-kB motif via IRF6's alias "LPS" colliding with the lipopolysaccharide
# treatment in "ThioMac-LPS-Expression". Aliases claimed by two predicted TFs are dropped.
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

motif_tokens <- function(name) {
  parts <- strsplit(name, "/")[[1]]
  common <- trimws(sub("\\(.*$", "", parts[1]))
  toks <- c(common, trimws(strsplit(common, "[:+]")[[1]]))
  if (length(parts) > 1) {
    src_toks <- trimws(strsplit(gsub("\\(.*?\\)", "", parts[2]), "[-_.]")[[1]])
    chip_at <- which(toupper(src_toks) == "CHIP")
    if (length(chip_at) > 0 && chip_at[1] > 1) {
      factor_tok <- src_toks[chip_at[1] - 1]
      toks <- c(toks, sub("\\.(FLAG|GFP|HA|V5|BIOTIN|MYC)$", "", factor_tok, ignore.case = TRUE))
    }
  }
  toks <- norm_sym(toks)
  setdiff(toks[nzchar(toks)],
          c("CHIP", "SEQ", "CHIPSEQ", "GSE", "HOMER", "EXPRESSION", "PROMOTERS",
            "GRO", "RNA", "DNASE", "FLAG", "GFP", "HA", "V5", "ENCODE"))
}

motif_label <- function(name) trimws(sub("\\(.*$", "", strsplit(name, "/")[[1]][1]))

motif_headers <- grep("^>", readLines(known_motifs), value = TRUE)
motif_index <- rbindlist(lapply(motif_headers, function(h) {
  f <- strsplit(sub("^>", "", h), "\t")[[1]]
  hits <- unique(unname(alias2tf[intersect(motif_tokens(f[2]), names(alias2tf))]))
  if (length(hits) == 0) return(NULL)
  data.table(Motif_Name = f[2], Symbol = hits)
}))
# rbindlist() of all-NULL is a column-less table: happens when none of the predicted TFs (e.g. a
# single-module test) has a HOMER-library motif
if (nrow(motif_index) == 0) motif_index <- data.table(Motif_Name = character(), Symbol = character())
motif_index[, `:=`(Relation = "direct", Motif_TF = Symbol, Motif_DB = "HOMER", Matrix_ID = NA_character_)]
cat("HOMER known motifs:", length(motif_headers), "| matched to",
    length(unique(motif_index$Symbol)), "of", length(predicted_tfs), "predicted TFs\n")

# Extended library: JASPAR/HOCOMOCO motifs are linked to predicted TFs through the mapping table written
# by build_extended_motif_library.py (never by alias guessing). Relation says whether the motif belongs to
# the TF itself ("direct") or to a same-family paralog standing in for it ("paralog_proxy").
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
# HOMER WRAPPERS (genome-based: regions are BED intervals, not gene lists)
# ------------------------------------------------------------------------------------------------------------------------------------------

run_findmotifs <- function(peak_file, bg_file, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  args <- c(peak_file, genome, out_dir, "-size", "given", "-nomotif",
            "-bg", bg_file, "-p", n_cpus)
  if (extended) args <- c(args, "-mknown", extended_lib)
  log <- file.path(out_dir, "homer.log")
  system2(file.path(homer_bin, "findMotifsGenome.pl"), args = as.character(args),
          stdout = log, stderr = log) == 0
}

# Per-enhancer attribution: rescan the module's own enhancers for ONE motif. -find writes one row
# per occurrence; column 1 is the peak ID, which maps back to a gene through peak2gene.
find_motif_genes <- function(peak_file, motif_file, out_dir, module_genes) {
  hits_file <- file.path(out_dir, paste0(tools::file_path_sans_ext(basename(motif_file)), ".hits.txt"))
  args <- c(peak_file, genome, file.path(out_dir, "find_tmp"), "-size", "given",
            "-find", motif_file, "-p", n_cpus)
  status <- system2(file.path(homer_bin, "findMotifsGenome.pl"), args = as.character(args),
                    stdout = hits_file, stderr = file.path(out_dir, "find.log"))
  if (status != 0 || !file.exists(hits_file)) return(character(0))
  hits <- tryCatch(fread(hits_file), error = function(e) NULL)
  if (is.null(hits) || nrow(hits) == 0) return(character(0))
  genes <- unique(unname(peak2gene[as.character(hits[[1]])]))
  intersect(genes[!is.na(genes)], module_genes)
}

parse_known <- function(out_dir) {
  f <- file.path(out_dir, "knownResults.txt")
  if (!file.exists(f)) return(NULL)
  kr <- tryCatch(fread(f), error = function(e) NULL)
  if (is.null(kr) || nrow(kr) == 0) return(NULL)
  setnames(kr, 1, "Motif_Name")
  qcol <- grep("^q-value", colnames(kr), value = TRUE)[1]
  pcol <- grep("^P-value$", colnames(kr), value = TRUE)[1]
  tcol <- grep("^% of Target", colnames(kr), value = TRUE)[1]
  bcol <- grep("^% of Background", colnames(kr), value = TRUE)[1]
  kr[, .(Motif_Name, P.value = get(pcol), q.value = get(qcol),
         Pct_target = get(tcol), Pct_bg = get(bcol))]
}

# ------------------------------------------------------------------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------------------------------------------------------------------

cat("\nTesting", length(module2genes), "modules on FANTOM5 enhancers (expression-correlation",
    "gene links, FDR < 1e-5), background = enhancers of other expressed genes\n\n")

validation <- list(); discovery <- list()

for (module in names(module2genes)) {
  genes <- unique(module2genes[[module]])
  tfs   <- module2tfs[[module]]
  mod_enh <- enh[Gene %in% genes]
  if (nrow(mod_enh) < MIN_ENHANCERS) {
    cat("Module", module, "- only", nrow(mod_enh), "enhancers, skipping\n"); next
  }
  bg_enh <- enh[!Gene %in% genes]

  out_dir <- file.path(runs_dir, paste0("module_", module))
  peak_file <- file.path(out_dir, "target.peaks")
  bg_file   <- file.path(out_dir, "background.peaks")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write_peaks(mod_enh, peak_file); write_peaks(bg_enh, bg_file)

  cat("Module", module, "-", nrow(mod_enh), "enhancers from",
      uniqueN(mod_enh$Gene), "genes | predicted TFs:", paste(tfs, collapse = ", "), "\n")

  if (!file.exists(file.path(out_dir, "knownResults.txt"))) {
    if (!run_findmotifs(peak_file, bg_file, out_dir)) {
      cat("  HOMER failed for module", module, "\n"); next
    }
  }
  kr <- parse_known(out_dir)
  if (is.null(kr)) { cat("  no knownResults for module", module, "\n"); next }

  # ---- tier 1a: validation of the module's own predicted TFs ----
  for (tf in tfs) {
    # A paralog proxy stands in only for a TF with no motif of its own; it never competes with a direct one.
    tf_map <- motif_index[Symbol == tf]
    if (any(tf_map$Relation == "direct")) tf_map <- tf_map[Relation == "direct"]
    tf_motifs <- unique(tf_map$Motif_Name)
    if (length(tf_motifs) == 0) {
      validation[[length(validation) + 1]] <- data.table(
        Module = module, TF = tf, Status = "no_motif_in_library",
        Relation = NA_character_, Motif_TF = NA_character_, Motif_DB = NA_character_,
        Matrix_ID = NA_character_, N_motifs_for_TF = 0L,
        Best_Motif = NA_character_, P.value = NA_real_, q.value = NA_real_,
        Pct_target = NA_character_, Pct_bg = NA_character_,
        N_genes = 0L, Genes = NA_character_)
      next
    }
    sub <- kr[Motif_Name %in% tf_motifs][order(q.value)]
    if (nrow(sub) == 0) next
    top <- sub[1]
    info <- tf_map[Motif_Name == top$Motif_Name][1]
    sig <- !is.na(top$q.value) && top$q.value <= q_threshold
    hit_genes <- character(0)
    if (sig) {
      mf <- file.path(out_dir, "knownResults", paste0("known", which(kr$Motif_Name == top$Motif_Name)[1], ".motif"))
      if (file.exists(mf)) hit_genes <- find_motif_genes(peak_file, mf, out_dir, genes)
    }
    validation[[length(validation) + 1]] <- data.table(
      Module = module, TF = tf,
      Status = if (sig) "tested_significant" else "tested_not_significant",
      Relation = info$Relation, Motif_TF = info$Motif_TF, Motif_DB = info$Motif_DB,
      Matrix_ID = info$Matrix_ID, N_motifs_for_TF = length(tf_motifs),
      Best_Motif = top$Motif_Name, P.value = top$P.value, q.value = top$q.value,
      Pct_target = top$Pct_target, Pct_bg = top$Pct_bg,
      N_genes = length(hit_genes),
      Genes = if (length(hit_genes)) paste(hit_genes, collapse = "|") else NA_character_)
    if (sig) cat("  -> VALIDATED:", tf, "via", motif_label(top$Motif_Name),
                 "(q =", top$q.value, ",", length(hit_genes), "genes)\n")
  }

  # ---- tier 1b: discovery - significant motifs NOT matching any predicted TF of this module ----
  predicted_motifs <- motif_index[Symbol %in% tfs, Motif_Name]
  disc <- kr[!Motif_Name %in% predicted_motifs & !is.na(q.value) & q.value <= q_threshold][order(q.value)]
  if (nrow(disc) > 0) {
    disc <- head(disc, N_DISCOVERY)
    for (i in seq_len(nrow(disc))) {
      mf <- file.path(out_dir, "knownResults", paste0("known", which(kr$Motif_Name == disc$Motif_Name[i])[1], ".motif"))
      hit_genes <- if (file.exists(mf)) find_motif_genes(peak_file, mf, out_dir, genes) else character(0)
      discovery[[length(discovery) + 1]] <- data.table(
        Module = module, Motif_Name = disc$Motif_Name[i], Motif_Label = motif_label(disc$Motif_Name[i]),
        P.value = disc$P.value[i], q.value = disc$q.value[i],
        Pct_target = disc$Pct_target[i], Pct_bg = disc$Pct_bg[i],
        N_genes = length(hit_genes),
        Genes = if (length(hit_genes)) paste(hit_genes, collapse = "|") else NA_character_)
    }
    cat("  -> discovery:", nrow(disc), "significant unpredicted motif(s):",
        paste(sapply(disc$Motif_Name, motif_label), collapse = ", "), "\n")
  }
}

# ------------------------------------------------------------------------------------------------------------------------------------------
# OUTPUT
# ------------------------------------------------------------------------------------------------------------------------------------------

val_df  <- if (length(validation)) rbindlist(validation) else data.table()
disc_df <- if (length(discovery))  rbindlist(discovery)  else data.table()

fwrite(val_df,  file.path(output_dir, paste0("HOMER_enhancer_enrichment_modules", lib_tag, ".csv")))
fwrite(disc_df, file.path(output_dir, paste0("HOMER_enhancer_discovery_modules", lib_tag, ".csv")))

write_mvf <- function(df, path, type, colour, label_col) {
  lines <- c(paste0("::TYPE=", type), paste0("::TITLE:", type), "::OBJECT=GENES",
             paste0("::COLOR=", colour))
  if (nrow(df) > 0) {
    keep <- df[!is.na(Genes) & nzchar(Genes)]
    if (nrow(keep) > 0) lines <- c(lines, apply(keep, 1, function(r)
      paste(r[["Module"]], r[["Genes"]], r[[label_col]], sep = "\t")))
  }
  writeLines(lines, path)
}
val_sig <- if (nrow(val_df)) val_df[Status == "tested_significant"] else val_df
# paralog-derived calls carry the paralog in the label so they stay distinguishable in the viewer
if (nrow(val_sig)) val_sig[, Label := ifelse(Relation == "paralog_proxy" & !is.na(Relation),
                                             sprintf("%s (via %s)", TF, Motif_TF), TF)]
write_mvf(val_sig, mvf_valid, paste0("HOMER_enhancer", lib_tag), "DARKCYAN", "Label")
write_mvf(disc_df, mvf_disc,  paste0("HOMER_enhancer_discovery", lib_tag), "CORAL", "Motif_Label")

cat("\n=== HOMER enhancer motif enrichment complete ===\n")
if (nrow(val_df)) {
  cat("Module-TF pairs evaluated:", nrow(val_df), "\n")
  for (s in c("no_motif_in_library", "tested_not_significant", "tested_significant"))
    cat("  ", s, ":", sum(val_df$Status == s), "\n")
  if (extended) {
    tested <- val_df[Status != "no_motif_in_library"]
    cat("  tested by relation:", paste(names(table(tested$Relation)), table(tested$Relation), collapse = ", "), "\n")
    cat("  significant via paralog proxy:", sum(val_df$Status == "tested_significant" & val_df$Relation == "paralog_proxy"), "\n")
  }
}
cat("Discovery:", nrow(disc_df), "significant unpredicted motif(s) across",
    if (nrow(disc_df)) uniqueN(disc_df$Module) else 0, "module(s)\n")
cat("Outputs in:", output_dir, "\n")
