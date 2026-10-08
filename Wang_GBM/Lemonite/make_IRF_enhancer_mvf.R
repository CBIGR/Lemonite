#!/usr/bin/Rscript

############################################################################################################################################
#### IRF-enhancer highlight track for the ModuleViewer heatmaps.
####
#### Reads the extended-library results of HOMER_enhancer_enrichment.R and keeps the module-TF pairs where
####   - the predicted TF is an IRF factor (IRF1..IRF9) and its motif was significantly enriched (q <= 0.05)
####     in the module's enhancers ("tested_significant"), and
####   - at least one enhancer linked to a module gene carries that motif.
#### Output: one .mvf row per such module, listing the GENES whose linked enhancer(s) carry the IRF motif.
#### The heatmap scripts draw it as a right-hand column next to the expression heatmap, the same way
#### metabolite-gene interactions from the PKN are drawn.
####
#### "Linked" = FANTOM5 expression-correlation enhancer-TSS association (FDR < 1e-5), i.e. a predicted, not
#### experimentally validated, enhancer-gene link; "carries the motif" = a PWM match above HOMER's detection
#### threshold, not evidence of IRF binding. Motif matches are for the IRF6 (or other IRF) matrix, so they
#### also reflect the near-identical binding preference of the whole IRF family.
####
#### Usage: Rscript make_IRF_enhancer_mvf.R [results.csv]
####   default results.csv: the full extended run if present, otherwise the single-module test run.
############################################################################################################################################

suppressMessages(library(data.table))

base_dir <- '/home/borisvdm/Documents/PhD/Lemonite/Wang_GBM/results/LemonTree/transcriptomics_clustering/'
full_csv   <- paste0(base_dir, 'Enrichment/HOMER_enhancer_extended/HOMER_enhancer_enrichment_modules_extended.csv')
subset_csv <- paste0(base_dir, 'Enrichment/HOMER_enhancer_extended_subset/HOMER_enhancer_enrichment_modules_extended_subset.csv')

args <- commandArgs(trailingOnly = TRUE)
in_csv <- if (length(args) >= 1) args[1] else if (file.exists(full_csv)) full_csv else subset_csv
out_mvf <- paste0(base_dir, 'ModuleViewer_files/HOMER_enhancer_IRF_interactions.mvf')

res <- fread(in_csv)
hits <- res[Status == "tested_significant" & grepl("^IRF[0-9]+$", TF) & !is.na(Genes) & nzchar(Genes)]

lines <- c("::TYPE=HOMER_enhancer_IRF", "::TITLE:HOMER_enhancer_IRF", "::OBJECT=GENES", "::COLOR=MAGENTA")
if (nrow(hits) > 0) {
  lines <- c(lines, hits[, paste(Module, Genes, TF, sep = "\t"), by = seq_len(nrow(hits))]$V1)
}
writeLines(lines, out_mvf)

cat("Input:  ", in_csv, "\n")
cat("Output: ", out_mvf, "\n")
if (nrow(hits) > 0) {
  for (i in seq_len(nrow(hits)))
    cat(sprintf("  module %s | %s | %s (q = %.2g) | %d genes\n", hits$Module[i], hits$TF[i],
                hits$Best_Motif[i], hits$q.value[i], hits$N_genes[i]))
} else cat("  no IRF-TF enhancer hits\n")
