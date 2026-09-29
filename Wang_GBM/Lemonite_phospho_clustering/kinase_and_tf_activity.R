#!/usr/bin/env Rscript
# Kinase activity (decoupleR + OmniPath KSN) and TF activity (decoupleR + CollecTRI) for the
# phospho-clustering variant. Mirrors nextflow/scripts/Preprocessing_TFA_Proteomics.R (the
# decouple() + run_consensus() pattern) but adds KINASE activity on the phospho matrix.
#
# Kinase activity: net = OmniPath KSN (source=kinase, target=GENE_residue). We re-key the
# phospho matrix rows GENE:RefSeq:Res[;Res2] -> GENE_Res (expanding multi-residue ids) so they
# match KSN targets, then decouple(mat_phospho, net=KSN).
# TF activity: net = CollecTRI, run on the same phospho matrix's parent PROTEIN layer is not
# available here, so TFA is run on the phospho matrix as a footprint proxy (documented caveat).
#
# Outputs (into --out-dir): kinase_activity.txt, TF_activity.txt  (feature x sample), plus
# LemonTree regulator files KinaseActivity.txt / kinaseactivity.txt for the regulator layer.
suppressMessages({library(decoupleR); library(data.table); library(dplyr); library(tidyr); library(tibble)})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default=NULL){ i<-which(args==flag); if(length(i)) args[i+1] else default }
expr_file <- getarg("--expr")        # LemonPreprocessed_expression.txt (phospho, symbol,id,samples)
out_dir   <- getarg("--out-dir")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading phospho matrix:", expr_file, "\n")
df <- fread(expr_file, data.table = FALSE)
rn <- df[[1]]                        # phosphosite id GENE:RefSeq:Res[;Res2]
mat <- as.matrix(df[, -(1:2), drop = FALSE]); storage.mode(mat) <- "double"
rownames(mat) <- rn

# ---- re-key rows to GENE_residue (expand multi-residue) for KSN matching ------------------
# GENE:RefSeq:S364;S368 -> c(GENE_S364, GENE_S368). Keep a long map row_site -> matrix row idx.
gene_of <- sub(":.*$", "", rn)
res_of  <- sub("^[^:]*:[^:]*:", "", rn)          # residue field (may be 'S364;S368' or 'S364.1')
res_of  <- sub("\\..*$", "", res_of)             # drop trailing .1 duplicates
keys <- list(); kidx <- integer(0)
for (i in seq_along(rn)) {
  for (res in strsplit(res_of[i], ";", fixed = TRUE)[[1]]) {
    keys[[length(keys)+1]] <- paste0(gene_of[i], "_", res); kidx <- c(kidx, i)
  }
}
site_key <- unlist(keys)
# Aggregate duplicate GENE_residue keys (different phospho ids can map to the same site key)
# to one row by mean, so decoupleR gets unique, matchable feature names.
mat_long <- mat[kidx, , drop = FALSE]
sums <- rowsum(mat_long, site_key)
cnts <- as.vector(table(site_key)[rownames(sums)])
mat_ksn <- sums / cnts

run_ulm_activity <- function(m, net, label, source_col="source", target_col="target"){
  # Univariate linear model (run_ulm) is the standard decoupleR method for kinase/TF activity
  # from footprints; unlike the full decouple() consensus (which includes mlm) it does NOT
  # choke on colinear regulators (common in KSNs where kinases share substrates).
  cat("Running decoupleR run_ulm for", label, "(", nrow(net), "edges )...\n")
  res <- tryCatch(run_ulm(mat = m, net = net, .source = source_col, .target = target_col,
                          minsize = 5),
                  error=function(e){cat("[WARN]",label,"run_ulm failed:",conditionMessage(e),"\n"); NULL})
  if (is.null(res)) return(NULL)
  wide <- res %>% dplyr::filter(statistic=="ulm") %>%
    pivot_wider(id_cols="condition", names_from="source", values_from="score") %>%
    column_to_rownames("condition") %>% t() %>% as.data.frame()
  cat("[OK]", label, ":", nrow(wide), "features x", ncol(wide), "samples\n")
  wide
}

# ---- KINASE activity: OmniPath KSN ---------------------------------------------------------
ksn <- tryCatch(decoupleR::get_ksn_omnipath(),
                error=function(e){cat("[WARN] KSN download failed:",conditionMessage(e),"\n"); NULL})
if (!is.null(ksn)) {
  kin <- run_ulm_activity(mat_ksn, as.data.frame(ksn), "kinase-activity")
  if (!is.null(kin)) {
    write.table(kin, file.path(out_dir, "kinase_activity.txt"), sep="\t", quote=FALSE, col.names=NA)
    # LemonTree regulator layer: rows = kinases, cols = samples (z-score per kinase)
    kz <- as.data.frame(t(scale(t(as.matrix(kin)))))
    reg <- data.frame(Gene_symbol=rownames(kz), Protein_id=rownames(kz), kz, check.names=FALSE)
    write.table(reg, file.path(out_dir, "LemonPreprocessed_kinaseactivity.txt"), sep="\t", quote=FALSE, row.names=FALSE)
    writeLines(rownames(kz), file.path(out_dir, "kinaseactivity.txt"))
    cat("[OK] wrote kinase_activity.txt + regulator files (", nrow(kz), "kinases )\n")
  }
}

# ---- TF activity: CollecTRI (footprint on phospho matrix; caveat: not expression) ----------
ct <- tryCatch(fread("/home/borisvdm/repo/LemonIte/nextflow/PKN/CollecTRI_network.txt"),
               error=function(e) NULL)
if (!is.null(ct)) {
  if(!all(c("source","target") %in% colnames(ct))) colnames(ct)[1:2] <- c("source","target")
  # TF targets are genes; collapse phospho matrix to gene level (max abs per gene) for TFA proxy
  gmat <- rowsum(abs(mat), gene_of); gmat <- gmat[rowSums(is.finite(gmat))>0,,drop=FALSE]
  tfa <- run_ulm_activity(gmat, as.data.frame(ct), "TF-activity")
  if (!is.null(tfa)) write.table(tfa, file.path(out_dir,"TF_activity.txt"), sep="\t", quote=FALSE, col.names=NA)
}
cat("Done.\n")
