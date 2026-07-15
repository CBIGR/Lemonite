#!/usr/bin/Rscript

setwd('/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/MOFA')

# Load name mapping for display of original metabolite names
if (!exists('name_lookup')) {
  name_mapping_file <- '../../data/name_mapping.tsv'
  if (file.exists(name_mapping_file)) {
    name_mapping_df <- read.table(name_mapping_file, sep='\t', header=TRUE, stringsAsFactors=FALSE)
    name_lookup <- setNames(name_mapping_df$original, name_mapping_df$cleaned)
    message(paste("Loaded name mapping:", length(name_lookup), "entries"))
  } else {
    name_lookup <- NULL
    warning("name_mapping.tsv not found")
  }
}

restore_names <- function(nms) {
  if (is.null(name_lookup)) return(gsub('_', ' ', nms))
  # First try direct lookup
  restored <- name_lookup[nms]
  # For names not found, try with colons replaced by underscores
  # (name_mapping keys use C16_1_CE but model features use C16:1_CE)
  not_found <- is.na(restored)
  if (any(not_found)) {
    normalized <- gsub(':', '_', nms[not_found])
    restored[not_found] <- name_lookup[normalized]
  }
  # Final fallback: replace underscores with spaces
  ifelse(is.na(restored), gsub('_', ' ', nms), restored)
}

# BiocManager::install("MOFA2")
library(DESeq2)
library(data.table)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(IMIFA)
library(MOFA2)
library(MultiAssayExperiment)
library(reticulate)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(tibble)
library(grid)

#use_python("/home/borisvdm/Software/Miniconda3/bin/python")

###########################################################################################
#### Read RNAseq data, filter for protein coding genes
###########################################################################################

RNAseq <- fread("../../data/GSE111889_host_tx_counts.tsv", header=TRUE, data.table=TRUE)
# Select for protein coding genes, use ensembl IDs as rownames
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = 'https://jan2024.archive.ensembl.org')
# # Retrieve all Ensembl Gene IDs, HGNC Symbols, and Gene Biotypes for all possible genes
# all_genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
#                    mart = ensembl)
# id_ensembl <- all_genes[all_genes$hgnc_symbol %in% RNAseq$count, ]

id_ensembl <- fread('/home/borisvdm/Documents/PhD/Lemonite/ensembl_mapping_jan2024.txt')


RNAseq <- merge(RNAseq, id_ensembl, by.x = 'count', by.y='hgnc_symbol')
RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$count), ]) # 5 non-unique ensembl ID's
rownames(RNAseq_coding) <- RNAseq_coding$count
RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
# Some genes have a different ensembl ID but the same gene symbol, they also all have the same counts -> group them into 1 entry
RNAseq_coding <- as.data.frame(RNAseq_coding %>% group_by(count) %>% filter(row_number()==1)) # Just keep the first entry, they have the same counts anyway
rownames(RNAseq_coding) <- RNAseq_coding$count; genes <- RNAseq_coding$count
RNAseq_coding$count <- NULL; RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$ensembl_gene_id <- NULL


###########################################################################################
#### Normalization, log transformation and scaling + selection for highly variable genes using DESeq2
#### Included samples:
#### - Ulcerative colitis and nonIBD
#### - Colon or rectum
#### - Transcriptomics loosely coupled to metabolomics (no exact week (moment of sampling) match required)
###########################################################################################
metadata <- fread('../../data/metadata_transriptomics_metabolomics_looseCoupled.txt', data.table =FALSE); metadata$V1 <- NULL
samples_rna <- metadata$`External ID`
RNAseq <- RNAseq_coding[, colnames(RNAseq_coding) %in% samples_rna ]

# Remove problematic samples upfront (MSM719ME and HSM5FZB5)
RNAseq$MSM719ME <- NULL; metadata <- metadata[metadata$`External ID` != 'MSM719ME', ]
RNAseq$HSM5FZB5 <- NULL; metadata <- metadata[metadata$`External ID` != 'HSM5FZB5', ]

DESeq_groups <- metadata[,c('diagnosis', 'sex', "biopsy_location", 'Age at diagnosis', 'Participant ID')]
rownames(DESeq_groups) <- metadata$`External ID`
DESeq_groups[] <- lapply(DESeq_groups, factor)

ord <- match(colnames(RNAseq), rownames(DESeq_groups))
DESeq_groups <- DESeq_groups[ord, ]

dds <- DESeqDataSetFromMatrix(countData = (RNAseq+1), colData = DESeq_groups, design = ~ biopsy_location + diagnosis + sex) 
keep <- rowSums(counts(dds) >=10) >= 3 # Dispersion plot looks better with some prefiltering
dds <- dds[keep, ]
dds <- DESeq(dds)
normcnt <- as.data.frame(counts(dds, normalized=TRUE))
log_normcnt <- log(normcnt)

# Select for highly variable genes
M <- log_normcnt
vars <- apply(M,1,var)
pdf('./variance_histogram.pdf')
hist(vars, 1100,xlim=c(0,2))
abline(v=0.35, col='red', lwd=3)
dev.off()
length(names(vars[vars>=0.35])) # 4429

variable_genes <- M[names(vars[vars>0.35]), ]
RNA_preprocessed <- as.data.frame(t(scale(t(variable_genes))))
#setwd('./4429genes')

pdf('./Normalized_expression.pdf')
RNA_preprocessed %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()

###########################################################################################
#### Preprocessing of metabolomics data
###########################################################################################
HILIC_pos <- fread('../../data/HILIC_pos.txt', fill=TRUE)
HILIC_neg <- fread('../../data/HILIC_neg.txt', fill=TRUE)
C18_pos <- fread('../../data/C18_pos.txt', fill=TRUE)
C18_neg <- fread('../../data/C18_neg.txt', fill=TRUE)

# Data has already been normalized
# the first row under the header contains factor information, we do not need this line
HILIC_pos <- HILIC_pos[-c(1), ]
HILIC_neg <- HILIC_pos[-c(1), ]
C18_pos <- C18_pos[-c(1), ]
C18_neg <- C18_neg[-c(1), ]

# HILIC_pos$Samples <- paste0(HILIC_pos$Samples, '_Hpos')
# HILIC_neg$Samples <- paste0(HILIC_neg$Samples, '_Hneg')
# C18_pos$Samples <- paste0(C18_pos$Samples, '_Cpos')
# C18_neg$Samples <- paste0(C18_neg$Samples, '_Cneg')

# Concatenate data into one dataframe
abundancies <- as.data.frame(rbind(HILIC_pos, HILIC_neg, C18_pos, C18_neg))
abundancies <-as.data.frame(abundancies %>% group_by(Samples) %>% filter(row_number()==1))
abundancies <- abundancies %>% mutate_all(na_if, "") # Change empty entries to NA
abundancies[,-c(1)] <- as.data.frame(sapply(abundancies[,-c(1)], as.numeric)) # Change column types to numeric
abundancies$Samples <- str_replace_all(abundancies$Samples, ' ', '_')
rownames(abundancies) <- abundancies$Samples; abundancies$Samples <- NULL

# Select for samples also present in the transcriptomics dataset
# abundancies <- abundancies[, colnames(abundancies) %in% metadata$metabolomics_id] # 46 different patients included in the transcriptomics also have metabolomics information
# colnames(abundancies) <- metadata$`Participant ID`[match(colnames(abundancies), metadata$metabolomics_id)]
abundancies[is.na(abundancies)] <- 0
abundancies <- log(abundancies + 1)
abundancies <- as.data.frame(t(pareto_scale(t(abundancies)))) # In principle normalization has already been done...

pdf('./Normalized_metabolomics_all_samples.pdf') # sample MSM719ME does not follow normalized data pattern
abundancies %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()

# Now create a complete dataset for all 76 samples remaining in RNA_preprocessed
metabolomics <- data.frame(matrix(ncol=length(colnames(RNA_preprocessed)), nrow = nrow(abundancies)))
rownames(metabolomics) <- rownames(abundancies); colnames(metabolomics) <- colnames(RNA_preprocessed)

RNA_preprocessed$HSM5FZB5 <- NULL; metabolomics$HSM5FZB5 <- NULL # This sample is not present in the metabolomics data

for (i in 1:length(colnames(RNA_preprocessed))){
  # i <- 2
  ext_id <- colnames(RNA_preprocessed)[i]
  patient <- metadata[metadata$`External ID` == ext_id, 'Participant ID']
  meta_id <- metadata[metadata$`External ID` == ext_id, 'metabolomics_id']
   #For some patients there are multiple entries, but this is how we choose to do the analysis
  metabolomics[,ext_id] <- abundancies[, meta_id]
}

pdf('./Normalized_metabolomics_matched_samples.pdf') # sample MSM719ME does not follow normalized data pattern
metabolomics %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()

###########################################################################################
#### Run MOFA
###########################################################################################
# metadata <- fread("../data/hmp2_metadata.csv", data.table = FALSE)[, c(1:5, 41, 71)]
# metadata <- metadata[(metadata$`External ID` %in% colnames(metabolomics) & metadata$data_type == 'host_transcriptomics'), ]
rownames(metadata) <- metadata$`External ID` 
coldata <- metadata[,c(3,4,7)]; rownames(coldata) <- metadata$`External ID`

coldata$sex[coldata$sex == 'Male'] <- 1; coldata$sex[coldata$sex == 'Female'] <- 2
coldata$biopsy_location[coldata$biopsy_location == 'Colon'] <- 1; coldata$biopsy_location[coldata$biopsy_location == 'Rectum'] <- 2
# coldata[is.na(coldata)] <- 0 # Some samples have NA for 'Age at diagnosis', MOFA cannot deal with this
coldata$sex <- as.numeric(coldata$sex); coldata$biopsy_location <- as.numeric(coldata$biopsy_location)

data <- list('Transcriptomics' = RNA_preprocessed, 'Metabolomics' = metabolomics)
exp <- MultiAssayExperiment(data, colData = coldata)
#MOFAobject_group <- create_mofa(exp, groups = 'diagnosis')

MOFAobject <- create_mofa_from_MultiAssayExperiment(exp, extract_metadata = TRUE)
MOFAobject <- set_covariates(MOFAobject, c('sex', 'biopsy_location'))

plot_data_overview(MOFAobject)
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)
model_opts <- get_default_model_options(MOFAobject)
model_opts
train_opts <- get_default_training_options(MOFAobject)
train_opts$verbose <- TRUE; train_opts$seed <- 42
train_opts$convergence_mode <- "slow"

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
#MOFA.trained <- run_mofa(MOFAobject, './MOFA1.hdf5', use_basilisk = TRUE)


MOFA.trained <- load_model('/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/MOFA/MOFA1.hdf5')
###########################################################################################
#### Downstream analysis
###########################################################################################
model <- MOFA.trained
plot_data_overview(model)
head(model@cache$variance_explained$r2_per_factor[[1]])# Check how much variance is explained per view

write.table(model@cache$variance_explained$r2_per_factor[[1]], file = './variance_explained_per_view_and_factor.txt', sep = '\t', row.names = TRUE, quote = FALSE)

plot_variance_explained(model, x="view", y="factor") # Factor 5!
ggsave('variance_explained.png')




plot_factor(model, 
            factor = 1:5,
            color_by = "diagnosis",
            shape_by = "biopsy_location"
)
ggsave('Projection_latent_factors.png')

p <- plot_factor(model, 
                 factors = c(1,2,3,4,5,6),
                 color_by = "diagnosis",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit

p <- p + 
  scale_color_manual(values=c("nonIBD"="black", "UC"="red")) +
  scale_fill_manual(values=c("nonIBD"="black", "UC"="red"))

print(p)
ggsave('Projection_latent_factors_grouped.png')

plot_factors(model, 
             factors = 1:15,
             color_by = "diagnosis"
)
ggsave('Factor_distribution_all_factors.png')

plot_factors(model, 
             factors = 1:8,
             color_by = "diagnosis"
)
ggsave('Factor_distribution_factors1-8.png')

plot_factors(model, 
             factors = c(1,2,3),
             color_by = "diagnosis"
)
ggsave('Factor_distribution_factors1-2-3.png')


# Kruskall Wallis test for diff sample distributions
Z <- get_factors(model)[[1]]              # samples x factors matrix
meta <- samples_metadata(model)
group <- factor(meta$diagnosis)

kw_res <- t(sapply(seq_len(ncol(Z)), function(i) {
  f <- Z[, i]
  valid <- complete.cases(f, group)
  
  kw <- kruskal.test(f[valid], group[valid])
  
  # eta-squared effect size (based on the H statistic)
  n <- sum(valid)
  eta2 <- (kw$statistic - length(levels(group)) + 1) / (n - length(levels(group)))
  
  c(KW_stat = unname(kw$statistic), p_value = kw$p.value, eta_squared = unname(eta2))
}))
rownames(kw_res) <- colnames(Z)
kw_res <- as.data.frame(kw_res)

kw_res$p_value_adj <- p.adjust(kw_res$p_value, method = "BH")

write.csv(kw_res, file = './factor_multiomic_original_kruskalwallis.csv', row.names = TRUE)
kw_res

# Correlation between sample distributions and diagnosis as covariate
Z <- get_factors(model)[[1]]              # samples x factors matrix
meta <- samples_metadata(model)
covariate <- meta$diagnosis

# replicate MOFA2's coercion of the categorical covariate to numeric
if (!is.numeric(covariate)) {
  covariate <- as.numeric(as.factor(covariate))
}

cor_res <- t(sapply(seq_len(ncol(Z)), function(i) {
  valid <- complete.cases(Z[, i], covariate)
  ct <- cor.test(Z[valid, i], covariate[valid])
  c(r = unname(ct$estimate), p_value = ct$p.value)
}))
rownames(cor_res) <- colnames(Z)
cor_res <- as.data.frame(cor_res)

cor_res$p_value_adj    <- p.adjust(cor_res$p_value, method = "BH")
cor_res$neg_log10_padj <- -log10(cor_res$p_value_adj)

write.csv(cor_res, file = './factor_diagnosis_original_correlation.csv', row.names = TRUE)

# Inspecting feature weights tells us more about the most informative features
# Look at features with top loading on factor 2 (dominated by metabolomics data)

plot_weights(model,
             view = "Metabolomics",
             factor = 5,
             nfeatures = 15,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
) + { if (!is.null(name_lookup)) scale_y_discrete(labels = restore_names) }
ggsave('feature_loading_factor5.png')




for (i in 1:15){
  
  
  
  p <- plot_top_weights(model,
                        view = "Metabolomics",
                        factor = i,
                        nfeatures = 10)
  
  if (!is.null(name_lookup)) p <- p + scale_y_discrete(labels = restore_names)

  # Customize axis label sizes
  p <- p + theme(
    axis.title.x = element_text(size = 16),  # X-axis label size
    #axis.title.y = element_text(size = 16),  # Y-axis label size
    axis.text.x  = element_text(size = 16),  # X-axis tick labels
    axis.text.y  = element_text(size = 20)   # Y-axis tick labels
  )
  
  #print(p)
  
  
  ggsave(paste0('./feature_weights/factor',i,'_feature_weights_metabolomics.png'), p)
  
  #pdf(paste0('factor',i,'_feature_weights_expression.pdf'))
  #i <- 2
  p <- plot_top_weights(model,
                        view = "Transcriptomics",
                        factor = i,
                        nfeatures = 10)
  
  # Customize axis label sizes
  p <- p + theme(
    axis.title.x = element_text(size = 16),  # X-axis label size
    #axis.title.y = element_text(size = 16),  # Y-axis label size
    axis.text.x  = element_text(size = 16),  # X-axis tick labels
    axis.text.y  = element_text(size = 20)   # Y-axis tick labels
  )
  
  #print(p)
  
  ggsave(paste0('./feature_weights/factor',i,'_feature_weights_expression.png'), p)
}

#pdf('./correlation_with_covariates.pdf')
print(correlate_factors_with_covariates(model, covariates = c("sex", "biopsy_location"), plot="log_pval"))
#dev.off()
ggsave('correlation_with_covariates.png')

#pdf('./Factor_correlation.pdf')
plot_factor_cor(MOFA.trained)
#dev.off()
ggsave('Factor_correlation.png')





#### Heatmaps with top weight features across factors

library(ComplexHeatmap)
library(circlize)

views <- c("Transcriptomics", "Metabolomics")
factors <- c(1,2,3)
top_n <- 3

output_file_expression <- "./heatmaps/combined_top3_factors1-2-3_expression_horizontal.pdf"
output_file_weights    <- "./heatmaps/combined_top3_factors1-2-3_weights_horizontal.pdf"

dir.create("./heatmaps", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# COLOUR SCALE
# ============================================================
col_fun_weights <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "#F7F7F7", "#B35806")
)

# ============================================================
# GET WEIGHTS
# ============================================================
all_weights <- get_weights(
  model,
  views = "all",
  factors = "all"
)

# ============================================================
# BUILD COMBINED WEIGHT MATRIX
# (same orientation as your 6-factor figure)
# Rows = selected factors
# Columns = top features grouped by view
# ============================================================
combined_weight_matrix <- NULL
column_split <- c()
column_labels <- c()

for (view in views) {
  
  cat("Processing:", view, "\n")
  
  weight_mat <- all_weights[[view]]
  
  # Check factor availability
  if (max(factors) > ncol(weight_mat)) {
    stop(
      paste0(
        "Requested factor ",
        max(factors),
        " but view ",
        view,
        " only contains ",
        ncol(weight_mat),
        " factors."
      )
    )
  }
  
  # ----------------------------------------------------------
  # Extract top features for requested factors
  # ----------------------------------------------------------
  top_features_per_factor <- lapply(
    factors,
    function(f) {
      names(
        sort(
          abs(weight_mat[, f]),
          decreasing = TRUE
        )
      )[1:top_n]
    }
  )
  
  top_features <- unique(unlist(top_features_per_factor))
  
  # ----------------------------------------------------------
  # Weight matrix:
  # rows = factors
  # cols = features
  # ----------------------------------------------------------
  view_weights <- t(
    weight_mat[
      top_features,
      factors,
      drop = FALSE
    ]
  )
  
  # Truncate long feature names
  colnames(view_weights) <- substr(
    colnames(view_weights),
    1,
    35
  )
  
  # ----------------------------------------------------------
  # Combine
  # ----------------------------------------------------------
  combined_weight_matrix <- cbind(
    combined_weight_matrix,
    view_weights
  )
  
  column_split <- c(
    column_split,
    rep(view, ncol(view_weights))
  )
  
  column_labels <- c(
    column_labels,
    colnames(view_weights)
  )
}

# ============================================================
# ROW NAMES
# ============================================================
rownames(combined_weight_matrix) <- paste0(
  "Factor",
  factors
)

colnames(combined_weight_matrix) <- column_labels

# ============================================================
# COLUMN ANNOTATION
# ============================================================
ha <- HeatmapAnnotation(
  View = column_split,
  show_annotation_name = FALSE
)

colnames(combined_weight_matrix) <- ifelse(
  nchar(colnames(combined_weight_matrix)) > 20,
  paste0(substr(colnames(combined_weight_matrix), 1, 20), "..."),
  colnames(combined_weight_matrix)
)

ht <- Heatmap(
  combined_weight_matrix,
  
  name = "Weight",
  
  col = col_fun_weights,
  
  cluster_rows = FALSE,
  
  # cluster all features together
  cluster_columns = TRUE,
  
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_side = "left",
  
  row_names_gp = gpar(
    fontsize = 12,
    fontface = "bold"
  ),
  
  column_names_gp = gpar(
    fontsize = 14
  ),
  
  column_names_rot = 45,
  
  column_title = paste(
    "Top", top_n,
    "features from Factors",
    paste(factors, collapse = ", ")
  ),
  
  column_title_gp = gpar(
    fontsize = 14,
    fontface = "bold"
  ),
  
  heatmap_legend_param = list(
    title = "Weight",
    at = c(-1, -0.5, 0, 0.5, 1)
  )
)

# ============================================================
# SAVE PDF
# ============================================================
pdf(
  output_file_weights,
  width = 16,
  height = 6
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()



###########################################################################################################################
#### GSEA Analysis for all factors (Supplementary Information)
###########################################################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(MOFAcellulaR)
library(enrichplot)

# ------------------------ VISUALIZATION FUNCTION ----------------------------
# Use shared GSEA/EnrichR plotting utilities
source('/home/borisvdm/repo/gsea_and_enrichr_plotting_utils.R')
cat('Sourced GSEA/EnrichR plotting utilities from /home/borisvdm/repo/gsea_and_enrichr_plotting_utils.R\n')

# Create GSEA results directory
gsea_dir <- './GSEA_all_factors'
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Get number of factors
n_factors <- model@dimensions$K

# Run GSEA for all factors
for (factor_num in 1:n_factors) {
  
  #factor_num <- 2
  
  cat(sprintf("\n=== Processing Factor %d ===\n", factor_num))
  
  # Create factor-specific directory
  factor_dir <- paste0(gsea_dir, "/Factor_", factor_num)
  dir.create(factor_dir, showWarnings = FALSE)
  
  # Extract gene weights for this factor
  gene_weights <- get_geneweights(model = model, factor = paste0("Factor", factor_num))
  
  # Filter for Transcriptomics view only
  gene_weights <- gene_weights[gene_weights$ctype == "Transcriptomics", ]
  
  if (nrow(gene_weights) == 0) {
    cat(sprintf("  No expression weights for factor %d, skipping...\n", factor_num))
    next
  }
  
  # Sort by absolute value and take top genes
  gene_weights <- gene_weights[order(abs(gene_weights$value), decreasing = TRUE), ]
  gene_weights <- gene_weights[1:min(10000, nrow(gene_weights)), ]
  gene_weights <- gene_weights[order(gene_weights$value, decreasing = TRUE), ]
  
  # Create named vector
  weights <- as.numeric(gene_weights$value)
  names(weights) <- gene_weights$feature
  weights <- weights[!is.na(weights)]
  
  if (length(weights) == 0) {
    cat(sprintf("  No valid weights for factor %d, skipping...\n", factor_num))
    next
  }
  
  # Convert gene symbols to Entrez IDs
  tryCatch({
    converted_entrez <- bitr(
      names(weights),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
    weights_entrez <- weights[converted_entrez$SYMBOL]
    names(weights_entrez) <- converted_entrez$ENTREZID
  }, error = function(e) {
    cat(sprintf("  Warning: Could not convert gene symbols to Entrez IDs for factor %d\n", factor_num))
    weights_entrez <- NULL
  })
  
  # GO Enrichment Analysis (include ALL, BP, MF, CC)
  for (db in c('ALL', 'BP', 'MF', 'CC')) {
    tryCatch({
      gse_result <- gseGO(
        geneList = weights,
        ont = db,
        keyType = "SYMBOL",
        nPermSimple = 10000,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.10,
        verbose = FALSE,
        OrgDb = org.Hs.eg.db,
        pAdjustMethod = "BH"
      )
      
      if (nrow(as.data.frame(gse_result)) > 0) {
        
        write.table(as.data.frame(gse_result),
                    file = paste0(factor_dir, "/Factor", factor_num, "_gseGO_", db, "_results.txt"),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        safe_plot_save(gse_result, paste0(factor_dir, "/Factor", factor_num, "_gseGO_", db, ".png"), db_name = paste0("Factor ", factor_num, " GO ", db))
        cat(sprintf("  Saved GO %s enrichment\n", db))
      }
    }, error = function(e) {
      cat(sprintf("  No enrichment found for GO %s\n", db))
    })
  }
}

## Aggregate per-factor GSEA results into combined files (GO, KEGG, Reactome)
agg_go <- data.frame()
agg_kegg <- data.frame()
agg_reactome <- data.frame()

for (factor_num in 1:n_factors) {
  factor_dir <- file.path(gsea_dir, paste0("Factor_", factor_num))

  # GO files (ALL, BP, MF, CC)
  for (db in c('ALL', 'BP', 'MF', 'CC')) {
    file_go <- file.path(factor_dir, paste0("Factor", factor_num, "_gseGO_", db, "_results.txt"))
    if (file.exists(file_go)) {
      res <- tryCatch(read.table(file_go, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        res$Factor <- factor_num
        res$Database <- db
        agg_go <- rbind(agg_go, res)
      }
    }
  }

  # KEGG file
  file_kegg <- file.path(factor_dir, paste0("Factor", factor_num, "_KEGG_results.txt"))
  if (file.exists(file_kegg)) {
    res <- tryCatch(read.table(file_kegg, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      res$Factor <- factor_num
      res$Database <- 'KEGG'
      agg_kegg <- rbind(agg_kegg, res)
    }
  }

  # Reactome file
  file_react <- file.path(factor_dir, paste0("Factor", factor_num, "_Reactome_results.txt"))
  if (file.exists(file_react)) {
    res <- tryCatch(read.table(file_react, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      res$Factor <- factor_num
      res$Database <- 'Reactome'
      agg_reactome <- rbind(agg_reactome, res)
    }
  }
}

# Put factor and database columns first
if (nrow(agg_go) > 0) {
  agg_go <- agg_go %>%
    dplyr::select(Factor, Database, everything())
}

if (nrow(agg_kegg) > 0) {
  agg_kegg <- agg_kegg %>%
    dplyr::select(Factor, Database, everything())
}

if (nrow(agg_reactome) > 0) {
  agg_reactome <- agg_reactome %>%
    dplyr::select(Factor, Database, everything())
}

# Write aggregated outputs
if (nrow(agg_go) > 0) {
  write.table(agg_go, file = file.path(gsea_dir, 'Combined_GSEA_GO_BP_MF_CC_all_factors.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}
if (nrow(agg_kegg) > 0) {
  write.table(agg_kegg, file = file.path(gsea_dir, 'Combined_GSEA_KEGG_all_factors.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}
if (nrow(agg_reactome) > 0) {
  write.table(agg_reactome, file = file.path(gsea_dir, 'Combined_GSEA_Reactome_all_factors.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}

cat(sprintf("\n=== GSEA Analysis Complete ===\n"))
cat(sprintf("Results saved in: %s\n", gsea_dir))


###########################################################################################
#### Extract and save feature weights for all factors (Supplementary Information)
###########################################################################################

# Get all factor weights for both views
all_weights <- get_weights(model, views = "all", factors = "all", as.data.frame = TRUE)

# Separate by view
expression_weights <- all_weights %>% 
  filter(view == "Transcriptomics") %>%
  dplyr::select(feature, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value, names_prefix = "Factor_")

metabolomics_weights <- all_weights %>%
  filter(view == "Metabolomics") %>%
  dplyr::select(feature, factor, value) %>%
  pivot_wider(names_from = factor, values_from = value, names_prefix = "Factor_")


# Replace NA with 0
expression_weights[is.na(expression_weights)] <- 0
metabolomics_weights[is.na(metabolomics_weights)] <- 0
lipidomics_weights[is.na(lipidomics_weights)] <- 0

# Rename feature column for clarity
colnames(expression_weights)[1] <- "Gene"
colnames(metabolomics_weights)[1] <- "Metabolite"
metabolomics_weights$Metabolite <- restore_names(metabolomics_weights$Metabolite)


# Save to tab-delimited files
write.table(expression_weights, 
            file = "./feature_weights/Supplementary_MOFA_expression_feature_weights_all_factors.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(metabolomics_weights,
            file = "./feature_weights/Supplementary_MOFA_metabolomics_feature_weights_all_factors.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
cat("\n=== MOFA Feature Weights Summary ===\n")
cat(sprintf("Total expression features: %d\n", nrow(expression_weights)))
cat(sprintf("Total metabolomics features: %d\n", nrow(metabolomics_weights)))
cat(sprintf("Total lipidomics features: %d\n", nrow(lipidomics_weights)))
cat(sprintf("Total factors: %d\n", ncol(expression_weights) - 1))
cat(sprintf("\nFiles saved:\n"))
cat(sprintf("  - %s\n", "./Supplementary_MOFA_expression_feature_weights_all_factors.txt"))
cat(sprintf("  - %s\n", "./Supplementary_MOFA_metabolomics_feature_weights_all_factors.txt"))


#######################################################################################################################
#### Follow-up analysis with COSMOS (for all factors)
#######################################################################################################################

setwd('./To_COSMOS')

# reduce_solution_network() has no default cutoff. The official MOON tutorial uses 1.5;
# the function's own doc example uses 0.4. The right value is data-specific: it must sit
# below the achievable MOON score range of the nodes you want to keep. Metabolites are the
# observed downstream nodes and reach only modest scores (~1.4 ceiling), so the tutorial's
# 1.5 prunes ALL of them. A cutoff of 1.0 retains the strongest metabolites.
solution_cutoff_rec_to_TFmetab <- 1.0   # run 1: Receptors -> TF + Metabolites
solution_cutoff_TF_lig         <- 1.5   # run 2: TF -> Ligands (gene-only, unchanged)

library(cosmosR)
library(liana)
#BiocManager::install("saezlab/decoupleR")
library(decoupleR)
#data manipulations
library(dplyr)
library(reshape2)
library(GSEABase)
library(tidyr)

#plotting
library(ggplot2)
library(ggfortify)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)
library(ggvenn)

weights <- get_weights(model, views = "all", factors = "all")
# Get the number of factors
n_factors <- model@dimensions$K
factor_cols <- colnames(weights$Transcriptomics)   # "Factor1" ... "Factor15"

# Metabolite name -> HMDB mapping table (shared across factors and reused in the MOON loop below)
# cwd here is results/MOFA/To_COSMOS, so go up 3 levels to reach the project's data/ dir
metab_to_hmdb_moon <- fread('../../../data/name_map (2).csv')

# Preload the COSMOS meta-network once for reuse across factors
data("meta_network")
meta_network_metab <- meta_network[grepl("HMDB", meta_network$source) | grepl("HMDB", meta_network$target), ]
meta_network_metab$source <- gsub("Metab__", "", meta_network_metab$source)
meta_network_metab$target <- gsub("Metab__", "", meta_network_metab$target)
meta_network_metab$source <- gsub("_.*", "", meta_network_metab$source)
meta_network_metab$target <- gsub("_.*", "", meta_network_metab$target)
meta_network_metabs <- unique(c(meta_network_metab$source, meta_network_metab$target))
meta_network_metabs[is.na(meta_network_metabs)] <- "No_HMDB_unknown"
# Summary table for COSMOS metabolite counts
cosmos_summary <- data.frame(
  Factor = integer(0),
  total_metabolites = integer(0),
  overlap_metabolites = integer(0)
)
for (fact in 1:n_factors) {
  
  cat(sprintf("\n=== Processing COSMOS for Factor %d ===\n", fact))
  
  # Create factor-specific directory
  fact_dir <- paste0("./Factor_", fact, "_COSMOS")
  dir.create(fact_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get metabolite weights for current factor
  metab_inputs <- weights$Metabolomics[, fact]

  # Match metabolites to HMDB via the shared Lloyd-Price name map (loaded once below).
  # Model feature names use hyphens (e.g. "21-deoxycortisol") while name_map Query uses
  # underscores (e.g. "21_deoxycortisol"); normalize before matching.
  metab_query <- gsub("-", "_", names(metab_inputs))
  metab_to_hmdb <- metab_to_hmdb_moon[metab_to_hmdb_moon$Query %in% metab_query, ]

  # Replace missing HMDB IDs with a placeholder using original metabolite name
  metab_to_hmdb$HMDB[is.na(metab_to_hmdb$HMDB)] <- paste0("No_HMDB_", metab_to_hmdb$Query[is.na(metab_to_hmdb$HMDB)])

  # Replace names in metab_inputs with HMDB/placeholder
  common_metabs <- intersect(metab_query, metab_to_hmdb$Query)
  names(metab_inputs)[match(common_metabs, metab_query)] <- metab_to_hmdb$HMDB[match(common_metabs, metab_to_hmdb$Query)]
  
  # Filter metabolites with weight > 0.2
  metab_inputs_toCosmos <- metab_inputs[abs(metab_inputs) > 0.1]
  total_metabs <- length(metab_inputs_toCosmos)
  overlap_metabs <- sum(names(metab_inputs_toCosmos) %in% meta_network_metabs)
  cosmos_summary <- rbind(cosmos_summary, data.frame(Factor = fact, total_metabolites = total_metabs, overlap_metabolites = overlap_metabs))
  
  # Prepare Venn diagram list
  venn_list <- list(
    "COSMOS Meta-network" = meta_network_metabs,
    "Metabolites with weight > 0.2" = names(metab_inputs_toCosmos),
    "Metabolites in dataset" = names(metab_inputs)
  )
  
  # Colorblind-friendly palette
  cb_colors <- c("#E69F00", "#56B4E9", "#009E73")
  
  # Plot Venn diagram
  p <- ggvenn(
    venn_list,
    fill_color = cb_colors,
    stroke_size = 0.5,
    set_name_size = 5,
    text_size = 4,
    show_percentage = FALSE
  )
  
  # Save plot
  ggsave(file.path(fact_dir, 'venn_diagram.png'), plot = p, bg = 'white')
  
  cat(sprintf("Saved Venn diagram for Factor %d in %s\n", fact, fact_dir))
}

# Create a summary barplot with stacked bars showing overlap vs non-overlap
# Create a summary barplot with stacked bars showing overlap vs non-overlap
cosmos_summary <- cosmos_summary %>%
  mutate(
    Factor = factor(
      paste0("Factor ", Factor),
      levels = paste0("Factor ", seq_len(n_factors))
    ),
    non_overlap_metabolites = total_metabolites - overlap_metabolites
  )

# Reshape data for stacked barplot
cosmos_summary_long <- cosmos_summary %>%
  pivot_longer(
    cols = c(overlap_metabolites, non_overlap_metabolites),
    names_to = "category",
    values_to = "count"
  ) %>%
  mutate(
    category = factor(
      category,
      levels = c("non_overlap_metabolites", "overlap_metabolites"),
      labels = c("Not in COSMOS+", "In COSMOS+")
    )
  )


# ==========================
# Plot ALL factors
# ==========================

cosmos_barplot_all <- ggplot(
  cosmos_summary_long,
  aes(x = Factor, y = count, fill = category)
) +
  geom_col(
    color = "black",
    width = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "Not in COSMOS+" = "#FEE0D2",
      "In COSMOS+" = "#A50F15"
    ),
    name = "Metabolite\nCategory"
  ) +
  
  geom_text(
    data = cosmos_summary,
    aes(
      x = Factor,
      y = total_metabolites,
      label = total_metabolites
    ),
    vjust = -0.5,
    size = 3.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "Number of metabolites with |weight| > 0.2"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 13)
  )


# ==========================
# Plot only Factors 1, 2, 3
# ==========================

cosmos_summary_long_123 <- cosmos_summary_long %>%
  filter(Factor %in% c("Factor 1", "Factor 2", 'Factor 3'))

cosmos_summary_123 <- cosmos_summary %>%
  filter(Factor %in% c("Factor 1", "Factor 2", 'Factor 3'))


cosmos_barplot_123 <- ggplot(
  cosmos_summary_long_123,
  aes(x = Factor, y = count, fill = category)
) +
  geom_col(
    color = "black",
    width = 0.7
  ) +
  
  scale_fill_manual(
    values = c(
      "Not in COSMOS+" = "#FEE0D2",
      "In COSMOS+" = "#A50F15"
    ),
    name = "Metabolite\nCategory"
  ) +
  
  geom_text(
    data = cosmos_summary_123,
    aes(
      x = Factor,
      y = total_metabolites,
      label = total_metabolites
    ),
    vjust = -0.5,
    size = 4.5,
    inherit.aes = FALSE
  ) +
  
  labs(
    title = NULL,
    x = NULL,
    y = "Number of metabolites with |weight| > 0.2"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 13),
    axis.title.y = element_text(size = 15),
    axis.text.y = element_text(size = 13)
  )


# Display plots
cosmos_barplot_all
cosmos_barplot_123


# Save plots
ggsave(
  file.path("./", "COSMOS_metabolite_counts_all_factors.png"),
  plot = cosmos_barplot_all,
  width = 8,
  height = 5,
  bg = "white"
)

ggsave(
  file.path("./", "COSMOS_metabolite_counts_by_factor_1-2-3.png"),
  plot = cosmos_barplot_123,
  width = 8,
  height = 5,
  bg = "white"
)


###########################################################################################
#### COSMOS+ / MOON analysis (per factor) + EnrichR functional enrichment
#### Ported from Wang_GBM/COSMOS_MOON_analysis.R, adapted for the 2-view (Transcriptomics +
#### Metabolomics) Lloyd-Price MOFA model.
###########################################################################################

library(stringr)
library(igraph)
library(ggraph)
library(ggrepel)
library(enrichR)

if (!exists("plot_enrichr")) source('/home/borisvdm/repo/gsea_and_enrichr_plotting_utils.R')

# ---- Shared setup (runs once) ----

# Build or load ligand-receptor resource
if (!file.exists("ligrec_ressource.RData")) {
  ligrec_ressource <- distinct(liana::decomplexify(liana::select_resource("Consensus")[[1]]))
  save(ligrec_ressource, file = "ligrec_ressource.RData")
} else {
  load(file = "ligrec_ressource.RData")
}
ligrec_geneset <- cosmosR::format_LR_ressource(ligrec_ressource)

# TF-target network from CollecTRI via decoupleR
if (!file.exists("dorothea_df.RData")) {
  dorothea_df <- decoupleR::get_collectri()
  save(dorothea_df, file = "dorothea_df.RData")
} else {
  load(file = "dorothea_df.RData")
}

# LIANA LR dataframe for node-type annotation
ligrec_df <- distinct(ligrec_ressource[, c("source_genesymbol", "target_genesymbol")])
names(ligrec_df) <- c("Node1", "Node2")
ligrec_df$Node1  <- gsub("-", "_", ligrec_df$Node1)
ligrec_df$Node2  <- gsub("-", "_", ligrec_df$Node2)
ligrec_df$Sign   <- 1
ligrec_df$Weight <- 1

# Genes expressed in the MOFA model
expressed_genes <- setNames(rep(1, nrow(weights$Transcriptomics)), rownames(weights$Transcriptomics))

# RNA_all: genes x factors matrix for cross-factor activity heatmaps
RNA_all <- as.data.frame(weights$Transcriptomics)

# HMDB ID mapper
data("HMDB_mapper_vec")

# Dorothea PKN in source/interaction/target format (reused each iteration)
dorothea_PKN_base           <- dorothea_df[, c(1, 3, 2)]
names(dorothea_PKN_base)[2] <- "interaction"


# ---- Cross-factor LR and TF activity heatmaps (computed once) ----

ligrec_factors_res <- decoupleR::run_ulm(
  mat = as.matrix(RNA_all), network = ligrec_geneset, .source = set, .target = gene, minsize = 2)
ligrec_factors_df <- cosmosR::wide_ulm_res(ligrec_factors_res)
ligrec_factors_df <- ligrec_factors_df[
  order(apply(ligrec_factors_df, 1, function(x) max(abs(x))), decreasing = TRUE)[
    seq_len(min(25, nrow(ligrec_factors_df)))], ]
palette_lr <- cosmosR::make_heatmap_color_palette(ligrec_factors_df)
pheatmap(ligrec_factors_df, show_rownames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
         color = palette_lr, angle_col = 315,
         filename = "mofa_top_ligrec_all_factors.pdf", width = 4, height = 4.3)

TF_factors_res <- decoupleR::run_ulm(
  mat = as.matrix(RNA_all), network = dorothea_df, minsize = 10)
TF_factors_df <- cosmosR::wide_ulm_res(TF_factors_res)
TF_factors_df <- TF_factors_df[
  order(apply(TF_factors_df, 1, function(x) max(abs(x))), decreasing = TRUE)[
    seq_len(min(25, nrow(TF_factors_df)))], ]
palette_tf <- cosmosR::make_heatmap_color_palette(TF_factors_df)
pheatmap(TF_factors_df, show_rownames = TRUE, cluster_cols = FALSE, cluster_rows = FALSE,
         color = palette_tf, angle_col = 315,
         filename = "mofa_top_TF_all_factors.pdf", width = 4, height = 4.3)


# ---- Per-factor COSMOS+ / MOON loop ----

for (selected_factor in 1:n_factors) {

  cat(sprintf("\n=== COSMOS+ / MOON: Factor %d / %d ===\n", selected_factor, n_factors))

  fact_dir <- paste0("./Factor_", selected_factor, "_MOON")
  dir.create(fact_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- 1. Extract weights for the current factor ----
  RNA <- data.frame(weights$Transcriptomics[, selected_factor])
  rownames(RNA) <- rownames(weights$Transcriptomics)

  # ---- 2. Compute LR and TF activities ----
  ligrec_high_vs_low <- decoupleR::run_ulm(
    mat = as.matrix(RNA), network = ligrec_geneset, .source = set, .target = gene, minsize = 2)
  ligrec_high_vs_low        <- ligrec_high_vs_low[ligrec_high_vs_low$statistic == "ulm", ]
  ligrec_high_vs_low_vector <- setNames(ligrec_high_vs_low$score, ligrec_high_vs_low$source)

  ligrec_top15 <- ligrec_high_vs_low[
    order(abs(ligrec_high_vs_low$score), decreasing = TRUE)[seq_len(min(15, nrow(ligrec_high_vs_low)))],
    c("source", "score")]
  ligrec_top15$source <- factor(ligrec_top15$source, levels = ligrec_top15$source)
  ggplot(ligrec_top15, aes(x = source, y = score)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 315, vjust = 0.5, hjust = 0)) +
    ggtitle(paste0("Top LR activities — Factor ", selected_factor))
  ggsave(file.path(fact_dir, "mofa_top_ligrec_barplot.pdf"), width = 8, height = 4)

  TF_high_vs_low <- decoupleR::run_ulm(
    mat = as.matrix(RNA), network = dorothea_df, minsize = 10)
  TF_high_vs_low        <- TF_high_vs_low[TF_high_vs_low$statistic == "ulm", ]
  TF_high_vs_low_vector <- setNames(TF_high_vs_low$score, TF_high_vs_low$source)
  TF_high_vs_low        <- as.data.frame(TF_high_vs_low[, c("source", "score")])
  rownames(TF_high_vs_low) <- TF_high_vs_low[, 1]

  TF_top10 <- TF_high_vs_low[
    order(abs(TF_high_vs_low$score), decreasing = TRUE)[seq_len(min(10, nrow(TF_high_vs_low)))], ]
  TF_top10$source <- factor(TF_top10$source, levels = TF_top10$source)
  ggplot(TF_top10, aes(x = source, y = score)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    ggtitle(paste0("Top TF activities — Factor ", selected_factor))
  ggsave(file.path(fact_dir, "mofa_top_TF_barplot.pdf"), width = 6, height = 4)

  ligrec_TF_moon_inputs <- list("ligrec" = ligrec_high_vs_low_vector, "TF" = TF_high_vs_low_vector)
  save(ligrec_TF_moon_inputs, file = file.path(fact_dir, "ligrec_TF_moon_inputs.Rdata"))

  # ---- 3. Prepare and filter MOON inputs ----
  signaling_input <- ligrec_TF_moon_inputs$TF
  ligrec_input    <- ligrec_TF_moon_inputs$ligrec

  RNA_input <- weights$Transcriptomics[, selected_factor]
  names(RNA_input) <- rownames(weights$Transcriptomics)

  pdf(file.path(fact_dir, "density_inputs.pdf"))
  {plot(density(RNA_input), main = paste0("RNA weights — Factor ", selected_factor))
   abline(v = -0.2, col = "red"); abline(v = 0.2, col = "red")}
  RNA_input <- ifelse(RNA_input > -0.2 & RNA_input < 0.2, 0, sign(RNA_input) * 10)

  {plot(density(signaling_input), main = paste0("TF activities — Factor ", selected_factor))
   abline(v = -0.5, col = "red"); abline(v = 3.5, col = "red")}
  TF_to_remove    <- signaling_input[signaling_input > -0.5 & signaling_input < 3.5]
  signaling_input <- signaling_input[signaling_input < -0.5 | signaling_input > 3.5]

  {plot(density(ligrec_input), main = paste0("LR activities — Factor ", selected_factor))
   abline(v = -0.5, col = "red"); abline(v = 2.5, col = "red")}
  dev.off()

  ligrec_input <- ligrec_input[ligrec_input > 2.5 | ligrec_input < -0.5]
  rec_inputs   <- ligrec_input
  names(rec_inputs) <- gsub(".+___", "", names(rec_inputs))
  rec_inputs    <- tapply(rec_inputs, names(rec_inputs), mean)
  rec_inputs    <- setNames(as.numeric(rec_inputs), names(rec_inputs))

  # ---- 4. Prepare metabolite inputs ----
  metab_raw <- setNames(weights$Metabolomics[, selected_factor], rownames(weights$Metabolomics))

  # Model feature names use hyphens; name_map Query uses underscores (see venn loop above)
  metab_query_fact <- gsub("-", "_", names(metab_raw))
  hmdb_map_fact <- metab_to_hmdb_moon[metab_to_hmdb_moon$Query %in% metab_query_fact, ]
  hmdb_map_fact$HMDB[is.na(hmdb_map_fact$HMDB)] <-
    paste0("No_HMDB_", hmdb_map_fact$Query[is.na(hmdb_map_fact$HMDB)])

  common_metabs_moon <- intersect(metab_query_fact, hmdb_map_fact$Query)
  metab_inputs_moon  <- metab_raw[match(common_metabs_moon, metab_query_fact)]
  metab_inputs_moon  <- as.numeric(scale(metab_inputs_moon, center = FALSE))
  names(metab_inputs_moon) <- hmdb_map_fact$HMDB[match(common_metabs_moon, hmdb_map_fact$Query)]
  metab_inputs_moon  <- cosmosR::prepare_metab_inputs(metab_inputs_moon, compartment_codes = c("m", "c"))
  metab_to_exclude   <- metab_inputs_moon[abs(metab_inputs_moon) < 0.2]
  metab_inputs_moon  <- metab_inputs_moon[abs(metab_inputs_moon) >= 0.2]

  # ---- 5. Filter the meta prior knowledge network ----
  data("meta_network")
  meta_network_base <- cosmosR::meta_network_cleanup(meta_network)
  meta_network_base <- meta_network_base[
    !(meta_network_base$source %in% names(TF_to_remove)) &
      !(meta_network_base$target %in% names(TF_to_remove)), ]
  meta_network_base <- meta_network_base[
    !(meta_network_base$source %in% names(metab_to_exclude)) &
      !(meta_network_base$target %in% names(metab_to_exclude)), ]

  meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(
    names(expressed_genes), meta_pkn = meta_network_base)

  TF_inputs <- scale(ligrec_TF_moon_inputs$TF, center = FALSE)
  TF_inputs  <- setNames(TF_inputs[, 1], rownames(TF_inputs))

  upstream_inputs   <- c(rec_inputs)
  downstream_inputs <- c(metab_inputs_moon, TF_inputs)

  upstream_inputs_filtered   <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs, meta_network_filtered)
  downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs, meta_network_filtered)

  n_steps <- 6
  meta_network_filtered <- cosmosR:::keep_controllable_neighbours(
    meta_network_filtered, n_steps, names(upstream_inputs_filtered))
  downstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    downstream_inputs_filtered, meta_network_filtered)
  meta_network_filtered <- cosmosR:::keep_observable_neighbours(
    meta_network_filtered, n_steps, names(downstream_inputs_filtered))
  upstream_inputs_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(
    upstream_inputs_filtered, meta_network_filtered)

  write.csv(meta_network_filtered,
            file = file.path(fact_dir, "meta_network_filtered.csv"), row.names = FALSE)

  meta_network_compressed_list <- cosmosR::compress_same_children(
    meta_network_filtered,
    sig_input   = upstream_inputs_filtered,
    metab_input = downstream_inputs_filtered)
  meta_network_compressed <- cosmosR::meta_network_cleanup(
    meta_network_compressed_list$compressed_network)

  # ---- 6. Build MOFA weight annotation table from the model's weight matrices ----
  MOFA_weights <- do.call(rbind, lapply(names(weights), function(v) {
    mat <- weights[[v]]
    data.frame(
      feature = rownames(mat),
      value   = mat[, factor_cols[selected_factor]],
      view    = v,
      stringsAsFactors = FALSE
    )
  }))
  MOFA_weights <- MOFA_weights %>%
    arrange(feature, view) %>%
    filter(!duplicated(feature))
  MOFA_weights <- MOFA_weights[, c("feature", "value")]
  colnames(MOFA_weights) <- c("Nodes", "mofa_weights")
  MOFA_weights$sign  <- sign(MOFA_weights$mofa_weights)
  MOFA_weights$mofa_weights <- abs(MOFA_weights$mofa_weights)
  MOFA_weights <- MOFA_weights %>%
    group_by(Nodes) %>%
    summarise(mofa_weights = max(mofa_weights, na.rm = TRUE),
              sign = max(sign, na.rm = TRUE), .groups = "drop") %>%
    as.data.frame()
  MOFA_weights$mofa_weights <- MOFA_weights$mofa_weights * MOFA_weights$sign
  MOFA_weights <- MOFA_weights[, c("Nodes", "mofa_weights")]

  MOFA_weights_metabaddon <- data.frame(
    mofa_weights = metab_inputs_moon,
    Nodes = sapply(names(metab_inputs_moon), function(x, HMDB_mapper_vec) {
      x       <- gsub("Metab__", "", x)
      suffixe <- stringr::str_extract(x, "_[a-z]$")
      x       <- gsub("_[a-z]$", "", x)
      if (x %in% names(HMDB_mapper_vec)) x <- paste0("Metab__", HMDB_mapper_vec[x])
      if (!is.na(suffixe)) x <- paste0(x, suffixe)
      x
    }, HMDB_mapper_vec = HMDB_mapper_vec)
  )
  MOFA_weights <- as.data.frame(rbind(MOFA_weights, MOFA_weights_metabaddon))

  lig_weights_ann <- ligrec_TF_moon_inputs$ligrec
  names(lig_weights_ann) <- gsub("___.+", "", names(lig_weights_ann))
  lig_weights_ann <- data.frame(
    Nodes           = names(tapply(lig_weights_ann, names(lig_weights_ann), mean)),
    feature_weights = as.numeric(tapply(lig_weights_ann, names(lig_weights_ann), mean)))

  rec_weights_ann <- ligrec_TF_moon_inputs$ligrec
  names(rec_weights_ann) <- gsub(".+___", "", names(rec_weights_ann))
  rec_weights_ann <- data.frame(
    Nodes           = names(tapply(rec_weights_ann, names(rec_weights_ann), mean)),
    feature_weights = as.numeric(tapply(rec_weights_ann, names(rec_weights_ann), mean)))

  TF_weights_ann <- ligrec_TF_moon_inputs$TF
  TF_weights_ann <- TF_weights_ann[!(names(TF_weights_ann) %in%
                                       c(rec_weights_ann$Nodes, lig_weights_ann$Nodes))]
  TF_weights_ann <- data.frame(Nodes = names(TF_weights_ann), feature_weights = TF_weights_ann)

  feature_weights <- as.data.frame(rbind(TF_weights_ann, rec_weights_ann, lig_weights_ann))

  # ---- 7. MOON run 1: Receptors → TF + Metabolites ----
  cat("  Running MOON: Receptors → TF + Metabolites\n")
  upstream_inputs_filtered   <- upstream_inputs_filtered[is.finite(upstream_inputs_filtered)]
  downstream_inputs_filtered <- downstream_inputs_filtered[is.finite(downstream_inputs_filtered)]
  moon_res_rec_to_TFmet <- data.frame()
  SIF_rec_to_TFmetab <- data.frame(source=character(), target=character(), sign=numeric(), Weight=numeric())
  ATT_rec_to_TFmetab <- data.frame(Nodes=character(), AvgAct=numeric(), mofa_weights=numeric(), feature_weights=numeric(), NodeType=numeric())
  tryCatch({
    meta_network_rec_to_TFmetab <- meta_network_compressed
    before <- 1; after <- 0; i <- 1
    while (before != after & i < 10) {
      before <- nrow(meta_network_rec_to_TFmetab)
      moon_res <- cosmosR::moon(
        upstream_input   = upstream_inputs_filtered,
        downstream_input = downstream_inputs_filtered,
        meta_network     = meta_network_rec_to_TFmetab,
        n_layers         = n_steps,
        statistic        = "ulm")
      meta_network_rec_to_TFmetab <- cosmosR::filter_incohrent_TF_target(
        moon_res, dorothea_df, meta_network_rec_to_TFmetab, RNA_input)
      after <- nrow(meta_network_rec_to_TFmetab)
      i <- i + 1
    }
    cat(if (i < 10) paste0("  Converged after ", i - 1, " iterations\n")
        else paste0("  Interrupted after ", i, " iterations. Convergence uncertain.\n"))

    moon_res <- cosmosR::decompress_moon_result(
      moon_res, meta_network_compressed_list, meta_network_rec_to_TFmetab)

    moon_res_rec_to_TFmet      <- moon_res[, c(4, 2, 3)]
    moon_res_rec_to_TFmet[, 1] <- cosmosR::translate_column_HMDB(moon_res_rec_to_TFmet[, 1], HMDB_mapper_vec)
    levels_rec     <- moon_res[, c(4, 3)]
    moon_res_score <- moon_res[, c(4, 2, 3)]; names(moon_res_score)[1] <- "source"

    pdf(file.path(fact_dir, "density_moon_scores_rec_to_TFmetab.pdf"))
    if (sum(!is.na(moon_res_score$score)) >= 2) {
      plot(density(moon_res_score$score, na.rm = TRUE), main = "MOON scores: Rec → TF/Metab")
      abline(v = 1); abline(v = -1)
    } else {
      plot.new(); title(main = "MOON scores: Rec → TF/Metab (insufficient data)")
    }
    dev.off()

    solution_network_1 <- cosmosR::reduce_solution_network(
      decoupleRnival_res = moon_res_score,
      meta_network       = meta_network_filtered,
      cutoff             = solution_cutoff_rec_to_TFmetab,
      upstream_input     = upstream_inputs_filtered,
      RNA_input          = RNA_input,
      n_steps            = n_steps)

    if (!is.null(solution_network_1$SIF) && nrow(solution_network_1$SIF) > 0) {
      SIF_rec_to_TFmetab <- solution_network_1$SIF; names(SIF_rec_to_TFmetab)[3] <- "sign"
      ATT_rec_to_TFmetab <- solution_network_1$ATT

      translated_1        <- cosmosR::translate_res(SIF_rec_to_TFmetab, ATT_rec_to_TFmetab, HMDB_mapper_vec)
      levels_translated_1 <- cosmosR::translate_res(SIF_rec_to_TFmetab, levels_rec, HMDB_mapper_vec)[[2]]
      SIF_rec_to_TFmetab  <- translated_1[[1]]
      ATT_rec_to_TFmetab  <- translated_1[[2]]

      ATT_rec_to_TFmetab <- merge(ATT_rec_to_TFmetab, MOFA_weights, all.x = TRUE)
      ATT_rec_to_TFmetab <- merge(ATT_rec_to_TFmetab, feature_weights, all.x = TRUE)
      names(ATT_rec_to_TFmetab)[2] <- "AvgAct"
      ATT_rec_to_TFmetab$NodeType <- ifelse(
        ATT_rec_to_TFmetab$Nodes %in% levels_translated_1[levels_translated_1$level == 0, 1], 1, 0)
      ATT_rec_to_TFmetab$NodeType <- ifelse(
        ATT_rec_to_TFmetab$Nodes %in% ligrec_df$Node1, 2,
        ifelse(ATT_rec_to_TFmetab$Nodes %in% ligrec_df$Node2, 3,
          ifelse(ATT_rec_to_TFmetab$Nodes %in% dorothea_df$source, 4,
                 ATT_rec_to_TFmetab$NodeType)))
      names(SIF_rec_to_TFmetab)[4] <- "Weight"
    } else {
      message("  [!] Empty Rec->TF/Metab network (Factor ", selected_factor, "); using empty placeholders.")
    }
  }, error = function(e) {
    message("  [!] MOON run 1 failed (Factor ", selected_factor, "): ", conditionMessage(e))
  })

  write.csv(SIF_rec_to_TFmetab,
            file = file.path(fact_dir, "SIF_rec_TFmetab.csv"), row.names = FALSE)
  write.csv(ATT_rec_to_TFmetab,
            file = file.path(fact_dir, "ATT_rec_TFmetab.csv"), row.names = FALSE)

  # ---- 8. MOON run 2: TF → Ligands ----
  cat("  Running MOON: TF → Ligands\n")
  dorothea_PKN_filtered <- cosmosR:::filter_pkn_expressed_genes(
    names(expressed_genes), meta_pkn = dorothea_PKN_base)

  upstream_inputs_2 <- setNames(TF_weights_ann$feature_weights, TF_weights_ann$Nodes)
  upstream_inputs_2 <- upstream_inputs_2[abs(upstream_inputs_2) > 2]

  lig_inputs_2 <- ligrec_TF_moon_inputs$ligrec
  names(lig_inputs_2) <- gsub("___.+", "", names(lig_inputs_2))
  lig_inputs_2 <- tapply(lig_inputs_2, names(lig_inputs_2), mean)
  lig_inputs_2 <- scale(as.numeric(lig_inputs_2), center = FALSE)
  lig_inputs_2 <- setNames(lig_inputs_2[, 1], names(tapply(
    ligrec_TF_moon_inputs$ligrec, gsub("___.+", "", names(ligrec_TF_moon_inputs$ligrec)), mean)))

  downstream_inputs_2 <- lig_inputs_2

  upstream_inputs_filtered_2   <- cosmosR:::filter_input_nodes_not_in_pkn(upstream_inputs_2, dorothea_PKN_filtered)
  downstream_inputs_filtered_2 <- cosmosR:::filter_input_nodes_not_in_pkn(downstream_inputs_2, dorothea_PKN_filtered)

  n_steps_2 <- 1
  dorothea_PKN_filtered <- cosmosR:::keep_controllable_neighbours(
    dorothea_PKN_filtered, n_steps_2, names(upstream_inputs_filtered_2))
  downstream_inputs_filtered_2 <- cosmosR:::filter_input_nodes_not_in_pkn(
    downstream_inputs_filtered_2, dorothea_PKN_filtered)
  dorothea_PKN_filtered <- cosmosR:::keep_observable_neighbours(
    dorothea_PKN_filtered, n_steps_2, names(downstream_inputs_filtered_2))
  upstream_inputs_filtered_2 <- cosmosR:::filter_input_nodes_not_in_pkn(
    upstream_inputs_filtered_2, dorothea_PKN_filtered)

  write.csv(dorothea_PKN_filtered,
            file = file.path(fact_dir, "meta_network_TF_lig.csv"), row.names = FALSE)

  moon_res_TF_lig <- data.frame()
  SIF_TF_lig <- data.frame(source=character(), target=character(), sign=numeric(), Weight=numeric())
  ATT_TF_lig <- data.frame(Nodes=character(), AvgAct=numeric(), mofa_weights=numeric(), feature_weights=numeric(), NodeType=numeric())
  tryCatch({
    upstream_inputs_filtered_2   <- upstream_inputs_filtered_2[is.finite(upstream_inputs_filtered_2)]
    downstream_inputs_filtered_2 <- downstream_inputs_filtered_2[is.finite(downstream_inputs_filtered_2)]
    meta_network_TF_lig <- dorothea_PKN_filtered
    before <- 1; after <- 0; i <- 1
    while (before != after & i < 10) {
      before <- nrow(meta_network_TF_lig)
      moon_res_2 <- cosmosR::moon(
        upstream_input   = upstream_inputs_filtered_2,
        downstream_input = downstream_inputs_filtered_2,
        meta_network     = meta_network_TF_lig,
        n_layers         = n_steps_2,
        statistic        = "ulm")
      meta_network_TF_lig <- cosmosR::filter_incohrent_TF_target(
        moon_res_2, dorothea_df, meta_network_TF_lig, RNA_input)
      after <- nrow(meta_network_TF_lig)
      i <- i + 1
    }
    cat(if (i < 10) paste0("  Converged after ", i - 1, " iterations\n")
        else paste0("  Interrupted after ", i, " iterations. Convergence uncertain.\n"))

    moon_res_TF_lig      <- moon_res_2
    moon_res_TF_lig[, 1] <- cosmosR::translate_column_HMDB(moon_res_TF_lig[, 1], HMDB_mapper_vec)
    levels_lig       <- moon_res_2[, c(1, 3)]
    moon_res_2_score <- moon_res_2[, c(1, 2, 3)]; names(moon_res_2_score)[1] <- "source"

    pdf(file.path(fact_dir, "density_moon_scores_TF_lig.pdf"))
    if (sum(!is.na(moon_res_2_score$score)) >= 2) {
      plot(density(moon_res_2_score$score, na.rm = TRUE), main = "MOON scores: TF → Ligands")
      abline(v = 1); abline(v = -1)
    } else {
      plot.new(); title(main = "MOON scores: TF → Ligands (insufficient data)")
    }
    dev.off()

    solution_network_2 <- cosmosR::reduce_solution_network(
      decoupleRnival_res = moon_res_2_score,
      meta_network       = as.data.frame(dorothea_PKN_filtered[, c(1, 3, 2)]),
      cutoff             = solution_cutoff_TF_lig,
      upstream_input     = upstream_inputs_filtered_2,
      RNA_input          = RNA_input,
      n_steps            = n_steps_2)

    if (!is.null(solution_network_2$SIF) && nrow(solution_network_2$SIF) > 0) {
      SIF_TF_lig <- solution_network_2$SIF; names(SIF_TF_lig)[3] <- "sign"
      ATT_TF_lig <- solution_network_2$ATT

      translated_2        <- cosmosR::translate_res(SIF_TF_lig, ATT_TF_lig, HMDB_mapper_vec)
      levels_translated_2 <- cosmosR::translate_res(SIF_TF_lig, levels_lig, HMDB_mapper_vec)[[2]]
      SIF_TF_lig <- translated_2[[1]]
      ATT_TF_lig <- translated_2[[2]]

      ATT_TF_lig <- merge(ATT_TF_lig, MOFA_weights, all.x = TRUE)
      ATT_TF_lig <- merge(ATT_TF_lig, feature_weights, all.x = TRUE)
      names(ATT_TF_lig)[2] <- "AvgAct"
      ATT_TF_lig$NodeType <- ifelse(
        ATT_TF_lig$Nodes %in% levels_translated_2[levels_translated_2$level == 0, 1], 1, 0)
      ATT_TF_lig$NodeType <- ifelse(
        ATT_TF_lig$Nodes %in% ligrec_df$Node1, 2,
        ifelse(ATT_TF_lig$Nodes %in% ligrec_df$Node2, 3,
          ifelse(ATT_TF_lig$Nodes %in% dorothea_df$source, 4, ATT_TF_lig$NodeType)))
      names(SIF_TF_lig)[4] <- "Weight"
    } else {
      message("  [!] Empty TF->Ligands network (Factor ", selected_factor, "); using empty placeholders.")
    }
  }, error = function(e) {
    message("  [!] MOON run 2 failed (Factor ", selected_factor, "): ", conditionMessage(e))
  })

  write.csv(SIF_TF_lig, file = file.path(fact_dir, "SIF_TF_lig.csv"), row.names = FALSE)
  write.csv(ATT_TF_lig, file = file.path(fact_dir, "ATT_TF_lig.csv"), row.names = FALSE)

  # ---- 9. Combine both MOON networks ----
  cat("  Combining MOON networks\n")
  combined_SIF_moon <- unique(as.data.frame(rbind(SIF_rec_to_TFmetab, SIF_TF_lig)))

  combined_ATT_moon <- as.data.frame(rbind(ATT_rec_to_TFmetab, ATT_TF_lig)) %>%
    group_by(Nodes) %>%
    summarise(across(everything(), \(x) mean(x, na.rm = TRUE)), .groups = "drop") %>%
    as.data.frame()

  ligrec_ressource_addon <- ligrec_ressource[
    ligrec_ressource$source_genesymbol %in% combined_SIF_moon$target &
      ligrec_ressource$target_genesymbol %in% combined_SIF_moon$source,
    c("source_genesymbol", "target_genesymbol")]
  ligrec_ressource_addon$sign   <- 1
  ligrec_ressource_addon$Weight <- 1
  names(ligrec_ressource_addon)[c(1, 2)] <- c("source", "target")
  combined_SIF_moon <- as.data.frame(rbind(combined_SIF_moon, unique(ligrec_ressource_addon)))

  write.csv(combined_SIF_moon,
            file = file.path(fact_dir, "combined_SIF_moon.csv"), row.names = FALSE)
  write.csv(combined_ATT_moon,
            file = file.path(fact_dir, "combined_ATT_moon.csv"), row.names = FALSE)
  write.csv(moon_res_rec_to_TFmet,
            file = file.path(fact_dir, "moon_res_rec_to_TFmet.csv"), row.names = FALSE)
  write.csv(moon_res_TF_lig,
            file = file.path(fact_dir, "moon_res_TF_lig.csv"), row.names = FALSE)

  # ---- 10. Network visualization ----
  cat("  Visualizing combined MOON network\n")
  if (nrow(combined_SIF_moon) > 0) {
    all_nodes <- union(combined_SIF_moon$source, combined_SIF_moon$target)
    att_full  <- merge(data.frame(Nodes = all_nodes, stringsAsFactors = FALSE),
                       combined_ATT_moon, by = "Nodes", all.x = TRUE)
    att_full$AvgAct[is.na(att_full$AvgAct)]   <- 0
    att_full$NodeType[is.na(att_full$NodeType)] <- 0

    g <- igraph::graph_from_data_frame(
      d        = combined_SIF_moon[, c("source", "target", "sign", "Weight")],
      directed = TRUE,
      vertices = att_full)

    node_type_map <- c("0" = "Other", "1" = "Level 0",
                       "2" = "Ligand", "3" = "Receptor", "4" = "TF")
    V(g)$NodeLabel <- node_type_map[as.character(V(g)$NodeType)]
    V(g)$NodeLabel[is.na(V(g)$NodeLabel)] <- "Other"

    p_net <- ggraph(g, layout = "stress") +
      geom_edge_link(aes(color = factor(sign), width = abs(Weight)),
                     arrow = arrow(length = unit(2, "mm"), type = "closed"),
                     end_cap = circle(3, "mm"), alpha = 0.6) +
      geom_node_point(aes(size  = pmax(1, abs(AvgAct), na.rm = TRUE),
                          color = NodeLabel), alpha = 0.9) +
      geom_node_text(aes(label = name), size = 2, repel = TRUE, max.overlaps = 30) +
      scale_edge_color_manual(
        values = c("-1" = "firebrick", "1" = "steelblue"), name = "Sign") +
      scale_edge_width(range = c(0.3, 1.5), guide = "none") +
      scale_color_manual(
        values = c("Other"    = "grey70",  "Level 0"  = "steelblue",
                   "Ligand"   = "coral",   "Receptor" = "forestgreen",
                   "TF"       = "gold"),
        name = "Node type") +
      scale_size(range = c(2, 8), name = "|Activity|") +
      theme_graph(base_family = "sans") +
      ggtitle(paste0("MOON network — Factor ", selected_factor))

    ggsave(file.path(fact_dir, "MOON_network.pdf"), p_net, width = 14, height = 10)
    cat(sprintf("  Saved MOON network plot in %s\n", fact_dir))
  } else {
    message("  [!] Empty combined network; skipping visualization (Factor ", selected_factor, ").")
  }

  # ---- 11. EnrichR functional enrichment ----
  cat("  Running EnrichR enrichment\n")
  enrich_genes <- combined_ATT_moon$Nodes[
    !grepl("^Metab__", combined_ATT_moon$Nodes) & !is.na(combined_ATT_moon$Nodes)]

  if (length(enrich_genes) >= 3) {
    enrich_dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human",
                    "Reactome_2022", "MSigDB_Hallmark_2020")
    enrich_res <- enrichR::enrichr(enrich_genes, databases = enrich_dbs)

    for (db in enrich_dbs) {
      df <- enrich_res[[db]]
      if (!is.null(df) && nrow(df) > 0) {
        write.csv(df,
          file = file.path(fact_dir,
            paste0("enrichr_", gsub("[^A-Za-z0-9]", "_", db), ".csv")),
          row.names = FALSE)
      }
    }

    for (db in c("GO_Biological_Process_2023", "KEGG_2021_Human")) {
      df <- enrich_res[[db]]
      if (!is.null(df) && nrow(df) > 0) {
        p_enr <- plot_enrichr(df, top_n = 20, numChar = 60,
                               title = paste0(db, " — Factor ", selected_factor))
        if (!is.null(p_enr)) {
          ggsave(
            file.path(fact_dir,
              paste0("enrichr_", gsub("[^A-Za-z0-9]", "_", db), ".pdf")),
            p_enr, width = 10, height = 6)
        }
      }
    }
    cat(sprintf("  Saved EnrichR results in %s\n", fact_dir))
  } else {
    message("  [!] Too few genes for enrichment (n=", length(enrich_genes),
            "); skipping (Factor ", selected_factor, ").")
  }

  cat(sprintf("  Factor %d complete. Results saved in: %s\n", selected_factor, fact_dir))

} # end per-factor loop

cat(sprintf("\nCOSMOS+ / MOON analysis complete for all %d factors.\n", n_factors))


###########################################################################################
#### Aggregate EnrichR results across all factors into one summary file
###########################################################################################

library(openxlsx)

enrich_dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human",
                "Reactome_2022", "MSigDB_Hallmark_2020")

enrichr_summary <- data.frame()
for (fact in 1:n_factors) {
  fact_dir <- paste0("./Factor_", fact, "_MOON")
  for (db in enrich_dbs) {
    f <- file.path(fact_dir, paste0("enrichr_", gsub("[^A-Za-z0-9]", "_", db), ".csv"))
    if (file.exists(f)) {
      df <- read.csv(f, stringsAsFactors = FALSE)
      if (nrow(df) > 0) {
        df$Factor   <- fact
        df$Database <- db
        enrichr_summary <- rbind(enrichr_summary, df)
      }
    }
  }
}

if (nrow(enrichr_summary) > 0) {
  enrichr_summary <- enrichr_summary %>% dplyr::select(Factor, Database, dplyr::everything())

  write.csv(enrichr_summary, file = "Combined_enrichr_all_factors.csv", row.names = FALSE)
  openxlsx::write.xlsx(enrichr_summary, file = "Combined_enrichr_all_factors.xlsx")

  cat(sprintf("\nSaved combined EnrichR summary (%d rows) to Combined_enrichr_all_factors.csv / .xlsx\n",
              nrow(enrichr_summary)))
} else {
  message("[!] No EnrichR result files found to aggregate.")
}

