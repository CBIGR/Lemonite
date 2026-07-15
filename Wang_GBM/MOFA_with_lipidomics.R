#!/usr/bin/Rscript 

setwd('/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/MOFA_with_lipidomics')
base_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/'
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
library(clusterProfiler)
#use_python("/home/boris/Software/Miniconda3/bin/python")

# %%

###########################################################################################
#### Input files
###########################################################################################

expression <- paste0(base_dir, 'data/fpkm_gene_expression.csv') 
metabolomics <- paste0(base_dir, 'data/metabolome.csv')
lipids_pos <- paste0(base_dir, 'data/lipidome_pos.csv')
lipids_neg <- paste0(base_dir, 'data/lipidome_neg.csv')

metadata <- paste0(base_dir, 'data/clinical_metadata.csv')
metadata2 <- paste0(base_dir, 'data/Additional_sample_annotations_suppl2.csv')



###########################################################################################
#### Select samples
###########################################################################################

samples_metabolomics <- colnames(fread(metabolomics, data.table = FALSE, header = TRUE)[,-c(1)])
samples_transcriptomics <- colnames(fread(expression, data.table = FALSE, header = TRUE)[,-c(1:4)])
samples_lipids_pos <- colnames(fread(lipids_pos, data.table = FALSE, header = TRUE)[,-c(1)])
samples_lipids_neg <- colnames(fread(lipids_neg, data.table = FALSE, header = TRUE)[,-c(1)])

metadata <- fread(metadata, data.table =FALSE); metadata$V1 <- NULL
metadata2 <- fread(metadata2, data.table =FALSE); metadata2$V1 <- NULL

# GTEx samples have a different format, remove the 'GTEX-' prefix and everything after the remaining '-'
#metadata$case_submitter_id <- gsub('^GTEX-(\\w+)-.*', 'PT-\\1', metadata$case_submitter_id) # Skip this line to skip control GTEx samples

# Add column 'Diagnosis' to metadata, this should be 'GBM' if 'case_submitter_id' starts with 'C3', else 'control'
metadata$diagnosis <- ifelse(grepl('^C3', metadata$case_submitter_id), 'GBM', 'control')

metadata[metadata == "'--"] <- 'Blank'


# Select samples that are present in all datasets
samples_all_omics <- Reduce(intersect, list(samples_metabolomics, samples_transcriptomics, samples_lipids_pos, samples_lipids_neg))
samples <- Reduce(intersect, list(samples_metabolomics, samples_transcriptomics, samples_lipids_pos, samples_lipids_neg, metadata$case_submitter_id))
metadata <- metadata[metadata$case_submitter_id %in% samples, ]



###########################################################################################
#### Read RNAseq data, filter for protein coding genes
###########################################################################################

RNAseq <- fread(expression, header=TRUE, data.table=TRUE)

# Some gene symbols got converted to datas (suppl data provided in Excel), convert these back
# Replace date-like strings with corresponding gene names
RNAseq <- RNAseq %>%
  mutate(gene_name = case_when(
    gene_name == "1-Dec"  ~ "DELEC1",
    gene_name == "10-Mar" ~ "MARCHF10",
    gene_name == "11-Mar" ~ "MARCHF11",
    gene_name == "12-Sep" ~ "SEPTIN12",
    gene_name == "14-Sep" ~ "SEPTIN14",
    gene_name == "3-Sep"  ~ "SEPTIN3",
    gene_name == "4-Mar"  ~ "MARCHF4",
    gene_name == "9-Mar"  ~ "MARCHF9",
    TRUE ~ gene_name  # Keep the original value if none of the conditions are met
  ))

id_ensembl = RNAseq[,c(2,1)]
write.table(id_ensembl[,c(2,1)], file = paste0(base_dir, 'results/MOFA_with_lipidomics/ensemble_mapping.txt'), quote=FALSE, sep = '\t', row.names = FALSE)
# RNAseq <- merge(RNAseq, id_ensembl, by.x = 'count', by.y='hgnc_symbol')

# Select protein coding genes
RNAseq <- RNAseq[RNAseq$gene_type == 'protein_coding', ]
RNAseq <- as.data.frame(RNAseq[!duplicated(RNAseq$gene_name), ]) # Only 58 duplicate gene names
rownames(RNAseq) <- RNAseq$gene_name
RNAseq <- RNAseq[, -c(1:4)]
RNAseq <- RNAseq[, colnames(RNAseq) %in% samples]
# Set NAs to zero
RNAseq[is.na(RNAseq)] <- 0

# For FPKM data, skip DESeq2 normalization and directly log-transform
# Ensure no zero or negative values before log transformation
RNAseq[RNAseq <= 0] <- 1e-6  # Small positive value for log transformation
M <- log(RNAseq + 1)

# Variance-based HVG selection with fixed threshold 0.7
vars <- apply(M, 1, var)
pdf('./expression_variance_histogram.pdf')
hist(vars, 1100, xlim = c(0, 7))
abline(v = 0.7, col = 'red', lwd = 3)
dev.off()

variable_genes <- M[vars > 0.7, , drop = FALSE]

# Z-score per gene (row-wise)
RNAseq <- as.data.frame(t(scale(t(variable_genes))))



###########################################################################################
#### Preprocessing of metabolomics data
###########################################################################################

abundancies <- fread(metabolomics, data.table = FALSE, header = TRUE)
rownames(abundancies) <- abundancies$Metabolite; abundancies$Metabolite <- NULL
#abundancies <- abundancies %>% mutate_all(na_if, "") # Change empty entries to NA

rownames(abundancies) <- str_replace_all(rownames(abundancies), ' ', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), '-', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), ':', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), '\\+', '_')

# Select samples
abundancies <- abundancies[, colnames(abundancies) %in% samples]

abundancies[is.na(abundancies)] <- 0
#abundancies <- log(abundancies + 1)
abundancies <- as.data.frame(t(pareto_scale(t(abundancies)))) # In principle normalization has already been done...

pdf(paste0(base_dir, 'results/MOFA_with_lipidomics/Normalized_metabolomics.pdf')) # sample MSM719ME does not follow normalized data pattern
abundancies %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()


###########################################################################################
#### Preprocessing lipidomics data
###########################################################################################

lipidomics_pos <- fread(lipids_pos, data.table = FALSE, header = TRUE)
lipidomics_neg <- fread(lipids_neg, data.table = FALSE, header = TRUE)

# Get lipid names
lipids_pos_names <- lipidomics_pos$Lipid
lipids_neg_names <- lipidomics_neg$Lipid

# Concatenate lipid names
lipids_names <- c(lipids_pos_names, lipids_neg_names)
# Replace weird characters
lipids_names <- str_replace_all(lipids_names, ' ', '_')
lipids_names <- str_replace_all(lipids_names, '-', '_')
lipids_names <- str_replace_all(lipids_names, ':', '_')
lipids_names <- str_replace_all(lipids_names, '\\+', '_')

# Select samples
lipidomics_pos <- lipidomics_pos[, colnames(lipidomics_pos) %in% samples]
lipidomics_neg <- lipidomics_neg[, colnames(lipidomics_neg) %in% samples]

lipidomics <- rbind(lipidomics_pos, lipidomics_neg, fill=FALSE)
# Remove the last row, no idea why this is added in the first place (248+334=582 and not 583)
lipidomics <- lipidomics[-nrow(lipidomics),]
lipidomics$symbol <- lipids_names; lipidomics$ensembl_gene_id <- lipids_names

# Rearrange df and place column 'lipid' first


lipidomics <- lipidomics[, c(ncol(lipidomics)-1, ncol(lipidomics),1:(ncol(lipidomics)-2))] # Move lipid column to first column

# Have a look at normalization

lipidomics[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))


# Set NAs to 0 and perform Pareto scaling
lipidomics[is.na(lipidomics)] <- 0
lipidomics[, -c(1,2)] <- as.data.frame(t(pareto_scale(t(lipidomics[, -c(1,2)]))))


pdf(paste0(base_dir, 'results/MOFA_with_lipidomics/Normalized_lipidomics.pdf'))
lipidomics[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()
# Write to file
write.table(lipidomics, paste0(base_dir, 'results/MOFA_with_lipidomics/LemonPreprocessed_lipidomics.txt'), sep = '\t', quote=FALSE, row.names=FALSE)

# Load name mapping for restoring original names in visualizations
name_mapping <- fread(paste0(base_dir, 'data/name_mapping.tsv'), data.table = FALSE)
name_lookup <- setNames(name_mapping$original, name_mapping$cleaned)


###########################################################################################
#### Run MOFA
###########################################################################################
# metadata <- fread("../data/hmp2_metadata.csv", data.table = FALSE)[, c(1:5, 41, 71)]
# metadata <- metadata[(metadata$`External ID` %in% colnames(metabolomics) & metadata$data_type == 'host_transcriptomics'), ]
rownames(metadata) <- metadata$case_submitter_id
rownames(metadata2) <- metadata2$case
# merge metadata on rownames
metadata_merge <- merge(metadata, metadata2, by = 'row.names'); rownames(metadata_merge) <- metadata_merge$Row.names; metadata_merge$Row.names <- NULL

cols_to_keep <- c('gender', 'multiomic', 'nmf_consensus', 'rna_wang_cancer_cell_2017')
coldata <- metadata_merge[,cols_to_keep]
# set columns as numerical factor
# coldata$gender <- as.numeric(as.factor(coldata$gender))
# coldata$multiomic <- as.numeric(as.factor(coldata$multiomic))
# coldata$nmf_consensus <- as.numeric(as.factor(coldata$nmf_consensus))
# coldata$rna_wang_cancer_cell_2017 <- as.numeric(as.factor(coldata$rna_wang_cancer_cell_2017))



# # Save original metadata values
coldata$gender_original <- coldata$gender
coldata$multiomic_original <- coldata$multiomic
coldata$nmf_consensus_original <- coldata$nmf_consensus
coldata$rna_wang_cancer_cell_2017_original <- coldata$rna_wang_cancer_cell_2017

# Convert columns to numerical factors
coldata$gender <- as.numeric(as.factor(coldata$gender))
coldata$multiomic <- as.numeric(as.factor(coldata$multiomic))
coldata$nmf_consensus <- as.numeric(as.factor(coldata$nmf_consensus))
coldata$rna_wang_cancer_cell_2017 <- as.numeric(as.factor(coldata$rna_wang_cancer_cell_2017))

# remove lipids with non-unique symbol
lipidomics <- lipidomics[!duplicated(lipidomics$symbol),]
rownames(lipidomics) <- lipidomics$symbol; lipidomics$symbol <- NULL; lipidomics$ensembl_gene_id <- NULL

# RNAseq: Convert to a SummarizedExperiment
RNAseq_matrix <- as.matrix(RNAseq)  # Ensure it's numeric
RNAseq_se <- SummarizedExperiment(assays = list(counts = RNAseq_matrix))

# Metabolomics: Convert to a SummarizedExperiment
abundancies_matrix <- as.matrix(abundancies)  # Ensure it's numeric
abundancies_se <- SummarizedExperiment(assays = list(counts = abundancies_matrix))

# Lipidomics: Convert to a SummarizedExperiment
lipidomics_matrix <- as.matrix(lipidomics)  # Ensure it's numeric
lipidomics_se <- SummarizedExperiment(assays = list(counts = lipidomics_matrix))

# Create list of experiments
data <- list(
  Transcriptomics = RNAseq_se,
  Metabolomics = abundancies_se,
  Lipidomics = lipidomics_se
)

# Create the MultiAssayExperiment
exp <- MultiAssayExperiment(data, colData = coldata)


MOFAobject <- create_mofa_from_MultiAssayExperiment(exp, extract_metadata = TRUE)
MOFAobject <- set_covariates(MOFAobject, c('gender', 'multiomic', 'rna_wang_cancer_cell_2017', 'nmf_consensus'))

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
# MOFA.trained <- run_mofa(MOFAobject, './MOFA1.hdf5', use_basilisk = TRUE)
# load model
MOFA.trained <- load_model('./MOFA1.hdf5')
# change MOFA view names in pretrained object to match current view names

# Rename metabolite/lipid features in MOFA model to original names for display
for (view in c('Metabolomics', 'Lipidomics')) {
  if (view %in% views_names(MOFA.trained)) {
    feat_names <- features_names(MOFA.trained)[[view]]
    original <- ifelse(feat_names %in% names(name_lookup), name_lookup[feat_names], feat_names)
    features_names(MOFA.trained)[[view]] <- original
  }
}


###########################################################################################
#### Downstream analysis
###########################################################################################
model <- MOFA.trained
plot_data_overview(model)
head(model@cache$variance_explained$r2_per_factor[[1]])# Check how much variance is explained per view
# write variance explained to file
variance_explained <- get_variance_explained(model)
variance_explained <- get_variance_explained(model)$r2_per_factor
# convert to df
variance_explained <- do.call(rbind, variance_explained)
write.table(variance_explained, file = './variance_explained_per_view_and_factor.txt', sep = '\t', quote = FALSE, row.names = TRUE)

plot <- plot_variance_explained(model, x="view", y="factor")
plot + theme(axis.text.x = element_text(size = 15))
ggsave('./variance_explained_per_factor.png')


plot <- plot_variance_explained(model, x="view", y="factor", plot_total = T)
plot + theme(axis.text.x = element_text(size = 15))
ggsave('./total_var_explained.png')

pdf('./Visualization_in_latent_space.pdf')
plot <- plot_factor(model, 
            factor = 1:6,
            color_by = "multiomic",
            shape_by = "gender"
)

# increase label size on x-axis
plot + theme(axis.text.x = element_text(size = 10))


dev.off()

p <- plot_factor(model, 
                 factors = c(1,2,3,4,5,6),
                 color_by = "multiomic",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = T,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

p <- p + 
  scale_color_manual(values=c("nonIBD"="black", "UC"="red")) +
  scale_fill_manual(values=c("nonIBD"="black", "UC"="red"))
pdf('./factor_loadings.pdf')
print(p)
dev.off()

plot_factors(model, 
             factors = 1:6,
             color_by = "multiomic_original",
)
ggsave('./Factor_distribution.png')

plot_factors(model, 
             factors = c(1,2, 3,4,6),
             color_by = "multiomic_original",
)
ggsave('./Factor_distribution_1-2-3-4-6.png')

Z <- get_factors(model)[[1]]              # samples x factors matrix
meta <- samples_metadata(model)
group <- factor(meta$multiomic_original)

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

# PCC between factors and multiomic subtype as covariate - not ideal since conversion to numeric needed
Z <- get_factors(model)[[1]]              # samples x factors matrix
meta <- samples_metadata(model)
covariate <- meta$multiomic_original

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

write.csv(cor_res, file = './factor_multiomic_original_correlation.csv', row.names = TRUE)



for (i in 1:15){
  #pdf(paste0('factor',i,'_feature_weights_metabolomics.pdf'))
  #i <- 4
  plot_top_weights(model,
                   view = "Metabolomics",
                   factor = i,
                   nfeatures = 10
  ) + theme(axis.text.y = element_text(size = 15))
  ggsave(paste0('./feature_weights/factor',i,'_feature_weights_metabolomics.png'))
  #dev.off()

  plot_top_weights(model,
                         view = "Transcriptomics",
                         factor = i,
                         nfeatures = 10
  ) + theme(axis.text.y = element_text(size = 15))
  #dev.off()
  ggsave(paste0('./feature_weights/factor',i,'_feature_weights_expression.png'))
  
  plot_top_weights(model, 
                         view = "Lipidomics",
                         factor = i,
                         nfeatures = 10
  ) + theme(axis.text.y = element_text(size = 15))
  
  ggsave(paste0('./feature_weights/factor',i,'_feature_weights_lipidomics.png'))
}


print(correlate_factors_with_covariates(model, covariates = c("gender", "multiomic", "rna_wang_cancer_cell_2017"), plot="log_pval"))
ggsave('./Factor_correlation_with_covariates.png')


# %%
library(MOFA2)
library(ComplexHeatmap)
library(circlize)
library(grid)

# ============================================================
# PARAMETERS
# ============================================================
views <- c("Metabolomics", "Lipidomics", "Transcriptomics")
factors <- c(1,2, 3, 4, 6)
top_n <- 3

output_file <- "./heatmaps/combined_top3_factors1-2-3-4-6_weights_horizontal.pdf"

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
  output_file,
  width = 20,
  height = 6
)

draw(
  ht,
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

dev.off()

cat("\nSaved:\n")
cat(output_file, "\n")



###########################################################################################
#### Factor exploration - feature weights and GSEA
###########################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
organism <- org.Hs.eg.db

set.seed(1234)

# Source GSEA plotting utilities
source('/home/borisvdm/repo/gsea_plotting_utils.R')

# Create GSEA results directory
gsea_dir <- paste0(base_dir, 'results/MOFA_with_lipidomics/GSEA_all_factors')
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Get the number of factors
n_factors <- model@dimensions$K

# Run GSEA for all factors
for (fact in 1:n_factors) {
  #fact <- 6
  cat(sprintf("\n=== Processing Factor %d ===\n", fact))
  
  # Create factor-specific directory
  fact_dir <- paste0(gsea_dir, "/Factor_", fact)
  dir.create(fact_dir, showWarnings = FALSE)
  
  # Extract expression weights for this factor
  weights <- get_weights(model, view = "Transcriptomics", as.data.frame = TRUE)
  fact_weights <- subset(weights, factor == paste0('Factor', fact))
  weights_vec <- fact_weights$value
  names(weights_vec) <- fact_weights$feature
  
  # Remove genes with NA weights or NA/empty names
  weights_vec <- weights_vec[!is.na(weights_vec) & !is.na(names(weights_vec)) & names(weights_vec) != ""]
  
  # Remove zero weights
  weights_vec <- weights_vec[weights_vec != 0]
  
  # Sort in decreasing order
  weights_vec <- sort(weights_vec, decreasing = TRUE)
  
  # GO Enrichment Analysis
  for (db in c('BP', 'MF', 'CC')) {
    tryCatch({
      gse_result <- gseGO(
        geneList = weights_vec,
        ont = db,
        keyType = "SYMBOL",
        nPermSimple = 10000,
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.10,
        verbose = FALSE,
        OrgDb = organism,
        pAdjustMethod = "BH",
        seed = TRUE
      )
      
      # Save results table
      if (nrow(gse_result) > 0) {
        write.table(as.data.frame(gse_result),
                   file = paste0(fact_dir, "/Factor", fact, "_gseGO_", db, "_results.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Create and save plot using new visualization
        db_title <- switch(db,
                           "BP" = "GO Biological Process",
                           "MF" = "GO Molecular Function",
                           "CC" = "GO Cellular Component",
                           db)
        safe_plot_save(gse_result, 
                      paste0(fact_dir, "/Factor", fact, "_gseGO_", db, ".png"),
                      db_name = db_title)
        cat(sprintf("  Saved GO %s enrichment\n", db))
      }
    }, error = function(e) {
      cat(sprintf("  No enrichment found for GO %s\n", db))
    })
  }

}

cat(sprintf("\n=== GSEA Analysis Complete ===\n"))
cat(sprintf("Results saved in: %s\n", gsea_dir))

## Aggregate per-factor GSEA results into combined GO file (BP, MF, CC only)
agg_go <- data.frame()

for (fact in 1:n_factors) {
  factor_dir <- file.path(gsea_dir, paste0("Factor_", fact))

  # GO files (BP, MF, CC) — skip 'ALL'
  for (db in c('BP', 'MF', 'CC')) {
    file_go <- file.path(factor_dir, paste0("Factor", fact, "_gseGO_", db, "_results.txt"))
    if (file.exists(file_go)) {
      res <- tryCatch(read.table(file_go, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        res$Factor <- fact
        res$Database <- db
        agg_go <- rbind(agg_go, res)
      }
    }
  }
}

# Write aggregated GO output
if (nrow(agg_go) > 0) {
  
  # Put columns factor and database first
  agg_go <- dplyr::select(agg_go, Factor, Database, dplyr::everything())
  
  write.table(agg_go, file = file.path(gsea_dir, 'Combined_GSEA_GO_BP_MF_CC_all_factors.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
}

###########################################################################################
#### Save feature weights for all factors (per-factor files AND one combined file)
###########################################################################################
# Get the number of factors
n_factors <- model@dimensions$K
# Create weights directory
weights_dir <- paste0(base_dir, 'results/MOFA_with_lipidomics/feature_weights_all_factors')
dir.create(weights_dir, showWarnings = FALSE, recursive = TRUE)

# Re-fetch weights as named list (per view) for per-factor file writing
weights <- get_weights(model, views = "all", factors = "all")

# Write per-factor files and plots (keep existing behaviour)
for (fact in 1:n_factors) {
  cat(sprintf("Saving weights for Factor %d\n", fact))
  
  # Expression weights (per-factor file)
  if ("Transcriptomics" %in% names(weights)) {
    expression_weights <- as.data.frame(weights$Transcriptomics[, fact])
    colnames(expression_weights) <- paste0("Factor_", fact)
    write.table(expression_weights, file = paste0(weights_dir, "/Factor", fact, "_expression_weights.txt"),
                sep = "\t", quote = FALSE)
  }
  
  # Metabolomics weights (per-factor file)
  if ("Metabolomics" %in% names(weights)) {
    metab_weights <- as.data.frame(weights$Metabolomics[, fact])
    colnames(metab_weights) <- paste0("Factor_", fact)
    write.table(metab_weights, file = paste0(weights_dir, "/Factor", fact, "_metabolomics_weights.txt"),
                sep = "\t", quote = FALSE)
  }
  
  # Lipidomics weights if present (per-factor file)
  if ("Lipidomics" %in% names(weights)) {
    lipid_weights <- as.data.frame(weights$Lipidomics[, fact])
    colnames(lipid_weights) <- paste0("Factor_", fact)
    write.table(lipid_weights, file = paste0(weights_dir, "/Factor", fact, "_lipidomics_weights.txt"),
                sep = "\t", quote = FALSE)
  }
  
  # Save top features plots (best-effort)
  tryCatch({
    pdf(paste0(weights_dir, "/Factor", fact, "_top_loadings.pdf"), width = 10, height = 8)
    plot_weights(model, view = "Transcriptomics", factor = fact, nfeatures = 15, 
                 text_size = 3, scale = TRUE)
    dev.off()
  }, error = function(e) {
    cat(sprintf("Could not save loadings plot for factor %d\n", fact))
  })
}

# Aggregate all weights across views and factors into a single combined file
cat(sprintf("Aggregating feature weights across all views and factors\n"))
weights <- get_weights(model, views = "all", factors = "all", as.data.frame = TRUE)

weights_expr <- subset(weights, view == "Transcriptomics")
weights_metab <- subset(weights, view == "Metabolomics")
weights_lipid <- subset(weights, view == "Lipidomics")

# convert to wide format on column 'factor'
weights_expr <- reshape2::dcast(weights_expr, feature ~ factor, value.var = "value")
weights_metab <- reshape2::dcast(weights_metab, feature ~ factor, value.var = "value")
weights_lipid <- reshape2::dcast(weights_lipid, feature ~ factor, value.var = "value")

# write to file
write.table(weights_expr, file = paste0(weights_dir, "/Combined_expression_weights_all_factors.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(weights_metab, file = paste0(weights_dir, "/Combined_metabolomics_weights_all_factors.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(weights_lipid, file = paste0(weights_dir, "/Combined_lipidomics_weights_all_factors.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Feature weights (per-factor files and combined) saved in: %s\n", weights_dir))



#######################################################################################################################
#### COSMOS+ / MOON analysis
#### Moved to a separate script: COSMOS_MOON_analysis.R
#### Run that script after this one has written the combined weight files under
####   results/MOFA_with_lipidomics/feature_weights_all_factors/
#######################################################################################################################
