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
library(dplyr)
library(tibble)
library(grid)

# ----------------------------
# Parameters
# ----------------------------
views <- c("Metabolomics", "Lipidomics", "Transcriptomics")
factors <- 1:6
top_n <- 3

output_file_expression <- "./heatmaps/combined_top3_6factors_expression.pdf"
output_file_weights    <- "./heatmaps/combined_top3_6factors_weights.pdf"

# ----------------------------
# Color-blind friendly palettes
# ----------------------------

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00",
  "#CC79A7", "#999999"
)

factor_labels <- paste0("F", factors)
factor_colors <- setNames(okabe_ito[seq_along(factor_labels)], factor_labels)

view_colors <- c(
  Transcriptomics = "#0072B2",
  Metabolomics   = "#E69F00",
  Lipidomics     = "#009E73"
)

col_fun_expression <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B35806")
)

col_fun_weights <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "#F7F7F7", "#B35806")
)

# ----------------------------
# Extract weights
# ----------------------------
all_weights <- get_weights(model, views = "all", factors = "all")

all_data <- list()
row_view_annotation   <- c()
row_factor_annotation <- c()

# ----------------------------
# Build heatmap data (NO suffixes)
# ----------------------------
for (view in views) {
  
  weight_mat <- all_weights[[view]][, factors, drop = FALSE]
  
  # Identify top features per factor
  top_features_per_factor <- lapply(factors, function(f) {
    names(sort(abs(weight_mat[, f]), decreasing = TRUE))[1:top_n]
  })
  
  names(top_features_per_factor) <- paste0("F", factors)
  
  # Union of all features
  all_features <- unique(unlist(top_features_per_factor))
  
  # Assign each feature to factor with strongest absolute weight
  assigned_factor <- sapply(all_features, function(feat) {
    f <- which.max(abs(weight_mat[feat, ]))
    paste0("F", factors[f])
  })
  
  # Extract expression data
  data_matrix <- get_data(model, views = view)[[view]][[1]]
  if (!is.matrix(data_matrix)) data_matrix <- as.matrix(data_matrix)
  
  heatmap_data <- data_matrix[all_features, , drop = FALSE]
  heatmap_data <- t(scale(t(heatmap_data)))
  
  all_data[[view]] <- heatmap_data
  row_view_annotation   <- c(row_view_annotation, rep(view, length(all_features)))
  row_factor_annotation <- c(row_factor_annotation, assigned_factor)
}

combined_data <- do.call(rbind, all_data)

# ----------------------------
# Row annotations
# ----------------------------
ha_row <- rowAnnotation(
  View   = row_view_annotation,
  Factor = row_factor_annotation,
  col = list(
    View   = view_colors,
    Factor = factor_colors
  )
)

# ----------------------------
# Column annotations
# ----------------------------
sample_metadata <- coldata[colnames(combined_data), ]

ha_col <- HeatmapAnnotation(
  Gender = sample_metadata$gender_original,
  Multiomic = sample_metadata$multiomic_original,
  RNA_subtype = sample_metadata$rna_wang_cancer_cell_2017_original,
  annotation_name_gp = gpar(fontsize = 10)
)

# ----------------------------
# Expression heatmap
# ----------------------------
pdf(output_file_expression, width = 10, height = 14)
ht_expression <- Heatmap(
  combined_data,
  name = "Z-score\n(Expression)",
  col = col_fun_expression,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  cluster_row_slices = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  left_annotation = ha_row,
  top_annotation  = ha_col,
  row_split = factor(row_view_annotation, levels = views),
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 8)
)
draw(ht_expression, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# ----------------------------
# Weights heatmap
# ----------------------------
weights_data <- matrix(
  nrow = nrow(combined_data),
  ncol = length(factors),
  dimnames = list(rownames(combined_data), paste0("Factor", factors))
)

for (i in seq_len(nrow(weights_data))) {
  feat <- rownames(weights_data)[i]
  view <- row_view_annotation[i]
  weights_data[i, ] <- all_weights[[view]][feat, factors]
}

# Row annotation for weights heatmap (only Factor, no View)
ha_row_weights <- rowAnnotation(
  Factor = row_factor_annotation,
  col = list(
    Factor = factor_colors
  )
)

pdf(output_file_weights, width = 10, height = 14)
ht_weights <- Heatmap(
  weights_data,
  name = "Feature\nWeight",
  col = col_fun_weights,
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  show_row_names = TRUE,
  left_annotation = ha_row_weights,
  row_split = factor(row_view_annotation, levels = views),
  row_gap = unit(2, "mm"),
  row_names_gp = gpar(fontsize = 10),
  column_names_rot = 45
)
draw(ht_weights, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

# %%




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

# Create weights directory
weights_dir <- paste0(base_dir, 'results/MOFA_with_lipidomics/feature_weights_all_factors')
dir.create(weights_dir, showWarnings = FALSE, recursive = TRUE)

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
#### Follow-up analysis with COSMOS (for all factors)
#######################################################################################################################

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

weights <- get_weights(model, views = "all", factors = "all")

# Get the number of factors
n_factors <- model@dimensions$K

setwd('./To_COSMOS')

for (fact in 1:n_factors) {
  
  cat(sprintf("\n=== Processing COSMOS for Factor %d ===\n", fact))
  
  # Create factor-specific directory
  fact_dir <- paste0("./Factor_", fact, "_COSMOS")
  dir.create(fact_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get metabolite weights for current factor
  metab_inputs <- weights$metab[, fact]
  
  # Read mapping table from metabolite to HMDB
  metab_to_hmdb <- fread('/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2/Preprocessing/name_map.csv')
  
  # Keep only metabolites present in the data
  metab_to_hmdb <- metab_to_hmdb[metab_to_hmdb$Query %in% names(metab_inputs), ]
  
  # Replace missing HMDB IDs with a placeholder using original metabolite name
  metab_to_hmdb$HMDB[is.na(metab_to_hmdb$HMDB)] <- paste0("No_HMDB_", metab_to_hmdb$Query[is.na(metab_to_hmdb$HMDB)])
  
  # Replace names in metab_inputs with HMDB/placeholder
  common_metabs <- intersect(names(metab_inputs), metab_to_hmdb$Query)
  names(metab_inputs)[match(common_metabs, names(metab_inputs))] <- metab_to_hmdb$HMDB[match(common_metabs, metab_to_hmdb$Query)]
  
  # Load meta-network
  data("meta_network")
  
  # Filter network to metabolites only
  meta_network_metab <- meta_network[grepl("HMDB", meta_network$source) | grepl("HMDB", meta_network$target), ]
  
  # Clean source and target names
  meta_network_metab$source <- gsub("Metab__", "", meta_network_metab$source)
  meta_network_metab$target <- gsub("Metab__", "", meta_network_metab$target)
  meta_network_metab$source <- gsub("_.*", "", meta_network_metab$source)
  meta_network_metab$target <- gsub("_.*", "", meta_network_metab$target)
  
  # Unique metabolites in network
  meta_network_metabs <- unique(c(meta_network_metab$source, meta_network_metab$target))
  
  # Replace NAs in network (if any) with placeholders
  meta_network_metabs[is.na(meta_network_metabs)] <- "No_HMDB_unknown"
  
  # Filter metabolites with weight > 0.2
  metab_inputs_toCosmos <- metab_inputs[abs(metab_inputs) > 0.2]
  
  # Prepare Venn diagram list
  venn_list <- list(
    "COSMOS Meta-network" = meta_network_metabs,
    "Metabolites with weight > 0.2" = names(metab_inputs_toCosmos),
    "Metabolites in dataset" = names(metab_inputs)
  )
  
  # Load ggvenn for plotting
  library(ggvenn)
  
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
