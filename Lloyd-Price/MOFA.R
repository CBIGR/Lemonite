#!/usr/bin/Rscript

setwd('/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/MOFA')

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

# Inspecting feature weights tells us more about the most informative features
# Look at features with top loading on factor 2 (dominated by metabolomics data)

plot_weights(model,
             view = "Metabolomics",
             factor = 5,
             nfeatures = 15,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)
ggsave('feature_loading_factor5.png')




for (i in 1:15){
  
  
  
  p <- plot_top_weights(model,
                        view = "Metabolomics",
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
factors <- 1:8
top_n <- 3

output_file_expression <- "./heatmaps/combined_top3_factors1-8_expression_horizontal.pdf"
output_file_weights    <- "./heatmaps/combined_top3_factors1-8_weights_horizontal.pdf"

dir.create("./heatmaps", showWarnings = FALSE, recursive = TRUE)

# ==============================
# Heatmap colors
# ==============================
col_fun_expression <- colorRamp2(
  c(-2, 0, 2),
  c("#2166AC", "#F7F7F7", "#B35806")
)

col_fun_weights <- colorRamp2(
  c(-1, 0, 1),
  c("#2166AC", "#F7F7F7", "#B35806")
)

# ==============================
# Extract weights
# ==============================
all_weights <- get_weights(model, views = "all", factors = "all")

all_data <- list()
feature_view_annotation <- c()

# ==============================
# Build heatmap data
# ==============================
for (view in views) {
  
  weight_mat <- all_weights[[view]][, factors, drop = FALSE]
  
  top_features_per_factor <- lapply(factors, function(f) {
    names(sort(abs(weight_mat[, f]), decreasing = TRUE))[1:top_n]
  })
  
  all_features <- unique(unlist(top_features_per_factor))
  
  data_matrix <- get_data(model, views = view)[[view]][[1]]
  data_matrix <- as.matrix(data_matrix)
  
  heatmap_data <- data_matrix[all_features, , drop = FALSE]
  heatmap_data <- t(scale(t(heatmap_data)))
  
  all_data[[view]] <- heatmap_data
  feature_view_annotation <- c(feature_view_annotation,
                               rep(view, length(all_features)))
}

combined_data <- do.call(rbind, all_data)

# =========================================================
# HORIZONTAL ORIENTATION → transpose matrices
# =========================================================
combined_data_t <- t(combined_data)

# ==============================
# Row annotation (samples ONLY)
# ==============================
sample_metadata <- DESeq_groups[rownames(combined_data_t), ]

ha_row_samples <- rowAnnotation(
  Sex       = sample_metadata$sex,
  Diagnosis = sample_metadata$diagnosis,
  Biopsy    = sample_metadata$biopsy_location,
  col = list(
    Sex       = c(Male = "#1F78B4", Female = "#E31A1C"),
    Diagnosis = c(UC = "#E5C494", nonIBD = "#B3CDE3"),
    Biopsy    = c(Colon = "#CAB2D6", Rectum = "#FFFF99")
  ),
  annotation_name_gp = gpar(fontsize = 10)
)

# ==============================
# Column split by view (optional grouping of features)
# ==============================
column_split_factor <- factor(feature_view_annotation, levels = views)

# ==============================
# Expression heatmap (HORIZONTAL)
# ==============================
pdf(output_file_expression, width = 14, height = 10)

ht_expression <- Heatmap(
  combined_data_t,
  name = "Z-score\n(Expression)",
  col = col_fun_expression,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",       # ← FORCE row names to the left
  show_column_names = TRUE,
  left_annotation = ha_row_samples,
  column_split = column_split_factor,
  column_gap = unit(2, "mm"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 8)
)

draw(ht_expression,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")

dev.off()
message("Expression heatmap saved to ", output_file_expression)

# ==============================
# Weights heatmap
# ==============================
weights_data <- matrix(
  nrow = nrow(combined_data),
  ncol = length(factors),
  dimnames = list(rownames(combined_data), paste0("Factor", factors))
)

for (i in seq_len(nrow(weights_data))) {
  feat <- rownames(weights_data)[i]
  view <- feature_view_annotation[i]
  weights_data[i, ] <- all_weights[[view]][feat, factors]
}

# transpose for horizontal layout
weights_data_t <- t(weights_data)

pdf(output_file_weights, width = 12, height = 6)

ht_weights <- Heatmap(
  weights_data_t,
  name = "Feature\nWeight",
  col = col_fun_weights,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_side = "left",        # ← FORCE row names to the left
  show_column_names = TRUE,
  column_split = column_split_factor,
  column_gap = unit(2, "mm"),
  column_names_rot = 45,
  column_names_gp = gpar(fontsize = 9)
)

draw(
  ht_weights,
  heatmap_legend_side = "bottom",
  padding = unit(c(5, 20, 5, 5), "mm")
)

dev.off()
message("Weights heatmap saved to ", output_file_weights)




###########################################################################################################################
#### GSEA Analysis for all factors (Supplementary Information)
###########################################################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(MOFAcellulaR)
library(enrichplot)

# ------------------------ VISUALIZATION FUNCTION ----------------------------
# Use shared GSEA plotting utilities
source('/home/borisvdm/repo/gsea_plotting_utils.R')
cat('Sourced GSEA plotting utilities from /home/borisvdm/repo/gsea_plotting_utils.R\n')

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
  
  # Read mapping table from metabolite to HMDB
  metab_to_hmdb <- fread('/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum/Preprocessing/name_map.csv')
  
  # Keep only metabolites present in the data
  metab_to_hmdb <- metab_to_hmdb[metab_to_hmdb$Query %in% names(metab_inputs), ]
  
  # Replace missing HMDB IDs with a placeholder using original metabolite name
  metab_to_hmdb$HMDB[is.na(metab_to_hmdb$HMDB)] <- paste0("No_HMDB_", metab_to_hmdb$Query[is.na(metab_to_hmdb$HMDB)])
  
  # Replace names in metab_inputs with HMDB/placeholder
  common_metabs <- intersect(names(metab_inputs), metab_to_hmdb$Query)
  names(metab_inputs)[match(common_metabs, names(metab_inputs))] <- metab_to_hmdb$HMDB[match(common_metabs, metab_to_hmdb$Query)]
  
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
cosmos_summary <- cosmos_summary %>%
  mutate(
    Factor = factor(Factor, levels = seq_len(n_factors)),
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
    category = factor(category, 
                      levels = c("non_overlap_metabolites", "overlap_metabolites"),
                      labels = c("Not in COSMOS+", "In COSMOS+"))
  )

cosmos_barplot <- ggplot(cosmos_summary_long, aes(x = Factor, y = count, fill = category)) +
  geom_col(color = "black", width = 0.7) +
  scale_fill_manual(
    values = c("Not in COSMOS+" = "#FEE0D2", "In COSMOS+" = "#A50F15"),
    name = "Metabolite\nCategory"
  ) +
  geom_text(
    data = cosmos_summary,
    aes(x = Factor, y = total_metabolites, label = total_metabolites, fill = NULL),
    vjust = -0.5, size = 3.5, inherit.aes = FALSE
  ) +
  labs(
    title = "Metabolites with |weight| > 0.2 per factor",
    x = "Factor",
    y = "Metabolites with |weight| > 0.2"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

cosmos_barplot

ggsave(file.path("./", "COSMOS_metabolite_counts_by_factor.png"), plot = cosmos_barplot, width = 8, height = 5)
cat("Saved COSMOS metabolite count summary plot in ./COSMOS_metabolite_counts_by_factor.png\n")
