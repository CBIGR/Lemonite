#!/usr/bin/Rscript

setwd('/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/LemonTree/')
base_dir <- '/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/'

# You will need the following (sub)directory structure:
# ./
# -- data
# -- results (in this run this is CollecTri_consensus)
#   ---- Preprocessing
#   ---- TFA
#   ---- LemonTree/Figures
#   ---- TF-metabolite_interactions

# BiocManager::install("MOFA2")
library(DESeq2)
library(data.table)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(IMIFA)
library(ggplot2)
library(decoupleR)
library(pheatmap)
#library(Org.Hs.eg.db)
# library(EnsDb.Hsapiens.v86)
# library(EnsDb.Hsapiens.v79)
library(AnnotationDbi)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

###########################################################################################
#### Input files
###########################################################################################

expression <- paste0(base_dir, 'data/GSE111889_host_tx_counts.tsv') 
metabolomics <- paste0(base_dir, 'data/Metabolomics_complete.txt')
metadata <- paste0(base_dir, 'data/metadata_transriptomics_metabolomics_looseCoupled.txt')
prior_network <- '/home/borisvdm/Documents/PhD/resources/networks/CollecTRI_network.txt'
perform_TFA <- TRUE
variability_threshold <- 0.35
DESeq_contrast <- c('diagnosis', 'UC', 'nonIBD')
metadata_to_include <- c('diagnosis', 'sex', "biopsy_location", 'Age at diagnosis', 'Participant ID') # Covariates to include
TFs <- 'data/lovering_TF_list.txt'

###########################################################################################
#### Read RNAseq data, filter for protein coding genes
###########################################################################################

RNAseq <- fread(expression, header=TRUE, data.table=TRUE)
# Select for protein coding genes, keep gene symbols as rownames
# ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl", host='https://www.ensembl.org')
# id_ensembl <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
#                     filters='hgnc_symbol', values=RNAseq$count, mart=ensembl)

# Use the Ensembl BioMart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = 'https://jan2024.archive.ensembl.org')
# Retrieve all Ensembl Gene IDs, HGNC Symbols, and Gene Biotypes for all possible genes
all_genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
                   mart = ensembl)

write.table(all_genes, file = '/home/borisvdm/Documents/PhD/Lemonite/ensembl_mapping_jan2024.txt', quote=FALSE, sep = '\t', row.names = FALSE)

id_ensembl <- all_genes[all_genes$hgnc_symbol %in% RNAseq$count, ]

# write.table(id_ensembl[,c(2,1)], file = '/home/boris/Documents/PhD/gut_brain/IBD/Lloyd-Price2019/results/LemonTree_CollecTri_consensus/Preprocessing/ensemble_mapping.txt', quote=FALSE, sep = '\t', row.names = FALSE)
#id_ensembl <- fread('/home/boris/Documents/PhD/gut_brain/IBD/Lloyd-Price2019/results/LemonTree_CollecTri_consensus/Preprocessing/ensemble_mapping.txt', header=TRUE, data.table=TRUE)

RNAseq <- merge(RNAseq, id_ensembl, by.x = 'count', by.y='hgnc_symbol')
RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$count), ]) # 5 non-unique gene symbols
rownames(RNAseq_coding) <- RNAseq_coding$count
RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
RNAseq_coding$count <- NULL; RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$ensembl_gene_id <- NULL

#####################################################################################################
#### Normalization, log transformation and scaling + selection for highly variable genes using DESeq2
#### Also do TFA inference using DecouplR
#####################################################################################################

metadata <- fread(metadata, data.table =FALSE); metadata$V1 <- NULL
DESeq_groups <- metadata[,metadata_to_include]
rownames(DESeq_groups) <- metadata$`External ID`; DESeq_groups[] <- lapply(DESeq_groups, factor)

RNAseq <- RNAseq_coding[, colnames(RNAseq_coding) %in% metadata$`External ID`] # Select samples
RNAseq$MSM719ME <- NULL; metadata <- metadata[metadata$`External ID` != 'MSM719ME', ]; DESeq_groups <- DESeq_groups[rownames(DESeq_groups) != 'MSM719ME', ]
RNAseq$HSM5FZB5 <- NULL; metadata <- metadata[metadata$`External ID` != 'HSM5FZB5', ]; DESeq_groups <- DESeq_groups[rownames(DESeq_groups) != 'HSM5FZB5', ]
ord <- match(colnames(RNAseq), rownames(DESeq_groups))
DESeq_groups <- DESeq_groups[ord, ]

dds <- DESeqDataSetFromMatrix(countData = (RNAseq+1), colData = DESeq_groups, design = ~ biopsy_location + diagnosis + sex) 
keep <- rowSums(counts(dds) >=10) >= 3 # Dispersion plot looks better with some prefiltering
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- data.frame(results(dds, contrast=DESeq_contrast))
#write.table(res, './DE_analysis/DE_analysis.txt', sep = '\t', quote=FALSE, row.names=TRUE, header=NA, sep = '\t')
up <- res[(res$log2FoldChange >=1 & res$padj <= 0.05), ]
down <- res[(res$log2FoldChange <=1 & res$padj <= 0.05), ]
normcnt <- as.data.frame(counts(dds, normalized=TRUE))
log_normcnt <- log(normcnt)

# Select for highly variable genes
M <- log_normcnt
vars <- apply(M,1,var)
pdf('./Preprocessing/expression_variance_histogram.pdf')
hist(vars, 1100,xlim=c(0,2))
abline(v=variability_threshold, col='red', lwd=3)
dev.off()
length(names(vars[vars>=variability_threshold])) # 4429
length(names(vars[vars>=0.35])) # 4429
variable_genes <- M[names(vars[vars>variability_threshold]), ]


# # Scale the complete dataframe
# RNA_preprocessed <- as.data.frame(t(scale(t(variable_genes))))
# pdf('./Preprocessing/Normalized_expression.pdf') 
# RNA_preprocessed_withTFA %>%
#   gather(Sample, Count) %>%
#   ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
# dev.off()

##### TFA inference on log transformed + normalized gene expression table (without HVG selection, we will just re-add these TFs later)
net <- fread(prior_network)
decoupled <- decouple(mat=log_normcnt, net=net, .source='source', .target='target')
consensus <- run_consensus(decoupled)
# Transform to wide matrix
TFA_df <- consensus %>%
  dplyr::filter(statistic == 'consensus') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  t %>%
  as.data.frame()

# Visualization
n_tfs <- 50

# Transform to wide matrix
sample_acts_mat <- consensus %>%
  dplyr::filter(statistic == 'consensus') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix()

# Get top tfs with more variable means across clusters
tfs <- consensus %>%
  group_by(source) %>%
  summarise(std = sd(score)) %>%
  arrange(-abs(std)) %>%
  head(n_tfs) %>%
  pull(source)
sample_acts_mat <- sample_acts_mat[,tfs]

# Scale per sample
sample_acts_mat <- scale(sample_acts_mat)

# Choose color palette
palette_length = 100
my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)

my_breaks <- c(seq(-3, 0, length.out=ceiling(palette_length/2) + 1),
               seq(0.05, 3, length.out=floor(palette_length/2)))

ann_colors = list(
  diagnosis = c(UC = "firebrick", nonIBD = "black"),
  sex = c(Female = "#1B9E77", Male = "#D95F02"),
  biopsy_location = c(Colon = "#7570B3", Rectum = "#E7298A")
)

plot <- pheatmap(t(sample_acts_mat), border_color = NA, color=my_color, breaks = my_breaks, 
                 annotation_col = DESeq_groups[,-c(4,5)], annotation_colors = ann_colors,
                 cluster_rows = T, show_colnames = FALSE, fontsize_row = 5)
save(TFA_df, file='./TFA/TFA_df.RData')
save_pheatmap_pdf(plot, paste0('./TFA/TFA_consensus_',n_tfs,'heatmap_annotated.pdf'))

load('./TFA/TFA_df.RData')
# write to file 
write.table(TFA_df, './TFA/TFA_consensus.txt', sep = '\t', quote=FALSE, row.names=TRUE)

# Read Lovering TF list
Lovering_TFs <- fread(paste0(base_dir, TFs), fill=TRUE)
Lovering_TF_list <- Lovering_TFs$HGNC_approved_gene_symbol

# Plot heatmap with TF expression instead
genes_to_plot <- colnames(sample_acts_mat)

# Transpose the matrix
transposed_matrix <- t(normcnt)

# Scale the transposed matrix to get Z-scores (mean = 0, standard deviation = 1)
scaled_matrix <- scale(transposed_matrix)

# Transpose back to the original orientation if required
z_score_matrix <- t(scaled_matrix)
# Select rows that are in genes_to_plot
z_score_matrix <- z_score_matrix[rownames(z_score_matrix) %in% genes_to_plot, ]

plot2 <- pheatmap(z_score_matrix, border_color = NA, color=my_color, breaks = my_breaks, 
                 annotation_col = DESeq_groups[,-c(4,5)], annotation_colors = ann_colors,
                 cluster_rows = T, show_colnames = FALSE, fontsize_row = 5)

save_pheatmap_pdf(plot2, paste0('./TFA/TFAgenes_ZscoredExpression_',n_tfs,'heatmap_annotated.pdf'))

## Create RNA_preprocessed_noTFA: HVGs + Lovering TFs (NO TFA, just expression)
# Identify Lovering TFs not in variable_genes that we need to add
lovering_tfs_not_in_hvg <- setdiff(Lovering_TF_list, rownames(variable_genes))
lovering_tfs_not_in_hvg <- intersect(lovering_tfs_not_in_hvg, rownames(log_normcnt))  # Must exist in log_normcnt

# Start with variable genes (will scale all)
RNA_preprocessed_noTFA <- variable_genes

# Add Lovering TFs not in HVGs (use expression from log_normcnt)
if (length(lovering_tfs_not_in_hvg) > 0) {
  RNA_preprocessed_noTFA <- rbind(RNA_preprocessed_noTFA, log_normcnt[lovering_tfs_not_in_hvg, colnames(RNA_preprocessed_noTFA), drop=FALSE])
}

# Scale all genes (HVGs + Lovering TFs)
RNA_preprocessed_noTFA <- as.data.frame(t(scale(t(RNA_preprocessed_noTFA))))

# Add Ensembl IDs and write to LemonPreprocessed_expression.txt
RNA_preprocessed_noTFA$symbol <- row.names(RNA_preprocessed_noTFA)
RNA_preprocessed_noTFA <- merge(RNA_preprocessed_noTFA, id_ensembl, by.x='symbol', by.y='hgnc_symbol')
RNA_preprocessed_noTFA <- as.data.frame(RNA_preprocessed_noTFA %>% group_by(symbol) %>% dplyr::filter(row_number()==1))
RNA_preprocessed_noTFA <- RNA_preprocessed_noTFA[, c(1,ncol(RNA_preprocessed_noTFA)-1,2:(ncol(RNA_preprocessed_noTFA)-2))]
names(RNA_preprocessed_noTFA)[2] <- 'ensembl_gene_id'

write.table(RNA_preprocessed_noTFA, './Preprocessing/LemonPreprocessed_expression.txt', sep = '\t', quote=FALSE, row.names=FALSE)


## Create RNA_preprocessed_withTFA: HVGs + Lovering TFs (with TFA where available, else expression)
# Filter TFA_df to only include Lovering TFs
TFA_df_lovering <- TFA_df[rownames(TFA_df) %in% Lovering_TF_list, , drop=FALSE]

# Start with variable genes
RNA_preprocessed_withTFA <- variable_genes

# Identify Lovering TFs that exist in variable_genes (will be replaced with TFA if available)
lovering_tfs_in_hvg <- intersect(rownames(variable_genes), rownames(TFA_df_lovering))

# Identify which of the non-HVG Lovering TFs have TFA scores
lovering_tfs_with_tfa <- intersect(lovering_tfs_not_in_hvg, rownames(TFA_df_lovering))

# Identify which of the non-HVG Lovering TFs don't have TFA (will use expression)
lovering_tfs_without_tfa <- setdiff(lovering_tfs_not_in_hvg, rownames(TFA_df_lovering))

# Scale non-Lovering TF genes in variable_genes
non_lovering_tf_genes <- setdiff(rownames(variable_genes), Lovering_TF_list)
if (length(non_lovering_tf_genes) > 0) {
  RNA_preprocessed_withTFA[non_lovering_tf_genes, ] <- t(scale(t(variable_genes[non_lovering_tf_genes, , drop=FALSE])))
}

# Replace Lovering TFs in HVGs with TFA scores (unscaled)
if (length(lovering_tfs_in_hvg) > 0) {
  RNA_preprocessed_withTFA[lovering_tfs_in_hvg, ] <- TFA_df_lovering[lovering_tfs_in_hvg, , drop=FALSE]
}

# Add Lovering TFs not in HVGs that have TFA scores (unscaled)
if (length(lovering_tfs_with_tfa) > 0) {
  RNA_preprocessed_withTFA <- rbind(RNA_preprocessed_withTFA, TFA_df_lovering[lovering_tfs_with_tfa, , drop=FALSE])
}

# Add Lovering TFs not in HVGs that don't have TFA (use scaled expression from log_normcnt)
if (length(lovering_tfs_without_tfa) > 0) {
  lovering_tfs_expr <- log_normcnt[lovering_tfs_without_tfa, colnames(RNA_preprocessed_withTFA), drop=FALSE]
  lovering_tfs_expr_scaled <- t(scale(t(lovering_tfs_expr)))
  RNA_preprocessed_withTFA <- rbind(RNA_preprocessed_withTFA, lovering_tfs_expr_scaled)
}

# Add Ensembl IDs for complete dataframe
RNA_preprocessed_withTFA_ids <- RNA_preprocessed_withTFA
RNA_preprocessed_withTFA_ids$symbol <- row.names(RNA_preprocessed_withTFA_ids)
RNA_preprocessed_withTFA_ids <- merge(RNA_preprocessed_withTFA_ids, id_ensembl, by.x='symbol', by.y='hgnc_symbol')
RNA_preprocessed_withTFA_ids <- as.data.frame(RNA_preprocessed_withTFA_ids %>% group_by(symbol) %>% dplyr::filter(row_number()==1))
RNA_preprocessed_withTFA_ids <- RNA_preprocessed_withTFA_ids[, c(1,ncol(RNA_preprocessed_withTFA_ids)-1,2:(ncol(RNA_preprocessed_withTFA_ids)-2))]
names(RNA_preprocessed_withTFA_ids)[2] <- 'ensembl_gene_id'

pdf('./Preprocessing/Normalized_expression.pdf')
RNA_preprocessed_noTFA[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()


###########################################################################################
#### Preprocessing of metabolomics data
###########################################################################################

abundancies <- fread(metabolomics, data.table = FALSE, header = TRUE)
rownames(abundancies) <- abundancies$V1; abundancies$V1 <- NULL
#abundancies <- abundancies %>% mutate_all(na_if, "") # Change empty entries to NA

rownames(abundancies) <- str_replace_all(rownames(abundancies), ' ', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), '-', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), ':', '_')
rownames(abundancies) <- str_replace_all(rownames(abundancies), '\\+', '_')
# Remove trailing '_' from rownames
rownames(abundancies) <- str_replace_all(rownames(abundancies), '_$', '')

abundancies[is.na(abundancies)] <- 0
abundancies <- log(abundancies + 1)
abundancies <- as.data.frame(t(pareto_scale(t(abundancies)))) # In principle normalization has already been done...

pdf('./Preprocessing/Normalized_metabolomics.pdf') # sample MSM719ME does not follow normalized data pattern
abundancies %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()

###########################################################################################
#### Create complete dataframe for LemonTree
###########################################################################################

# Add a second column with metabolomics ID's to abundancy table
abundancies$ensembl_gene_id <- row.names(abundancies); abundancies$symbol <- row.names(abundancies)
abundancies <- abundancies[, c(ncol(abundancies),(ncol(abundancies)-1),1:(ncol(abundancies)-2))]
write.table(abundancies, './Preprocessing/LemonPreprocessed_metabolomics.txt', sep = '\t', quote=FALSE, row.names=FALSE)

# Create complete dataframe with TFA
complete_df <- rbind(RNA_preprocessed_withTFA_ids, abundancies, fill=TRUE)
write.table(complete_df, './Preprocessing/LemonPreprocessed_complete.txt', sep = '\t', quote=FALSE, row.names=FALSE)
write.table(rownames(abundancies), './Preprocessing/metabolites.txt', quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(DESeq_groups, './Preprocessing/DESeq_groups.txt', quote=FALSE, sep='\t', row.names = TRUE)

# Write Lovering TFs that are present in the final dataset
genes_in_output <- RNA_preprocessed_withTFA_ids$symbol
Lovering_TFs_in_output <- Lovering_TFs %>% dplyr::filter(HGNC_approved_gene_symbol %in% genes_in_output)
write.table(Lovering_TFs_in_output$HGNC_approved_gene_symbol, file = './Preprocessing/lovering_TFs.txt', row.names=FALSE, quote=FALSE, col.names=FALSE)
