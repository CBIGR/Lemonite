#!/usr/bin/Rscript

setwd('/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results')
base_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/'
DE_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/DE_analysis/'
TFA_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/TFA/'

'/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree_bulk/no_controls/Variability_0.7_5817genes_old/'

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

expression <- paste0(base_dir, 'data/fpkm_gene_expression.csv') 
metabolomics <- paste0(base_dir, 'data/metabolome.csv')
lipids_pos <- paste0(base_dir, 'data/lipidome_pos.csv')
lipids_neg <- paste0(base_dir, 'data/lipidome_neg.csv')

metadata <- paste0(base_dir, 'data/clinical_metadata.csv')
metadata2 <- paste0(base_dir, 'data/Additional_sample_annotations_suppl2.csv')

prior_network <- '/home/borisvdm/Documents/PhD/resources/networks/CollecTRI_network.txt'
perform_TFA <- TRUE
variability_threshold <- 0.7
# HVGs <- 6000

TFs <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/data/lovering_TF_list.txt'


###########################################################################################
#### Build directory structure
##########################################################################################

if (!dir.exists(DE_dir)) {dir.create(DE_dir)}
if (!dir.exists(TFA_dir)) {dir.create(TFA_dir)}
if (!dir.exists(file.path(base_dir, 'results/LemonTree'))) {dir.create(file.path(base_dir, 'results/LemonTree'))}

run_id <- paste0('Variability_', variability_threshold)
base_dir <- paste0(base_dir, 'results/LemonTree/', run_id, '/')

if (!dir.exists(base_dir)) {dir.create(base_dir)}
#if (!dir.exists(file.path(base_dir, 'results'))) {dir.create(file.path(base_dir, 'results'))}
#if (!dir.exists(file.path(base_dir, 'LemonTree'))) {dir.create(file.path(base_dir, 'LemonTree'))}
if (!dir.exists(file.path(base_dir, 'Preprocessing'))) {dir.create(file.path(base_dir, 'Preprocessing'))}
if (!dir.exists(file.path(base_dir, 'Figures'))) {dir.create(file.path(base_dir, 'Figures'))}
if (!dir.exists(file.path(base_dir, 'Lemon_out'))) {dir.create(file.path(base_dir, 'Lemon_out'))}
#if (!dir.exists(file.path(base_dir, 'results/TF-metabolite_interactions'))) {dir.create(file.path(base_dir, 'results/TF-metabolite_interactions'))}

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
write.table(id_ensembl[,c(2,1)], file = paste0(base_dir, '/Preprocessing/ensemble_mapping.txt'), quote=FALSE, sep = '\t', row.names = FALSE)
# RNAseq <- merge(RNAseq, id_ensembl, by.x = 'count', by.y='hgnc_symbol')

# Select protein coding genes
RNAseq <- RNAseq[RNAseq$gene_type == 'protein_coding', ]
RNAseq <- as.data.frame(RNAseq[!duplicated(RNAseq$gene_name), ]) # Only 58 duplicate gene names
rownames(RNAseq) <- RNAseq$gene_name
RNAseq <- RNAseq[, -c(1:4)]
RNAseq <- RNAseq[, colnames(RNAseq) %in% samples]
# Set NAs to zero
RNAseq[is.na(RNAseq)] <- 0


#####################################################################################################
#### Normalization, log transformation and scaling + selection for highly variable genes using DESeq2
#### Also do TFA inference using DecouplR
#####################################################################################################

# Metadata file: https://portal.gdc.cancer.gov/repository?filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22content%22%3A%7B%22field%22%3A%22cases.case_id%22%2C%22value%22%3A%5B%22set_id%3A6i0q7YsBB7vVQnmuUESr%22%5D%7D%2C%22op%22%3A%22IN%22%7D%5D%7D&searchTableTab=cases
# - click 'clinical' and then select tsv format

intersect(colnames(RNAseq), metadata$case_submitter_id)
metadata <- metadata[metadata$case_submitter_id %in% colnames(RNAseq), ]

metadata <- metadata[,c(2,6,12,21,115, 119, 159)]
metadata$age_at_diagnosis <- round(as.numeric(metadata$age_at_diagnosis)/365) # Days are pretty useless

# Combine metadata files
metadata2 <- metadata2[,c(1,2,3,6)]
metadata <- merge(metadata, metadata2, by.x = 'case_submitter_id', by.y = 'case', all.x = TRUE)

DESeq_groups <- metadata
rownames(DESeq_groups) <- metadata$case_submitter_id
DESeq_groups$case_submitter_id <- NULL
DESeq_groups[, -c(3)] <- lapply(DESeq_groups[, -c(3)], factor)

RNAseq <- RNAseq[, colnames(RNAseq) %in% rownames(DESeq_groups)] # Select samples
ord <- match(colnames(RNAseq), rownames(DESeq_groups))
DESeq_groups <- DESeq_groups[ord, ]

write.table(round(RNAseq), '/home/borisvdm/Documents/PhD/Test_LemonIte_nextflow_Wang/data/Counts.tsv', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

# These counts are already fpkm-normalized, but check first
pdf(paste0(base_dir, 'Preprocessing/Normalized_fkpm_expression.pdf'))
log(RNAseq+1) %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()


# Normalization with DESeq2

dds <- DESeqDataSetFromMatrix(countData = (round(RNAseq)+1), colData = DESeq_groups, design = ~ 1) # Round because this should have already been fpkm normalized, but doesn't seem very good
keep <- rowSums(counts(dds) >=10) >= 3 # Dispersion plot looks better with some prefiltering
dds <- dds[keep, ]
dds <- DESeq(dds)
#res <- data.frame(results(dds, contrast=DESeq_contrast))
#write.table(res, './DE_analysis/DE_analysis.txt', sep = '\t', quote=FALSE, row.names=TRUE, header=NA, sep = '\t')
#up <- res[(res$log2FoldChange >=1 & res$padj <= 0.05), ]
#down <- res[(res$log2FoldChange <=1 & res$padj <= 0.05), ]
normcnt <- as.data.frame(counts(dds, normalized=TRUE))

# Log transformation
RNAseq <- log(normcnt)
#RNAseq <- log(RNAseq+1)

# This data should already be normalized, let's have a look
pdf(paste0(base_dir, 'Preprocessing/Normalized_expression.pdf'))
RNAseq %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()


# Select for highly variable genes
M <- RNAseq
vars <- apply(M,1,var)
pdf(paste0(base_dir, 'Preprocessing/expression_variance_histogram.pdf'))
hist(vars, 1100,xlim=c(0,7))
abline(v=variability_threshold, col='red', lwd=3)
dev.off()
length(names(vars[vars>=variability_threshold])) # 5817 with threshold == 0.7
variable_genes <- M[names(vars[vars>variability_threshold]), ]

RNA_preprocessed_noTFA <- as.data.frame(t(scale(t(variable_genes))))

#Add a column with Ensembl ID to RNAseq
RNA_preprocessed_noTFA$symbol <- row.names(RNA_preprocessed_noTFA)
RNA_preprocessed_noTFA <- merge(RNA_preprocessed_noTFA, id_ensembl, by.x='symbol', by.y='gene_name')

RNA_preprocessed_noTFA <- as.data.frame(RNA_preprocessed_noTFA %>% group_by(symbol) %>% filter(row_number()==1))
RNA_preprocessed_noTFA <- RNA_preprocessed_noTFA[, c(1,ncol(RNA_preprocessed_noTFA),2:(ncol(RNA_preprocessed_noTFA)-1))]

write.table(RNA_preprocessed_noTFA, paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression.txt'), sep = '\t', quote=FALSE, row.names=FALSE)


# # Scale the complete dataframe
# RNA_preprocessed <- as.data.frame(t(scale(t(variable_genes))))
# pdf('./Preprocessing/Normalized_expression.pdf') 
# RNA_preprocessed_withTFA %>%
#   gather(Sample, Count) %>%
#   ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
# dev.off()AL356413.1


#####################################################################################################
#### Also do TFA inference using DecouplR and CollecTri
#####################################################################################################

##### TFA inference on log transformed + normalized gene expression table (without HVG selection, we will just re-add these TFs later)
net <- fread(prior_network)
#decoupled <- decouple(mat=RNAseq, net=net, .source='source', .target='target') # RNAseq is the log-transformed df with all protein coding genes
#consensus <- run_consensus(decoupled)
# Transform to wide matrix
TFA_df <- consensus %>%
  filter(statistic == 'consensus') %>%
  pivot_wider(id_cols = 'condition', names_from = 'source',
              values_from = 'score') %>%
  column_to_rownames('condition') %>%
  t %>%
  as.data.frame()

# Visualization
n_tfs <- 50

# Transform to wide matrix
sample_acts_mat <- consensus %>%
  filter(statistic == 'consensus') %>%
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
  diagnosis = c(GBM = "red", control = 'green'),
  multiomic = c(nmf1 = "black", nmf2 = "pink", nmf3 = "brown", `IDH mutant` = "grey"),
  progression_or_recurrence = c(yes = "blue", `not reported` = "black", Blank = 'white', no = 'yellow'),
  gender = c(female = "#1B9E77", male = "#D95F02")
)

plot <- pheatmap(t(sample_acts_mat), border_color = NA, color=my_color, breaks = my_breaks, 
                 annotation_col = DESeq_groups[,c(2,4,6,8)], annotation_colors = ann_colors,
                 cluster_rows = T, show_colnames = FALSE, fontsize_row = 2)
#save(TFA_df, file='./TFA/TFA_df.RData')
save_pheatmap_pdf(plot, paste0(TFA_dir, 'TFA_consensus_',n_tfs,'heatmap_annotated.pdf'))
# 

TFA_df[,-c(1,2)] %>% 
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4)) # These values are already scaled! SO no need for second scaling?


# # Scale the complete dataframe (so with TFA)
#RNA_preprocessed <- as.data.frame(t(scale(t(RNA_preprocessed))))
pdf(paste0(base_dir, 'Preprocessing/Normalized_expression.pdf'))
RNA_preprocessed[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()
write.table(RNA_preprocessed, paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression_TFA.txt'), sep = '\t', quote=FALSE, row.names=FALSE)



## Replace expression data by TFA scores
ind_a <- which(RNA_preprocessed_noTFA$symbol %in% rownames(TFA_df)) # These indices will need to be replaced by TFA
ind_b <- which(rownames(TFA_df) %in% RNA_preprocessed_noTFA$symbol) # THese indices indicate were to find TFA values

## Replace expression data by TFA scores
matched_rows <- match(rownames(variable_genes), rownames(TFA_df))
variable_genes[ind_a, ] <- TFA_df[ind_b, ]
new_genes <- setdiff(rownames(TFA_df), rownames(variable_genes))
variable_genes <- rbind(variable_genes, TFA_df[new_genes,])

# # Scale the complete dataframe
RNA_preprocessed <- as.data.frame(t(scale(t(variable_genes)))) # Variable genes contains unscaled data
write.table(RNA_preprocessed, './LemonTree/LemonPreprocessed_expression_TFA.txt', sep = '\t', quote=FALSE, row.names=FALSE)


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

write.table(abundancies, '/home/borisvdm/Documents/PhD/Test_LemonIte_nextflow_Wang/data/Metabolomics.txt', sep = '\t', quote=FALSE, row.names=TRUE, col.names=NA)

abundancies[is.na(abundancies)] <- 0
#abundancies <- log(abundancies + 1)
abundancies <- as.data.frame(t(pareto_scale(t(abundancies)))) # In principle normalization has already been done...

pdf(paste0(base_dir, 'Preprocessing/Normalized_metabolomics.pdf')) # sample MSM719ME does not follow normalized data pattern
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


pdf(paste0(base_dir, 'Preprocessing/Normalized_lipidomics.pdf'))
lipidomics[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()
# Write to file
write.table(lipidomics, paste0(base_dir, 'Preprocessing/LemonPreprocessed_lipidomics.txt'), sep = '\t', quote=FALSE, row.names=FALSE)

###########################################################################################
#### Create complete dataframe for LemonTree
###########################################################################################

# Add a column with Ensembl ID to RNAseq
RNA_preprocessed$symbol <- row.names(RNA_preprocessed)
RNA_preprocessed <- merge(RNA_preprocessed, id_ensembl, by.x='symbol', by.y='gene_name')
RNA_preprocessed <- as.data.frame(RNA_preprocessed %>% group_by(symbol) %>% filter(row_number()==1))
RNA_preprocessed <- RNA_preprocessed[, c(1,ncol(RNA_preprocessed),2:(ncol(RNA_preprocessed)-1))]
names(RNA_preprocessed)[2] <- 'ensembl_gene_id'
#write.table(RNA_preprocessed, paste0(base_dir, 'Preprocessing/LemonPreprocessed_expression_TFA.txt'), sep = '\t', quote=FALSE, row.names=FALSE)

# Add a second column with metabolomics ID's to abundancy table
abundancies$ensembl_gene_id <- row.names(abundancies); abundancies$symbol <- row.names(abundancies)
abundancies <- abundancies[, c(ncol(abundancies),(ncol(abundancies)-1),1:(ncol(abundancies)-2))]
write.table(abundancies, paste0(base_dir, 'Preprocessing/LemonPreprocessed_metabolomics.txt'), sep = '\t', quote=FALSE, row.names=FALSE)



complete_df <- rbind(RNA_preprocessed, abundancies, lipidomics, fill=TRUE)
# Remove the last row
complete_df <- complete_df[-nrow(complete_df),]

write.table(complete_df, paste0(base_dir, 'Preprocessing/LemonPreprocessed_complete.txt'), sep = '\t', quote=FALSE, row.names=FALSE)
write.table(lipidomics$symbol, paste0(base_dir, 'Preprocessing/lipids.txt'), quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(abundancies$symbol, paste0(base_dir, 'Preprocessing/metabolites.txt'), quote = FALSE, row.names = FALSE, col.names=FALSE)
write.table(DESeq_groups, paste0(base_dir, 'Preprocessing/DESeq_groups.txt'), quote=FALSE, sep='\t', row.names = TRUE)


genes <- RNA_preprocessed$symbol
TFs <- fread(TFs, fill=TRUE)
TFs <- TFs %>% filter(TFs$HGNC_approved_gene_symbol %in% genes)
TFs <- as.data.frame(TFs$HGNC_approved_gene_symbol)
write.table(TFs, file = paste0(base_dir, 'Preprocessing/lovering_TFs.txt'), row.names=FALSE, quote=FALSE, col.names=FALSE)


