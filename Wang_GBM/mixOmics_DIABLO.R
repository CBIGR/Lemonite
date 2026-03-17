#!/usr/bin/Rscript

library(mixOmics)
library(data.table)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(IMIFA)
library(ggplot2)

 
setwd('/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/MixOmics')
base_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/'

########################################################################################
# Load datasets and do some exploration
######################################################################################### 

expression <- paste0(base_dir, 'data/fpkm_gene_expression.csv') 
metabolomics <- paste0(base_dir, 'data/metabolome.csv')
lipids_pos <- paste0(base_dir, 'data/lipidome_pos.csv')
lipids_neg <- paste0(base_dir, 'data/lipidome_neg.csv')

metadata <- paste0(base_dir, 'data/clinical_metadata.csv')
metadata2 <- paste0(base_dir, 'data/Additional_sample_annotations_suppl2.csv')



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

# Standardized preprocessing (align with Preprocessing_and_TFA_noControls_with_Proteomics.R):
# 1) DESeq2 normalization (SKIP for FPKM data - already normalized)
# 2) log transform
# 3) High-variance genes with variance > 0.7
# 4) Z-score per gene

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

pdf(paste0(base_dir, 'results/MixOmics/Normalized_metabolomics.pdf')) # sample MSM719ME does not follow normalized data pattern
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


pdf(paste0(base_dir, 'results/MixOmics/Normalized_lipidomics.pdf'))
lipidomics[, -c(1,2)] %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4))
dev.off()
# Write to file
write.table(lipidomics, paste0(base_dir, './LemonPreprocessed_lipidomics.txt'), sep = '\t', quote=FALSE, row.names=FALSE)

# Drop non-unique symbols
lipidomics <- lipidomics[!duplicated(lipidomics$symbol),]
rownames(lipidomics) <- lipidomics$symbol; lipidomics$symbol <- NULL; lipidomics$ensembl_gene_id <- NULL

###########################################################################################
#### Initial exploratory analyses
###########################################################################################

RNAseq <- as.data.frame(t(RNAseq))
abundancies <- as.data.frame(t(abundancies))
lipidomics <- as.data.frame(t(lipidomics))

rownames(metadata) <- metadata$case_submitter_id
rownames(metadata2) <- metadata2$case
# merge metadata on rownames
metadata_merge <- merge(metadata, metadata2, by = 'row.names'); rownames(metadata_merge) <- metadata_merge$Row.names; metadata_merge$Row.names <- NULL

cols_to_keep <- c('gender', 'multiomic', 'nmf_consensus', 'rna_wang_cancer_cell_2017')
coldata <- metadata_merge[,cols_to_keep]


pca.gene <- pca(RNAseq, ncomp = 10, center = TRUE, scale = TRUE)
pca.mets <- pca(abundancies, ncomp = 10, center = TRUE, scale = TRUE)
pca.lipids <- pca(lipidomics[, -c(1,2)], ncomp = 10, center = TRUE, scale = TRUE)

# 2 components seems to be quite ok
plot(pca.gene)
plot(pca.mets)
plot(pca.lipids)


plotIndiv(pca.gene, comp = c(1, 2), 
          group = coldata$multiomic, 
          ind.names = rownames(coldata), 
          legend = TRUE, title = 'Gene expression PCA comp 1 - 2')
ggsave(paste0(base_dir, 'results/MixOmics/PCA_gene_expression.png'))

plotIndiv(pca.mets, comp = c(1, 2), 
          group = coldata$multiomic, 
          ind.names = rownames(coldata), 
          legend = TRUE, title = 'Metabolomics PCA comp 1 - 2')
ggsave(paste0(base_dir, 'results/MixOmics/PCA_metabolomics.png'))


plotIndiv(pca.lipids, comp = c(1, 2), 
          group = coldata$multiomic, 
          ind.names = rownames(coldata), 
          legend = TRUE, title = 'Lipidomics PCA comp 1 - 2')

ggsave(paste0(base_dir, 'results/MixOmics/PCA_lipidomics.png'))


# Set multi-omics subtypes as Y variable (response)
Y <- as.factor(coldata$multiomic)
summary(Y)

###########################################################################################
#### Initial analysis - pairwise PLS comparisons
###########################################################################################

list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)


# set a list of all the X dataframes
data = list(Transcriptomics = RNAseq,
            Lipidomics = lipidomics, 
            Metabolomics = abundancies)

# Ensure all data is properly scaled for mixOmics (center and scale)
data$Transcriptomics <- scale(data$Transcriptomics, center = TRUE, scale = TRUE)
data$Metabolomics <- scale(data$Metabolomics, center = TRUE, scale = TRUE)  
data$Lipidomics <- scale(data$Lipidomics, center = TRUE, scale = TRUE)

# Convert back to data frames
data$Transcriptomics <- as.data.frame(data$Transcriptomics)
data$Metabolomics <- as.data.frame(data$Metabolomics)
data$Lipidomics <- as.data.frame(data$Lipidomics)

# generate three pairwise PLS models
pls1 <- spls(data[["Lipidomics"]], data[["Transcriptomics"]], 
             keepX = list.keepX, keepY = list.keepY) 
pls2 <- spls(data[["Lipidomics"]], data[["Metabolomics"]], 
             keepX = list.keepX, keepY = list.keepY)
pls3 <- spls(data[["Transcriptomics"]], data[["Metabolomics"]], 
             keepX = list.keepX, keepY = list.keepY)

# plot features of first PLS
plotVar(pls1, cutoff = 0.5, title = "(a) lipidomics vs mRNA", 
        legend = c("lipidomics", "mRNA"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of second PLS
plotVar(pls2, cutoff = 0.5, title = "(b) lipidomics vs metabolomics", 
        legend = c("lipidomics", "metabolomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) mRNA vs metabolomics", 
        legend = c("mRNA", "metabolomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# calculate correlation of miRNA and mRNA
cor(pls1$variates$X, pls1$variates$Y) 


# for square matrix filled with 0.1s
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design

# form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 10, design = design)


# remove features with zero variance from all datasets
library(caret)
nzv_mrna <- caret::nearZeroVar(data$Transcriptomics)
nzv_metab <- caret::nearZeroVar(data$Metabolomics)
nzv_lipid <- caret::nearZeroVar(data$Lipidomics)

if(length(nzv_mrna) > 0) data$Transcriptomics <- data$Transcriptomics[, -nzv_mrna]
if(length(nzv_metab) > 0) data$Metabolomics <- data$Metabolomics[, -nzv_metab]
if(length(nzv_lipid) > 0) data$Lipidomics <- data$Lipidomics[, -nzv_lipid]

cat("After near-zero variance removal:\n")
cat("Transcriptomics features:", ncol(data$Transcriptomics), "\n")
cat("Metabolomics features:", ncol(data$Metabolomics), "\n")
cat("Lipidomics features:", ncol(data$Lipidomics), "\n")

#### Tune the number of features to keep
# set grid of values for each component to test (comprehensive range from low to high using step sizes)
# Generate sequences up to full dataset size with appropriate step sizes
test.keepX = list(
  Lipidomics = seq(5, ncol(data$Lipidomics), by = 50),     # Step size 50 for lipidomics
  Transcriptomics = seq(10, ncol(data$Transcriptomics), by = 500),              # Step size 500 for mRNA
  Metabolomics = seq(5, ncol(data$Metabolomics), by = 20)  # Step size 20 for metabolomics
)

# run component number tuning with repeated CV (or load if already exists)
perf_file <- './perf_diablo.rds'
if (file.exists(perf_file)) {
  cat("Loading existing performance results from", perf_file, "\n")
  perf.diablo <- readRDS(perf_file)
} else {
  cat("Running component number tuning with repeated CV (this may take a while)...\n")
  perf.diablo = perf(basic.diablo.model, validation = 'Mfold',
                     folds = 10, nrepeat = 10, cpus = 6)
  # Save the results
  saveRDS(perf.diablo, perf_file)
  cat("Performance results saved to", perf_file, "\n")
}


pdf(paste0(base_dir, 'results/MixOmics/Diablo_tuning__elbow_plot.pdf'))
plot(perf.diablo) # plot output of tuning
dev.off()


# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 


# # run the feature selection tuning
# library(BiocParallel)
# BPPARAM <- MulticoreParam(workers = 6)  
# 
# # Try tuning with error handling
# tune.diablo <- tryCatch({
#   tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
#                    test.keepX = test.keepX, design = design,
#                    validation = 'Mfold', folds = 5, nrepeat = 1,
#                    dist = "centroids.dist", BPPARAM = BPPARAM)
# }, error = function(e) {
#   cat("Error in tune.block.splsda:", e$message, "\n")
#   cat("Trying with fewer workers...\n")
#   BPPARAM_simple <- SerialParam()  # Fall back to serial processing
#   tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
#                    test.keepX = test.keepX, design = design,
#                    validation = 'Mfold', folds = 5, nrepeat = 1,
#                    dist = "centroids.dist", BPPARAM = BPPARAM_simple)
# })
# list.keepX = tune.diablo$choice.keepX # set the optimal values of features to retain
# list.keepX

# # save tune.diablo
# saveRDS(tune.diablo, file = paste0(base_dir, 'results/MixOmics/tune_diablo.rds'))

# load model
tune.diablo <- readRDS(paste0(base_dir, 'results/MixOmics/tune_diablo.rds'))
list.keepX = tune.diablo$choice.keepX # set the optimal values of features to retain
list.keepX

names(list.keepX)[names(list.keepX) == "lipidomics"]   <- "Lipidomics"
names(list.keepX)[names(list.keepX) == "mRNA"]         <- "Transcriptomics"
names(list.keepX)[names(list.keepX) == "metabolomics"] <- "Metabolomics"


# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)
final.diablo.model$design # design matrix for the final model



for (comp in 1:ncomp) {
  png(paste0(base_dir, 'results/MixOmics/Diablo_plot_comp', comp, '.png'))
  plotDiablo(final.diablo.model, ncomp = comp, title = paste('DIABLO plot comp', comp))
  dev.off()
}


# Sample projection into latent components
png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotIndiv_comp1-3.png'))
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots', comp = c(1,3))
dev.off()

png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotIndiv_comp1-2.png'))
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots', comp = c(1,2))
dev.off()

png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotIndiv_comp2,-3.png'))
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots', comp = c(2,3))
dev.off()




# Correlation between variables and components
png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotVar_com1-2.png'))
plotVar(final.diablo.model, var.names = FALSE, comp = c(1,2),
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'))
dev.off()

# Correlation between variables and components
png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotVar_comp2-3.png'))
plotVar(final.diablo.model, var.names = FALSE, comp = c(2,3),
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'))
dev.off()

# Correlation between variables and components
png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotVar_comp1-3.png'))
plotVar(final.diablo.model, var.names = FALSE, comp = c(1,3),
        style = 'graphics', legend = TRUE,
        pch = c(16, 17, 15), cex = c(2,2,2), 
        col = c('darkorchid', 'brown1', 'lightgreen'))
dev.off()




png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp1.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(1), margin=c(8,20), legend.position = "right", size.legend=0.7)

dev.off()

# component 2
png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp2.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(2), margin=c(8,20), legend.position = "right", size.legend=0.7)
dev.off()

# component 3
png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp3.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(3), margin=c(8,20), legend.position = "right", size.legend=0.7)
dev.off()

# component 4
png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp4.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(4), margin=c(8,20), legend.position = "right", size.legend=0.7)
dev.off()

# component 5
png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp5.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(5), margin=c(8,20), legend.position = "right", size.legend=0.7)
dev.off()

# component 6
png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp6.png'), width = 12, height = 7, units = "in", res = 300)
cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1', 'lightgreen'),
          comp = c(6), margin=c(8,20), legend.position = "right", size.legend=0.7)
dev.off()



# create circosplot for each component
for (comp in 1:ncomp) {
  png(paste0(base_dir, 'results/MixOmics/Diablo_final_circosPlot_comp', comp, '.png'), width = 9, height = 9, units = 'in', res = 600)
  circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
             color.blocks= c('darkorchid', 'brown1', 'lightgreen'),
             color.cor = c("chocolate3","grey20"), size.labels = 1.5, comp = comp, size.variables = 0.6)
  dev.off()
}


save(final.diablo.model, file = './final.diablo.model.RData')

final.diablo.model <- get(load('./final.diablo.model.RData'))


###########################################################################################
#### Save feature weights for all components
###########################################################################################

# Create weights directory
weights_dir <- paste0(base_dir, 'results/MixOmics/feature_weights_all_components')
dir.create(weights_dir, showWarnings = FALSE, recursive = TRUE)

# Get the number of components
n_components <- final.diablo.model$ncomp
for (comp in 1:n_components) {
  cat(sprintf("Saving weights for Component %d\n", comp))
  
  # Transcriptomics weights
  transcriptomics_weights <- as.data.frame(final.diablo.model$loadings$Transcriptomics[, comp])
  colnames(transcriptomics_weights) <- paste0("Component_", comp)
  write.table(transcriptomics_weights, file = paste0(weights_dir, "/Component", comp, "_Transcriptomics_weights.txt"),
              sep = "\t", quote = FALSE)
  
  # Metabolomics weights
  if ("Metabolomics" %in% names(final.diablo.model$loadings)) {
    metab_weights <- as.data.frame(final.diablo.model$loadings$Metabolomics[, comp])
    colnames(metab_weights) <- paste0("Component_", comp)
    write.table(metab_weights, file = paste0(weights_dir, "/Component", comp, "_Metabolomics_weights.txt"),
                sep = "\t", quote = FALSE)
  }
  
  # Lipidomics weights if present
  if ("Lipidomics" %in% names(final.diablo.model$loadings)) {
    lipid_weights <- as.data.frame(final.diablo.model$loadings$Lipidomics[, comp])
    colnames(lipid_weights) <- paste0("Component_", comp)
    write.table(lipid_weights, file = paste0(weights_dir, "/Component", comp, "_Lipidomics_weights.txt"),
                sep = "\t", quote = FALSE)
  }
  
  # Save top features plots
  tryCatch({
    pdf(paste0(weights_dir, "/Component", comp, "_top_loadings.pdf"), width = 10, height = 8)
    plotLoadings(final.diablo.model, comp = comp, contrib = 'max', method = 'median',
                 size.name = 0.7, ndisplay = 5, title = paste('DIABLO Loadings comp', comp))
    dev.off()
  }, error = function(e) {
    cat(sprintf("Could not save loadings plot for component %d\n", comp))
  })
}

cat(sprintf("Feature weights saved in: %s\n", weights_dir))

#########

# Now create barplots that show feature weights on the different components, grouped per omics type

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

# --- 1. Extract loadings ---
loadings_list <- lapply(final.diablo.model$loadings, function(x) x[, 1:4, drop = FALSE])
blocks <- rep(names(loadings_list), sapply(loadings_list, nrow))
feature_names <- unlist(lapply(loadings_list, rownames))
loadings_matrix <- do.call(rbind, loadings_list)

# --- 2. Convert to long format ---
df_long <- as.data.frame(loadings_matrix) %>%
  mutate(
    Feature = feature_names,
    Block = blocks
  ) %>%
  pivot_longer(
    cols = starts_with("Comp"),
    names_to = "Component",
    values_to = "Loading"
  )

# remove block Y
df_long <- df_long %>% filter(Block != "Y")

# --- 3. Keep only features selected for their component ---
df_long_nonzero <- df_long %>% filter(Loading != 0)

# --- 4. Select top 5 features per block per component ---
topN <- 5
df_top <- df_long_nonzero %>%
  group_by(Component, Block) %>%
  slice_max(order_by = abs(Loading), n = topN) %>%
  ungroup()

# --- 5. Reorder features independently per component ---
df_top <- df_top %>%
  group_by(Component) %>%
  mutate(
    Feature_plot = factor(Feature, levels = Feature[order(abs(Loading))])
  ) %>%
  ungroup()

# --- 6. Create plots per Component (same strategy as before) ---
plots <- lapply(unique(df_top$Component), function(comp) {
  
  df_sub <- df_top %>% filter(Component == comp)
  comp_num <- gsub("Comp|comp", "", comp)
  
  ggplot(df_sub, aes(x = Feature_plot, y = Loading, fill = Block)) +
    geom_col(width = 0.8, color = "black", size = 0.5) +
    scale_x_discrete(
      labels = function(x)
        ifelse(nchar(x) > 20, paste0(substr(x, 1, 20), "..."), x)
    ) +
    scale_fill_manual(
      values = c("Lipidomics" = "#F8766D", "Transcriptomics" = "#00BA38", "Metabolomics" = "#619CFF"),
      name = "Omics Type"
    ) +
    theme_bw() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 11,
        margin = margin(t = 10)
      ),
      plot.margin = margin(10, 20, 80, 20),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      y = "Feature Loading",
      title = paste("Component", comp_num)
    )
})

# --- 7. Arrange side-by-side with shared legend ---
final_plot <- wrap_plots(
  plots,
  ncol = length(plots),
  guides = "collect"
) +
  plot_annotation(title = "Top 5 DIABLO Features per Component and Block")

# --- 8. Display ---
final_plot

ggsave(paste0(base_dir, 'results/MixOmics/Top_DIABLO_features_per_component.png'),
       width = 18, height = 7)


# Commented out: vertical version that was overwriting horizontal version
# (Keeping the horizontal orientation as per user request with 45-degree x-axis labels)
# 
# #######
# 
# # Now create a heatmap that shows feature weights on the different components, grouped per omics type
# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(patchwork) 
# 
# # --- 1. Extract loadings ---
# loadings_list <- lapply(final.diablo.model$loadings, function(x) x[, 1:4, drop = FALSE])
# blocks <- rep(names(loadings_list), sapply(loadings_list, nrow))
# feature_names <- unlist(lapply(loadings_list, rownames))
# loadings_matrix <- do.call(rbind, loadings_list)
# 
# # --- 2. Convert to long format ---
# df_long <- as.data.frame(loadings_matrix) %>%
#   mutate(Feature = feature_names,
#          Block = blocks) %>%
#   pivot_longer(
#     cols = starts_with("Comp"),
#     names_to = "Component",
#     values_to = "Loading"
#   )
# 
# # remove block Y
# df_long <- df_long %>% filter(Block != "Y")
# 
# # --- 3. Keep only features selected for their component ---
# df_long_nonzero <- df_long %>% filter(Loading != 0)
# 
# # --- 4. Select top 5 features per block per component ---
# topN <- 5
# df_top <- df_long_nonzero %>%
#   group_by(Component, Block) %>%
#   slice_max(order_by = abs(Loading), n = topN) %>%
#   ungroup()
# 
# # --- 5. Reorder features independently per Component ---
# df_top <- df_top %>%
#   group_by(Component) %>%
#   mutate(Feature_plot = factor(Feature, levels = Feature[order(abs(Loading))])) %>%
#   ungroup()
# 
# # --- 6. Create separate plots per Component ---
# plots <- lapply(unique(df_top$Component), function(comp) {
#   df_sub <- df_top %>% filter(Component == comp)
#   
#   # Extract the component number from 'comp1', 'comp2', etc.
#   comp_num <- gsub("comp", "", comp)
#   
#   ggplot(df_sub, aes(x = Feature_plot, y = Loading, fill = Block)) +
#     geom_col(width = 0.8) +  # thinner bars
#     coord_flip() +
#     scale_x_discrete(position = "top") +   # move y-axis labels to right (top after coord_flip)
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(size = 8, margin = margin(l = 15)),  # adjust margin
#       axis.title.x = element_text(margin = margin(l = 15)),          
#       plot.margin = margin(5, 20, 5, 40),
#       legend.position = "left"  # move legend to left
#     ) +
#     labs(
#       y = "Feature",
#       x = "Feature Loading",
#       title = paste("Component", comp_num)  # Nice title: Component 1, 2, ...
#     )
# })
# 
# 
# 
# # Arrange plots in a single column
# final_plot <- wrap_plots(plots, ncol = 1) + 
#   plot_annotation(title = "Top 5 DIABLO Features per Component and Block")
# 
# # Display plot
# final_plot
# 
# ggsave(paste0(base_dir, 'results/MixOmics/Top_DIABLO_features_per_component.png'),
#        width = 8, height = 12)
# 
# # --- Create version with legend and consistent colors ---
# plots_with_legend <- lapply(unique(df_top$Component), function(comp) {
#   df_sub <- df_top %>% filter(Component == comp)
#   
#   # Extract the component number from 'comp1', 'comp2', etc.
#   comp_num <- gsub("comp", "", comp)
#   
#   ggplot(df_sub, aes(x = Feature_plot, y = Loading, fill = Block)) +
#     geom_col(width = 0.8, color = "black", size = 0.5) +  # thinner bars with black outline
#     coord_flip() +
#     scale_x_discrete(position = "top") +   # move y-axis labels to right (top after coord_flip)
#     scale_fill_manual(
#       values = c("Lipidomics" = "#F8766D", "Transcriptomics" = "#00BA38", "Metabolomics" = "#619CFF"),
#       name = "Omics Type"
#     ) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(size = 8, margin = margin(l = 15)),
#       axis.title.x = element_text(margin = margin(l = 15)),          
#       plot.margin = margin(5, 20, 5, 40),
#       legend.position = "right",
#       legend.title = element_text(size = 10, face = "bold"),
#       legend.text = element_text(size = 9)
#     ) +
#     labs(
#       y = "Feature",
#       x = "Feature Loading",
#       title = paste("Component", comp_num)
#     )
# })
# 
# final_plot_with_legend <- wrap_plots(plots_with_legend, ncol = 1, guides = "collect") + 
#   plot_annotation(title = "Top 5 DIABLO Features per Component and Block") &
#   theme(legend.position = "right")
# 
# ggsave(paste0(base_dir, 'results/MixOmics/Top_DIABLO_features_per_component_withlegend.png'),
#        plot = final_plot_with_legend,
#        width = 10, height = 12)
# 
# cat("Saved both versions:\n")
# cat("  - Top_DIABLO_features_per_component.png (original)\n")
# cat("  - Top_DIABLO_features_per_component_withlegend.png (with legend)\n")




###########################################################################################
#### GSEA Analysis for all components (Supplementary Information)
###########################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
organism <- org.Hs.eg.db

# Source GSEA plotting utilities
source('/home/borisvdm/repo/gsea_plotting_utils.R')

# Create GSEA results directory
gsea_dir <- paste0(base_dir, 'results/MixOmics/GSEA_all_components')
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Get the number of components
n_components <- final.diablo.model$ncomp

# Run GSEA for all components
for (comp in 1:n_components) {
  #comp <- 2
  cat(sprintf("\n=== Processing Component %d ===\n", comp))
  
  # Create component-specific directory
  comp_dir <- paste0(gsea_dir, "/Component_", comp)
  dir.create(comp_dir, showWarnings = FALSE)
  
  # Extract Transcriptomics weights for this component
  comp_weights <- as.data.frame(final.diablo.model$loadings$Transcriptomics[, comp])
  weights <- comp_weights[, 1]
  names(weights) <- rownames(comp_weights)
  
  # Remove zero weights and sort
  weights <- weights[weights != 0]
  
  if (length(weights) == 0) {
    cat(sprintf("  No non-zero weights for component %d, skipping...\n", comp))
    next
  }
  
  weights <- sort(weights, decreasing = TRUE)
  
  # Convert gene symbols to Entrez IDs for KEGG and Reactome
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
    cat(sprintf("  Warning: Could not convert gene symbols to Entrez IDs for component %d\n", comp))
    weights_entrez <- NULL
  })
  
  # GO Enrichment Analysis
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
        OrgDb = organism,
        pAdjustMethod = "BH"
      )
      
      # Save results table
      if (nrow(as.data.frame(gse_result)) > 0) {
        write.table(as.data.frame(gse_result),
                   file = paste0(comp_dir, "/Component", comp, "_gseGO_", db, "_results.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Create and save plot using new visualization
        db_title <- switch(db,
                           "ALL" = "GO All Ontologies",
                           "BP" = "GO Biological Process",
                           "MF" = "GO Molecular Function",
                           "CC" = "GO Cellular Component",
                           db)
        safe_plot_save(gse_result, 
                      paste0(comp_dir, "/Component", comp, "_gseGO_", db, ".png"),
                      db_name = db_title)
        cat(sprintf("  Saved GO %s enrichment\n", db))
      }
    }, error = function(e) {
      cat(sprintf("  No enrichment found for GO %s\n", db))
    })
  }
  
  # KEGG Enrichment Analysis
  if (!is.null(weights_entrez) && length(weights_entrez) > 0) {
    tryCatch({
      kegg_result <- gseKEGG(
        geneList = weights_entrez,
        organism = "hsa",
        minGSSize = 3,
        maxGSSize = 800,
        pvalueCutoff = 0.10,
        verbose = FALSE,
        pAdjustMethod = "BH"
      )
      
      if (nrow(as.data.frame(kegg_result)) > 0) {
        write.table(as.data.frame(kegg_result),
                   file = paste0(comp_dir, "/Component", comp, "_KEGG_results.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        p <- dotplot(kegg_result, showCategory = 10, split = ".sign") +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 10))
        ggsave(paste0(comp_dir, "/Component", comp, "_KEGG.png"),
               plot = p, width = 12, height = 8)
        cat(sprintf("  Saved KEGG enrichment\n"))
      }
    }, error = function(e) {
      cat(sprintf("  No KEGG enrichment found\n"))
    })
    
    # Reactome Enrichment Analysis
    tryCatch({
      reactome_result <- gsePathway(
        geneList = weights_entrez,
        organism = "human",
        nPerm = 1000,
        minGSSize = 10,
        pvalueCutoff = 0.10,
        verbose = FALSE,
        pAdjustMethod = 'BH'
      )
      
      if (nrow(as.data.frame(reactome_result)) > 0) {
        write.table(as.data.frame(reactome_result),
                   file = paste0(comp_dir, "/Component", comp, "_Reactome_results.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        p <- dotplot(reactome_result, showCategory = 10, split = ".sign") +
          facet_grid(. ~ .sign) +
          theme(axis.text.y = element_text(size = 10))
        ggsave(paste0(comp_dir, "/Component", comp, "_Reactome.png"),
               plot = p, width = 12, height = 8)
        cat(sprintf("  Saved Reactome enrichment\n"))
      }
    }, error = function(e) {
      cat(sprintf("  No Reactome enrichment found\n"))
    })
  }
}

cat(sprintf("\n=== GSEA Analysis Complete ===\n"))
cat(sprintf("Results saved in: %s\n", gsea_dir))

# Aggregate per-component GSEA results into a single combined file
agg_list <- list()
idx <- 1
for (comp in 1:n_components) {
  comp_dir <- file.path(gsea_dir, paste0("Component_", comp))

  # GO (BP, MF, CC) - skip 'ALL'
  for (db in c('BP', 'MF', 'CC')) {
    file_go <- file.path(comp_dir, paste0("Component", comp, "_gseGO_", db, "_results.txt"))
    if (file.exists(file_go)) {
      res <- tryCatch(read.table(file_go, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        res <- cbind(Component = comp, Database = db, res)
        agg_list[[idx]] <- res
        idx <- idx + 1
      }
    }
  }

  # KEGG
  file_kegg <- file.path(comp_dir, paste0("Component", comp, "_KEGG_results.txt"))
  if (file.exists(file_kegg)) {
    res <- tryCatch(read.table(file_kegg, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      res <- cbind(Component = comp, Database = 'KEGG', res)
      agg_list[[idx]] <- res
      idx <- idx + 1
    }
  }

  # Reactome
  file_react <- file.path(comp_dir, paste0("Component", comp, "_Reactome_results.txt"))
  if (file.exists(file_react)) {
    res <- tryCatch(read.table(file_react, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(res) && nrow(res) > 0) {
      res <- cbind(Component = comp, Database = 'Reactome', res)
      agg_list[[idx]] <- res
      idx <- idx + 1
    }
  }
}

if (length(agg_list) > 0) {
  agg_all <- do.call(rbind, agg_list)
  cols <- colnames(agg_all)
  rest <- setdiff(cols, c('Component', 'Database'))
  agg_all <- agg_all[, c('Component', 'Database', rest)]
  out_file <- file.path(gsea_dir, 'Combined_GSEA_all_components.txt')
  write.table(agg_all, file = out_file, sep = '\t', row.names = FALSE, quote = FALSE)
  cat(sprintf("Combined GSEA results written to: %s\n", out_file))
} else {
  cat("No per-component GSEA result files found to aggregate.\n")
}



