#!/usr/bin/Rscript

library(mixOmics)
library(data.table)
library(DESeq2)
library(tidyverse)
library(dplyr)
library(biomaRt)
library(IMIFA)
library(ggplot2)
library(caret)

setwd('/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results')
base_dir <- "/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/"

###########################################################################################
#### Input files and metadata
###########################################################################################

expression <- paste0(base_dir, 'data/GSE111889_host_tx_counts.tsv') 
metabolomics <- paste0(base_dir, 'data/Metabolomics_complete.txt')
metadata <- paste0(base_dir, 'data/metadata_transriptomics_metabolomics_looseCoupled.txt')
metadata <- fread(metadata, data.table = FALSE)
# remove samples MSM719ME and HSM5FZB5
metadata <- metadata[metadata$`External ID` != 'MSM719ME', ]; metadata <- metadata[metadata$`External ID` != 'HSM5FZB5', ]
metadata$V1 <- NULL; rownames(metadata) <- metadata$`External ID`; metadata$`External ID` <- NULL
# remove columns 'week_num' and 'metabolomics_id' from metadata
metadata_to_include <- c('diagnosis', 'sex', "biopsy_location", 'Age at diagnosis', 'Participant ID') # Covariates to include
metadata <- metadata[, metadata_to_include]

###########################################################################################
#### Read RNAseq data, filter for protein coding genes
###########################################################################################

RNAseq <- fread(expression, header=TRUE, data.table=TRUE)
# Select for protein coding genes, keep gene symbols as rownames
# Use the Ensembl BioMart
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = 'https://jan2024.archive.ensembl.org')
# # Retrieve all Ensembl Gene IDs, HGNC Symbols, and Gene Biotypes for all possible genes
# all_genes <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id', 'gene_biotype'),
#                    mart = ensembl)
# id_ensembl <- all_genes[all_genes$hgnc_symbol %in% RNAseq$count, ]

id_ensembl <- fread('/home/borisvdm/Documents/PhD/Lemonite/ensembl_mapping_jan2024.txt')

RNAseq <- merge(RNAseq, id_ensembl, by.x = 'count', by.y='hgnc_symbol')
RNAseq_coding <- as.data.frame(RNAseq[!duplicated(RNAseq$count), ]) # 5 non-unique gene symbols
rownames(RNAseq_coding) <- RNAseq_coding$count
RNAseq_coding <- RNAseq_coding[RNAseq_coding$gene_biotype == 'protein_coding', ]
RNAseq_coding$count <- NULL; RNAseq_coding$gene_biotype <- NULL; RNAseq_coding$ensembl_gene_id <- NULL

###########################################################################################
#### Normalization, log transformation and scaling + selection for highly variable genes using DESeq2
###########################################################################################

# Select samples
RNAseq <- RNAseq_coding[, colnames(RNAseq_coding) %in% rownames(metadata)]

# Order metadata to match RNAseq columns
ord <- match(colnames(RNAseq), rownames(metadata))
metadata <- metadata[ord, ]

dds <- DESeqDataSetFromMatrix(countData = (RNAseq+1), colData = metadata, design = ~ biopsy_location + diagnosis + sex) 
keep <- rowSums(counts(dds) >=10) >= 3 # Dispersion plot looks better with some prefiltering
dds <- dds[keep, ]
dds <- DESeq(dds)
normcnt <- as.data.frame(counts(dds, normalized=TRUE))
log_normcnt <- log(normcnt)

# Select for highly variable genes
M <- log_normcnt
vars <- apply(M,1,var)
#pdf('./Preprocessing/expression_variance_histogram.pdf')
hist(vars, 1100,xlim=c(0,2))
abline(v=0.35, col='red', lwd=3)
#dev.off()
length(names(vars[vars>=0.35])) # 4429
variable_genes <- M[names(vars[vars>0.35]), ]

RNAseq <- as.data.frame(t(scale(t(variable_genes))))


###########################################################################################
#### Read metabolomics data, scale dataset
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

abundancies %>%
  gather(Sample, Count) %>%
  ggplot(aes(Sample, Count)) + geom_boxplot() + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size = 4)) #  sample MSM719ME does not follow normalized data pattern -> So it was removed from the analysis


###########################################################################################
#### Initial exploratory analyses
###########################################################################################

RNAseq <- as.data.frame(t(RNAseq))
abundancies <- as.data.frame(t(abundancies))

# remove columns that are completely Na
RNAseq <- RNAseq[, colSums(is.na(RNAseq)) != nrow(RNAseq)]
abundancies <- abundancies[, colSums(is.na(abundancies)) != nrow(abundancies)]

coldata <- metadata


pca.gene <- pca(RNAseq, ncomp = 10, center = TRUE, scale = TRUE)
pca.mets <- pca(abundancies, ncomp = 10, center = TRUE, scale = TRUE)

# 2 components seems to be quite ok
plot(pca.gene)
plot(pca.mets)


plotIndiv(pca.gene, comp = c(1, 2), 
          group = coldata$diagnosis, 
          ind.names = rownames(coldata), 
          legend = TRUE, title = 'Gene expression PCA comp 1 - 2')
ggsave(paste0(base_dir, 'results/MixOmics/PCA_gene_expression.png'))

plotIndiv(pca.mets, comp = c(1, 2), 
          group = coldata$diagnosis, 
          ind.names = rownames(coldata), 
          legend = TRUE, title = 'Metabolomics PCA comp 1 - 2')
ggsave(paste0(base_dir, 'results/MixOmics/PCA_metabolomics.png'))


# Set multi-omics subtypes as Y variable (response)
Y <- as.factor(coldata$diagnosis)
summary(Y)



###########################################################################################
#### Initial analysis - pairwise PLS comparisons
###########################################################################################

list.keepX = c(25, 25) # select arbitrary values of features to keep
list.keepY = c(25, 25)


# set a list of all the X dataframes
data = list(Transcriptomics = RNAseq,
            Metabolomics = abundancies)

# Identify near-zero variance features in the Transcriptomics block
nzv_mRNA <- nearZeroVar(data$Transcriptomics, saveMetrics = TRUE)

# Check features with zero variance
zero_variance_features <- rownames(nzv_mRNA)[nzv_mRNA$zeroVar]
print(zero_variance_features)
# Remove features with zero variance
data$Transcriptomics <- data$Transcriptomics[, !colnames(data$Transcriptomics) %in% zero_variance_features]

# Check dimensions of the filtered Transcriptomics block
dim(data$Transcriptomics)


pls3 <- spls(data[["Transcriptomics"]], data[["Metabolomics"]], 
             keepX = list.keepX, keepY = list.keepY)

# plot features of third PLS
plotVar(pls3, cutoff = 0.5, title = "(c) Transcriptomics vs Metabolomics", 
        legend = c("Transcriptomics", "Metabolomics"), 
        var.names = FALSE, style = 'graphics', 
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'lightgreen'))

# calculate correlation of miRNA and mRNA
cor(pls3$variates$X, pls3$variates$Y) 


# for square matrix filled with 0.1s
design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0 # set diagonal to 0s

design


# form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 10, design = design)

# run component number tuning with repeated CV (or load if already exists)
perf_file <- '/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/MixOmics/perf_diablo.rds'
if (file.exists(perf_file)) {
  cat("Loading existing performance results from", perf_file, "\n")
  perf.diablo <- readRDS(perf_file)
} else {
  cat("Running component number tuning with repeated CV (this may take a while)...\n")
  perf.diablo = perf(basic.diablo.model, validation = 'Mfold',
                     folds = 10, nrepeat = 10, cpus = 8)
  # Save the results
  saveRDS(perf.diablo, perf_file)
  cat("Performance results saved to", perf_file, "\n")
}
#perf.diablo <- readRDS(perf_file)

plot(perf.diablo) # plot output of tuning

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote


pdf(paste0(base_dir, 'results/MixOmics/Diablo_tuning__elbow_plot.pdf'))
plot(perf.diablo) # plot output of tuning
dev.off()


#### Tune the number of features to keep
# set grid of values for each component to test
test.keepX = list (Transcriptomics = c(seq(100, ncol(RNAseq), 500)),
                   Metabolomics = c(seq(5, ncol(abundancies), 5)))

# run the feature selection tuning (or load if already exists)
tune_file <- '/home/borisvdm/Documents/PhD/Lemonite/Lloyd-Price_IBD/results/MixOmics/tune_diablo.rds'
if (file.exists(tune_file)) {
  cat("Loading existing tuning results from", tune_file, "\n")
  tune.diablo <- readRDS(tune_file)
} else {
  cat("Running feature selection tuning (this may take a while)...\n")
  
  # Set up parallel processing
  library(BiocParallel)
  BPPARAM <- MulticoreParam(workers = 8)  # Adjust number of workers as needed
  
  tune.diablo = tune.block.splsda(X = data, Y = Y, ncomp = ncomp,
                                  test.keepX = test.keepX, design = design,
                                  validation = 'Mfold', folds = 5, nrepeat = 1,
                                  dist = "centroids.dist",
                                  BPPARAM = BPPARAM)
  # Save the results
  saveRDS(tune.diablo, tune_file)
  cat("Tuning results saved to", tune_file, "\n")
}

names(tune.diablo$choice.keepX) <- c("Transcriptomics", "Metabolomics")


list.keepX = tune.diablo$choice.keepX # set the optimal values of features to retain
list.keepX


# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp,
                                  keepX = list.keepX, design = design)
final.diablo.model$design # design matrix for the final model

# the features selected to form the first component
selectVar(final.diablo.model, block = 'Transcriptomics', comp = 1)$Transcriptomics$name


plotDiablo(final.diablo.model, ncomp = 1)

plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE,
          title = 'DIABLO Sample Plots')


plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE,
          title = 'DIABLO')

circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)

# save session
save(final.diablo.model, file = './final.diablo.model.RData')


# read in the final model
load('./final.diablo.model.RData')


# the features selected to form the first component
selectVar(final.diablo.model, block = 'Transcriptomics', comp = 1)$Transcriptomics$name 
#selectVar(final.diablo.model, block = 'lipidomics', comp = 1)$lipidomics$name
selectVar(final.diablo.model, block = 'Metabolomics', comp = 1)$Metabolomics$name



for (comp in 1:final.diablo.model$ncomp) {
  pdf(paste0(base_dir, 'results/MixOmics/Diablo_final_plotDiablo_comp', comp, '.pdf'), width = 8, height = 6)
  plotDiablo(final.diablo.model, ncomp = comp)
  dev.off()
}



png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotIndiv.png'))
plotIndiv(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO Sample Plots', comp = c(1,2))
dev.off()


# difficult to interpret
plotArrow(final.diablo.model, ind.names = FALSE, legend = TRUE, 
          title = 'DIABLO')


png(paste0(base_dir, 'results/MixOmics/Diablo_final_plotVar.png'))
plotVar(final.diablo.model, var.names = FALSE, 
        style = 'graphics', legend = TRUE,
        pch = c(16, 17), cex = c(2,2), 
        col = c('darkorchid', 'brown1'))
dev.off()

png(paste0(base_dir, 'results/MixOmics/Diablo_final_circosPlot.png'))
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= c('darkorchid', 'brown1'),
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
dev.off()

# Create cimDIABLO plots for each component
for (comp in 1:final.diablo.model$ncomp) {
  png(paste0(base_dir, 'results/MixOmics/Diablo_final_cimDiablo_comp', comp, '.png'), width = 1200, height = 700)
  cimDiablo(final.diablo.model, color.blocks = c('darkorchid', 'brown1'),
            comp = comp, margin=c(8,20), legend.position = "right", size.legend=0.7)
  dev.off()
}

# Create circosPlot for each component
for (comp in 1:final.diablo.model$ncomp) {
  png(paste0(base_dir, 'results/MixOmics/Diablo_final_circosPlot_comp', comp, '.png'), width = 600, height = 600)
  circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
             color.blocks= c('darkorchid', 'brown1'),
             color.cor = c("chocolate3","grey20"), size.labels = 1.5, comp = comp, size.variables = 0.6)
  dev.off()
}

# network(final.diablo.model, blocks = c(1,2),
#         color.node = c('darkorchid', 'brown1'), cutoff = 0.4)

for (comp in 1:final.diablo.model$ncomp) {
  pdf(paste0(base_dir, 'results/MixOmics/Diablo_loadings_comp', comp, '.pdf'), width = 10, height = 10)
  plotLoadings(final.diablo.model, comp = comp, contrib = 'max', method = 'median', 
               size.name = 1,  # Reduced font size for feature names
               ndisplay = 10, 
               title = paste0('DIABLO Loadings comp ', comp), 
               size.subtitle = 0.8, 
               legend = FALSE,
               size.legend = 0.2,
               name.var = NULL)  # Use default shortened names if names are too long
  dev.off()
}



###########################################################################################
#### GSEA Analysis for all components (Supplementary Information)
###########################################################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
organism <- org.Hs.eg.db

# ------------------------ VISUALIZATION FUNCTION ----------------------------
# Use shared GSEA plotting utilities
source('/home/borisvdm/repo/gsea_plotting_utils.R')
cat('Sourced GSEA plotting utilities from /home/borisvdm/repo/gsea_plotting_utils.R\n')

# Create GSEA results directory
gsea_dir <- paste0(base_dir, 'results/MixOmics/GSEA_all_components')
dir.create(gsea_dir, showWarnings = FALSE, recursive = TRUE)

# Get the number of components
n_components <- final.diablo.model$ncomp

# Run GSEA for all components
for (comp in 1:n_components) {
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
        
        safe_plot_save(kegg_result, paste0(comp_dir, "/Component", comp, "_KEGG.png"), db_name = "KEGG")
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
        
        safe_plot_save(reactome_result, paste0(comp_dir, "/Component", comp, "_Reactome.png"), db_name = "Reactome")
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

  # GO files (ALL, BP, MF, CC)
  for (db in c('ALL', 'BP', 'MF', 'CC')) {
    file_go <- file.path(comp_dir, paste0("Component", comp, "_gseGO_", db, "_results.txt"))
    if (file.exists(file_go)) {
      res <- tryCatch(read.table(file_go, header = TRUE, sep = "\t", stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(res) && nrow(res) > 0) {
        # ensure Component and Database columns are present and populated
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
  # ensure Component is first column and Database second
  cols <- colnames(agg_all)
  rest <- setdiff(cols, c('Component', 'Database'))
  agg_all <- agg_all[, c('Component', 'Database', rest)]
  # write combined file
  out_file <- file.path(gsea_dir, 'Combined_GSEA_all_components.txt')
  write.table(agg_all, file = out_file, sep = '\t', row.names = FALSE, quote = FALSE)
  cat(sprintf("Combined GSEA results written to: %s\n", out_file))
} else {
  cat("No per-component GSEA result files found to aggregate.\n")
}


###########################################################################################
#### Extract and save feature weights for all components (Supplementary Information)
###########################################################################################

# Get the number of components
n_components <- final.diablo.model$ncomp

# Initialize lists to store all feature weights
all_mrna_weights <- list()
all_metabolomics_weights <- list()

# Extract weights for each component
for (comp in 1:n_components) {
  # Transcriptomics weights
  transcriptomics_weights <- as.data.frame(final.diablo.model$loadings$Transcriptomics[, comp])
  colnames(transcriptomics_weights) <- paste0("Component_", comp, "_weight")
  transcriptomics_weights$Feature <- rownames(transcriptomics_weights)
  # Keep only non-zero weights
  transcriptomics_weights <- transcriptomics_weights[transcriptomics_weights[, 1] != 0, ]
  transcriptomics_weights <- transcriptomics_weights[order(abs(transcriptomics_weights[, 1]), decreasing = TRUE), ]
  all_mrna_weights[[comp]] <- transcriptomics_weights
  
  # Metabolomics weights
  metab_weights <- as.data.frame(final.diablo.model$loadings$Metabolomics[, comp])
  colnames(metab_weights) <- paste0("Component_", comp, "_weight")
  metab_weights$Feature <- rownames(metab_weights)
  # Keep only non-zero weights
  metab_weights <- metab_weights[metab_weights[, 1] != 0, ]
  metab_weights <- metab_weights[order(abs(metab_weights[, 1]), decreasing = TRUE), ]
  all_metabolomics_weights[[comp]] <- metab_weights
}

# Create comprehensive data frames with all components
# For mRNA features
mrna_all_components <- data.frame(Feature = character(), stringsAsFactors = FALSE)
for (comp in 1:n_components) {
  comp_data <- all_mrna_weights[[comp]]
  if (nrow(comp_data) > 0) {
    colnames(comp_data) <- c(paste0("Component_", comp), "Feature")
    if (comp == 1) {
      mrna_all_components <- comp_data
    } else {
      mrna_all_components <- merge(mrna_all_components, comp_data, by = "Feature", all = TRUE)
    }
  }
}
# Reorder columns to have Feature first
mrna_all_components <- mrna_all_components[, c("Feature", grep("Component", colnames(mrna_all_components), value = TRUE))]
# Replace NA with 0
mrna_all_components[is.na(mrna_all_components)] <- 0

# For metabolomics features
metab_all_components <- data.frame(Feature = character(), stringsAsFactors = FALSE)
for (comp in 1:n_components) {
  comp_data <- all_metabolomics_weights[[comp]]
  if (nrow(comp_data) > 0) {
    colnames(comp_data) <- c(paste0("Component_", comp), "Feature")
    if (comp == 1) {
      metab_all_components <- comp_data
    } else {
      metab_all_components <- merge(metab_all_components, comp_data, by = "Feature", all = TRUE)
    }
  }
}
# Reorder columns to have Feature first
metab_all_components <- metab_all_components[, c("Feature", grep("Component", colnames(metab_all_components), value = TRUE))]
# Replace NA with 0
metab_all_components[is.na(metab_all_components)] <- 0

# Save to files
write.table(mrna_all_components, 
            file = paste0(base_dir, "results/MixOmics/Supplementary_mRNA_feature_weights_all_components.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(metab_all_components,
            file = paste0(base_dir, "results/MixOmics/Supplementary_metabolomics_feature_weights_all_components.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Print summary
cat("\n=== Feature Weights Summary ===\n")
cat(sprintf("Total mRNA features with non-zero weights: %d\n", nrow(mrna_all_components)))
cat(sprintf("Total metabolomics features with non-zero weights: %d\n", nrow(metab_all_components)))
cat(sprintf("\nFiles saved:\n"))
cat(sprintf("  - %s\n", paste0(base_dir, "results/MixOmics/Supplementary_mRNA_feature_weights_all_components.txt")))
cat(sprintf("  - %s\n", paste0(base_dir, "results/MixOmics/Supplementary_metabolomics_feature_weights_all_components.txt")))


###########################################################################################
#### Create Top DIABLO Features visualization 
###########################################################################################

library(patchwork)

# --- 1. Extract loadings ---
loadings_list <- lapply(final.diablo.model$loadings, function(x) x[, 1:n_components, drop = FALSE])
blocks <- rep(names(loadings_list), sapply(loadings_list, nrow))
feature_names <- unlist(lapply(loadings_list, rownames))
loadings_matrix <- do.call(rbind, loadings_list)

# --- 2. Convert to long format ---
df_long <- as.data.frame(loadings_matrix) %>%
  mutate(Feature = feature_names,
         Block = blocks) %>%
  pivot_longer(
    cols = starts_with("comp"),
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

# --- 5. Reorder features independently per Component ---
df_top <- df_top %>%
  group_by(Component) %>%
  mutate(Feature_plot = factor(Feature, levels = Feature[order(abs(Loading))])) %>%
  ungroup()

# --- 6. Create separate plots per Component (matching Wang_GBM colors) ---
plots <- lapply(unique(df_top$Component), function(comp) {
  df_sub <- df_top %>% filter(Component == comp)
  
  # Extract the component number from 'comp1', 'comp2', etc.
  comp_num <- gsub("comp", "", comp)
  
  ggplot(df_sub, aes(x = Feature_plot, y = Loading, fill = Block)) +
    geom_col(width = 0.8) +
    coord_flip() +
    scale_x_discrete(position = "top") +
    scale_fill_manual(
      values = c("Transcriptomics" = "#00BA38", "Metabolomics" = "#619CFF"),
      name = "Omics Type"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = 8, margin = margin(l = 15)),
      axis.title.x = element_text(margin = margin(l = 15)),
      plot.margin = margin(5, 20, 5, 40),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    ) +
    labs(
      y = "Feature",
      x = "Feature Loading",
      title = paste("Component", comp_num)
    )
})

# Arrange plots in a single column with shared legend
final_plot <- wrap_plots(plots, ncol = 1, guides = "collect") + 
  plot_annotation(title = "Top 5 DIABLO Features per Component and Block") &
  theme(legend.position = "right")

# Display plot
final_plot

ggsave(paste0(base_dir, 'results/MixOmics/Top_DIABLO_features_per_component_withlegend.png'),
       plot = final_plot,
       width = 10, height = 12)

cat("\n=== Top DIABLO Features Visualization Created ===\n")
cat(sprintf("Saved: %s\n", paste0(base_dir, 'results/MixOmics/Top_DIABLO_features_per_component_withlegend.png')))
cat("Colors: Transcriptomics (green), Metabolomics (blue)\n")


