#!/usr/bin/Rscript

###########################################################################################
#### COSMOS+ / MOON Analysis — standalone script
####
#### Prerequisites: run MOFA_with_lipidomics.R first so that the combined feature-weight
#### files exist under results/MOFA_with_lipidomics/feature_weights_all_factors/:
####   Combined_expression_weights_all_factors.txt
####   Combined_metabolomics_weights_all_factors.txt
####   Combined_lipidomics_weights_all_factors.txt
###########################################################################################

base_dir    <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/'
weights_dir <- paste0(base_dir, 'results/MOFA_with_lipidomics/feature_weights_all_factors/')
cosmos_dir  <- paste0(base_dir, 'results/MOFA_with_lipidomics/To_COSMOS/')

dir.create(cosmos_dir, showWarnings = FALSE, recursive = TRUE)
setwd(cosmos_dir)

library(cosmosR)
library(liana)
library(decoupleR)
library(dplyr)
library(reshape2)
library(GSEABase)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggvenn)
library(stringr)
library(data.table)
library(igraph)
library(ggraph)
library(ggrepel)
library(enrichR)

source('/home/borisvdm/repo/gsea_and_enrichr_plotting_utils.R')

###########################################################################################
#### Tunable parameters
####
#### reduce_solution_network() has no default cutoff. The official MOON tutorial uses 1.5;
#### the function's own doc example uses 0.4. The right value is data-specific: it must sit
#### below the achievable MOON score range of the nodes you want to keep. In this dataset
#### metabolite MOON scores top out at ~1.43, so the tutorial's 1.5 prunes ALL metabolites.
#### A cutoff of 1.0 retains the strongest metabolites while staying reasonably stringent.
###########################################################################################
solution_cutoff_rec_to_TFmetab <- 1.0   # run 1: Receptors -> TF + Metabolites
solution_cutoff_TF_lig         <- 1.5   # run 2: TF -> Ligands (gene-only, unchanged)


###########################################################################################
#### Load feature weights from files
###########################################################################################

expr_wide  <- fread(paste0(weights_dir, 'Combined_expression_weights_all_factors.txt'),
                    data.table = FALSE)
metab_wide <- fread(paste0(weights_dir, 'Combined_metabolomics_weights_all_factors.txt'),
                    data.table = FALSE)
lipid_wide <- fread(paste0(weights_dir, 'Combined_lipidomics_weights_all_factors.txt'),
                    data.table = FALSE)

# Convert to feature × factor matrices (mirrors get_weights() list output)
.to_matrix <- function(df) {
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- df[[1]]
  mat
}

weights <- list(
  Transcriptomics = .to_matrix(expr_wide),
  Metabolomics    = .to_matrix(metab_wide),
  Lipidomics      = .to_matrix(lipid_wide)
)

n_factors    <- ncol(weights$Transcriptomics)
factor_cols  <- colnames(weights$Transcriptomics)   # e.g. "Factor1" ... "Factor15"

cat(sprintf("Loaded weights: %d factors, %d transcriptomic features, %d metabolites, %d lipids\n",
            n_factors,
            nrow(weights$Transcriptomics),
            nrow(weights$Metabolomics),
            nrow(weights$Lipidomics)))


###########################################################################################
#### Shared setup (runs once)
###########################################################################################

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

# Genes expressed in the MOFA model (clean feature names, no _Transcriptomics suffix)
expressed_genes <- setNames(
  rep(1, nrow(weights$Transcriptomics)),
  gsub("_Transcriptomics$", "", rownames(weights$Transcriptomics))
)

# RNA_all: genes × factors matrix for cross-factor activity heatmaps
RNA_all <- as.data.frame(weights$Transcriptomics)
rownames(RNA_all) <- gsub("_Transcriptomics$", "", rownames(RNA_all))

# HMDB ID mapper
data("HMDB_mapper_vec")

# Metabolite name → HMDB mapping table
metab_to_hmdb_moon <- fread(
  paste0(base_dir, 'results/LemonTree/noProteomics_percentile2_divide_by_sum/Preprocessing/name_map.csv')
)

# Dorothea PKN in source/interaction/target format (reused each iteration)
dorothea_PKN_base        <- dorothea_df[, c(1, 3, 2)]
names(dorothea_PKN_base)[2] <- "interaction"


###########################################################################################
#### COSMOS meta-network overlap / Venn diagrams (per factor)
###########################################################################################

data("meta_network")
meta_network_metab <- meta_network[
  grepl("HMDB", meta_network$source) | grepl("HMDB", meta_network$target), ]
meta_network_metab$source <- gsub("_.*", "", gsub("Metab__", "", meta_network_metab$source))
meta_network_metab$target <- gsub("_.*", "", gsub("Metab__", "", meta_network_metab$target))
meta_network_metabs <- unique(c(meta_network_metab$source, meta_network_metab$target))
meta_network_metabs[is.na(meta_network_metabs)] <- "No_HMDB_unknown"

cosmos_summary <- data.frame(
  Factor              = integer(0),
  total_metabolites   = integer(0),
  overlap_metabolites = integer(0)
)

for (fact in 1:n_factors) {

  cat(sprintf("\n=== Processing COSMOS for Factor %d ===\n", fact))

  fact_dir <- paste0("./Factor_", fact, "_COSMOS")
  dir.create(fact_dir, showWarnings = FALSE, recursive = TRUE)

  metab_inputs <- weights$Metabolomics[, fact]

  metab_to_hmdb <- metab_to_hmdb_moon[metab_to_hmdb_moon$Query %in% names(metab_inputs), ]
  metab_to_hmdb$HMDB[is.na(metab_to_hmdb$HMDB)] <-
    paste0("No_HMDB_", metab_to_hmdb$Query[is.na(metab_to_hmdb$HMDB)])

  common_metabs <- intersect(names(metab_inputs), metab_to_hmdb$Query)
  names(metab_inputs)[match(common_metabs, names(metab_inputs))] <-
    metab_to_hmdb$HMDB[match(common_metabs, metab_to_hmdb$Query)]

  metab_inputs_toCosmos <- metab_inputs[abs(metab_inputs) > 0.2]
  total_metabs    <- length(metab_inputs_toCosmos)
  overlap_metabs  <- sum(names(metab_inputs_toCosmos) %in% meta_network_metabs)
  cosmos_summary  <- rbind(cosmos_summary,
                           data.frame(Factor = fact,
                                      total_metabolites   = total_metabs,
                                      overlap_metabolites = overlap_metabs))

  venn_list <- list(
    "COSMOS Meta-network"           = meta_network_metabs,
    "Metabolites with weight > 0.2" = names(metab_inputs_toCosmos),
    "Metabolites in dataset"        = names(metab_inputs)
  )
  cb_colors <- c("#E69F00", "#56B4E9", "#009E73")
  p <- ggvenn(venn_list, fill_color = cb_colors, stroke_size = 0.5,
              set_name_size = 5, text_size = 4, show_percentage = FALSE)
  ggsave(file.path(fact_dir, 'venn_diagram.png'), plot = p, bg = 'white')
  cat(sprintf("  Saved Venn diagram in %s\n", fact_dir))
}

# NOTE: `weights` and `n_factors` are already defined above from the loaded weight
# files; this standalone script has no live `model` object, so we reuse them as-is.

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
  filter(Factor %in% c("Factor 1",'Factor 2', "Factor 3", "Factor 4", "Factor 6"))

cosmos_summary_123 <- cosmos_summary %>%
  filter(Factor %in% c("Factor 1", 'Factor 2', "Factor 3", "Factor 4", "Factor 6"))


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
  file.path("./", "COSMOS_metabolite_counts_by_factor_1-2-3-4-6.png"),
  plot = cosmos_barplot_123,
  width = 8,
  height = 5,
  bg = "white"
)


###########################################################################################
#### Cross-factor LR and TF activity heatmaps (computed once)
###########################################################################################

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


###########################################################################################
#### Per-factor COSMOS+ / MOON loop
###########################################################################################

for (selected_factor in 1:n_factors) {

  cat(sprintf("\n=== COSMOS+ / MOON: Factor %d / %d ===\n", selected_factor, n_factors))

  fact_dir <- paste0("./Factor_", selected_factor, "_MOON")
  dir.create(fact_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- 1. Extract weights for the current factor ----
  RNA <- data.frame(weights$Transcriptomics[, selected_factor])
  rownames(RNA) <- gsub("_Transcriptomics$", "", rownames(RNA))

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
  names(RNA_input) <- gsub("_Transcriptomics$", "", names(RNA_input))

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

  hmdb_map_fact <- metab_to_hmdb_moon[metab_to_hmdb_moon$Query %in% names(metab_raw), ]
  hmdb_map_fact$HMDB[is.na(hmdb_map_fact$HMDB)] <-
    paste0("No_HMDB_", hmdb_map_fact$Query[is.na(hmdb_map_fact$HMDB)])

  common_metabs_moon <- intersect(names(metab_raw), hmdb_map_fact$Query)
  metab_inputs_moon  <- metab_raw[common_metabs_moon]
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

  # ---- 6. Build MOFA weight annotation table from loaded weight files ----
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
  MOFA_weights$Nodes <- gsub("_Transcriptomics$", "", MOFA_weights$Nodes)
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
      plot(density(moon_res_score$score, na.rm = TRUE), main = "MOON scores: Rec \u2192 TF/Metab")
      abline(v = 1); abline(v = -1)
    } else {
      plot.new(); title(main = "MOON scores: Rec \u2192 TF/Metab (insufficient data)")
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
      plot(density(moon_res_2_score$score, na.rm = TRUE), main = "MOON scores: TF \u2192 Ligands")
      abline(v = 1); abline(v = -1)
    } else {
      plot.new(); title(main = "MOON scores: TF \u2192 Ligands (insufficient data)")
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
    # Ensure every node referenced in edges has an ATT row
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
      ggtitle(paste0("MOON network \u2014 Factor ", selected_factor))

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

    # Dotplots for GO BP and KEGG (top 20 by adjusted p-value), shared plotting style
    for (db in c("GO_Biological_Process_2023", "KEGG_2021_Human")) {
      df <- enrich_res[[db]]
      if (!is.null(df) && nrow(df) > 0) {
        p_enr <- plot_enrichr(df, top_n = 20, numChar = 60,
                               title = paste0(db, " \u2014 Factor ", selected_factor))
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
