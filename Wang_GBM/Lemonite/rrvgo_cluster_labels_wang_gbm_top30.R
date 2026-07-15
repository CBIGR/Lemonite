source("/home/borisvdm/repo/LemonIte/Wang_GBM/Lemonite/rrvgo_cluster_labels.R")

enrichment_path <- "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum/Enrichment/Modules_gsea/module_members/Gsea_all_enriched_pathways_up_per_module.csv"
cluster_path <- "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum/Networks/megaGO_exploration/cluster_assignments_comparison.csv"
output_dir <- "/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/LemonTree/noProteomics_percentile2_divide_by_sum/Networks/megaGO_exploration/top_30"

cluster_column <- "top_30"
ontology <- "BP"
top_n <- 30L
rrvgo_threshold <- 0.7
rrvgo_method <- "Rel"
worker_cores <- 8L

results <- run_rrvgo_cluster_labels(
  enrichment_path = enrichment_path,
  cluster_path = cluster_path,
  output_dir = output_dir,
  cluster_column = cluster_column,
  ontology = ontology,
  top_n = top_n,
  rrvgo_threshold = rrvgo_threshold,
  rrvgo_method = rrvgo_method,
  worker_cores = worker_cores
)

results$cluster_labels
head(results$module_labels)
head(results$reduced_terms_audit)