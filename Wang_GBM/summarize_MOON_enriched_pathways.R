#!/usr/bin/Rscript

###########################################################################################
#### Summarize MOON results per MOFA factor: enriched-pathway tables + HTML overview
####
#### Reads the per-factor outputs written by COSMOS_MOON_analysis.R and produces:
####   MOON_enriched_pathways_per_factor.csv   - tidy long table, significant terms only
####   MOON_enriched_pathways_per_factor.xlsx  - one sheet per factor (+ "All_factors")
####   MOON_overview.html                      - single self-contained report with all
####                                             plots (base64-embedded) + pathway tables
####
#### Self-contained HTML: PDF plots are rasterized to PNG with `pdftoppm` (poppler) and
#### embedded as base64, so the .html can be opened/shared without any side-car files.
###########################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
  library(base64enc)
})

# The To_COSMOS output directory can be passed as the first command-line argument;
# it defaults to the Wang_GBM analysis. Any dataset with the same per-factor layout
# (Factor_N_MOON / Factor_N_COSMOS) works without further changes.
cosmos_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(cosmos_dir) || !nzchar(cosmos_dir)) {
  cosmos_dir <- '/home/borisvdm/Documents/PhD/thesis_Mirte/Wang2021/results/MOFA_with_lipidomics/To_COSMOS/'
}
cosmos_dir <- paste0(sub('/*$', '', cosmos_dir), '/')   # ensure single trailing slash
if (!dir.exists(cosmos_dir)) stop("cosmos_dir does not exist: ", cosmos_dir)
cat(sprintf("Using cosmos_dir: %s\n", cosmos_dir))

padj_cutoff <- 0.05   # "enriched" = Adjusted.P.value below this
top_n_html  <- 25     # max pathway rows shown per (factor, database) block in the HTML

enrich_dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human",
                "Reactome_2022", "MSigDB_Hallmark_2020")

# Discover factor directories present on disk
fact_dirs <- list.dirs(cosmos_dir, recursive = FALSE, full.names = FALSE)
fact_dirs <- fact_dirs[grepl("^Factor_[0-9]+_MOON$", fact_dirs)]
factors   <- sort(unique(as.integer(gsub("Factor_([0-9]+)_MOON", "\\1", fact_dirs))))
cat(sprintf("Found %d factor MOON directories.\n", length(factors)))


###########################################################################################
#### 1. Aggregate significantly enriched pathways
###########################################################################################

summary_all <- data.frame()
for (fact in factors) {
  fact_dir <- file.path(cosmos_dir, paste0("Factor_", fact, "_MOON"))
  for (db in enrich_dbs) {
    f <- file.path(fact_dir, paste0("enrichr_", gsub("[^A-Za-z0-9]", "_", db), ".csv"))
    if (!file.exists(f)) next
    df <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(df) == 0) next
    df <- df[!is.na(df$Adjusted.P.value) & df$Adjusted.P.value < padj_cutoff, ]
    if (nrow(df) == 0) next
    df$Factor   <- fact
    df$Database <- db
    df$nGenes   <- ifelse(nzchar(df$Genes), lengths(strsplit(df$Genes, ";")), 0L)
    keep <- c("Factor", "Database", "Term", "Overlap", "nGenes",
              "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
    summary_all <- rbind(summary_all, df[, intersect(keep, colnames(df))])
  }
}

if (nrow(summary_all) == 0) stop("No significantly enriched pathways found across any factor.")

summary_all <- summary_all %>%
  arrange(Factor, Database, Adjusted.P.value, desc(Combined.Score))

csv_out <- file.path(cosmos_dir, "MOON_enriched_pathways_per_factor.csv")
write.csv(summary_all, csv_out, row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "All_factors")
writeData(wb, "All_factors", summary_all, withFilter = TRUE)
freezePane(wb, "All_factors", firstRow = TRUE)
for (fact in sort(unique(summary_all$Factor))) {
  sht <- paste0("Factor_", fact)
  addWorksheet(wb, sht)
  writeData(wb, sht, summary_all[summary_all$Factor == fact, ], withFilter = TRUE)
  freezePane(wb, sht, firstRow = TRUE)
}
xlsx_out <- file.path(cosmos_dir, "MOON_enriched_pathways_per_factor.xlsx")
saveWorkbook(wb, xlsx_out, overwrite = TRUE)


###########################################################################################
#### 2. Build the single self-contained HTML overview
###########################################################################################

has_pdftoppm <- nzchar(Sys.which("pdftoppm"))
if (!has_pdftoppm)
  message("[!] 'pdftoppm' not found; PDF plots will be linked instead of embedded.")

# Convert a PDF (first page) to a base64-embedded <img>, or fall back to a link.
embed_plot <- function(path, alt = "", max_width = "900px") {
  if (!file.exists(path)) return("")
  ext <- tolower(tools::file_ext(path))
  png_path <- NA_character_
  if (ext == "png") {
    png_path <- path
  } else if (ext == "pdf" && has_pdftoppm) {
    tmp <- tempfile()
    # -r 130 dpi, first page only; pdftoppm appends "-1.png"
    system2("pdftoppm", c("-png", "-r", "130", "-f", "1", "-l", "1",
                          shQuote(path), shQuote(tmp)), stdout = FALSE, stderr = FALSE)
    cand <- paste0(tmp, c("-1.png", "-01.png", ".png"))
    cand <- cand[file.exists(cand)]
    if (length(cand)) png_path <- cand[1]
  }
  if (is.na(png_path) || !file.exists(png_path)) {
    return(sprintf('<p><a href="%s">%s (open plot)</a></p>',
                   normalizePath(path), htmlesc(alt)))
  }
  uri <- dataURI(file = png_path, mime = "image/png")
  sprintf('<figure><img src="%s" alt="%s" style="max-width:%s;width:100%%;height:auto;border:1px solid #ddd;border-radius:6px"/><figcaption>%s</figcaption></figure>',
          uri, htmlesc(alt), max_width, htmlesc(alt))
}

htmlesc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

# Render a data.frame as an HTML table
df_to_table <- function(df, numeric_signif = TRUE) {
  if (is.null(df) || nrow(df) == 0) return("<p><em>No rows.</em></p>")
  fmt <- function(v) {
    if (is.numeric(v) && numeric_signif) signif(v, 3) else v
  }
  hdr <- paste0("<tr>", paste0("<th>", htmlesc(colnames(df)), "</th>", collapse = ""), "</tr>")
  rows <- apply(df, 1, function(r) {
    paste0("<tr>", paste0("<td>", htmlesc(sapply(r, function(x) x)), "</td>", collapse = ""), "</tr>")
  })
  paste0("<table>", hdr, paste0(rows, collapse = ""), "</table>")
}

# How many metabolite nodes survived into each factor's final network (the headline metric)
metab_count <- function(fact) {
  f <- file.path(cosmos_dir, paste0("Factor_", fact, "_MOON"), "combined_ATT_moon.csv")
  if (!file.exists(f)) return(NA_integer_)
  att <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(att) || !"Nodes" %in% colnames(att)) return(0L)
  sum(grepl("^Metab__", att$Nodes))
}
node_count <- function(fact) {
  f <- file.path(cosmos_dir, paste0("Factor_", fact, "_MOON"), "combined_ATT_moon.csv")
  if (!file.exists(f)) return(0L)
  att <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(att)) 0L else nrow(att)
}

# ---- assemble per-factor sections ----
sections <- character()
toc      <- character()

for (fact in factors) {
  fact_dir <- file.path(cosmos_dir, paste0("Factor_", fact, "_MOON"))
  cosmos_subdir <- file.path(cosmos_dir, paste0("Factor_", fact, "_COSMOS"))
  if (!dir.exists(fact_dir)) next

  n_nodes  <- node_count(fact)
  n_metab  <- metab_count(fact)
  n_terms  <- sum(summary_all$Factor == fact)

  toc <- c(toc, sprintf(
    '<li><a href="#factor-%d">Factor %d</a> &mdash; %d network nodes (%d metabolites), %d enriched pathways</li>',
    fact, fact, n_nodes, n_metab, n_terms))

  plots <- c(
    embed_plot(file.path(fact_dir, "MOON_network.pdf"),
               sprintf("Factor %d — combined MOON network", fact)),
    embed_plot(file.path(cosmos_subdir, "venn_diagram.png"),
               sprintf("Factor %d — metabolite / meta-network overlap", fact)),
    embed_plot(file.path(fact_dir, "density_inputs.pdf"),
               sprintf("Factor %d — input distributions (RNA / TF / LR)", fact)),
    embed_plot(file.path(fact_dir, "density_moon_scores_rec_to_TFmetab.pdf"),
               sprintf("Factor %d — MOON scores: Rec → TF/Metab", fact)),
    embed_plot(file.path(fact_dir, "density_moon_scores_TF_lig.pdf"),
               sprintf("Factor %d — MOON scores: TF → Ligands", fact)),
    embed_plot(file.path(fact_dir, "mofa_top_TF_barplot.pdf"),
               sprintf("Factor %d — top TF activities", fact)),
    embed_plot(file.path(fact_dir, "mofa_top_ligrec_barplot.pdf"),
               sprintf("Factor %d — top ligand-receptor activities", fact)),
    embed_plot(file.path(fact_dir, "enrichr_GO_Biological_Process_2023.pdf"),
               sprintf("Factor %d — GO BP enrichment", fact)),
    embed_plot(file.path(fact_dir, "enrichr_KEGG_2021_Human.pdf"),
               sprintf("Factor %d — KEGG enrichment", fact))
  )
  plots <- plots[nzchar(plots)]

  # enriched-pathway tables for this factor (top N per database)
  fdat <- summary_all[summary_all$Factor == fact, ]
  tabs <- ""
  if (nrow(fdat) > 0) {
    for (db in enrich_dbs) {
      d <- fdat[fdat$Database == db, ]
      if (nrow(d) == 0) next
      d <- head(d, top_n_html)
      d <- d[, c("Term", "Overlap", "nGenes", "Adjusted.P.value", "Combined.Score")]
      tabs <- paste0(tabs,
        sprintf("<h4>%s (top %d of %d)</h4>", htmlesc(db), nrow(d),
                sum(fdat$Database == db)),
        df_to_table(d))
    }
  } else {
    tabs <- "<p><em>No significantly enriched pathways for this factor.</em></p>"
  }

  sections <- c(sections, sprintf(
    '<section id="factor-%d"><h2>Factor %d</h2>
     <p class="meta"><b>%d</b> nodes in combined network &middot; <b>%d</b> metabolite nodes &middot; <b>%d</b> enriched pathways (adj. p &lt; %.2f)</p>
     <div class="plots">%s</div>
     <h3>Enriched pathways</h3>%s</section>',
    fact, fact, n_nodes, n_metab, n_terms, padj_cutoff,
    paste0(plots, collapse = "\n"), tabs))
}

# ---- global plots ----
global_plots <- c(
  embed_plot(file.path(cosmos_dir, "COSMOS_metabolite_counts_all_factors.png"),
             "Metabolites with |weight| > 0.2 per factor (in vs not in COSMOS+)"),
  embed_plot(file.path(cosmos_dir, "mofa_top_TF_all_factors.pdf"),
             "Top TF activities across all factors"),
  embed_plot(file.path(cosmos_dir, "mofa_top_ligrec_all_factors.pdf"),
             "Top ligand-receptor activities across all factors")
)
global_plots <- global_plots[nzchar(global_plots)]

# ---- overview metabolite table ----
overview_tbl <- data.frame(
  Factor             = factors,
  Network_nodes      = sapply(factors, node_count),
  Metabolite_nodes   = sapply(factors, metab_count),
  Enriched_pathways  = sapply(factors, function(f) sum(summary_all$Factor == f))
)

css <- "
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;margin:0;color:#222;line-height:1.5}
header{background:#1a3a5c;color:#fff;padding:24px 32px}
header h1{margin:0 0 6px 0;font-size:24px}
header p{margin:0;opacity:.85;font-size:14px}
main{max-width:1000px;margin:0 auto;padding:24px 32px}
h2{border-bottom:2px solid #1a3a5c;padding-top:28px;margin-top:36px}
h3{color:#1a3a5c;margin-top:24px}
h4{color:#555;margin:18px 0 6px 0}
.meta{color:#666;font-size:14px}
table{border-collapse:collapse;width:100%;font-size:13px;margin:8px 0 18px 0}
th,td{border:1px solid #ddd;padding:5px 8px;text-align:left;vertical-align:top}
th{background:#eef3f8}
tr:nth-child(even){background:#fafafa}
figure{margin:14px 0}
figcaption{font-size:13px;color:#666;margin-top:4px}
.plots{display:flex;flex-direction:column;gap:8px}
ul.toc{columns:2;font-size:14px}
.callout{background:#fff6e5;border-left:4px solid #e0a106;padding:12px 16px;border-radius:4px;margin:16px 0;font-size:14px}
a{color:#1a3a5c}
"

html <- paste0(
  "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>MOON / COSMOS+ overview</title><style>", css, "</style></head><body>",
  "<header><h1>MOON / COSMOS+ overview</h1>",
  sprintf("<p>Per-factor MOON networks and enriched pathways &middot; generated %s</p></header>",
          format(Sys.time(), "%Y-%m-%d %H:%M")),
  "<main>",
  "<div class='callout'><b>Note on metabolites:</b> metabolites are retained in a factor's network only if their MOON score exceeds the run-1 <code>reduce_solution_network</code> cutoff and they remain reachable from an upstream receptor. With the previous cutoff of 1.5 (above the ~1.43 metabolite score ceiling) no metabolites survived; the cutoff is now configurable in <code>COSMOS_MOON_analysis.R</code>.</div>",
  "<h2 style='border:none;margin-top:0'>Summary</h2>",
  df_to_table(overview_tbl, numeric_signif = FALSE),
  "<h3>Contents</h3><ul class='toc'>", paste0(toc, collapse = ""), "</ul>",
  if (length(global_plots)) paste0("<h3>Cross-factor plots</h3><div class='plots'>",
                                   paste0(global_plots, collapse = "\n"), "</div>") else "",
  paste0(sections, collapse = "\n"),
  "</main></body></html>"
)

html_out <- file.path(cosmos_dir, "MOON_overview.html")
writeLines(html, html_out, useBytes = TRUE)


###########################################################################################
#### 3. Console report
###########################################################################################

per_factor <- summary_all %>%
  group_by(Factor) %>%
  summarise(n_terms = n(), n_databases = dplyr::n_distinct(Database), .groups = "drop")

cat(sprintf("\nSignificantly enriched pathways (Adjusted.P.value < %.2f):\n", padj_cutoff))
print(as.data.frame(per_factor), row.names = FALSE)
cat(sprintf("\nTotal: %d enriched terms across %d factors.\n",
            nrow(summary_all), length(unique(summary_all$Factor))))
cat("Saved:\n")
cat(sprintf("  %s\n  %s\n  %s\n", csv_out, xlsx_out, html_out))
