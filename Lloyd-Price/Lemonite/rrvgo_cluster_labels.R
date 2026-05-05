suppressPackageStartupMessages({
  options(stringsAsFactors = FALSE)
})

parse_args <- function(args) {
  parsed <- list()

  for (arg in args) {
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unsupported argument format: %s", arg), call. = FALSE)
    }

    key_value <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- key_value[[1]]
    value <- if (length(key_value) > 1) paste(key_value[-1], collapse = "=") else TRUE
    parsed[[key]] <- value
  }

  parsed
}

extract_go_id <- function(term) {
  if (is.na(term) || !nzchar(term)) {
    return(NA_character_)
  }

  matched <- regmatches(term, regexpr("GO:[0-9]+", term))
  if (!length(matched) || identical(matched, character(0)) || matched[[1]] == "") {
    return(NA_character_)
  }

  matched[[1]]
}

ensure_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      sprintf(
        "Missing required R packages: %s. Install them before running this helper.",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

score_terms <- function(p_adjust_values) {
  safe_values <- pmax(as.numeric(p_adjust_values), .Machine$double.xmin)
  sum(-log10(safe_values))
}

default_worker_cores <- function() {
  detected_cores <- suppressWarnings(parallel::detectCores(logical = FALSE))

  if (is.na(detected_cores) || detected_cores < 1L) {
    return(1L)
  }

  min(8L, as.integer(detected_cores))
}

select_fallback_label <- function(cluster_terms, reason) {
  if (nrow(cluster_terms) == 0) {
    return(list(
      representative_go_id = NA_character_,
      representative_term = "No BP enrichment",
      label_source = "fallback_no_bp_terms",
      fallback_reason = reason,
      module_support = 0L,
      audit = data.frame(
        go = NA_character_,
        term = "No BP enrichment",
        aggregate_score = NA_real_,
        module_support = 0L,
        parent = NA_character_,
        parentTerm = "No BP enrichment",
        reduction_source = reason,
        stringsAsFactors = FALSE
      )
    ))
  }

  fallback_row <- cluster_terms[order(cluster_terms$p.adjust, cluster_terms$Term), ][1, , drop = FALSE]

  list(
    representative_go_id = fallback_row$GO_ID[[1]],
    representative_term = fallback_row$Term[[1]],
    label_source = paste0("fallback_", reason),
    fallback_reason = reason,
    module_support = 1L,
    audit = data.frame(
      go = fallback_row$GO_ID[[1]],
      term = fallback_row$Term[[1]],
      aggregate_score = score_terms(fallback_row$p.adjust),
      module_support = 1L,
      parent = fallback_row$GO_ID[[1]],
      parentTerm = fallback_row$Term[[1]],
      reduction_source = reason,
      stringsAsFactors = FALSE
    )
  )
}

label_cluster <- function(cluster_name,
                          members,
                          cluster_terms,
                          semdata,
                          ontology,
                          rrvgo_method,
                          rrvgo_threshold,
                          score_method) {
  aggregated_terms <- cluster_terms %>%
    filter(!is.na(GO_ID)) %>%
    group_by(GO_ID) %>%
    summarise(
      representative_term = Term[which.min(p.adjust)],
      aggregate_score = score_terms(p.adjust),
      module_support = n_distinct(Module),
      min_p_adjust = min(p.adjust),
      .groups = "drop"
    ) %>%
    arrange(desc(aggregate_score), desc(module_support), min_p_adjust, representative_term)

  fallback <- NULL
  representative_go_id <- NA_character_
  representative_term <- NA_character_
  label_source <- NA_character_
  fallback_reason <- NA_character_
  representative_support <- NA_integer_
  audit_df <- data.frame()

  if (nrow(cluster_terms) == 0) {
    fallback <- select_fallback_label(cluster_terms, "no_bp_terms")
  } else if (nrow(aggregated_terms) == 0) {
    fallback <- select_fallback_label(cluster_terms, "no_valid_go_id")
  } else if (nrow(aggregated_terms) == 1) {
    representative_go_id <- aggregated_terms$GO_ID[[1]]
    representative_term <- aggregated_terms$representative_term[[1]]
    label_source <- "single_go_term"
    fallback_reason <- NA_character_
    representative_support <- aggregated_terms$module_support[[1]]
    audit_df <- aggregated_terms %>%
      transmute(
        go = GO_ID,
        term = representative_term,
        aggregate_score,
        module_support,
        parent = GO_ID,
        parentTerm = representative_term,
        reduction_source = label_source
      )
  } else {
    scores <- aggregated_terms$aggregate_score
    names(scores) <- aggregated_terms$GO_ID

    sim_matrix <- tryCatch(
      rrvgo::calculateSimMatrix(
        x = names(scores),
        orgdb = "org.Hs.eg.db",
        semdata = semdata,
        ont = ontology,
        method = rrvgo_method
      ),
      error = function(err) err
    )

    if (inherits(sim_matrix, "error") || !is.matrix(sim_matrix) || nrow(sim_matrix) < 2) {
      fallback <- select_fallback_label(cluster_terms, "rrvgo_sim_matrix_failed")
    } else {
      valid_rows <- rownames(sim_matrix)[rowSums(!is.na(sim_matrix)) > 0]
      if (length(valid_rows) < 2) {
        fallback <- select_fallback_label(cluster_terms, "rrvgo_insufficient_terms_after_filter")
      } else {
        sim_matrix <- sim_matrix[valid_rows, valid_rows, drop = FALSE]
        scores <- scores[valid_rows]

        reduced_terms <- tryCatch(
          rrvgo::reduceSimMatrix(
            simMatrix = sim_matrix,
            scores = scores,
            threshold = rrvgo_threshold,
            orgdb = "org.Hs.eg.db"
          ),
          error = function(err) err
        )

        if (inherits(reduced_terms, "error") || nrow(reduced_terms) == 0) {
          fallback <- select_fallback_label(cluster_terms, "rrvgo_reduce_failed")
        } else {
          reduced_terms <- reduced_terms %>%
            as.data.frame(stringsAsFactors = FALSE) %>%
            left_join(
              aggregated_terms %>%
                rename(
                  go = GO_ID,
                  aggregate_score_total = aggregate_score,
                  module_support_total = module_support,
                  best_term = representative_term
                ),
              by = "go"
            )

          representative_parent <- reduced_terms %>%
            group_by(parent, parentTerm) %>%
            summarise(
              total_score = sum(score, na.rm = TRUE),
              total_module_support = sum(module_support_total, na.rm = TRUE),
              n_terms = dplyr::n(),
              .groups = "drop"
            ) %>%
            arrange(desc(total_score), desc(total_module_support), desc(n_terms), parentTerm) %>%
            slice_head(n = 1)

          representative_go_id <- representative_parent$parent[[1]]
          representative_term <- representative_parent$parentTerm[[1]]
          label_source <- "rrvgo_parentTerm"
          fallback_reason <- NA_character_
          representative_support <- aggregated_terms %>%
            filter(GO_ID == representative_go_id) %>%
            pull(module_support) %>%
            {
              if (length(.) == 0) NA_integer_ else as.integer(max(., na.rm = TRUE))
            }

          audit_df <- reduced_terms %>%
            transmute(
              go,
              term = dplyr::coalesce(best_term, term),
              aggregate_score = aggregate_score_total,
              module_support = module_support_total,
              parent,
              parentTerm,
              reduction_source = label_source
            )
        }
      }
    }
  }

  if (!is.null(fallback)) {
    representative_go_id <- fallback$representative_go_id
    representative_term <- fallback$representative_term
    label_source <- fallback$label_source
    fallback_reason <- fallback$fallback_reason
    representative_support <- fallback$module_support
    audit_df <- fallback$audit
  }

  list(
    cluster_label = data.frame(
      MegaGO_cluster = cluster_name,
      representative_go_id = representative_go_id,
      MegaGO_label = representative_term,
      n_modules = nrow(members),
      n_terms = nrow(cluster_terms),
      module_support = representative_support,
      score_method = score_method,
      rrvgo_threshold = rrvgo_threshold,
      rrvgo_method = rrvgo_method,
      label_source = label_source,
      fallback_reason = fallback_reason,
      stringsAsFactors = FALSE
    ),
    reduced_terms = audit_df %>%
      mutate(
        MegaGO_cluster = cluster_name,
        representative_go_id = representative_go_id,
        representative_label = representative_term,
        stringsAsFactors = FALSE
      )
  )
}

run_rrvgo_cluster_labels <- function(enrichment_path,
                                     cluster_path,
                                     output_dir,
                                     cluster_column = "top_30",
                                     ontology = "BP",
                                     top_n = 30L,
                                     rrvgo_threshold = 0.7,
                                     rrvgo_method = "Rel",
                                     worker_cores = default_worker_cores()) {
  ensure_packages(c("dplyr", "rrvgo", "org.Hs.eg.db", "GOSemSim"))

  suppressPackageStartupMessages({
    library(dplyr)
    library(rrvgo)
  })

  enrichment_path <- normalizePath(enrichment_path, mustWork = TRUE)
  cluster_path <- normalizePath(cluster_path, mustWork = TRUE)
  top_n <- as.integer(top_n)
  rrvgo_threshold <- as.numeric(rrvgo_threshold)
  worker_cores <- as.integer(worker_cores)
  score_method <- "sum(-log10(p.adjust))"

  if (is.na(worker_cores) || worker_cores < 1L) {
    stop("Argument 'cores' must be a positive integer.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  enrichment <- read.csv(enrichment_path, check.names = FALSE, stringsAsFactors = FALSE)
  required_columns <- c("Module", "Database", "Term", "p.adjust")
  missing_columns <- required_columns[!required_columns %in% colnames(enrichment)]
  if (length(missing_columns) > 0) {
    stop(
      sprintf(
        "Enrichment file is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  cluster_assignments <- read.csv(cluster_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  if (!cluster_column %in% colnames(cluster_assignments)) {
    stop(
      sprintf("Cluster column '%s' not present in %s", cluster_column, cluster_path),
      call. = FALSE
    )
  }

  module_clusters <- data.frame(
    Module = rownames(cluster_assignments),
    MegaGO_cluster = as.character(cluster_assignments[[cluster_column]]),
    stringsAsFactors = FALSE
  )

  bp_terms <- enrichment %>%
    mutate(
      Module = as.character(Module),
      Database = as.character(Database),
      Term = as.character(Term),
      p.adjust = as.numeric(p.adjust),
      GO_ID = vapply(Term, extract_go_id, character(1))
    ) %>%
    filter(Database == ontology, !is.na(p.adjust)) %>%
    group_by(Module) %>%
    arrange(p.adjust, Term, .by_group = TRUE) %>%
    {
      if (is.na(top_n) || top_n <= 0) {
        .
      } else {
        slice_head(., n = top_n)
      }
    } %>%
    ungroup()

  assigned_modules <- module_clusters %>%
    filter(!is.na(MegaGO_cluster), MegaGO_cluster != "N/A")

  clusters_to_label <- assigned_modules %>%
    filter(MegaGO_cluster != "Unassigned") %>%
    pull(MegaGO_cluster) %>%
    unique() %>%
    sort()

  cluster_mask <- assigned_modules$MegaGO_cluster %in% clusters_to_label
  cluster_members_by_name <- split(
    assigned_modules[cluster_mask, , drop = FALSE],
    assigned_modules$MegaGO_cluster[cluster_mask]
  )

  cluster_terms_input <- bp_terms %>%
    inner_join(assigned_modules[cluster_mask, , drop = FALSE], by = "Module")

  empty_cluster_terms <- cluster_terms_input[0, , drop = FALSE]
  cluster_terms_by_name <- if (nrow(cluster_terms_input) > 0) {
    split(cluster_terms_input, cluster_terms_input$MegaGO_cluster)
  } else {
    list()
  }

  cluster_results <- list()

  if (length(clusters_to_label) > 0) {
    message(sprintf("Preparing GO semantic data once for ontology %s...", ontology))
    semdata <- GOSemSim::godata("org.Hs.eg.db", ont = ontology, keytype = "ENTREZID")
    message(sprintf("Labeling %d MegaGO clusters using %d core(s)...", length(clusters_to_label), worker_cores))

    process_cluster <- function(cluster_name) {
      members <- cluster_members_by_name[[cluster_name]]
      cluster_terms <- cluster_terms_by_name[[cluster_name]]

      if (is.null(cluster_terms)) {
        cluster_terms <- empty_cluster_terms
      }

      label_cluster(
        cluster_name = cluster_name,
        members = members,
        cluster_terms = cluster_terms,
        semdata = semdata,
        ontology = ontology,
        rrvgo_method = rrvgo_method,
        rrvgo_threshold = rrvgo_threshold,
        score_method = score_method
      )
    }

    cluster_results <- if (worker_cores > 1L && .Platform$OS.type == "unix") {
      parallel::mclapply(
        clusters_to_label,
        process_cluster,
        mc.cores = worker_cores,
        mc.preschedule = FALSE
      )
    } else {
      if (worker_cores > 1L && .Platform$OS.type != "unix") {
        message("Parallel execution is only enabled on Unix-like systems; falling back to serial execution.")
      }

      lapply(clusters_to_label, process_cluster)
    }
  }

  cluster_label_rows <- lapply(cluster_results, function(result) result$cluster_label)
  reduced_term_rows <- lapply(cluster_results, function(result) result$reduced_terms)

  cluster_labels <- if (length(cluster_label_rows) > 0) {
    bind_rows(cluster_label_rows) %>%
      arrange(MegaGO_cluster)
  } else {
    data.frame(
      MegaGO_cluster = character(),
      representative_go_id = character(),
      MegaGO_label = character(),
      n_modules = integer(),
      n_terms = integer(),
      module_support = integer(),
      score_method = character(),
      rrvgo_threshold = numeric(),
      rrvgo_method = character(),
      label_source = character(),
      fallback_reason = character(),
      stringsAsFactors = FALSE
    )
  }

  reduced_terms_audit <- if (length(reduced_term_rows) > 0) {
    bind_rows(reduced_term_rows) %>%
      arrange(MegaGO_cluster, desc(aggregate_score), term)
  } else {
    data.frame(
      go = character(),
      term = character(),
      aggregate_score = numeric(),
      module_support = integer(),
      parent = character(),
      parentTerm = character(),
      reduction_source = character(),
      MegaGO_cluster = character(),
      representative_go_id = character(),
      representative_label = character(),
      stringsAsFactors = FALSE
    )
  }

  module_labels <- module_clusters %>%
    left_join(
      cluster_labels %>%
        select(MegaGO_cluster, MegaGO_label, representative_go_id, label_source),
      by = "MegaGO_cluster"
    ) %>%
    mutate(
      label_source = ifelse(MegaGO_cluster == "Unassigned", "unassigned", label_source),
      MegaGO_label = ifelse(MegaGO_cluster == "Unassigned", NA_character_, MegaGO_label),
      representative_go_id = ifelse(MegaGO_cluster == "Unassigned", NA_character_, representative_go_id)
    ) %>%
    arrange(suppressWarnings(as.numeric(Module)), Module)

  cluster_labels_path <- file.path(output_dir, paste0("rrvgo_cluster_labels_", cluster_column, ".csv"))
  module_labels_path <- file.path(output_dir, paste0("rrvgo_module_labels_", cluster_column, ".csv"))
  reduced_terms_path <- file.path(output_dir, paste0("rrvgo_reduced_terms_", cluster_column, ".csv"))

  write.csv(cluster_labels, cluster_labels_path, row.names = FALSE, na = "")
  write.csv(module_labels, module_labels_path, row.names = FALSE, na = "")
  write.csv(reduced_terms_audit, reduced_terms_path, row.names = FALSE, na = "")

  message(sprintf("Wrote cluster labels to %s", cluster_labels_path))
  message(sprintf("Wrote module labels to %s", module_labels_path))
  message(sprintf("Wrote reduced-term audit to %s", reduced_terms_path))

  invisible(list(
    cluster_labels = cluster_labels,
    module_labels = module_labels,
    reduced_terms_audit = reduced_terms_audit,
    cluster_labels_path = cluster_labels_path,
    module_labels_path = module_labels_path,
    reduced_terms_path = reduced_terms_path
  ))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  required_args <- c("enrichment", "clusters", "out-dir")
  missing_args <- required_args[!required_args %in% names(args)]
  if (length(missing_args) > 0) {
    stop(
      sprintf("Missing required arguments: %s", paste(missing_args, collapse = ", ")),
      call. = FALSE
    )
  }

  run_rrvgo_cluster_labels(
    enrichment_path = args[["enrichment"]],
    cluster_path = args[["clusters"]],
    output_dir = args[["out-dir"]],
    cluster_column = if (!is.null(args[["cluster-column"]])) args[["cluster-column"]] else "top_30",
    ontology = if (!is.null(args[["ontology"]])) args[["ontology"]] else "BP",
    top_n = if (!is.null(args[["top-n"]])) as.integer(args[["top-n"]]) else 30L,
    rrvgo_threshold = if (!is.null(args[["threshold"]])) as.numeric(args[["threshold"]]) else 0.7,
    rrvgo_method = if (!is.null(args[["method"]])) args[["method"]] else "Rel",
    worker_cores = if (!is.null(args[["cores"]])) as.integer(args[["cores"]]) else default_worker_cores()
  )
}

if (sys.nframe() == 0L) {
  main()
}