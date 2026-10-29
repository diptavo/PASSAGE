#!/usr/bin/env Rscript
# From the repository: Rscript inst/examples/run_kidney.R
# In R: source(system.file("examples", "run_kidney.R", package = "PASSAGE"))

run_kidney <- function(mode = "all", B = 199L, output = paste0("results/kidney_", mode),
                       method = "pca", K = 3L, n_features = 64L,
                       descriptive_gp = FALSE, seed = 1L, input = NULL) {
  if (!requireNamespace("PASSAGE", quietly = TRUE)) {
    stop("Install PASSAGE first: remotes::install_github('diptavo/PASSAGE')")
  }
  mode <- match.arg(mode, c("all", "gp", "competitive", "fixed"))
  method <- match.arg(method, c("pca", "factor_analysis"))
  if (length(B) != 1L || !is.finite(B) || B < 1 || B != as.integer(B)) {
    stop("B must be a positive integer")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max - 3L || seed != as.integer(seed)) {
    stop("seed must be a nonnegative integer below the integer limit")
  }
  if (method == "factor_analysis" && is.null(K)) stop("Factor analysis needs a fixed K")
  if (mode %in% c("all", "competitive") &&
      !requireNamespace("irlba", quietly = TRUE) &&
      !requireNamespace("RSpectra", quietly = TRUE)) {
    stop("Install irlba for the full-background example: install.packages('irlba')")
  }
  if (descriptive_gp && (!requireNamespace("GpGp", quietly = TRUE) ||
                         !requireNamespace("fields", quietly = TRUE))) {
    stop("Descriptive fits need install.packages(c('GpGp', 'fields'))")
  }
  if (is.null(input)) {
    input <- if (file.exists("inst/extdata/kidney_KC2.rds")) {
      "inst/extdata/kidney_KC2.rds"
    } else system.file("extdata", "kidney_KC2.rds", package = "PASSAGE")
  }
  if (!nzchar(input) || !file.exists(input)) stop("KC2 example data not found; update PASSAGE")
  if (file.exists(output) && !dir.exists(output)) stop("Output is a file: ", output)
  if (dir.exists(output) && length(list.files(output, all.files = TRUE, no.. = TRUE))) {
    stop("Choose a new or empty output directory: ", output)
  }
  b <- PASSAGE::passage_validate_bundle(input)
  pathways <- c("HALLMARK_HYPOXIA", "HALLMARK_ANGIOGENESIS",
                "HALLMARK_INFLAMMATORY_RESPONSE")
  if (!all(pathways %in% names(b$pathway_sets))) stop("Missing example pathway definitions")
  if (!dir.create(output, recursive = TRUE, showWarnings = FALSE) && !dir.exists(output)) {
    stop("Cannot create output directory: ", output)
  }
  started <- proc.time()[[3L]]
  coverage <- data.frame(pathway = pathways, n_input = lengths(b$pathway_sets[pathways]),
    n_observed = vapply(b$pathway_sets[pathways], function(g) {
      as.integer(sum(unique(g) %in% colnames(b$Y)))
    }, integer(1)))
  write.csv(coverage, file.path(output, "pathway_coverage.csv"), row.names = FALSE)
  offset <- apply(b$coords, 2L, min)
  xy <- sweep(b$coords, 2L, offset, "-")
  coordinate_scale <- max(xy)
  xy <- xy / coordinate_scale
  ranges <- c(0.05, 0.1, 0.2, 0.4)
  configuration <- list(
    PASSAGE_version = as.character(utils::packageVersion("PASSAGE")),
    input_md5 = unname(tools::md5sum(input)), sample = b$sample,
    n_spots = nrow(b$Y), n_background_genes = ncol(b$Y),
    covariates = colnames(b$X), pathways = pathways, mode = mode, B = B, seed = seed,
    pathway_method = method, pathway_K = K, pathway_variance_threshold = 0.95,
    coordinate_offset = offset, coordinate_scale = coordinate_scale,
    gp_ranges = ranges, gp_n_features = n_features, fixed_kernel_range = 0.25,
    competitive_K = 6L, competitive_neighbors = 20L,
    competitive_ranges = ranges, descriptive_gp = descriptive_gp,
    primary_self_contained_test = "gp", fixed_test_role = "sensitivity analysis",
    interpretation = "Research example; study-specific calibration and replication required"
  )
  saveRDS(configuration, file.path(output, "configuration.rds"))
  capture.output(utils::sessionInfo(), file = file.path(output, "sessionInfo.txt"))
  results <- list()
  add_result <- function(pathway, test, n_genes, K, p) {
    data.frame(pathway = pathway, test = test, n_genes = n_genes, K = K, p = p)
  }
  message("KC2: ", nrow(b$Y), " spots, ", ncol(b$Y), " background genes; B = ", B)
  if (mode != "competitive") {
    cache <- if (mode %in% c("all", "gp")) PASSAGE::passage_prepare_gp_kernel(
      xy, X = b$X, ranges = ranges, n_features = n_features) else NULL
    for (i in seq_along(pathways)) {
      pathway <- pathways[i]
      message("Pathway estimation and self-contained testing: ", pathway)
      fit <- PASSAGE::passage_fit_pathway(b$Y, b$pathway_sets[[pathway]], X = b$X,
        method = method, K = K, variance_threshold = 0.95, seed = seed + i)
      out <- file.path(output, pathway)
      dir.create(out)
      saveRDS(fit, file.path(out, "pathway_fit.rds"))
      write.csv(data.frame(spot = rownames(fit$scores), xy, fit$scores),
                file.path(out, "component_scores.csv"), row.names = FALSE)
      write.csv(data.frame(gene = rownames(fit$loadings), fit$loadings),
                file.path(out, "component_loadings.csv"), row.names = FALSE)
      if (mode %in% c("all", "gp")) {
        gp <- PASSAGE::passage_test_pathway_gp(fit, kernel_cache = cache, B = B, seed = seed + i)
        saveRDS(gp, file.path(out, "gp_test.rds"))
        write.csv(gp$component_fits, file.path(out, "gp_component_fits.csv"), row.names = FALSE)
        results[[length(results) + 1L]] <- add_result(pathway, "gp", length(fit$gene_names), fit$K, gp$p)
      }
      if (mode %in% c("all", "fixed")) {
        fixed <- PASSAGE::passage_test_pathway_fit(fit, coords = xy,
          kernel_range = 0.25, B = B, seed = seed + i)
        saveRDS(fixed, file.path(out, "fixed_test.rds"))
        results[[length(results) + 1L]] <- add_result(pathway, "fixed_omnibus",
          length(fit$gene_names), fit$K, fixed$p_omnibus_bonferroni)
      }
      if (descriptive_gp) {
        descriptive <- PASSAGE::passage_fit_pathway_gp(fit, xy, seed = seed + i)
        saveRDS(descriptive, file.path(out, "descriptive_gp.rds"))
        write.csv(descriptive$summary, file.path(out, "descriptive_gp.csv"), row.names = FALSE)
      }
    }
  }
  if (mode %in% c("all", "competitive")) {
    message("Competitive testing against the full gene background")
    competitive <- PASSAGE::passage_analyze_bundle(b, pathways = b$pathway_sets[pathways],
      B = B, seed = seed, K = 6L, m = 20L, range_grid = ranges)
    saveRDS(competitive, file.path(output, "competitive_result.rds"))
    write.csv(competitive$summary, file.path(output, "competitive_summary.csv"), row.names = FALSE)
    for (i in seq_len(nrow(competitive$summary))) {
      row <- competitive$summary[i, ]
      results[[length(results) + 1L]] <- add_result(row$pathway, "competitive",
        row$n_genes, 6L, row$p_competitive)
    }
  }
  summary <- do.call(rbind, results)
  rownames(summary) <- NULL
  summary$q_BH_within_test <- ave(summary$p, summary$test, FUN = function(p) stats::p.adjust(p, "BH"))
  summary$q_BY_within_test <- ave(summary$p, summary$test, FUN = function(p) stats::p.adjust(p, "BY"))
  summary$q_BH_all_reported <- stats::p.adjust(summary$p, "BH")
  summary$q_BY_all_reported <- stats::p.adjust(summary$p, "BY")
  summary$mc_floor <- ifelse(summary$test == "fixed_omnibus",
                            min(1, 2 / (B + 1)), 1 / (B + 1))
  write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
  configuration$elapsed_seconds <- proc.time()[[3L]] - started
  saveRDS(configuration, file.path(output, "configuration.rds"))
  writeLines("PASS", file.path(output, "analysis.complete"))
  print(summary, row.names = FALSE)
  message("Finished in ", round(configuration$elapsed_seconds, 1), " seconds. Results: ", output)
  invisible(summary)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 5L || any(args %in% c("--help", "-h"))) {
    cat("Usage: Rscript inst/examples/run_kidney.R [MODE] [B] [OUTPUT] [METHOD] [K]\n",
        "Defaults: all 199 results/kidney_all pca 3\n",
        "MODE: all, gp, competitive, fixed; METHOD: pca, factor_analysis; K: integer or auto (PCA only)\n")
    quit(status = if (length(args) > 5L) 1L else 0L)
  }
  mode <- if (length(args) >= 1L) args[1L] else "all"
  B <- if (length(args) >= 2L) as.numeric(args[2L]) else 199L
  output <- if (length(args) >= 3L) args[3L] else paste0("results/kidney_", mode)
  method <- if (length(args) >= 4L) args[4L] else "pca"
  K <- if (length(args) >= 5L) {
    if (args[5L] == "auto") NULL else as.numeric(args[5L])
  } else 3L
  run_kidney(mode = mode, B = B, output = output, method = method, K = K)
}
