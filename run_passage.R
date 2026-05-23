# ============================================================================
# run_passage.R
#
# PASSAGE: top-level driver.  One call:
#   - Fits the engine on (Y, locs)
#   - Precomputes all per-factor eigenvalues for the requested hypotheses
#   - Runs H1/H2/H3 (+ optionally H4) on each pathway in `pathways`
#   - Computes the PVE estimators per pathway
#   - Returns a tidy data.frame plus the engine + per-pathway result objects
#
# Typical usage:
#
#   results <- run_passage(
#     Y = log_counts,           # N x G
#     locs = coords,            # N x 2
#     pathways = gobp_pathways, # named list of gene-name vectors
#     gene_names = rownames(...),
#     Z_CT = cell_type_props,   # N x n_celltypes  (optional, enables H2)
#     K = 6,
#     hypotheses = c("H1", "H2", "H3"),
#     adjust_method = "BH"
#   )
#   head(results$summary_table)
#
# Dependencies: all PASSAGE files (engine_pca, h1, h2, h3, h4, pve, omnibus)
# ============================================================================


# ----------------------------------------------------------------------------
# Main entry point.
# ----------------------------------------------------------------------------
run_passage <- function(
    Y, locs, pathways,
    gene_names = NULL,
    Z_CT = NULL,
    K = 6L,
    hypotheses = c("H1", "H2", "H3"),
    engine_args = list(),
    h3_K_BG = NULL,
    h3_bg_method = "size_matched",
    pve_compute = c("cca", "loo", "range", "meangene"),
    weight_schemes = c("equal", "var", "range"),
    adjust_method = "BH",
    min_pathway_size = 5L,
    max_pathway_size = 500L,
    n_cores = 1L,
    verbose = TRUE
) {

  if (is.null(gene_names)) {
    gene_names <- colnames(Y)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene_", seq_len(ncol(Y)))
      colnames(Y) <- gene_names
    }
  }
  if (verbose) {
    cat(sprintf("run_passage: N=%d locs, G=%d genes, K=%d factors, %d pathways requested\n",
                nrow(Y), ncol(Y), K, length(pathways)))
  }

  # ---------------------------------------------------------------------
  # 1. Filter pathways by size
  # ---------------------------------------------------------------------
  if (is.null(names(pathways))) names(pathways) <- paste0("pw_", seq_along(pathways))

  pw_sizes_raw <- vapply(pathways, function(p) length(intersect(p, gene_names)),
                         integer(1))
  keep <- pw_sizes_raw >= min_pathway_size & pw_sizes_raw <= max_pathway_size
  pathways <- pathways[keep]
  if (verbose) {
    cat(sprintf("  retained %d/%d pathways with %d <= size <= %d\n",
                length(pathways), length(keep), min_pathway_size,
                max_pathway_size))
  }
  if (length(pathways) == 0L) {
    stop("No pathways passed the size filter.")
  }

  # ---------------------------------------------------------------------
  # 2. Fit the engine
  # ---------------------------------------------------------------------
  if (verbose) cat("  fitting engine (PCA two-stage) ...\n")
  default_engine_args <- list(K = K, m = 20L, ordering = "maxmin",
                              smoothness = 1.5, D_orthogonalize = TRUE,
                              verbose = FALSE)
  engine_args_full <- modifyList(default_engine_args, engine_args)
  engine <- do.call(fit_engine_pca,
                    c(list(Y = Y, locs = locs), engine_args_full))

  # ---------------------------------------------------------------------
  # 3. Precompute eigenvalues per hypothesis
  # ---------------------------------------------------------------------
  precomp <- list()
  if ("H1" %in% hypotheses) {
    if (verbose) cat("  H1 precompute ...\n")
    precomp$H1 <- passage_h1_precompute(engine, verbose = FALSE)
  }
  if ("H2" %in% hypotheses) {
    if (is.null(Z_CT)) {
      warning("H2 requested but Z_CT not provided; skipping H2.")
      hypotheses <- setdiff(hypotheses, "H2")
    } else {
      if (verbose) cat("  H2 precompute ...\n")
      precomp$H2 <- passage_h2_precompute(engine, Z_CT, verbose = FALSE)
    }
  }
  if ("H3" %in% hypotheses) {
    if (verbose) {
      cat("  H3: a background engine is fit per pathway (size-matched).\n")
      cat("       Precompute reuses across pathways of the same size when possible.\n")
    }
    # H3 precompute happens per-pathway since the background depends on the
    # pathway's gene set. We cache by size bin to amortise across pathways.
    precomp$H3_cache <- new.env(hash = TRUE)
  }

  # ---------------------------------------------------------------------
  # 4. Run hypothesis tests per pathway
  # ---------------------------------------------------------------------
  if (verbose) cat("  running per-pathway tests ...\n")

  one_pathway <- function(pw_name) {
    pw_genes <- pathways[[pw_name]]
    pw_genes <- intersect(pw_genes, gene_names)
    # H3 precomp: cache by background gene-set size for amortisation
    precomp_h3_pw <- NULL
    if ("H3" %in% hypotheses) {
      cache_key <- as.character(length(pw_genes))
      if (!is.null(precomp$H3_cache[[cache_key]])) {
        precomp_h3_pw <- precomp$H3_cache[[cache_key]]
      } else {
        # Fit a background engine matched to this pathway size, precompute
        bg <- .fit_background_engine(
          engine, Y, pw_genes,
          K_BG = h3_K_BG, bg_method = h3_bg_method,
          gene_names = gene_names, verbose = FALSE)
        precomp_h3_pw <- passage_h3_precompute(
          engine, V_BG = bg$V_BG, Z_CT = Z_CT, verbose = FALSE)
        precomp$H3_cache[[cache_key]] <- precomp_h3_pw
      }
    }
    res <- tryCatch(
      passage_test_pathway(
        engine, Y, pw_genes,
        Z_CT = Z_CT,
        precomp_h1 = precomp$H1,
        precomp_h2 = precomp$H2,
        precomp_h3 = precomp_h3_pw,
        hypotheses = hypotheses,
        compute_pve = TRUE,
        pve_compute = pve_compute,
        weight_schemes = weight_schemes,
        gene_names = gene_names,
        verbose = FALSE),
      error = function(e) {
        message(sprintf("  pathway %s failed: %s", pw_name, conditionMessage(e)))
        NULL
      })
    res
  }

  per_pathway_results <- if (n_cores > 1L && requireNamespace("parallel", quietly = TRUE)) {
    parallel::mclapply(names(pathways), one_pathway, mc.cores = n_cores)
  } else {
    lapply(names(pathways), one_pathway)
  }
  names(per_pathway_results) <- names(pathways)

  # ---------------------------------------------------------------------
  # 5. Build tidy summary data.frame
  # ---------------------------------------------------------------------
  if (verbose) cat("  building summary table ...\n")
  rows <- list()
  for (nm in names(per_pathway_results)) {
    r <- per_pathway_results[[nm]]
    if (is.null(r)) next
    rows[[nm]] <- as.data.frame(r, pathway_name = nm)
  }
  summary_table <- do.call(rbind, rows)
  rownames(summary_table) <- NULL

  # FDR adjustments
  if (!is.null(adjust_method) && nrow(summary_table) > 0L) {
    for (col in c("p_H1", "p_H2", "p_H3")) {
      if (col %in% colnames(summary_table)) {
        summary_table[[paste0(col, "_adj")]] <- p.adjust(
          summary_table[[col]], method = adjust_method)
      }
    }
  }

  # Sort by adjusted H2 (if available) else H1
  sort_col <- if ("p_H2_adj" %in% names(summary_table)) "p_H2_adj"
              else if ("p_H1_adj" %in% names(summary_table)) "p_H1_adj"
              else "p_H1"
  if (sort_col %in% names(summary_table)) {
    summary_table <- summary_table[order(summary_table[[sort_col]]), , drop = FALSE]
  }

  out <- list(
    summary_table = summary_table,
    per_pathway = per_pathway_results,
    engine = engine,
    precomp = precomp[c("H1", "H2")],   # drop the per-size H3 cache
    pathways_input = pathways,
    hypotheses_run = hypotheses,
    adjust_method = adjust_method
  )
  class(out) <- c("passage_run", "list")
  out
}


# ----------------------------------------------------------------------------
# print method
# ----------------------------------------------------------------------------
print.passage_run <- function(x, ...) {
  cat("PASSAGE run\n")
  cat(sprintf("  hypotheses tested : %s\n",
              paste(x$hypotheses_run, collapse = ", ")))
  cat(sprintf("  pathways scored   : %d\n", nrow(x$summary_table)))
  cat("  Top 5 pathways (by sort column):\n")
  show_cols <- c("pathway", "pathway_size",
                 grep("^p_", names(x$summary_table), value = TRUE),
                 "pathway_specific_share", "R2_cca", "PSVS_range")
  show_cols <- intersect(show_cols, names(x$summary_table))
  print(head(x$summary_table[, show_cols, drop = FALSE], 5L))
  invisible(x)
}
