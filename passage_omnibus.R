# ============================================================================
# passage_omnibus.R
#
# PASSAGE: Layer 4 / Cross-hypothesis decomposition and unified reporting.
#
# Combines the per-hypothesis test outputs (H1, H2, H3) into the variance-
# share decomposition that is the actual scientific deliverable:
#
#   total spatial signal     = H1 statistic
#   cell-type-explained      = (H1 - H2) / H1
#   background-explained     = (H2 - H3) / H2
#   pathway-specific spatial = H3 / H1
#
# Per pathway, this three-way breakdown is more informative than three
# separate p-values; it tells the biologist *why* a pathway is spatial.
#
# Also provides:
#   - A unified `passage_test_pathway()` driver that runs H1, H2, H3 and
#     returns the full result with decomposition.
#   - Tidy data.frame conversion for downstream aggregation across pathways.
#
# Dependencies: passage_score_h1.R, _h2.R, _h3.R
# ============================================================================


# ----------------------------------------------------------------------------
# Compute decomposition shares from H1/H2/H3 result objects.
#
# Uses the joint Q-statistics (sum across factors, equal weighting) as the
# canonical effect-size proxy.  Each Q has dimension proportional to spatial
# variance after the relevant covariate adjustment.
# ----------------------------------------------------------------------------
passage_decomposition <- function(h1, h2 = NULL, h3 = NULL) {
  # Q values from the equal-weighting joint test
  Q1 <- if (!is.null(h1)) h1$joint$equal$Q_joint else NA_real_
  Q2 <- if (!is.null(h2)) h2$joint$equal$Q_joint else NA_real_
  Q3 <- if (!is.null(h3)) h3$joint$equal$Q_joint else NA_real_

  Q1_safe <- max(Q1, 1e-12)
  Q2_safe <- max(Q2, 1e-12)

  cell_type_share <- if (!is.na(Q1) && !is.na(Q2)) {
    pmax(pmin((Q1 - Q2) / Q1_safe, 1), 0)
  } else NA_real_
  background_share <- if (!is.na(Q2) && !is.na(Q3)) {
    pmax(pmin((Q2 - Q3) / Q2_safe, 1), 0)
  } else NA_real_
  specific_share <- if (!is.na(Q1) && !is.na(Q3)) {
    pmax(pmin(Q3 / Q1_safe, 1), 0)
  } else NA_real_

  list(
    Q_H1 = Q1, Q_H2 = Q2, Q_H3 = Q3,
    cell_type_share      = cell_type_share,
    background_share     = background_share,
    pathway_specific_share = specific_share,
    p_H1 = if (!is.null(h1)) h1$p_omnibus else NA_real_,
    p_H2 = if (!is.null(h2)) h2$p_omnibus else NA_real_,
    p_H3 = if (!is.null(h3)) h3$p_omnibus else NA_real_
  )
}


# ----------------------------------------------------------------------------
# All-in-one pathway tester: runs H1, optionally H2, optionally H3 + PVE,
# returns a structured `passage_pathway_result`.
# ----------------------------------------------------------------------------
passage_test_pathway <- function(engine, Y, pathway,
                                 Z_CT = NULL, V_BG = NULL,
                                 precomp_h1 = NULL, precomp_h2 = NULL,
                                 precomp_h3 = NULL,
                                 hypotheses = c("H1", "H2", "H3"),
                                 K_BG = NULL,
                                 bg_method = "size_matched",
                                 compute_pve = TRUE,
                                 pve_compute = c("cca", "loo", "range", "meangene"),
                                 weight_schemes = c("equal", "var", "range"),
                                 gene_names = NULL,
                                 verbose = FALSE) {

  results <- list()

  # H1
  if ("H1" %in% hypotheses) {
    if (verbose) cat("  H1 (any spatial signal)\n")
    results$H1 <- passage_h1(engine, Y, pathway, precomp = precomp_h1,
                             weight_schemes = weight_schemes,
                             gene_names = gene_names, verbose = FALSE)
  }

  # H2
  if ("H2" %in% hypotheses) {
    if (is.null(Z_CT) && is.null(precomp_h2)) {
      warning("H2 requested but Z_CT not provided; skipping.")
    } else {
      if (verbose) cat("  H2 (beyond cell type)\n")
      results$H2 <- passage_h2(engine, Y, pathway, Z_CT = Z_CT,
                               precomp = precomp_h2,
                               weight_schemes = weight_schemes,
                               gene_names = gene_names, verbose = FALSE)
    }
  }

  # H3
  if ("H3" %in% hypotheses) {
    if (verbose) cat("  H3 (beyond background)\n")
    results$H3 <- passage_h3(engine, Y, pathway, precomp = precomp_h3,
                             Z_CT = Z_CT, V_BG = V_BG,
                             K_BG = K_BG, bg_method = bg_method,
                             weight_schemes = weight_schemes,
                             gene_names = gene_names, verbose = FALSE)
  }

  # Decomposition
  decomp <- passage_decomposition(
    h1 = results$H1, h2 = results$H2, h3 = results$H3)

  # PVE
  pve <- if (compute_pve) {
    if (verbose) cat("  PVE estimators\n")
    passage_pve(engine, Y, pathway, gene_names = gene_names,
                compute = pve_compute)
  } else NULL

  out <- list(
    hypotheses = results,
    decomposition = decomp,
    pve = pve,
    pathway_size = if (!is.null(results$H1)) results$H1$pathway_size else NA_integer_
  )
  class(out) <- c("passage_pathway_result", "list")
  out
}


# ----------------------------------------------------------------------------
# Tidy conversion: turn a passage_pathway_result into a one-row data.frame
# suitable for rbind across many pathways.
# ----------------------------------------------------------------------------
as.data.frame.passage_pathway_result <- function(x, ..., pathway_name = NA) {
  data.frame(
    pathway       = pathway_name,
    pathway_size  = x$pathway_size,
    p_H1          = x$decomposition$p_H1,
    p_H2          = x$decomposition$p_H2,
    p_H3          = x$decomposition$p_H3,
    Q_H1          = x$decomposition$Q_H1,
    Q_H2          = x$decomposition$Q_H2,
    Q_H3          = x$decomposition$Q_H3,
    cell_type_share        = x$decomposition$cell_type_share,
    background_share       = x$decomposition$background_share,
    pathway_specific_share = x$decomposition$pathway_specific_share,
    R2_cca       = if (!is.null(x$pve)) x$pve$summary["R2_cca"]      else NA_real_,
    R2_loo       = if (!is.null(x$pve)) x$pve$summary["R2_loo"]      else NA_real_,
    PSVS_range   = if (!is.null(x$pve)) x$pve$summary["PSVS_range"]  else NA_real_,
    mean_propSV  = if (!is.null(x$pve)) x$pve$summary["mean_propSV"] else NA_real_,
    spasset_subset_size_H1 = if (!is.null(x$hypotheses$H1$spasset))
      length(x$hypotheses$H1$spasset$best_subset) else NA_integer_,
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# print method
# ----------------------------------------------------------------------------
print.passage_pathway_result <- function(x, ...) {
  cat("PASSAGE pathway result\n")
  cat(sprintf("  Pathway size : %d\n", x$pathway_size))
  cat("  Hypothesis p-values:\n")
  for (h in names(x$hypotheses)) {
    cat(sprintf("    %s  : %s\n", h,
                format.pval(x$hypotheses[[h]]$p_omnibus, digits = 3)))
  }
  cat("  Decomposition:\n")
  cat(sprintf("    cell-type explained share        : %.3f\n",
              x$decomposition$cell_type_share %|% NA))
  cat(sprintf("    background-spatial explained share: %.3f\n",
              x$decomposition$background_share %|% NA))
  cat(sprintf("    pathway-specific spatial share   : %.3f\n",
              x$decomposition$pathway_specific_share %|% NA))
  if (!is.null(x$pve)) {
    cat("  Effect sizes:\n")
    print(round(x$pve$summary, 4))
  }
  invisible(x)
}

# small util for printing NAs nicely
`%|%` <- function(a, b) if (is.na(a)) b else a
