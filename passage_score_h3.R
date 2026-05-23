# ============================================================================
# passage_score_h3.R
#
# PASSAGE: Layer 2 / H3 score test (beyond background spatial pattern)
#
# Tests
#   H_0^(3):  sigma_k^2 = 0  for all k,
#             AFTER adjusting for cell-type proportions AND for the dominant
#             spatial patterns of the surrounding (non-pathway) transcriptome.
#
# Mean-model design:
#   X^(3) = [1_N, Z_CT, V_hat_BG],
# where V_hat_BG is the N x K_BG matrix of factor scores from a "background"
# engine fit on a size-matched random sample of non-pathway genes.
#
# Interpretation of H_0 rejection:
#   "This pathway carries a spatial pattern that is not expressible in the
#    background spatial basis of the tissue."  This is the most stringent
#    test in the H1-H3 hierarchy and identifies pathways with characteristic
#    spatial structure.
#
# Computational architecture:
#   - For each pathway, fit a "background engine" on a size-matched random
#     sample of non-pathway genes (or use a precomputed shared background
#     for speed).
#   - Construct X^(3) using V_hat_BG from that background fit.
#   - Run the H2-style score machinery with M^(3) = I - X^(3) (X^(3)'X^(3))^{-1} X^(3)'.
#
# Calibration concern:
#   The background basis may be over-rich if K_BG is too large; H3 would
#   never reject.  Size-matching mitigates this.  Sensitivity analysis over
#   K_BG and over which non-pathway genes are used should be standard.
#
# Dependencies: Matrix, CompQuadForm, engine_pca.R, passage_score_h1.R
# ============================================================================


# ----------------------------------------------------------------------------
# Internal: fit a background engine on a size-matched random sample of
# non-pathway genes.  Returns just the V_scores (N x K_BG) at the original
# location indexing -- that's all H3 needs from the background.
# ----------------------------------------------------------------------------
.fit_background_engine <- function(engine, Y, pathway,
                                   K_BG = NULL,
                                   bg_method = c("size_matched", "all_nonpath",
                                                 "random_fixed"),
                                   bg_size = NULL,
                                   m = NULL,
                                   smoothness = NULL,
                                   gene_names = NULL,
                                   seed = 20260522,
                                   verbose = FALSE) {
  bg_method <- match.arg(bg_method)
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)
  if (is.null(K_BG)) K_BG <- engine$K
  if (is.null(m)) m <- engine$vecchia$m
  if (is.null(smoothness)) smoothness <- engine$fit_diagnostics$smoothness

  # Choose background gene set
  non_P <- setdiff(seq_len(engine$G), P)
  n_target <- switch(bg_method,
    "size_matched"  = length(P),
    "all_nonpath"   = length(non_P),
    "random_fixed"  = if (is.null(bg_size)) 500L else as.integer(bg_size)
  )
  n_target <- min(n_target, length(non_P))

  set.seed(seed)
  bg_idx <- if (n_target == length(non_P)) {
    non_P
  } else {
    sort(sample(non_P, n_target))
  }

  if (verbose) cat(sprintf(
    "  fitting background engine on %d genes (K_BG=%d)\n", n_target, K_BG))

  Y_bg <- Y[, bg_idx, drop = FALSE]
  bg_engine <- fit_engine_pca(
    Y = Y_bg, locs = engine$vecchia$locs_ord[engine$vecchia$inv_ord, , drop = FALSE],
    K = K_BG, m = m, ordering = "maxmin", smoothness = smoothness,
    D_orthogonalize = FALSE, verbose = FALSE)
  list(V_BG = bg_engine$V_scores, K_BG = K_BG, bg_idx = bg_idx)
}


# ----------------------------------------------------------------------------
# Precompute per-factor eigenvalues for H3.
#
# Builds X^(3) = [1, Z_CT, V_BG] (or [1, V_BG] if no Z_CT) and computes
# the per-factor eigenvalues of M^(3) Qtilde_k M^(3).
#
# Requires a `background_fit` from .fit_background_engine() or pass a
# pre-fitted V_BG matrix (N x K_BG) directly.
# ----------------------------------------------------------------------------
passage_h3_precompute <- function(engine, V_BG, Z_CT = NULL,
                                  verbose = TRUE, max_dense_N = 15000) {
  K <- engine$K; N <- engine$N
  stopifnot(nrow(V_BG) == N)
  if (!is.null(Z_CT)) stopifnot(nrow(Z_CT) == N)
  if (N > max_dense_N) {
    warning(sprintf("N = %d > max_dense_N = %d; H3 precompute may run out of memory.",
                    N, max_dense_N))
  }

  ord <- engine$vecchia$ord
  V_BG_ord <- V_BG[ord, , drop = FALSE]

  X_ord <- if (is.null(Z_CT)) {
    cbind(1.0, V_BG_ord)
  } else {
    cbind(1.0, Z_CT[ord, , drop = FALSE], V_BG_ord)
  }
  qx <- qr(X_ord)
  H <- qr.Q(qx)
  H <- H %*% t(H)
  M_ord <- diag(N) - H

  if (verbose) cat(sprintf(
    "passage_h3_precompute: K=%d, K_BG=%d, q_covar=%d, N=%d\n",
    K, ncol(V_BG_ord), ncol(X_ord) - 1L, N))

  eigvals_per_factor <- vector("list", K)
  for (k in seq_len(K)) {
    if (verbose) cat(sprintf("  factor %d/%d ... ", k, K))
    t0 <- Sys.time()
    Qk <- engine$Qtilde_list[[k]]
    MQ <- M_ord %*% Qk
    MQM <- MQ %*% M_ord
    rm(MQ); gc()
    MQM <- (MQM + t(MQM)) / 2
    eig <- eigen(MQM, symmetric = TRUE, only.values = TRUE)
    eigvals_per_factor[[k]] <- pmax(eig$values, 0)
    rm(MQM)
    if (verbose) {
      cat(sprintf("done (%.1fs)\n",
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
  }

  list(
    eigvals_per_factor = eigvals_per_factor,
    K = K, N = N, M_ord = M_ord, X_ord = X_ord,
    V_BG = V_BG, Z_CT = Z_CT,
    K_BG = ncol(V_BG), q_covar = ncol(X_ord) - 1L,
    hypothesis = "H3"
  )
}


# ----------------------------------------------------------------------------
# H3 test for a single pathway.  Reuses the H2 score-test internals with
# the H3 precomp (which contains M^(3) and the H3 eigenvalues).
# ----------------------------------------------------------------------------
passage_h3 <- function(engine, Y, pathway, precomp = NULL,
                      Z_CT = NULL, V_BG = NULL,
                      K_BG = NULL,
                      bg_method = "size_matched",
                      bg_size = NULL,
                      weight_schemes = c("equal", "var", "range"),
                      run_burden = TRUE, run_spasset = TRUE,
                      spasset_thresholds = c(2, 5, 10, 25, 50, 100, 250, 500),
                      gene_names = NULL, pvalue_method = "davies",
                      verbose = FALSE) {
  # If no precomp, must build one (background fit + design + eigenvalues)
  if (is.null(precomp)) {
    if (is.null(V_BG)) {
      if (verbose) cat("passage_h3: fitting background engine ...\n")
      bgfit <- .fit_background_engine(
        engine, Y, pathway, K_BG = K_BG,
        bg_method = bg_method, bg_size = bg_size,
        gene_names = gene_names, verbose = verbose)
      V_BG <- bgfit$V_BG
    }
    if (verbose) cat("passage_h3: precomputing H3 eigenvalues ...\n")
    precomp <- passage_h3_precompute(engine, V_BG, Z_CT = Z_CT,
                                     verbose = verbose)
  }
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)

  # H3 tests share machinery with H2: the only difference is M_ord and the
  # precomputed eigenvalues, both of which live in `precomp`.  Re-use H2
  # score-test functions by passing the H3 precomp.

  joint_results <- list()
  for (ws in weight_schemes) {
    joint_results[[ws]] <- passage_h2_joint(
      engine, Y, P, precomp = precomp,
      weight_scheme = ws, pvalue_method = pvalue_method)
  }
  burden_result <- if (run_burden) {
    passage_h2_burden(engine, Y, P, precomp = precomp,
                      pvalue_method = pvalue_method)
  } else NULL
  spasset_result <- if (run_spasset) {
    passage_h2_spasset(engine, Y, P, precomp = precomp,
                       thresholds = spasset_thresholds,
                       pvalue_method = pvalue_method)
  } else NULL

  comp_p <- sapply(joint_results, function(x) x$p_acat)
  names(comp_p) <- paste0("joint_", names(joint_results))
  if (!is.null(burden_result)) comp_p <- c(comp_p, burden = burden_result$p)
  if (!is.null(spasset_result)) comp_p <- c(comp_p, spasset = spasset_result$p_acat)
  p_omnibus <- acat_combine(comp_p)

  out <- list(
    p_omnibus = p_omnibus, p_components = comp_p,
    joint = joint_results, burden = burden_result, spasset = spasset_result,
    pathway_size = length(P), K = engine$K, N = engine$N,
    K_BG = precomp$K_BG, q_covar = precomp$q_covar,
    hypothesis = "H3"
  )
  class(out) <- c("passage_h3_result", "passage_h1_result", "list")
  out
}
