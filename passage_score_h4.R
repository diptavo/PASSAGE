# ============================================================================
# passage_score_h4.R
#
# PASSAGE: Layer 2 / H4 score test (regional differential)
#
# Tests
#   H_0^(4):  sigma_{k,r1}^2 = sigma_{k,r2}^2  for all k,
#             i.e. the pathway has the same spatial variance components across
#             two pre-specified regions r1 and r2.
#
# Mean-model design X^(4) = X^(2) (intercept + cell-type proportions).
#
# The spatial kernel is region-stratified.  Concretely, instead of one
# Qtilde_k per factor, we use a block-diagonal stratified version:
#
#   Qtilde_k^(strat) = blockdiag( Qtilde_k restricted to r1 locations,
#                                 Qtilde_k restricted to r2 locations )
#
# and test whether the contribution of factor k differs between the two
# blocks.  Operationally, this becomes a SKAT-style test on the contrast
# kernel:
#
#   K_contrast = w1 * Q_r1 - w2 * Q_r2
#
# with weights to balance the two regions' contributions.
#
# Practically simpler implementation (used here):
#   Compute the H1 score statistic separately on region r1 and region r2,
#   and test whether the difference is significant via a permutation / 
#   variance-stabilised contrast.
#
# This file implements the simpler region-stratified approach.  The
# permutation calibration is the primary calibration method given the
# subtlety of the analytical null.
#
# Dependencies: Matrix, CompQuadForm, passage_score_h1.R, passage_score_h2.R
# ============================================================================


# ----------------------------------------------------------------------------
# Internal: compute a per-region score statistic using a restricted Qtilde.
#
# For region indicator `region_mask` (logical vector of length N), compute
#   Q_region_k = (M_region Y_region a_k)^T  Qtilde_k_region  (M_region Y_region a_k)
# where Qtilde_k_region is the Qtilde restricted to the region's submatrix.
# The restriction preserves only the rows/cols of locations within the region
# but inherits the conditional structure as if those locations alone existed.
#
# Returns the per-factor Q statistics for one region.
# ----------------------------------------------------------------------------
.q_per_factor_in_region <- function(engine, Y_P_ord_centered, A_P, region_mask_ord,
                                    eigvals_per_factor_region = NULL) {
  K <- engine$K
  Q_k_vec <- numeric(K)
  for (k in seq_len(K)) {
    Q_full <- engine$Qtilde_list[[k]]
    Q_sub <- Q_full[region_mask_ord, region_mask_ord, drop = FALSE]
    Y_sub <- Y_P_ord_centered[region_mask_ord, , drop = FALSE]
    z_k <- as.numeric(Y_sub %*% A_P[, k])
    Q_k_vec[k] <- as.numeric(crossprod(z_k, Q_sub %*% z_k))
  }
  Q_k_vec
}


# ----------------------------------------------------------------------------
# Permutation-based calibration of region contrast.
#
# Permutes the region labels among locations B times.  For each permutation,
# recomputes the contrast statistic.  Reports an empirical p-value.
#
# This is the conservative, slow-but-honest calibration.  For large studies,
# we could derive an analytical mixed-chi^2 reference under exchangeability,
# but the analytical version is left for v2.
# ----------------------------------------------------------------------------
passage_h4 <- function(engine, Y, pathway, region, Z_CT = NULL,
                      n_perm = 1000L,
                      contrast = c("difference", "ratio"),
                      weight_scheme = c("equal", "var", "range"),
                      gene_names = NULL,
                      verbose = FALSE) {
  contrast <- match.arg(contrast)
  weight_scheme <- match.arg(weight_scheme)
  region <- as.factor(region)
  stopifnot(length(region) == engine$N)
  lvls <- levels(region)
  if (length(lvls) != 2L) {
    stop("H4 currently supports two-region contrasts only (got ",
         length(lvls), " levels)")
  }

  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)
  K <- engine$K
  N <- engine$N
  A_P <- engine$A_hat[P, , drop = FALSE]

  # Reorder Y, region, Z_CT to Vecchia ordering
  ord <- engine$vecchia$ord
  Y_P_ord <- Y[ord, P, drop = FALSE]
  region_ord <- region[ord]

  # Residualise by intercept + Z_CT if provided (H4 baseline is H2-style)
  if (is.null(Z_CT)) {
    Y_P_ord_centered <- sweep(Y_P_ord, 2L, colMeans(Y_P_ord), FUN = "-")
  } else {
    X_ord <- cbind(1.0, Z_CT[ord, , drop = FALSE])
    qx <- qr(X_ord)
    fitted <- qr.fitted(qx, Y_P_ord)
    Y_P_ord_centered <- Y_P_ord - fitted
  }

  # Per-region per-factor Q
  mask1 <- region_ord == lvls[1]
  mask2 <- region_ord == lvls[2]
  n1 <- sum(mask1); n2 <- sum(mask2)
  if (min(n1, n2) < 10L) {
    warning(sprintf("Region '%s' has only %d locations; H4 calibration unreliable.",
                    if (n1 < n2) lvls[1] else lvls[2], min(n1, n2)))
  }

  Q1_vec <- .q_per_factor_in_region(engine, Y_P_ord_centered, A_P, mask1)
  Q2_vec <- .q_per_factor_in_region(engine, Y_P_ord_centered, A_P, mask2)

  # Normalise by region size so that statistic doesn't trivially depend on n_r
  Q1_norm <- Q1_vec / n1
  Q2_norm <- Q2_vec / n2

  # Weights
  w_k <- switch(weight_scheme,
    equal  = rep(1, K),
    var    = engine$theta_hat[, "sigma2"],
    range  = engine$theta_hat[, "phi"])

  T_obs <- switch(contrast,
    "difference" = sum(w_k * abs(Q1_norm - Q2_norm)),
    "ratio"      = sum(w_k * abs(log((Q1_norm + 1e-6) / (Q2_norm + 1e-6))))
  )

  # Permutation: shuffle region labels, recompute T_perm
  if (verbose) cat(sprintf("H4 permutation calibration: %d perms, n1=%d, n2=%d\n",
                           n_perm, n1, n2))
  T_perm_vec <- numeric(n_perm)
  set.seed(20260522)
  for (b in seq_len(n_perm)) {
    perm_idx <- sample.int(N)
    region_perm <- region_ord[perm_idx]
    m1 <- region_perm == lvls[1]
    m2 <- region_perm == lvls[2]
    Q1b <- .q_per_factor_in_region(engine, Y_P_ord_centered, A_P, m1) / sum(m1)
    Q2b <- .q_per_factor_in_region(engine, Y_P_ord_centered, A_P, m2) / sum(m2)
    T_perm_vec[b] <- switch(contrast,
      "difference" = sum(w_k * abs(Q1b - Q2b)),
      "ratio"      = sum(w_k * abs(log((Q1b + 1e-6) / (Q2b + 1e-6))))
    )
  }

  # Empirical p-value with continuity correction
  p_emp <- (1 + sum(T_perm_vec >= T_obs)) / (1 + n_perm)

  out <- list(
    p = p_emp,
    T_obs = T_obs,
    T_perm_quantiles = quantile(T_perm_vec, c(0.05, 0.5, 0.95, 0.99)),
    Q_per_region = list(region1 = Q1_vec, region2 = Q2_vec),
    Q_norm_per_region = list(region1 = Q1_norm, region2 = Q2_norm),
    region_levels = lvls,
    region_sizes = c(n1, n2),
    n_perm = n_perm,
    contrast = contrast,
    weight_scheme = weight_scheme,
    pathway_size = length(P),
    hypothesis = "H4"
  )
  class(out) <- c("passage_h4_result", "list")
  out
}


print.passage_h4_result <- function(x, ...) {
  cat("PASSAGE H4 test (regional differential)\n")
  cat(sprintf("  Pathway size : %d\n", x$pathway_size))
  cat(sprintf("  Regions      : %s (n=%d), %s (n=%d)\n",
              x$region_levels[1], x$region_sizes[1],
              x$region_levels[2], x$region_sizes[2]))
  cat(sprintf("  Contrast     : %s, weighting: %s\n",
              x$contrast, x$weight_scheme))
  cat(sprintf("  T_obs        : %.4g\n", x$T_obs))
  cat(sprintf("  Permutation p: %s  (B = %d)\n",
              format.pval(x$p, digits = 3), x$n_perm))
  invisible(x)
}
