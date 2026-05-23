# ============================================================================
# passage_score_h2.R
#
# PASSAGE: Layer 2 / H2 score test (cell-type adjusted)
#
# Tests
#   H_0^(2):  sigma_k^2 = 0  for all k = 1, ..., K,
#             AFTER adjusting for cell-type proportions in the mean model.
#
# Mean-model design: X^(2) = [1_N, Z_CT], where Z_CT is the N x q matrix of
# (deconvolved or measured) cell-type proportions.
#
# The score statistic is identical in form to H1, with the centring projection
# M^(1) = I - 11^T/N replaced by M^(2) = I - X^(2) (X^(2)^T X^(2))^{-1} X^(2)^T.
#
# Per-factor eigenvalues of M^(2) Qtilde_k M^(2) are precomputed once per
# (engine, cell-type-design) pair and reused across pathway tests.
#
# Sensitivity to deconvolution choice should be assessed by running with two
# tools (CARD, cell2location) and reporting both H2 p-values.
#
# Dependencies: Matrix, CompQuadForm, passage_score_h1.R (sourced)
# ============================================================================

# Assumes passage_score_h1.R has been sourced for acat_combine, mixed_chi2_pvalue,
# .resolve_pathway helpers.

# ----------------------------------------------------------------------------
# Build the residual projection M^(2) from a covariate matrix.
# Assumes Z_CT is in ORIGINAL location indexing (matches input Y / locs);
# this function reorders to Vecchia ordering to match Qtilde.
# Returns:
#   M_ord   : N x N residual projection in Vecchia ordering (dense)
#   X_ord   : N x (q+1) design matrix [1, Z_CT] in Vecchia ordering
#   df      : N - rank(X) for the residual degrees of freedom
# ----------------------------------------------------------------------------
.build_M_h2 <- function(Z_CT, engine) {
  N <- engine$N
  Z_CT <- as.matrix(Z_CT)
  stopifnot(nrow(Z_CT) == N)
  ord <- engine$vecchia$ord
  Z_ord <- Z_CT[ord, , drop = FALSE]
  X_ord <- cbind(1.0, Z_ord)
  qx <- qr(X_ord)
  # M = I - X (X'X)^{-1} X'
  H <- qr.Q(qx)
  H <- H %*% t(H)
  M_ord <- diag(N) - H
  list(M_ord = M_ord, X_ord = X_ord, df = N - qx$rank,
       Z_ord = Z_ord)
}

# ----------------------------------------------------------------------------
# Precompute per-factor eigenvalues of M^(2) Qtilde_k M^(2).
#
# This is the analogue of passage_h1_precompute() with M^(2) in place of M^(1).
# Cost: same order as H1 (K dense N x N eigendecompositions).
# Must be recomputed if cell-type design changes.
# ----------------------------------------------------------------------------
passage_h2_precompute <- function(engine, Z_CT, verbose = TRUE,
                                  max_dense_N = 15000) {
  K <- engine$K; N <- engine$N
  if (N > max_dense_N) {
    warning(sprintf("N = %d > max_dense_N = %d; eigendecomposition may run out of memory.",
                    N, max_dense_N))
  }
  m2 <- .build_M_h2(Z_CT, engine)
  M_ord <- m2$M_ord

  if (verbose) {
    cat(sprintf("passage_h2_precompute: K=%d, N=%d, q_covar=%d\n",
                K, N, ncol(m2$X_ord) - 1L))
  }

  eigvals_per_factor <- vector("list", K)
  for (k in seq_len(K)) {
    if (verbose) cat(sprintf("  factor %d/%d ... ", k, K))
    t0 <- Sys.time()
    # M Qk M is dense; form it via two sparse-aware multiplications
    Qk <- engine$Qtilde_list[[k]]
    # M %*% Qk: dense N x N
    MQ <- M_ord %*% Qk
    # MQ %*% M: dense N x N
    MQM <- MQ %*% M_ord
    rm(MQ); gc()
    MQM <- (MQM + t(MQM)) / 2
    eig <- eigen(MQM, symmetric = TRUE, only.values = TRUE)
    eigvals_per_factor[[k]] <- pmax(eig$values, 0)
    rm(MQM)
    if (verbose) {
      cat(sprintf("done (%.1fs, max eig %.3e)\n",
                  as.numeric(difftime(Sys.time(), t0, units = "secs")),
                  max(eigvals_per_factor[[k]])))
    }
  }

  list(
    eigvals_per_factor = eigvals_per_factor,
    K = K, N = N,
    M_ord = M_ord,            # cache for residual computation in tests
    X_ord = m2$X_ord,
    Z_CT  = Z_CT,
    q_covar = ncol(m2$X_ord) - 1L,
    hypothesis = "H2"
  )
}


# ----------------------------------------------------------------------------
# H2 joint test: identical structure to passage_h1_joint, but residualisation
# uses M^(2) projection (the dense M_ord cached in precomp).
# ----------------------------------------------------------------------------
passage_h2_joint <- function(engine, Y, pathway, precomp,
                             weight_scheme = c("equal", "var", "invvar", "range"),
                             gene_names = NULL,
                             pvalue_method = "davies") {
  weight_scheme <- match.arg(weight_scheme)
  stopifnot(!is.null(precomp$M_ord),
            identical(precomp$hypothesis, "H2"))
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)
  K <- engine$K; N <- engine$N
  A_P <- engine$A_hat[P, , drop = FALSE]
  D_P <- engine$D_hat[P]

  ord <- engine$vecchia$ord
  Y_P_ord <- Y[ord, P, drop = FALSE]
  # Residualise via M^(2): R = M_ord %*% Y_P_ord
  R_ord <- precomp$M_ord %*% Y_P_ord

  Q_k_vec <- numeric(K); c_k_vec <- numeric(K); p_k_vec <- numeric(K)
  for (k in seq_len(K)) {
    a_k <- A_P[, k]
    z_k <- as.numeric(R_ord %*% a_k)
    Qk_mat <- engine$Qtilde_list[[k]]
    Q_k <- as.numeric(crossprod(z_k, Qk_mat %*% z_k))
    c_k <- sum(a_k * a_k * D_P)
    p_k <- mixed_chi2_pvalue(Q_k, c_k * precomp$eigvals_per_factor[[k]],
                             method = pvalue_method)
    Q_k_vec[k] <- Q_k; c_k_vec[k] <- c_k; p_k_vec[k] <- p_k
  }

  w_k <- switch(weight_scheme,
    equal  = rep(1, K),
    var    = engine$theta_hat[, "sigma2"],
    invvar = 1 / pmax(engine$theta_hat[, "sigma2"], 1e-6),
    range  = engine$theta_hat[, "phi"])

  Q_joint <- sum(w_k * Q_k_vec)
  combined_eigvals <- unlist(lapply(seq_len(K), function(k) {
    w_k[k] * c_k_vec[k] * precomp$eigvals_per_factor[[k]]
  }))
  p_pooled <- mixed_chi2_pvalue(Q_joint, combined_eigvals, method = pvalue_method)
  p_acat <- acat_combine(p_k_vec, weights = w_k)

  list(
    Q_joint = Q_joint, p_pooled = p_pooled, p_acat = p_acat,
    Q_per_factor = Q_k_vec, p_per_factor = p_k_vec,
    c_per_factor = c_k_vec, w_per_factor = w_k,
    weight_scheme = weight_scheme, pathway_size = length(P)
  )
}


# ----------------------------------------------------------------------------
# H2 burden test: 1_p / sqrt(p) direction, M^(2)-residualised data.
# ----------------------------------------------------------------------------
passage_h2_burden <- function(engine, Y, pathway, precomp,
                              kernel_factor = 1L,
                              gene_names = NULL,
                              pvalue_method = "davies") {
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)
  p_size <- length(P)
  D_P <- engine$D_hat[P]
  ord <- engine$vecchia$ord
  Y_P_ord <- Y[ord, P, drop = FALSE]
  R_ord <- precomp$M_ord %*% Y_P_ord

  a_burd <- rep(1 / sqrt(p_size), p_size)
  z_burd <- as.numeric(R_ord %*% a_burd)
  Qk_mat <- engine$Qtilde_list[[kernel_factor]]
  Q_burd <- as.numeric(crossprod(z_burd, Qk_mat %*% z_burd))
  c_burd <- mean(D_P)
  p_burd <- mixed_chi2_pvalue(Q_burd,
                              c_burd * precomp$eigvals_per_factor[[kernel_factor]],
                              method = pvalue_method)
  list(Q = Q_burd, p = p_burd, c = c_burd,
       kernel_factor = kernel_factor, pathway_size = p_size)
}


# ----------------------------------------------------------------------------
# H2 SpASSET subset search. Mirrors H1 with M^(2)-residualised data.
# ----------------------------------------------------------------------------
passage_h2_spasset <- function(engine, Y, pathway, precomp,
                               thresholds = c(2, 5, 10, 25, 50, 100, 250, 500),
                               weight_scheme = "equal",
                               scoring_kernel_factor = 1L,
                               gene_names = NULL,
                               pvalue_method = "davies") {
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)
  p_size <- length(P)
  D_P <- engine$D_hat[P]
  ord <- engine$vecchia$ord
  Y_P_ord <- Y[ord, P, drop = FALSE]
  R_ord <- precomp$M_ord %*% Y_P_ord

  Qk_score <- engine$Qtilde_list[[scoring_kernel_factor]]
  U_per_gene <- vapply(seq_len(p_size), function(j) {
    yj <- R_ord[, j]
    as.numeric(crossprod(yj, Qk_score %*% yj)) / max(D_P[j], 1e-12)
  }, numeric(1))

  gene_order <- order(U_per_gene, decreasing = TRUE)
  thresholds <- unique(c(thresholds[thresholds <= p_size], p_size))
  thresholds <- sort(thresholds[thresholds >= 2L])

  p_values <- numeric(length(thresholds))
  Q_values <- numeric(length(thresholds))
  for (i in seq_along(thresholds)) {
    r <- thresholds[i]
    sub_P <- P[gene_order[seq_len(r)]]
    jr <- passage_h2_joint(engine, Y, sub_P, precomp = precomp,
                           weight_scheme = weight_scheme,
                           pvalue_method = pvalue_method)
    p_values[i] <- jr$p_acat
    Q_values[i] <- jr$Q_joint
  }
  best_idx <- which.min(p_values)
  best_r <- thresholds[best_idx]
  best_subset <- P[gene_order[seq_len(best_r)]]

  list(
    Q_per_threshold = Q_values,
    p_per_threshold = setNames(p_values, paste0("r=", thresholds)),
    thresholds = thresholds, best_threshold = best_r,
    best_subset = best_subset, best_p = p_values[best_idx],
    p_acat = acat_combine(p_values), gene_scores = U_per_gene,
    gene_ranking = gene_order, pathway_size = p_size
  )
}


# ----------------------------------------------------------------------------
# Main H2 entry point.
# ----------------------------------------------------------------------------
passage_h2 <- function(engine, Y, pathway, Z_CT = NULL, precomp = NULL,
                       weight_schemes = c("equal", "var", "range"),
                       run_burden = TRUE, run_spasset = TRUE,
                       spasset_thresholds = c(2, 5, 10, 25, 50, 100, 250, 500),
                       gene_names = NULL, pvalue_method = "davies",
                       verbose = FALSE) {
  if (is.null(precomp)) {
    if (is.null(Z_CT)) {
      stop("Either precomp or Z_CT (cell-type proportions) must be provided")
    }
    if (verbose) cat("passage_h2: precomputing eigenvalues for H2 design...\n")
    precomp <- passage_h2_precompute(engine, Z_CT, verbose = verbose)
  }
  P <- .resolve_pathway(pathway, engine_G = engine$G, gene_names = gene_names)

  joint_results <- list()
  for (ws in weight_schemes) {
    joint_results[[ws]] <- passage_h2_joint(engine, Y, P, precomp,
                                            weight_scheme = ws,
                                            pvalue_method = pvalue_method)
  }
  burden_result <- if (run_burden) {
    passage_h2_burden(engine, Y, P, precomp, pvalue_method = pvalue_method)
  } else NULL
  spasset_result <- if (run_spasset) {
    passage_h2_spasset(engine, Y, P, precomp,
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
    hypothesis = "H2"
  )
  class(out) <- c("passage_h2_result", "passage_h1_result", "list")
  out
}
