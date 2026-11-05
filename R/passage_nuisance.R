# Measured covariates are preferred. Additional factors may be estimated ONLY
# from prespecified negative-control genes, never selected by pathway results.
passage_prepare_nuisance <- function(Y, X = NULL, control_genes = NULL,
                                      n_factors = 0L) {
  Y <- passage_check_y(Y)
  n_factors <- .passage_gp_integer(n_factors, "n_factors", 0L)
  if (!is.null(X) && !is.null(rownames(X)) && !is.null(rownames(Y)) &&
      !identical(rownames(X), rownames(Y))) stop("X spot IDs must match Y")
  X <- passage_prepare_design(X, nrow(Y))
  rownames(X) <- rownames(Y)
  qx <- qr(X)
  if (qx$rank != ncol(X) || nrow(Y) - qx$rank - n_factors < 3L)
    stop("design must have full rank and leave at least three residual degrees of freedom")
  controls <- integer()
  if (n_factors > 0L) {
    if (is.null(control_genes)) stop("n_factors > 0 requires prespecified negative-control genes")
    if (is.character(control_genes) && any(!control_genes %in% colnames(Y)))
      stop("unknown control genes")
    if (is.numeric(control_genes) && (anyNA(control_genes) ||
        any(!is.finite(control_genes)) || any(control_genes != as.integer(control_genes)) ||
        any(control_genes < 1L | control_genes > ncol(Y)))) stop("invalid control indices")
    controls <- passage_resolve_pathway(control_genes, colnames(Y))
    if (length(controls) < n_factors + 1L) stop("need more control genes than nuisance factors")
    residuals <- qr.resid(qx, Y[, controls, drop = FALSE])
    scales <- sqrt(colSums(residuals^2) / (nrow(Y) - qx$rank))
    keep <- is.finite(scales) & scales > sqrt(.Machine$double.eps)
    if (sum(keep) <= n_factors) stop("too few variable control genes")
    residuals <- sweep(residuals[, keep, drop = FALSE], 2L, scales[keep], "/")
    # Deterministic SVD makes the adjustment independent of randomization draws.
    s <- svd(residuals, nu = n_factors, nv = 0L)
    if (sum(s$d > max(s$d) * 1e-10) < n_factors) stop("control matrix rank is too small")
    factors <- s$u[, seq_len(n_factors), drop = FALSE]
    colnames(factors) <- paste0("control_PC", seq_len(n_factors))
    X <- cbind(X, factors)
    qx <- qr(X)
  } else if (!is.null(control_genes)) stop("control_genes supplied with n_factors=0")
  leverage <- rowSums(qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]^2)
  list(X = X, control_genes = colnames(Y)[controls],
    analysis_genes = colnames(Y)[setdiff(seq_len(ncol(Y)), controls)],
    diagnostics = list(design_rank = qx$rank, residual_df = nrow(Y) - qx$rank,
      max_leverage = max(leverage), leverage = leverage, n_factors = n_factors),
    assumptions = paste("Measured covariates must represent the intended adjustment.",
      "Control genes must share nuisance variation but lack the biological signal of interest;",
      "exclude them from the analysis universe. Estimated control factors are not known fixed covariates.",
      "Adjustment can remove real biology when controls or covariates are inappropriate."))
}
