if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) source(f)
} else {
  library(PASSAGE)
}

if (requireNamespace("GpGp", quietly = TRUE)) {
  set.seed(10292027)
  n <- 48L
  coords <- cbind(x = runif(n), y = runif(n))
  ids <- paste0("gp_spot_", seq_len(n))
  rownames(coords) <- ids
  patterns <- cbind(sin(4 * coords[, 1L]), cos(5 * coords[, 2L]))
  loadings <- rbind(c(1, 0.1), c(0.8, -0.2), c(1.2, 0.2), c(0.9, 0),
                    c(0.1, 1), c(-0.2, 0.8), c(0.2, 1.2), c(0, 0.9))
  Y <- patterns %*% t(loadings) + matrix(rnorm(n * 8L, sd = 0.15), n, 8L)
  dimnames(Y) <- list(ids, paste0("gp_gene_", seq_len(ncol(Y))))
  fit <- passage_fit_pathway(Y, colnames(Y), K = 2L, seed = 17L)
  saved_fit <- fit
  saved_seed <- .Random.seed
  gp <- passage_fit_pathway_gp(fit, coords, m_seq = c(5L, 10L),
                               max_iter = 5L, seed = 19L)
  stopifnot(identical(fit, saved_fit), identical(.Random.seed, saved_seed),
            nrow(gp$summary) == 2L, length(gp$gp_fits) == 2L,
            all(is.finite(gp$covariance_parameters$estimate)),
            all(is.finite(gp$summary$noise_variance)),
            identical(gp$scores_source, "unsmoothed_expression"),
            is.null(gp[["p", exact = TRUE]]),
            is.null(gp[["predictions", exact = TRUE]]))
  for (k in seq_len(fit$K)) {
    # GpGp sees the estimated unsmoothed component, not a reconstructed GP mean.
    stopifnot(identical(unname(gp$gp_fits[[k]]$y), unname(fit$scores[, k])),
              identical(gp$gp_fits[[k]]$X, fit$X))
    if (is.logical(gp$gp_fits[[k]]$conv) && length(gp$gp_fits[[k]]$conv) == 1L) {
      stopifnot(identical(gp$summary$converged[k], gp$gp_fits[[k]]$conv))
    }
  }

  # Check default covariance semantics against the actual installed GpGp kernel.
  # Its nugget parameter is relative; its alpha range is rho/sqrt(3) in PASSAGE.
  probe <- matrix(c(0, 0, 0.2, 0), ncol = 2L, byrow = TRUE)
  covparms <- c(2, 0.4, 0.3)
  covmat <- GpGp::matern15_isotropic(covparms, probe)
  rho <- sqrt(3) * covparms[2L]
  distance <- sqrt(sum((probe[1L, ] - probe[2L, ])^2))
  z <- sqrt(3) * distance / rho
  stopifnot(abs(covmat[1L, 1L] - 2.6) < 1e-10,
            abs(covmat[1L, 2L] - covparms[1L] * (1 + z) * exp(-z)) < 1e-10,
            max(abs(gp$summary$range_passage_rho - sqrt(3) * gp$summary$range_gp_alpha)) < 1e-12,
            max(abs(gp$summary$noise_variance - gp$summary$spatial_variance *
                    gp$summary$relative_nugget)) < 1e-12)

  # A row-order mismatch is an error, never an implicit reassignment of spots.
  stopifnot(inherits(try(passage_fit_pathway_gp(fit,
    coords[n:1L, , drop = FALSE]), silent = TRUE), "try-error"))

  # Tiny optimization budget may not converge; preserve the package's flag.
  # Do not assume optimizer convergence for a test to pass.
  one <- passage_fit_pathway(Y, colnames(Y), K = 1L, seed = 23L)
  rm(".Random.seed", envir = .GlobalEnv)
  capped <- passage_fit_pathway_gp(one, coords, m_seq = c(100L, 200L),
                                   max_iter = 1L, seed = 29L)
  stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
            identical(capped$configuration$m_seq_used, n - 1L),
            all(is.finite(capped$covariance_parameters$estimate)))
  conv <- capped$gp_fits[[1L]]$conv
  if (is.logical(conv) && length(conv) == 1L) {
    stopifnot(identical(capped$summary$converged[1L], conv))
    if (identical(conv, FALSE)) stopifnot(!capped$summary$converged[1L])
  }
  cat("Optional GpGp pathway-fit regression tests passed; one-iteration convergence flag:",
      capped$summary$converged[1L], "\n")
} else {
  cat("Optional GpGp pathway-fit regression tests skipped: GpGp unavailable\n")
}
