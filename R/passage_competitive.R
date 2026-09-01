# Competitive and covariance-aware pathway spatial variability metrics.

passage_rank_bins <- function(x, n_bins = 10L) {
  x <- as.numeric(x)
  n_bins <- as.integer(n_bins)
  if (n_bins <= 1L || length(unique(x[is.finite(x)])) <= 1L) {
    return(rep(1L, length(x)))
  }
  out <- rep(NA_integer_, length(x))
  ok <- is.finite(x)
  r <- rank(x[ok], ties.method = "average")
  out[ok] <- pmax(1L, pmin(n_bins, ceiling(r / max(r) * n_bins)))
  out[!ok] <- 1L
  out
}

passage_make_gene_bins <- function(Y,
                                   gene_names = NULL,
                                   n_mean_bins = 10L,
                                   n_detect_bins = 10L,
                                   n_var_bins = 0L,
                                   detect_threshold = 0) {
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  mean_expr <- colMeans(Y)
  detection <- colMeans(Y > detect_threshold)
  variance <- apply(Y, 2L, stats::var)
  mean_bin <- passage_rank_bins(mean_expr, n_mean_bins)
  detect_bin <- passage_rank_bins(detection, n_detect_bins)
  if (n_var_bins > 1L) {
    var_bin <- passage_rank_bins(variance, n_var_bins)
    bin <- paste(mean_bin, detect_bin, var_bin, sep = ":")
  } else {
    var_bin <- rep(1L, ncol(Y))
    bin <- paste(mean_bin, detect_bin, sep = ":")
  }
  data.frame(
    gene = gene_names,
    index = seq_along(gene_names),
    mean_expr = mean_expr,
    detection = detection,
    variance = variance,
    mean_bin = mean_bin,
    detect_bin = detect_bin,
    var_bin = var_bin,
    bin = bin,
    stringsAsFactors = FALSE
  )
}

passage_matrix_cov <- function(X, center = TRUE) {
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) < 2L) {
    return(matrix(0, ncol(X), ncol(X)))
  }
  if (center) {
    X <- sweep(X, 2L, colMeans(X), "-")
  }
  crossprod(X) / max(1L, nrow(X) - 1L)
}

passage_shrink_cov <- function(E,
                               method = c("shrink", "sample", "diagonal"),
                               shrinkage = NULL) {
  method <- match.arg(method)
  E <- as.matrix(E)
  storage.mode(E) <- "double"
  m <- ncol(E)
  if (m == 0L) return(matrix(0, 0L, 0L))
  if (m == 1L) {
    return(matrix(stats::var(as.numeric(E)), 1L, 1L))
  }
  if (method == "shrink" && requireNamespace("corpcor", quietly = TRUE)) {
    return(as.matrix(corpcor::cov.shrink(E, verbose = FALSE)))
  }
  S <- passage_matrix_cov(E)
  D <- diag(diag(S), nrow = m, ncol = m)
  if (method == "sample") return(S)
  if (method == "diagonal") return(D)
  if (is.null(shrinkage)) {
    shrinkage <- min(0.95, max(0.02, m / max(1, nrow(E) + m)))
  }
  (1 - shrinkage) * S + shrinkage * D
}

passage_cov_eigen <- function(S) {
  S <- as.matrix(S)
  if (nrow(S) == 0L) return(numeric(0))
  vals <- tryCatch(eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) rep(NA_real_, nrow(S)))
  vals <- as.numeric(vals)
  vals[!is.finite(vals) | vals < 0] <- 0
  vals
}

passage_effective_rank_from_values <- function(lambda) {
  lambda <- as.numeric(lambda)
  lambda <- lambda[is.finite(lambda) & lambda > 0]
  if (length(lambda) == 0L) return(NA_real_)
  sum(lambda)^2 / sum(lambda^2)
}

passage_generalized_pve_values <- function(S_spatial, S_residual, ridge = 1e-6) {
  S_spatial <- as.matrix(S_spatial)
  S_total <- S_spatial + as.matrix(S_residual)
  m <- ncol(S_total)
  if (m == 0L) return(numeric(0))
  sc <- mean(diag(S_total))
  if (!is.finite(sc) || sc <= 0) sc <- 1
  S_total <- (S_total + t(S_total)) / 2 + diag(ridge * sc, m)
  S_spatial <- (S_spatial + t(S_spatial)) / 2
  vals <- tryCatch(eigen(qr.solve(S_total, S_spatial), only.values = TRUE)$values,
                   error = function(e) rep(NA_real_, m))
  vals <- as.numeric(Re(vals))
  vals[!is.finite(vals)] <- NA_real_
  pmin(pmax(vals, 0), 1)
}

passage_mean_abs_cor <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) < 2L) return(NA_real_)
  s <- apply(X, 2L, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (sum(keep) < 2L) return(NA_real_)
  C <- suppressWarnings(stats::cor(X[, keep, drop = FALSE]))
  vals <- abs(C[upper.tri(C)])
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) NA_real_ else mean(vals)
}

passage_pc1_fraction <- function(X) {
  X <- as.matrix(X)
  if (ncol(X) < 2L) return(NA_real_)
  s <- apply(X, 2L, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (sum(keep) < 2L) return(NA_real_)
  vals <- passage_cov_eigen(passage_matrix_cov(X[, keep, drop = FALSE]))
  if (sum(vals) <= 0) NA_real_ else max(vals) / sum(vals)
}

passage_make_expression_coherence_balance <- function(Y,
                                                      precomp = NULL,
                                                      X = NULL,
                                                      metrics = c("expr_mean_abs_cor", "expr_pc1_fraction")) {
  Y <- passage_check_y(Y)
  qrx <- if (!is.null(precomp)) precomp$qr else qr(passage_prepare_design(X, nrow(Y), intercept = TRUE))
  R <- passage_residualize_with_qr(Y, qrx)
  metrics <- match.arg(metrics, several.ok = TRUE)
  function(idx) {
    idx <- as.integer(idx)
    out <- numeric(length(metrics))
    names(out) <- metrics
    if ("expr_mean_abs_cor" %in% metrics) {
      out[["expr_mean_abs_cor"]] <- passage_mean_abs_cor(R[, idx, drop = FALSE])
    }
    if ("expr_pc1_fraction" %in% metrics) {
      out[["expr_pc1_fraction"]] <- passage_pc1_fraction(R[, idx, drop = FALSE])
    }
    out
  }
}

passage_factor_coherence <- function(a, method = c("signed", "absolute", "none")) {
  method <- match.arg(method)
  a <- as.numeric(a)
  den <- length(a) * sum(a^2)
  if (!is.finite(den) || den <= 0) return(0)
  if (method == "none") return(1)
  num <- if (method == "absolute") sum(abs(a))^2 else sum(a)^2
  pmin(pmax(num / den, 0), 1)
}

passage_pathway_covariance_metrics <- function(engine,
                                               Y,
                                               pathway,
                                               precomp = NULL,
                                               X = NULL,
                                               gene_names = NULL,
                                               covariance = c("shrink", "sample", "diagonal"),
                                               shrinkage = NULL,
                                               scale_floor_frac = 0.05,
                                               coherence = c("signed", "absolute", "none"),
                                               ridge = 1e-6) {
  covariance <- match.arg(covariance)
  coherence <- match.arg(coherence)
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  qrx <- if (!is.null(precomp)) precomp$qr else qr(passage_prepare_design(X, nrow(Y), intercept = TRUE))
  R_obs <- passage_residualize_with_qr(Y[, P, drop = FALSE], qrx)
  F_hat <- engine$V %*% t(engine$A[P, , drop = FALSE])
  F_res <- passage_residualize_with_qr(F_hat, qrx)
  E_res <- R_obs - F_res
  S_spatial <- passage_matrix_cov(F_res)
  S_residual <- passage_shrink_cov(E_res, method = covariance, shrinkage = shrinkage)
  tr_sp <- sum(diag(S_spatial))
  tr_res <- sum(diag(S_residual))
  den <- max(tr_sp + tr_res, .Machine$double.eps)
  lambda_sp <- passage_cov_eigen(S_spatial)
  lambda_res <- passage_cov_eigen(S_residual)
  gpve <- passage_generalized_pve_values(S_spatial, S_residual, ridge = ridge)

  df <- max(1L, nrow(Y) - qrx$rank)
  V_res <- passage_residualize_with_qr(engine$V, qrx)
  factor_var <- colSums(sweep(V_res, 2L, colMeans(V_res), "-")^2) / df
  A_P <- engine$A[P, , drop = FALSE]
  raw_factor <- as.numeric(factor_var * colSums(A_P^2))
  bbox <- apply(engine$coords, 2L, range)
  tissue_diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  rel_range <- engine$theta$effective_range / max(tissue_diameter, .Machine$double.eps)
  scale_weight <- rel_range / pmax(rel_range + scale_floor_frac, .Machine$double.eps)
  scale_weight <- pmin(pmax(scale_weight, 0), 1)
  coh <- vapply(seq_len(engine$K), function(k) {
    passage_factor_coherence(A_P[, k], method = coherence)
  }, numeric(1))
  epsv_num <- sum(raw_factor * scale_weight, na.rm = TRUE)
  cepsv_num <- sum(raw_factor * scale_weight * coh, na.rm = TRUE)

  sp_diag <- pmax(diag(S_spatial), 0)
  res_diag <- pmax(diag(S_residual), 0)
  prop_sv <- sp_diag / pmax(sp_diag + res_diag, .Machine$double.eps)
  rho_bar <- NA_real_
  if (ncol(E_res) >= 2L) {
    C_res <- suppressWarnings(stats::cor(E_res))
    rv <- C_res[upper.tri(C_res)]
    rv <- rv[is.finite(rv)]
    if (length(rv) > 0L) rho_bar <- mean(rv)
  }
  rho_pos <- if (is.finite(rho_bar)) max(rho_bar, 0) else NA_real_
  m <- length(P)
  residual_eff_genes_vif <- if (is.finite(rho_pos)) {
    pmin(m, pmax(1, m / (1 + (m - 1) * rho_pos)))
  } else {
    NA_real_
  }
  pc1 <- if (sum(lambda_sp) > 0) max(lambda_sp) / sum(lambda_sp) else NA_real_

  summary <- c(
    pathway_size = m,
    mean_propSV_conditional = mean(prop_sv),
    median_propSV_conditional = stats::median(prop_sv),
    cwPVE_trace = tr_sp / den,
    cwPVE_top = if (length(gpve) > 0L) max(gpve, na.rm = TRUE) else NA_real_,
    cwPVE_mean = if (length(gpve) > 0L) mean(gpve, na.rm = TRUE) else NA_real_,
    ePSV = epsv_num / den,
    cEPSV = cepsv_num / den,
    trace_spatial = tr_sp,
    trace_residual = tr_res,
    spatial_eff_rank = passage_effective_rank_from_values(lambda_sp),
    residual_eff_rank = passage_effective_rank_from_values(lambda_res),
    residual_mean_cor = rho_bar,
    residual_eff_genes_vif = residual_eff_genes_vif,
    pc1_spatial_fraction = pc1,
    mean_abs_spatial_cor = passage_mean_abs_cor(F_res),
    mean_factor_coherence = mean(coh),
    tissue_diameter = tissue_diameter
  )
  list(
    summary = summary,
    propSV = stats::setNames(prop_sv, gene_names[P]),
    generalized_pve = gpve,
    spatial_eigen = lambda_sp,
    residual_eigen = lambda_res,
    factor_table = data.frame(
      factor = seq_len(engine$K),
      factor_var = factor_var,
      raw_spatial = raw_factor,
      effective_range = engine$theta$effective_range,
      relative_range = rel_range,
      scale_weight = scale_weight,
      coherence = coh,
      ePSV_component = raw_factor * scale_weight,
      cEPSV_component = raw_factor * scale_weight * coh
    ),
    covariance = covariance,
    coherence = coherence
  )
}

passage_gene_spatial_stats <- function(engine,
                                       Y,
                                       precomp = NULL,
                                       X = NULL,
                                       gene_names = NULL) {
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  qrx <- if (!is.null(precomp)) precomp$qr else qr(passage_prepare_design(X, nrow(Y), intercept = TRUE))
  R_obs <- passage_residualize_with_qr(Y, qrx)
  F_hat <- engine$V %*% t(engine$A)
  F_res <- passage_residualize_with_qr(F_hat, qrx)
  E_res <- R_obs - F_res
  df <- max(1L, nrow(Y) - qrx$rank)
  center_ss <- function(M) colSums(sweep(M, 2L, colMeans(M), "-")^2) / df
  spatial_var <- center_ss(F_res)
  residual_var <- center_ss(E_res)
  propSV <- spatial_var / pmax(spatial_var + residual_var, .Machine$double.eps)
  data.frame(
    gene = gene_names,
    index = seq_along(gene_names),
    spatial_var = spatial_var,
    residual_var = residual_var,
    propSV = propSV,
    stringsAsFactors = FALSE
  )
}

passage_competitive_gene_stat_z <- function(pathway,
                                            gene_stat,
                                            gene_bins,
                                            gene_names = NULL,
                                            finite_population = FALSE) {
  if (is.null(gene_names)) gene_names <- gene_bins$gene
  P <- passage_resolve_pathway(pathway, gene_names)
  stat <- as.numeric(gene_stat)
  if (length(stat) != nrow(gene_bins)) stop("gene_stat must align with gene_bins")
  P <- P[is.finite(stat[P])]
  if (length(P) == 0L) stop("pathway has no finite gene statistics")
  obs <- mean(stat[P])
  m <- length(P)
  bins <- gene_bins$bin
  expected <- 0
  var0 <- 0
  bin_rows <- list()
  for (b in unique(bins[P])) {
    in_p <- P[bins[P] == b]
    mh <- length(in_p)
    bg <- which(bins == b & is.finite(stat))
    bg_ex <- setdiff(bg, P)
    if (length(bg_ex) >= 2L) bg <- bg_ex
    vals <- stat[bg]
    mu <- mean(vals)
    vv <- if (length(vals) >= 2L) stats::var(vals) else 0
    Nh <- length(vals)
    fpc <- if (finite_population && Nh > 1L) max(0, 1 - mh / Nh) else 1
    expected <- expected + (mh / m) * mu
    var0 <- var0 + (mh / m)^2 * vv / mh * fpc
    bin_rows[[b]] <- data.frame(bin = b, m_h = mh, N_h = Nh, mean_h = mu, var_h = vv)
  }
  z <- if (var0 > 0) (obs - expected) / sqrt(var0) else if (obs > expected) Inf else -Inf
  p <- stats::pnorm(z, lower.tail = FALSE)
  list(
    statistic = obs,
    expected = expected,
    variance = var0,
    z = z,
    p = p,
    pathway_size = m,
    bins = do.call(rbind, bin_rows)
  )
}

passage_gene_score_z <- function(engine,
                                 Y,
                                 precomp,
                                 gene_names = NULL,
                                 weight_scheme = c("equal", "var", "range", "invvar")) {
  weight_scheme <- match.arg(weight_scheme)
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  if (length(gene_names) != ncol(Y)) stop("gene_names must align with columns of Y")
  if (is.null(precomp) || is.null(precomp$qr)) stop("precomp from passage_h_precompute() is required")

  R <- passage_residualize_with_qr(Y, precomp$qr)
  w <- passage_factor_weights(engine, weight_scheme)
  q_joint <- mean_joint <- var_joint <- rep(0, ncol(Y))
  for (k in seq_len(engine$K)) {
    vc <- engine$vecchia[[k]]
    R_ord <- R[vc$ord, , drop = FALSE]
    q_base <- colSums(R_ord * as.matrix(engine$K_score[[k]] %*% R_ord))
    a2 <- engine$A[, k]^2
    ck <- a2 * engine$D
    qk <- a2 * q_base
    mean_k <- ck * precomp$moments[[k]]$trace_mk
    var_k <- 2 * ck^2 * precomp$moments[[k]]$trace_mkmk
    q_joint <- q_joint + w[[k]] * qk
    mean_joint <- mean_joint + w[[k]] * mean_k
    var_joint <- var_joint + (w[[k]]^2) * var_k
  }
  z <- (q_joint - mean_joint) / sqrt(pmax(var_joint, .Machine$double.eps))
  data.frame(
    gene = gene_names,
    index = seq_along(gene_names),
    score_Q_gene = q_joint,
    score_mean_gene = mean_joint,
    score_var_gene = var_joint,
    score_z_gene = z,
    stringsAsFactors = FALSE
  )
}

passage_robust_pathway_score_stat <- function(engine,
                                              Y,
                                              pathway,
                                              precomp,
                                              gene_names = NULL,
                                              weight_scheme = c("equal", "var", "range", "invvar"),
                                              gene_scores = NULL,
                                              residual_matrix = NULL,
                                              winsor_q = 0.01,
                                              corr = c("abs", "positive", "none"),
                                              output = c("statistic", "details")) {
  weight_scheme <- match.arg(weight_scheme)
  corr <- match.arg(corr)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  if (is.null(gene_scores)) {
    gene_scores <- passage_gene_score_z(
      engine, Y, precomp = precomp, gene_names = gene_names,
      weight_scheme = weight_scheme
    )
  }
  z <- as.numeric(gene_scores$score_z_gene[P])
  ok <- is.finite(z)
  P <- P[ok]
  z <- z[ok]
  if (!length(z)) stop("pathway has no finite gene-level score z values")
  if (length(z) >= 3L && is.finite(winsor_q) && winsor_q > 0) {
    qq <- stats::quantile(z, probs = c(winsor_q, 1 - winsor_q), na.rm = TRUE, names = FALSE)
    z <- pmin(pmax(z, qq[[1L]]), qq[[2L]])
  }
  m <- length(z)
  rho_bar <- 0
  if (m >= 2L && corr != "none") {
    R <- if (!is.null(residual_matrix)) {
      residual_matrix[, P, drop = FALSE]
    } else {
      passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
    }
    s <- apply(R, 2L, stats::sd)
    keep <- is.finite(s) & s > 1e-10
    if (sum(keep) >= 2L) {
      C <- suppressWarnings(stats::cor(R[, keep, drop = FALSE]))
      rv <- C[upper.tri(C)]
      rv <- rv[is.finite(rv)]
      if (length(rv)) {
        rho_bar <- if (corr == "abs") mean(abs(rv)) else max(mean(rv), 0)
      }
    }
  }
  rho_bar <- pmin(pmax(rho_bar, 0), 0.99)
  m_eff <- if (m > 1L && corr != "none") m / (1 + (m - 1) * rho_bar) else m
  m_eff <- pmin(m, pmax(1, m_eff))
  stat <- mean(z) * sqrt(m_eff)
  details <- list(
    statistic = stat,
    mean_gene_z = mean(z),
    median_gene_z = stats::median(z),
    max_gene_z = max(z),
    pathway_size = m,
    effective_size = m_eff,
    residual_mean_correlation = rho_bar,
    gene_table = data.frame(
      gene = gene_names[P],
      index = P,
      score_z_gene = z,
      stringsAsFactors = FALSE
    )
  )
  if (output == "statistic") stat else details
}

passage_standardize_gene_stat_by_bin <- function(stat,
                                                 gene_bins,
                                                 method = c("robust", "rank", "active"),
                                                 center = c("median", "mean"),
                                                 active_tail = 0.90,
                                                 min_bin_n = 8L) {
  method <- match.arg(method)
  center <- match.arg(center)
  stat <- as.numeric(stat)
  if (length(stat) != nrow(gene_bins)) stop("stat must align with gene_bins")
  if (!"bin" %in% names(gene_bins)) stop("gene_bins must contain a bin column")
  z <- rep(NA_real_, length(stat))
  ok_all <- is.finite(stat)
  rank_score <- function(v) {
    r <- rank(v, ties.method = "average")
    stats::qnorm(pmin(pmax((r - 0.5) / length(v), 1e-6), 1 - 1e-6))
  }
  active_score <- function(v) {
    p0 <- 1 - active_tail
    thr <- stats::quantile(v, probs = active_tail, na.rm = TRUE, names = FALSE, type = 8)
    (as.numeric(v >= thr) - p0) / sqrt(pmax(p0 * (1 - p0), .Machine$double.eps))
  }
  if (method == "rank") {
    global_rank <- rank_score(stat[ok_all])
    global_z <- rep(NA_real_, length(stat))
    global_z[ok_all] <- global_rank
    for (bb in unique(gene_bins$bin)) {
      idx <- which(gene_bins$bin == bb & ok_all)
      z[idx] <- if (length(idx) >= min_bin_n) rank_score(stat[idx]) else global_z[idx]
    }
    return(z)
  }
  if (method == "active") {
    global_active <- active_score(stat[ok_all])
    global_z <- rep(NA_real_, length(stat))
    global_z[ok_all] <- global_active
    for (bb in unique(gene_bins$bin)) {
      idx <- which(gene_bins$bin == bb & ok_all)
      z[idx] <- if (length(idx) >= min_bin_n) active_score(stat[idx]) else global_z[idx]
    }
    return(z)
  }
  global_center <- if (center == "median") stats::median(stat[ok_all]) else mean(stat[ok_all])
  global_scale <- if (center == "median") {
    1.4826 * stats::mad(stat[ok_all], constant = 1)
  } else {
    stats::sd(stat[ok_all])
  }
  if (!is.finite(global_scale) || global_scale <= 0) global_scale <- 1
  for (bb in unique(gene_bins$bin)) {
    idx <- which(gene_bins$bin == bb & ok_all)
    if (length(idx) >= min_bin_n) {
      cc <- if (center == "median") stats::median(stat[idx]) else mean(stat[idx])
      ss <- if (center == "median") {
        1.4826 * stats::mad(stat[idx], constant = 1)
      } else {
        stats::sd(stat[idx])
      }
      if (!is.finite(ss) || ss <= 0) {
        cc <- global_center
        ss <- global_scale
      }
    } else {
      cc <- global_center
      ss <- global_scale
    }
    in_bin <- which(gene_bins$bin == bb & ok_all)
    z[in_bin] <- (stat[in_bin] - cc) / ss
  }
  z
}

passage_gene_pERSA <- function(engine,
                               Y,
                               precomp,
                               gene_names = NULL,
                               gene_bins = NULL,
                               scale_floor_frac = 0.05,
                               reliability = c("snr", "none"),
                               transform = c("logit", "identity"),
                               standardization = c("robust", "rank", "active"),
                               active_tail = 0.90,
                               ridge_frac = 0.05) {
  reliability <- match.arg(reliability)
  transform <- match.arg(transform)
  standardization <- match.arg(standardization)
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  if (length(gene_names) != ncol(Y)) stop("gene_names must align with columns of Y")
  if (is.null(precomp) || is.null(precomp$qr)) stop("precomp from passage_h_precompute() is required")

  V_res <- passage_residualize_with_qr(engine$V, precomp$qr)
  F_res <- V_res %*% t(engine$A)
  R_obs <- passage_residualize_with_qr(Y, precomp$qr)
  E_res <- R_obs - F_res
  df <- max(1L, nrow(Y) - precomp$qr$rank)
  center_ss <- function(M) colSums(sweep(M, 2L, colMeans(M), "-")^2) / df
  C_V <- passage_matrix_cov(V_res)
  factor_var <- pmax(diag(C_V), 0)
  raw_factor <- sweep(engine$A^2, 2L, factor_var, "*")
  spatial_var <- center_ss(F_res)
  residual_var <- center_ss(E_res)
  total_var <- pmax(spatial_var + residual_var, 0)
  ridge <- ridge_frac * stats::median(total_var[is.finite(total_var) & total_var > 0])
  if (!is.finite(ridge) || ridge <= 0) ridge <- .Machine$double.eps

  bbox <- apply(engine$coords, 2L, range)
  tissue_diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  rel_range <- engine$theta$effective_range / max(tissue_diameter, .Machine$double.eps)
  scale_weight <- rel_range / pmax(rel_range + scale_floor_frac, .Machine$double.eps)
  scale_weight <- pmin(pmax(scale_weight, 0), 1)

  rel_weight <- matrix(1, nrow = nrow(raw_factor), ncol = ncol(raw_factor))
  if (reliability == "snr") {
    local_noise <- pmax(residual_var, 0) / max(1L, engine$K)
    rel_weight <- raw_factor / pmax(raw_factor + local_noise + ridge / max(1L, engine$K),
                                    .Machine$double.eps)
  }
  numerator <- rowSums(sweep(raw_factor * rel_weight, 2L, scale_weight, "*"), na.rm = TRUE)
  pERSA <- numerator / pmax(total_var + ridge, .Machine$double.eps)
  pERSA[!is.finite(pERSA)] <- NA_real_
  pERSA <- pmin(pmax(pERSA, 0), 1)
  transformed <- if (transform == "logit") {
    p <- pmin(pmax(pERSA, 1e-6), 1 - 1e-6)
    stats::qlogis(p)
  } else {
    pERSA
  }
  pERSA_z <- if (!is.null(gene_bins)) {
    passage_standardize_gene_stat_by_bin(
      transformed, gene_bins,
      method = standardization,
      active_tail = active_tail
    )
  } else {
    ok <- is.finite(transformed)
    if (standardization == "rank") {
      r <- rank(transformed[ok], ties.method = "average")
      out <- rep(NA_real_, length(transformed))
      out[ok] <- stats::qnorm(pmin(pmax((r - 0.5) / length(r), 1e-6), 1 - 1e-6))
      out
    } else if (standardization == "active") {
      p0 <- 1 - active_tail
      thr <- stats::quantile(transformed[ok], probs = active_tail, na.rm = TRUE, names = FALSE, type = 8)
      (as.numeric(transformed >= thr) - p0) / sqrt(pmax(p0 * (1 - p0), .Machine$double.eps))
    } else {
      cc <- stats::median(transformed[ok])
      ss <- 1.4826 * stats::mad(transformed[ok], constant = 1)
      if (!is.finite(ss) || ss <= 0) ss <- stats::sd(transformed[ok])
      if (!is.finite(ss) || ss <= 0) ss <- 1
      (transformed - cc) / ss
    }
  }
  data.frame(
    gene = gene_names,
    index = seq_along(gene_names),
    pERSA = pERSA,
    pERSA_transformed = transformed,
    pERSA_z = pERSA_z,
    pERSA_numerator = numerator,
    spatial_var = spatial_var,
    residual_var = residual_var,
    total_var = total_var,
    stringsAsFactors = FALSE
  )
}

passage_pERSA_pathway_score_stat <- function(engine,
                                             Y,
                                             pathway,
                                             precomp,
                                             gene_names = NULL,
                                             gene_activity = NULL,
                                             gene_bins = NULL,
                                             residual_matrix = NULL,
                                             winsor_q = 0.01,
                                             corr = c("abs", "positive", "none"),
                                             output = c("statistic", "details")) {
  corr <- match.arg(corr)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  if (is.null(gene_activity)) {
    gene_activity <- passage_gene_pERSA(
      engine, Y, precomp = precomp, gene_names = gene_names,
      gene_bins = gene_bins
    )
  }
  z <- as.numeric(gene_activity$pERSA_z[P])
  ok <- is.finite(z)
  P <- P[ok]
  z <- z[ok]
  if (!length(z)) stop("pathway has no finite pERSA z values")
  if (length(z) >= 3L && is.finite(winsor_q) && winsor_q > 0) {
    qq <- stats::quantile(z, probs = c(winsor_q, 1 - winsor_q), na.rm = TRUE, names = FALSE)
    z <- pmin(pmax(z, qq[[1L]]), qq[[2L]])
  }
  m <- length(z)
  rho_bar <- 0
  if (m >= 2L && corr != "none") {
    R <- if (!is.null(residual_matrix)) {
      residual_matrix[, P, drop = FALSE]
    } else {
      passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
    }
    s <- apply(R, 2L, stats::sd)
    keep <- is.finite(s) & s > 1e-10
    if (sum(keep) >= 2L) {
      C <- suppressWarnings(stats::cor(R[, keep, drop = FALSE]))
      rv <- C[upper.tri(C)]
      rv <- rv[is.finite(rv)]
      if (length(rv)) rho_bar <- if (corr == "abs") mean(abs(rv)) else max(mean(rv), 0)
    }
  }
  rho_bar <- pmin(pmax(rho_bar, 0), 0.99)
  m_eff <- if (m > 1L && corr != "none") m / (1 + (m - 1) * rho_bar) else m
  m_eff <- pmin(m, pmax(1, m_eff))
  stat <- mean(z) * sqrt(m_eff)
  positive <- pmax(z, 0)
  driver_weight <- if (sum(positive) > 0) positive / sum(positive) else rep(1 / m, m)
  details <- list(
    statistic = stat,
    mean_pERSA_z = mean(z),
    median_pERSA_z = stats::median(z),
    max_pERSA_z = max(z),
    pathway_size = m,
    effective_size = m_eff,
    residual_mean_correlation = rho_bar,
    gene_table = data.frame(
      gene = gene_names[P],
      index = P,
      pERSA = gene_activity$pERSA[P],
      pERSA_z = z,
      driver_weight = driver_weight,
      stringsAsFactors = FALSE
    )
  )
  if (output == "statistic") stat else details
}

passage_effective_pathway_size <- function(Y,
                                           pathway,
                                           precomp,
                                           gene_names = NULL,
                                           residual_matrix = NULL,
                                           corr = c("abs", "positive", "none")) {
  corr <- match.arg(corr)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) <= 1L || corr == "none") return(length(P))
  R <- if (!is.null(residual_matrix)) {
    residual_matrix[, P, drop = FALSE]
  } else {
    passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
  }
  s <- apply(R, 2L, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (sum(keep) < 2L) return(length(P))
  C <- suppressWarnings(stats::cor(R[, keep, drop = FALSE]))
  rv <- C[upper.tri(C)]
  rv <- rv[is.finite(rv)]
  if (!length(rv)) return(length(P))
  rho_bar <- if (corr == "abs") mean(abs(rv)) else max(mean(rv), 0)
  rho_bar <- pmin(pmax(rho_bar, 0), 0.99)
  m <- sum(keep)
  pmin(m, pmax(1, m / (1 + (m - 1) * rho_bar)))
}

passage_gene_activity_pathway_score_stat <- function(activity,
                                                     pathway,
                                                     Y,
                                                     precomp,
                                                     gene_names = NULL,
                                                     residual_matrix = NULL,
                                                     winsor_q = 0.01,
                                                     corr = c("abs", "positive", "none"),
                                                     output = c("statistic", "details")) {
  corr <- match.arg(corr)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(activity) != length(gene_names)) stop("activity must align with gene_names")
  z <- as.numeric(activity[P])
  ok <- is.finite(z)
  P <- P[ok]
  z <- z[ok]
  if (!length(z)) stop("pathway has no finite gene activity values")
  if (length(z) >= 3L && is.finite(winsor_q) && winsor_q > 0) {
    qq <- stats::quantile(z, probs = c(winsor_q, 1 - winsor_q), na.rm = TRUE, names = FALSE)
    z <- pmin(pmax(z, qq[[1L]]), qq[[2L]])
  }
  m_eff <- passage_effective_pathway_size(
    Y, P, precomp = precomp, gene_names = gene_names,
    residual_matrix = residual_matrix, corr = corr
  )
  stat <- mean(z) * sqrt(m_eff)
  positive <- pmax(z, 0)
  driver_weight <- if (sum(positive) > 0) positive / sum(positive) else rep(1 / length(z), length(z))
  details <- list(
    statistic = stat,
    mean_activity = mean(z),
    median_activity = stats::median(z),
    max_activity = max(z),
    pathway_size = length(z),
    effective_size = m_eff,
    gene_table = data.frame(
      gene = gene_names[P],
      index = P,
      activity = z,
      driver_weight = driver_weight,
      stringsAsFactors = FALSE
    )
  )
  if (output == "statistic") stat else details
}

passage_sparse_topk_pathway_score_stat <- function(engine,
                                                   Y,
                                                   pathway,
                                                   precomp,
                                                   gene_names = NULL,
                                                   gene_scores = NULL,
                                                   residual_matrix = NULL,
                                                   topk_grid = c(3L, 5L, 10L, 15L, 25L, 50L),
                                                   weight_scheme = c("equal", "var", "range", "invvar"),
                                                   output = c("statistic", "details")) {
  weight_scheme <- match.arg(weight_scheme)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  if (is.null(gene_scores)) {
    gene_scores <- passage_gene_score_z(
      engine, Y, precomp = precomp, gene_names = gene_names,
      weight_scheme = weight_scheme
    )
  }
  z_all <- as.numeric(gene_scores$score_z_gene[P])
  ok <- is.finite(z_all)
  P <- P[ok]
  z_all <- z_all[ok]
  if (!length(z_all)) stop("pathway has no finite gene score values")
  ord <- order(z_all, decreasing = TRUE)
  grid <- sort(unique(pmin(length(P), as.integer(topk_grid))))
  grid <- grid[grid >= 1L]
  vals <- vapply(grid, function(k) {
    keep <- P[ord[seq_len(k)]]
    passage_robust_pathway_score_stat(
      engine, Y, keep, precomp = precomp, gene_names = gene_names,
      gene_scores = gene_scores, residual_matrix = residual_matrix,
      weight_scheme = weight_scheme
    )
  }, numeric(1))
  best <- which.max(vals)
  details <- list(
    statistic = vals[[best]],
    best_k = grid[[best]],
    topk_grid = grid,
    score_by_k = stats::setNames(vals, paste0("top", grid)),
    selected_genes = gene_names[P[ord[seq_len(grid[[best]])]]]
  )
  if (output == "statistic") details$statistic else details
}

passage_gene_lse <- function(gene_score_z,
                             gene_bins = NULL,
                             max_iter = 80L,
                             tol = 1e-6) {
  z <- as.numeric(gene_score_z)
  ok <- is.finite(z)
  zw <- z
  if (sum(ok) >= 10L) {
    qq <- stats::quantile(z[ok], probs = c(0.005, 0.995), na.rm = TRUE, names = FALSE)
    zw[ok] <- pmin(pmax(z[ok], qq[[1L]]), qq[[2L]])
  }
  z_ok <- zw[ok]
  mu0 <- stats::median(z_ok)
  sd0 <- 1.4826 * stats::mad(z_ok, constant = 1)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- stats::sd(z_ok)
  if (!is.finite(sd0) || sd0 <= 0) sd0 <- 1
  mu1 <- as.numeric(stats::quantile(z_ok, 0.90, names = FALSE))
  if (!is.finite(mu1) || mu1 <= mu0) mu1 <- mu0 + sd0
  sd1 <- stats::sd(z_ok[z_ok >= mu1])
  if (!is.finite(sd1) || sd1 <= 0) sd1 <- sd0
  pi1 <- 0.10
  ll_old <- -Inf
  for (ii in seq_len(max_iter)) {
    d0 <- (1 - pi1) * stats::dnorm(z_ok, mean = mu0, sd = sd0)
    d1 <- pi1 * stats::dnorm(z_ok, mean = mu1, sd = sd1)
    den <- pmax(d0 + d1, .Machine$double.eps)
    post <- d1 / den
    pi1 <- pmin(pmax(mean(post), 0.01), 0.75)
    mu0 <- sum((1 - post) * z_ok) / pmax(sum(1 - post), .Machine$double.eps)
    mu1 <- sum(post * z_ok) / pmax(sum(post), .Machine$double.eps)
    if (mu1 < mu0 + 0.05) mu1 <- mu0 + 0.05
    sd0 <- sqrt(sum((1 - post) * (z_ok - mu0)^2) / pmax(sum(1 - post), .Machine$double.eps))
    sd1 <- sqrt(sum(post * (z_ok - mu1)^2) / pmax(sum(post), .Machine$double.eps))
    sd0 <- pmax(sd0, 0.25)
    sd1 <- pmax(sd1, 0.25)
    ll <- sum(log(den))
    if (is.finite(ll_old) && abs(ll - ll_old) < tol * (1 + abs(ll_old))) break
    ll_old <- ll
  }
  d0 <- (1 - pi1) * stats::dnorm(zw, mean = mu0, sd = sd0)
  d1 <- pi1 * stats::dnorm(zw, mean = mu1, sd = sd1)
  posterior <- d1 / pmax(d0 + d1, .Machine$double.eps)
  logit_post <- stats::qlogis(pmin(pmax(posterior, 1e-6), 1 - 1e-6))
  lse_z <- if (!is.null(gene_bins)) {
    passage_standardize_gene_stat_by_bin(logit_post, gene_bins, method = "rank")
  } else {
    r <- rank(logit_post[ok], ties.method = "average")
    out <- rep(NA_real_, length(z))
    out[ok] <- stats::qnorm(pmin(pmax((r - 0.5) / length(r), 1e-6), 1 - 1e-6))
    out
  }
  data.frame(
    posterior_spatial = posterior,
    logit_posterior_spatial = logit_post,
    lse_z = lse_z,
    stringsAsFactors = FALSE
  )
}

passage_pathway_factor_hsic_stat <- function(engine,
                                             Y,
                                             pathway,
                                             precomp,
                                             gene_names = NULL,
                                             residual_matrix = NULL,
                                             spatial_features = NULL,
                                             n_pc = 5L,
                                             output = c("statistic", "details")) {
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  R <- if (!is.null(residual_matrix)) residual_matrix[, P, drop = FALSE] else passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
  if (is.null(spatial_features)) {
    spatial_features <- passage_residualize_with_qr(engine$V, precomp$qr)
  }
  sx <- apply(R, 2L, stats::sd)
  R <- R[, is.finite(sx) & sx > 1e-10, drop = FALSE]
  if (ncol(R) == 0L) return(if (output == "statistic") NA_real_ else list(statistic = NA_real_))
  pc <- tryCatch({
    sv <- svd(scale(R, center = TRUE, scale = FALSE), nu = min(n_pc, ncol(R), nrow(R) - 1L), nv = 0)
    sv$u %*% diag(sv$d[seq_len(ncol(sv$u))], ncol(sv$u))
  }, error = function(e) scale(R[, seq_len(min(n_pc, ncol(R))), drop = FALSE], center = TRUE, scale = FALSE))
  X <- scale(pc)
  S <- scale(spatial_features)
  X[!is.finite(X)] <- 0
  S[!is.finite(S)] <- 0
  C <- crossprod(X, S) / max(1L, nrow(X) - 1L)
  stat <- sum(C^2) / pmax(sum(stats::var(X)) * sum(stats::var(S)), .Machine$double.eps)
  if (output == "statistic") stat else list(statistic = stat, cross_cov = C, n_pc = ncol(X))
}

passage_pathway_activity_hotspot_stat <- function(engine,
                                                  Y,
                                                  pathway,
                                                  precomp,
                                                  gene_names = NULL,
                                                  residual_matrix = NULL,
                                                  weight_scheme = c("equal", "var", "range", "invvar"),
                                                  output = c("statistic", "details")) {
  weight_scheme <- match.arg(weight_scheme)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  R <- if (!is.null(residual_matrix)) residual_matrix[, P, drop = FALSE] else passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
  s <- apply(R, 2L, stats::sd)
  keep <- is.finite(s) & s > 1e-10
  if (sum(keep) == 0L) return(if (output == "statistic") NA_real_ else list(statistic = NA_real_))
  Z <- scale(R[, keep, drop = FALSE])
  Z[!is.finite(Z)] <- 0
  activity <- rowMeans(Z)
  activity <- as.numeric(scale(activity))
  activity[!is.finite(activity)] <- 0
  w <- passage_factor_weights(engine, weight_scheme)
  q <- 0
  for (k in seq_len(engine$K)) {
    vc <- engine$vecchia[[k]]
    a <- activity[vc$ord]
    q <- q + w[[k]] * as.numeric(crossprod(a, engine$K_score[[k]] %*% a)) / length(a)
  }
  if (output == "statistic") q else list(statistic = q, activity = activity)
}

passage_sample_matched_indices <- function(pathway,
                                           gene_bins,
                                           gene_names = NULL,
                                           B = 999L,
                                           seed = NULL,
                                           exclude_pathway = TRUE,
                                           replace = TRUE,
                                           sampler = c("independent", "module", "mcmc", "importance", "knockoff"),
                                           module_col = "module",
                                           mcmc_burnin = 200L,
                                           mcmc_thin = 20L,
                                           importance_oversample = 20L,
                                           importance_temperature = 1,
                                           balance_fun = NULL,
                                           balance_target = NULL,
                                           balance_tol = NULL,
                                           balance_max_tries = 1L) {
  sampler <- match.arg(sampler)
  if (is.null(gene_names)) gene_names <- gene_bins$gene
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in gene_bins")
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  B <- as.integer(B)
  bins <- gene_bins$bin
  modules <- if (module_col %in% names(gene_bins)) gene_bins[[module_col]] else bins
  balance_max_tries <- as.integer(balance_max_tries)
  use_balance <- !is.null(balance_fun) && !is.null(balance_target) &&
    !is.null(balance_tol) && balance_max_tries > 1L
  if (use_balance) {
    balance_names <- names(balance_target)
    balance_target <- as.numeric(balance_target)
    names(balance_target) <- balance_names
    balance_tol <- as.numeric(balance_tol)
    if (length(balance_tol) == 1L) balance_tol <- rep(balance_tol, length(balance_target))
    if (is.null(names(balance_target))) stop("balance_target must be a named numeric vector")
    if (is.null(names(balance_tol))) names(balance_tol) <- names(balance_target)
    balance_tol <- balance_tol[names(balance_target)]
    if (any(!is.finite(balance_tol) | balance_tol <= 0)) {
      stop("balance_tol must contain positive finite values")
    }
  }
  sample_one <- function() {
    idx <- integer(length(P))
    for (jj in seq_along(P)) {
      cand <- which(bins == bins[P[jj]])
      if (sampler %in% c("module", "mcmc", "importance", "knockoff")) {
        cand_mod <- cand[modules[cand] == modules[P[jj]]]
        if (length(cand_mod) > 0L) cand <- cand_mod
      }
      if (exclude_pathway) {
        cand_ex <- setdiff(cand, P)
        if (length(cand_ex) > 0L) cand <- cand_ex
      }
      if (length(cand) == 0L) cand <- seq_len(nrow(gene_bins))
      idx[jj] <- sample(cand, size = 1L, replace = replace)
    }
    idx
  }
  balance_distance <- function(value) {
    value <- as.numeric(value[names(balance_target)])
    bad <- !is.finite(value)
    if (any(bad)) value[bad] <- Inf
    excess <- pmax(abs(value - balance_target) - balance_tol, 0)
    max(excess / balance_tol)
  }
  sample_candidate <- function() {
    if (!use_balance) return(list(idx = sample_one(), dist = NA_real_, value = NULL))
    best_idx <- NULL
    best_val <- NULL
    best_dist <- Inf
    best_tt <- 1L
    for (tt in seq_len(balance_max_tries)) {
      idx <- sample_one()
      val <- balance_fun(idx)
      dist <- balance_distance(val)
      if (is.finite(dist) && dist < best_dist) {
        best_dist <- dist
        best_idx <- idx
        best_val <- as.numeric(val[names(balance_target)])
        best_tt <- tt
      }
      if (is.finite(dist) && dist <= 0) break
    }
    list(idx = best_idx, dist = best_dist, value = best_val, attempts = best_tt)
  }
  if (sampler == "importance") {
    M <- max(B, as.integer(B) * max(1L, as.integer(importance_oversample)))
    candidates <- vector("list", M)
    dist <- numeric(M)
    values <- if (use_balance) {
      matrix(NA_real_, nrow = M, ncol = length(balance_target),
             dimnames = list(NULL, names(balance_target)))
    } else {
      NULL
    }
    for (mm in seq_len(M)) {
      cc <- sample_candidate()
      candidates[[mm]] <- cc$idx
      dist[mm] <- if (use_balance && is.finite(cc$dist)) cc$dist else 0
      if (use_balance) values[mm, ] <- cc$value
    }
    temp <- max(as.numeric(importance_temperature), .Machine$double.eps)
    weights <- exp(-dist / temp)
    weights[!is.finite(weights)] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(weights))
    keep <- sample(seq_along(candidates), size = B, replace = TRUE, prob = weights)
    out <- candidates[keep]
    attr(out, "weights") <- weights[keep] / mean(weights[keep])
    if (use_balance) {
      attr(out, "balance") <- list(
        target = balance_target,
        tolerance = balance_tol,
        accepted = dist[keep] <= 0,
        attempts = rep(NA_integer_, B),
        values = values[keep, , drop = FALSE],
        acceptance_rate = mean(dist[keep] <= 0),
        sampler = sampler
      )
    }
    return(out)
  }
  if (sampler == "mcmc") {
    current <- sample_candidate()
    current_idx <- current$idx
    current_dist <- if (use_balance && is.finite(current$dist)) current$dist else 0
    total_iter <- as.integer(mcmc_burnin) + B * max(1L, as.integer(mcmc_thin))
    out <- vector("list", B)
    accepted <- logical(total_iter)
    saved <- 0L
    for (it in seq_len(total_iter)) {
      prop_idx <- current_idx
      jj <- sample.int(length(P), 1L)
      cand <- which(bins == bins[P[jj]] & modules == modules[P[jj]])
      if (length(cand) == 0L) cand <- which(bins == bins[P[jj]])
      if (exclude_pathway) {
        cand_ex <- setdiff(cand, P)
        if (length(cand_ex) > 0L) cand <- cand_ex
      }
      if (length(cand) == 0L) cand <- seq_len(nrow(gene_bins))
      prop_idx[jj] <- sample(cand, size = 1L, replace = replace)
      prop_dist <- if (use_balance) balance_distance(balance_fun(prop_idx)) else 0
      log_acc <- -(prop_dist - current_dist) / max(as.numeric(importance_temperature), .Machine$double.eps)
      if (!is.finite(log_acc) || log(stats::runif(1)) <= min(0, log_acc)) {
        current_idx <- prop_idx
        current_dist <- prop_dist
        accepted[it] <- TRUE
      }
      if (it > as.integer(mcmc_burnin) && ((it - as.integer(mcmc_burnin)) %% max(1L, as.integer(mcmc_thin)) == 0L)) {
        saved <- saved + 1L
        out[[saved]] <- current_idx
      }
    }
    attr(out, "balance") <- list(
      target = if (use_balance) balance_target else NULL,
      tolerance = if (use_balance) balance_tol else NULL,
      accepted = accepted,
      attempts = rep(NA_integer_, B),
      values = NULL,
      acceptance_rate = mean(accepted),
      sampler = sampler
    )
    return(out)
  }
  out <- vector("list", B)
  accepted <- logical(B)
  attempts <- integer(B)
  balance_values <- if (use_balance) {
    matrix(NA_real_, nrow = B, ncol = length(balance_target),
           dimnames = list(NULL, names(balance_target)))
  } else {
    NULL
  }
  for (bb in seq_len(B)) {
    if (!use_balance) {
      if (sampler == "knockoff") {
        out[[bb]] <- sample_one()
      } else {
        out[[bb]] <- sample_one()
      }
      attempts[[bb]] <- 1L
      accepted[[bb]] <- TRUE
    } else {
      cc <- sample_candidate()
      out[[bb]] <- cc$idx
      attempts[[bb]] <- cc$attempts
      accepted[[bb]] <- is.finite(cc$dist) && cc$dist <= 0
      balance_values[bb, ] <- cc$value
    }
  }
  if (use_balance) {
    attr(out, "balance") <- list(
      target = balance_target,
      tolerance = balance_tol,
      accepted = accepted,
      attempts = attempts,
      values = balance_values,
      acceptance_rate = mean(accepted)
    )
  }
  out
}

passage_competitive_permutation_test <- function(pathway,
                                                 score_fun,
                                                 gene_bins,
                                                 gene_names = NULL,
                                                 B = 999L,
                                                 seed = NULL,
                                                 observed = NULL,
                                                 alternative = c("greater", "less", "two.sided"),
                                                 sampler = c("independent", "module", "mcmc", "importance", "knockoff"),
                                                 module_col = "module",
                                                 mcmc_burnin = 200L,
                                                 mcmc_thin = 20L,
                                                 importance_oversample = 20L,
                                                 importance_temperature = 1,
                                                 balance_fun = NULL,
                                                 balance_target = NULL,
                                                 balance_tol = NULL,
                                                 balance_max_tries = 1L) {
  alternative <- match.arg(alternative)
  sampler <- match.arg(sampler)
  if (is.null(gene_names)) gene_names <- gene_bins$gene
  P <- passage_resolve_pathway(pathway, gene_names)
  if (is.null(observed)) observed <- score_fun(P)
  matched <- passage_sample_matched_indices(
    pathway, gene_bins, gene_names, B = B, seed = seed,
    sampler = sampler, module_col = module_col,
    mcmc_burnin = mcmc_burnin, mcmc_thin = mcmc_thin,
    importance_oversample = importance_oversample,
    importance_temperature = importance_temperature,
    balance_fun = balance_fun, balance_target = balance_target,
    balance_tol = balance_tol, balance_max_tries = balance_max_tries
  )
  balance <- attr(matched, "balance")
  weights <- attr(matched, "weights")
  null <- vapply(matched, score_fun, numeric(1))
  ok <- is.finite(null)
  null <- null[ok]
  if (!is.null(weights)) {
    weights <- as.numeric(weights)[ok]
    weights[!is.finite(weights) | weights < 0] <- 0
    if (sum(weights) <= 0) weights <- rep(1, length(null))
  }
  if (length(null) == 0L) {
    return(list(statistic = observed, p = NA_real_, null = null, B = B, balance = balance))
  }
  if (is.null(weights)) {
    p <- switch(
      alternative,
      greater = (1 + sum(null >= observed)) / (1 + length(null)),
      less = (1 + sum(null <= observed)) / (1 + length(null)),
      two.sided = {
        mu <- mean(null)
        (1 + sum(abs(null - mu) >= abs(observed - mu))) / (1 + length(null))
      }
    )
    null_mean <- mean(null)
    null_sd <- stats::sd(null)
  } else {
    wsum <- sum(weights)
    null_mean <- sum(weights * null) / wsum
    null_sd <- sqrt(sum(weights * (null - null_mean)^2) / wsum)
    tail <- switch(
      alternative,
      greater = sum(weights[null >= observed]) / wsum,
      less = sum(weights[null <= observed]) / wsum,
      two.sided = sum(weights[abs(null - null_mean) >= abs(observed - null_mean)]) / wsum
    )
    p <- (1 + length(null) * tail) / (1 + length(null))
  }
  list(
    statistic = observed,
    p = p,
    null = null,
    null_weights = weights,
    null_mean = null_mean,
    null_sd = null_sd,
    enrichment = observed / pmax(null_mean, .Machine$double.eps),
    B = length(null),
    alternative = alternative,
    sampler = sampler,
    balance = balance
  )
}

passage_pathway_score_stat <- function(engine,
                                       Y,
                                       pathway,
                                       precomp,
                                       gene_names = NULL,
                                       weight_scheme = c("equal", "var", "range", "invvar"),
                                       output = c("z", "Q")) {
  weight_scheme <- match.arg(weight_scheme)
  output <- match.arg(output)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  h <- passage_score_test(
    engine = engine,
    Y = Y,
    pathway = pathway,
    precomp = precomp,
    weight_schemes = weight_scheme,
    run_burden = FALSE,
    run_spasset = FALSE,
    gene_names = gene_names,
    calibration = "moment"
  )
  joint <- h$joint[[weight_scheme]]
  if (is.null(joint)) {
    stop("weight_scheme not found in score result: ", weight_scheme)
  }
  q <- as.numeric(joint$Q)
  if (output == "Q") return(q)
  (q - as.numeric(joint$mean)) / sqrt(pmax(as.numeric(joint$var), .Machine$double.eps))
}

passage_conditional_competitive_test <- function(engine,
                                                Y,
                                                pathway,
                                                gene_bins = NULL,
                                                precomp = NULL,
                                                X = NULL,
                                                gene_names = NULL,
                                                statistic = c("score_sparse_topk_z", "score_lse_z", "score_factor_hsic", "score_coherence", "score_activity_hotspot", "score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z", "score_robust_z", "score_z", "score_Q", "mean_propSV_conditional", "cEPSV", "ePSV", "cwPVE_trace", "cwPVE_top"),
                                                B = 999L,
                                                seed = NULL,
                                                covariance = c("shrink", "sample", "diagonal"),
                                                coherence = c("signed", "absolute", "none"),
                                                score_weight_scheme = c("equal", "var", "range", "invvar"),
                                                sampler = c("independent", "module", "mcmc", "importance", "knockoff"),
                                                module_col = "module",
                                                mcmc_burnin = 200L,
                                                mcmc_thin = 20L,
                                                importance_oversample = 20L,
                                                importance_temperature = 1,
                                                balance_fun = NULL,
                                                balance_target = NULL,
                                                balance_tol = NULL,
                                                balance_max_tries = 1L) {
  statistic <- match.arg(statistic)
  covariance <- match.arg(covariance)
  coherence <- match.arg(coherence)
  score_weight_scheme <- match.arg(score_weight_scheme)
  sampler <- match.arg(sampler)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  if (is.null(gene_bins)) gene_bins <- passage_make_gene_bins(Y, gene_names = gene_names)
  if (statistic %in% c("score_sparse_topk_z", "score_lse_z", "score_factor_hsic", "score_coherence", "score_activity_hotspot", "score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z", "score_robust_z", "score_z", "score_Q")) {
    if (is.null(precomp)) {
      precomp <- passage_h_precompute(engine, X = X)
    }
    residual_matrix <- passage_residualize_with_qr(Y, precomp$qr)
    if (statistic == "score_sparse_topk_z") {
      gene_scores <- passage_gene_score_z(
        engine, Y, precomp = precomp, gene_names = gene_names,
        weight_scheme = score_weight_scheme
      )
      score_fun <- function(idx) {
        passage_sparse_topk_pathway_score_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          gene_scores = gene_scores, residual_matrix = residual_matrix,
          weight_scheme = score_weight_scheme
        )
      }
    } else if (statistic == "score_lse_z") {
      gene_scores <- passage_gene_score_z(
        engine, Y, precomp = precomp, gene_names = gene_names,
        weight_scheme = score_weight_scheme
      )
      lse <- passage_gene_lse(gene_scores$score_z_gene, gene_bins = gene_bins)
      score_fun <- function(idx) {
        passage_gene_activity_pathway_score_stat(
          lse$lse_z, idx, Y, precomp = precomp, gene_names = gene_names,
          residual_matrix = residual_matrix
        )
      }
    } else if (statistic == "score_factor_hsic") {
      spatial_features <- passage_residualize_with_qr(engine$V, precomp$qr)
      score_fun <- function(idx) {
        passage_pathway_factor_hsic_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          residual_matrix = residual_matrix, spatial_features = spatial_features
        )
      }
    } else if (statistic == "score_coherence") {
      fast_ctx <- passage_competitive_fast_context(engine, Y, precomp = precomp, gene_names = gene_names)
      score_fun <- function(idx) passage_fast_pc1_spatial_fraction(idx, fast_ctx)
    } else if (statistic == "score_activity_hotspot") {
      score_fun <- function(idx) {
        passage_pathway_activity_hotspot_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          residual_matrix = residual_matrix, weight_scheme = score_weight_scheme
        )
      }
    } else if (statistic %in% c("score_pERSA_z", "score_pERSA_rank_z", "score_ePSA_rank_z", "score_pERSA_active_z")) {
      persa_reliability <- if (statistic == "score_ePSA_rank_z") "none" else "snr"
      persa_standardization <- switch(
        statistic,
        score_pERSA_z = "robust",
        score_pERSA_active_z = "active",
        "rank"
      )
      gene_activity <- passage_gene_pERSA(
        engine, Y, precomp = precomp, gene_names = gene_names,
        gene_bins = gene_bins,
        reliability = persa_reliability,
        standardization = persa_standardization
      )
      residual_matrix <- passage_residualize_with_qr(Y, precomp$qr)
      score_fun <- function(idx) {
        passage_pERSA_pathway_score_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          gene_activity = gene_activity, gene_bins = gene_bins,
          residual_matrix = residual_matrix
        )
      }
    } else if (statistic == "score_robust_z") {
      gene_scores <- passage_gene_score_z(
        engine, Y, precomp = precomp, gene_names = gene_names,
        weight_scheme = score_weight_scheme
      )
      score_fun <- function(idx) {
      passage_robust_pathway_score_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          weight_scheme = score_weight_scheme, gene_scores = gene_scores
        )
      }
    } else {
      score_output <- if (statistic == "score_z") "z" else "Q"
      score_fun <- function(idx) {
        passage_pathway_score_stat(
          engine, Y, idx, precomp = precomp, gene_names = gene_names,
          weight_scheme = score_weight_scheme, output = score_output
        )
      }
    }
    obs <- score_fun(passage_resolve_pathway(pathway, gene_names))
    perm <- passage_competitive_permutation_test(pathway, score_fun, gene_bins, gene_names,
                                                 B = B, seed = seed, observed = obs,
                                                 sampler = sampler, module_col = module_col,
                                                 mcmc_burnin = mcmc_burnin, mcmc_thin = mcmc_thin,
                                                 importance_oversample = importance_oversample,
                                                 importance_temperature = importance_temperature,
                                                 balance_fun = balance_fun,
                                                 balance_target = balance_target,
                                                 balance_tol = balance_tol,
                                                 balance_max_tries = balance_max_tries)
    return(list(statistic = statistic, analytic = NULL, permutation = perm))
  }
  if (statistic == "mean_propSV_conditional") {
    gs <- passage_gene_spatial_stats(engine, Y, precomp = precomp, X = X, gene_names = gene_names)
    z <- passage_competitive_gene_stat_z(pathway, gs$propSV, gene_bins, gene_names)
    score_fun <- function(idx) mean(gs$propSV[idx], na.rm = TRUE)
    perm <- passage_competitive_permutation_test(pathway, score_fun, gene_bins, gene_names,
                                                 B = B, seed = seed, observed = z$statistic,
                                                 sampler = sampler, module_col = module_col,
                                                 mcmc_burnin = mcmc_burnin, mcmc_thin = mcmc_thin,
                                                 importance_oversample = importance_oversample,
                                                 importance_temperature = importance_temperature,
                                                 balance_fun = balance_fun,
                                                 balance_target = balance_target,
                                                 balance_tol = balance_tol,
                                                 balance_max_tries = balance_max_tries)
    return(list(statistic = statistic, analytic = z, permutation = perm))
  }
  metric_fun <- function(idx) {
    passage_pathway_covariance_metrics(
      engine, Y, idx, precomp = precomp, X = X, gene_names = gene_names,
      covariance = covariance, coherence = coherence
    )$summary[[statistic]]
  }
  obs <- metric_fun(passage_resolve_pathway(pathway, gene_names))
  perm <- passage_competitive_permutation_test(pathway, metric_fun, gene_bins, gene_names,
                                               B = B, seed = seed, observed = obs,
                                               sampler = sampler, module_col = module_col,
                                               mcmc_burnin = mcmc_burnin, mcmc_thin = mcmc_thin,
                                               importance_oversample = importance_oversample,
                                               importance_temperature = importance_temperature,
                                               balance_fun = balance_fun,
                                               balance_target = balance_target,
                                               balance_tol = balance_tol,
                                               balance_max_tries = balance_max_tries)
  list(statistic = statistic, analytic = NULL, permutation = perm)
}

passage_pathway_coherence_test <- function(engine,
                                          Y,
                                          pathway,
                                          gene_bins = NULL,
                                          precomp = NULL,
                                          X = NULL,
                                          gene_names = NULL,
                                          statistic = c("pc1_spatial_fraction", "mean_abs_spatial_cor"),
                                          B = 999L,
                                          seed = NULL,
                                          covariance = c("shrink", "sample", "diagonal"),
                                          coherence = c("signed", "absolute", "none"),
                                          balance_fun = NULL,
                                          balance_target = NULL,
                                          balance_tol = NULL,
                                          balance_max_tries = 1L) {
  statistic <- match.arg(statistic)
  covariance <- match.arg(covariance)
  coherence <- match.arg(coherence)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  if (is.null(gene_bins)) gene_bins <- passage_make_gene_bins(Y, gene_names = gene_names)
  metric_fun <- function(idx) {
    passage_pathway_covariance_metrics(
      engine, Y, idx, precomp = precomp, X = X, gene_names = gene_names,
      covariance = covariance, coherence = coherence
    )$summary[[statistic]]
  }
  obs <- metric_fun(passage_resolve_pathway(pathway, gene_names))
  perm <- passage_competitive_permutation_test(pathway, metric_fun, gene_bins, gene_names,
                                               B = B, seed = seed, observed = obs,
                                               balance_fun = balance_fun,
                                               balance_target = balance_target,
                                               balance_tol = balance_tol,
                                               balance_max_tries = balance_max_tries)
  list(statistic = statistic, permutation = perm)
}

passage_region_enrichment_test <- function(engine,
                                           Y,
                                           pathway,
                                           region,
                                           precomp,
                                           gene_names = NULL,
                                           n_perm = 999L,
                                           seed = NULL,
                                           weight_scheme = c("equal", "var", "range", "invvar")) {
  weight_scheme <- match.arg(weight_scheme)
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  if (length(region) != nrow(Y)) stop("length(region) must equal nrow(Y)")
  region <- factor(region)
  if (nlevels(region) != 2L) stop("region must have exactly two levels")
  if (is.null(precomp) || is.null(precomp$qr)) stop("precomp from passage_h_precompute() is required")

  R_P <- passage_residualize_with_qr(Y[, P, drop = FALSE], precomp$qr)
  A_P <- engine$A[P, , drop = FALSE]
  w <- passage_factor_weights(engine, weight_scheme)
  q_region <- function(region_factor) {
    out <- matrix(0, nrow = engine$K, ncol = 2L)
    colnames(out) <- levels(region_factor)
    for (k in seq_len(engine$K)) {
      vc <- engine$vecchia[[k]]
      region_ord <- region_factor[vc$ord]
      z <- as.numeric(R_P[vc$ord, , drop = FALSE] %*% A_P[, k])
      for (rr in seq_along(levels(region_factor))) {
        mask <- region_ord == levels(region_factor)[rr]
        z_sub <- z[mask]
        K_sub <- engine$K_score[[k]][mask, mask, drop = FALSE]
        out[k, rr] <- as.numeric(crossprod(z_sub, K_sub %*% z_sub)) / max(sum(mask), 1L)
      }
    }
    out
  }
  Q_obs <- q_region(region)
  T_obs <- sum(w * abs(Q_obs[, 2L] - Q_obs[, 1L]))
  if (!is.null(seed)) set.seed(seed)
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("n_perm must be at least 1")
  T_perm <- numeric(n_perm)
  for (bb in seq_len(n_perm)) {
    region_perm <- factor(sample(region), levels = levels(region))
    Qb <- q_region(region_perm)
    T_perm[bb] <- sum(w * abs(Qb[, 2L] - Qb[, 1L]))
  }
  list(
    p = (1 + sum(T_perm >= T_obs)) / (1 + n_perm),
    T_obs = T_obs,
    null = T_perm,
    Q_region = Q_obs,
    region_levels = levels(region),
    region_sizes = as.integer(table(region)),
    n_perm = n_perm,
    weight_scheme = weight_scheme
  )
}

passage_competitive_fast_context <- function(engine,
                                             Y,
                                             precomp = NULL,
                                             X = NULL,
                                             gene_names = NULL,
                                             scale_floor_frac = 0.05,
                                             coherence = c("signed", "absolute", "none")) {
  coherence <- match.arg(coherence)
  if (!inherits(engine, "passage_engine")) stop("engine must be a PASSAGE engine")
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  qrx <- if (!is.null(precomp)) precomp$qr else qr(passage_prepare_design(X, nrow(Y), intercept = TRUE))
  R_obs <- passage_residualize_with_qr(Y, qrx)
  V_res <- passage_residualize_with_qr(engine$V, qrx)
  F_res <- V_res %*% t(engine$A)
  E_res <- R_obs - F_res
  df <- max(1L, nrow(Y) - qrx$rank)
  center_ss <- function(M) colSums(sweep(M, 2L, colMeans(M), "-")^2) / df
  C_V <- passage_matrix_cov(V_res)
  factor_var <- diag(C_V)
  spatial_var <- center_ss(F_res)
  residual_var <- center_ss(E_res)
  bbox <- apply(engine$coords, 2L, range)
  tissue_diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  rel_range <- engine$theta$effective_range / max(tissue_diameter, .Machine$double.eps)
  scale_weight <- rel_range / pmax(rel_range + scale_floor_frac, .Machine$double.eps)
  scale_weight <- pmin(pmax(scale_weight, 0), 1)
  ev <- eigen((C_V + t(C_V)) / 2, symmetric = TRUE)
  vals <- pmax(ev$values, 0)
  C_half <- ev$vectors %*% (sqrt(vals) * t(ev$vectors))
  list(
    A = engine$A,
    C_V = C_V,
    C_half = C_half,
    factor_var = factor_var,
    spatial_var = spatial_var,
    residual_var = residual_var,
    scale_weight = scale_weight,
    coherence = coherence,
    gene_names = gene_names
  )
}

passage_fast_cEPSV <- function(idx, ctx) {
  idx <- as.integer(idx)
  A_P <- ctx$A[idx, , drop = FALSE]
  raw_factor <- as.numeric(ctx$factor_var * colSums(A_P^2))
  coh <- vapply(seq_len(ncol(A_P)), function(k) {
    passage_factor_coherence(A_P[, k], method = ctx$coherence)
  }, numeric(1))
  num <- sum(raw_factor * ctx$scale_weight * coh, na.rm = TRUE)
  den <- sum(ctx$spatial_var[idx] + ctx$residual_var[idx], na.rm = TRUE)
  num / pmax(den, .Machine$double.eps)
}

passage_fast_ePSV <- function(idx, ctx) {
  idx <- as.integer(idx)
  A_P <- ctx$A[idx, , drop = FALSE]
  raw_factor <- as.numeric(ctx$factor_var * colSums(A_P^2))
  num <- sum(raw_factor * ctx$scale_weight, na.rm = TRUE)
  den <- sum(ctx$spatial_var[idx] + ctx$residual_var[idx], na.rm = TRUE)
  num / pmax(den, .Machine$double.eps)
}

passage_fast_pc1_spatial_fraction <- function(idx, ctx) {
  idx <- as.integer(idx)
  A_P <- ctx$A[idx, , drop = FALSE]
  if (nrow(A_P) == 0L) return(NA_real_)
  small <- ctx$C_half %*% crossprod(A_P) %*% ctx$C_half
  vals <- passage_cov_eigen(small)
  if (sum(vals) <= 0) return(NA_real_)
  max(vals) / sum(vals)
}

passage_fast_mean_factor_coherence <- function(idx, ctx) {
  idx <- as.integer(idx)
  A_P <- ctx$A[idx, , drop = FALSE]
  if (nrow(A_P) == 0L) return(NA_real_)
  coh <- vapply(seq_len(ncol(A_P)), function(k) {
    passage_factor_coherence(A_P[, k], method = ctx$coherence)
  }, numeric(1))
  mean(coh, na.rm = TRUE)
}

passage_make_factor_coherence_balance <- function(ctx,
                                                  metrics = c("loading_pc1_fraction",
                                                              "loading_mean_factor_coherence")) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  function(idx) {
    idx <- as.integer(idx)
    out <- numeric(length(metrics))
    names(out) <- metrics
    if ("loading_pc1_fraction" %in% metrics) {
      out[["loading_pc1_fraction"]] <- passage_fast_pc1_spatial_fraction(idx, ctx)
    }
    if ("loading_mean_factor_coherence" %in% metrics) {
      out[["loading_mean_factor_coherence"]] <- passage_fast_mean_factor_coherence(idx, ctx)
    }
    out
  }
}

passage_combine_balance_functions <- function(...) {
  funs <- list(...)
  funs <- funs[!vapply(funs, is.null, logical(1))]
  if (length(funs) == 0L) stop("at least one balance function is required")
  function(idx) {
    unlist(lapply(funs, function(f) f(idx)), use.names = TRUE)
  }
}
