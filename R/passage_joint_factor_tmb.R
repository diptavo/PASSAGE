# Joint latent spatial factor engine for
#   Y_i = X_i B + A u(s_i) + e_i.
#
# This is the production-oriented wrapper around a TMB implementation with a
# Vecchia prior for independent latent spatial factors.

passage_joint_kernel_code <- function(kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian")) {
  kernel <- match.arg(kernel)
  switch(kernel,
    matern12 = 1L,
    exponential = 1L,
    matern32 = 2L,
    matern52 = 3L,
    gaussian = 4L
  )
}

passage_joint_default_cpp <- function(cpp_file = NULL) {
  if (!is.null(cpp_file)) {
    if (!file.exists(cpp_file)) stop("cpp_file does not exist: ", cpp_file)
    return(normalizePath(cpp_file, mustWork = TRUE))
  }
  sys <- system.file("tmb", "passage_joint_factor_tmb.cpp", package = "PASSAGE")
  if (nzchar(sys) && file.exists(sys)) return(sys)
  candidates <- c(
    file.path(getwd(), "inst", "tmb", "passage_joint_factor_tmb.cpp"),
    file.path(getwd(), "inst", "tmb", "passage_joint_factor_tmb.cpp"),
    file.path(getwd(), "PASSAGE", "inst", "tmb", "passage_joint_factor_tmb.cpp"),
    file.path(dirname(getwd()), "PASSAGE", "inst", "tmb", "passage_joint_factor_tmb.cpp")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) return(normalizePath(hit[[1L]], mustWork = TRUE))
  stop("Could not locate passage_joint_factor_tmb.cpp")
}

passage_compile_joint_factor_tmb <- function(cpp_file = NULL,
                                             build_dir = file.path(tempdir(), "SpaPath_tmb"),
                                             force = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("TMB is not installed. Install TMB in the R environment used for joint factor fitting.")
  }
  src <- passage_joint_default_cpp(cpp_file)
  dir.create(build_dir, recursive = TRUE, showWarnings = FALSE)
  dst <- file.path(build_dir, basename(src))
  if (force || !file.exists(dst) || file.info(src)$mtime > file.info(dst)$mtime) {
    file.copy(src, dst, overwrite = TRUE)
  }

  dll <- tools::file_path_sans_ext(basename(dst))
  so_file <- file.path(build_dir, paste0(dll, .Platform$dynlib.ext))
  if (force || !file.exists(so_file) || file.info(dst)$mtime > file.info(so_file)$mtime) {
    oldwd <- getwd()
    on.exit(setwd(oldwd), add = TRUE)
    setwd(build_dir)
    TMB::compile(basename(dst))
  }
  loaded_paths <- vapply(getLoadedDLLs(), function(x) x[["path"]], character(1))
  if (!normalizePath(so_file, mustWork = TRUE) %in% normalizePath(loaded_paths, mustWork = FALSE)) {
    dyn.load(so_file)
  }
  dll
}

passage_joint_site_order_maxmin <- function(coords) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if (n <= 2L) return(seq_len(n))
  center <- colMeans(coords)
  first <- which.min(rowSums((coords - matrix(center, n, ncol(coords), byrow = TRUE))^2))
  ord <- integer(n)
  ord[[1L]] <- first
  min_dist <- rep(Inf, n)
  for (tt in 2:n) {
    last <- ord[[tt - 1L]]
    d <- sqrt(rowSums((coords - matrix(coords[last, ], n, ncol(coords), byrow = TRUE))^2))
    min_dist <- pmin(min_dist, d)
    min_dist[ord[seq_len(tt - 1L)]] <- -Inf
    ord[[tt]] <- which.max(min_dist)
  }
  ord
}

passage_joint_vecchia_neighbors <- function(coords, m = 30L) {
  coords <- as.matrix(coords)
  n <- nrow(coords)
  m <- min(as.integer(m), max(0L, n - 1L))
  if (m < 1L) stop("m must be at least 1 and smaller than nrow(coords)")
  NN <- matrix(-1L, nrow = n, ncol = m)
  for (i in seq_len(n)) {
    if (i == 1L) next
    prev <- seq_len(i - 1L)
    d <- sqrt(rowSums((coords[prev, , drop = FALSE] -
      matrix(coords[i, ], length(prev), ncol(coords), byrow = TRUE))^2))
    keep <- order(d, seq_along(d))[seq_len(min(m, length(prev)))]
    NN[i, seq_along(keep)] <- prev[keep] - 1L
  }
  NN
}

passage_joint_near_psd <- function(S, eps = 1e-8) {
  S <- (as.matrix(S) + t(as.matrix(S))) / 2
  ee <- eigen(S, symmetric = TRUE)
  vals <- pmax(ee$values, eps)
  out <- ee$vectors %*% (vals * t(ee$vectors))
  (out + t(out)) / 2
}

passage_joint_lower_triangular_from_pca <- function(R, K) {
  p <- ncol(R)
  K <- min(as.integer(K), p)
  S <- crossprod(R) / max(1L, nrow(R))
  eig <- eigen(passage_joint_near_psd(S), symmetric = TRUE)
  lam <- pmax(eig$values[seq_len(K)], 1e-8)
  A_pca <- eig$vectors[, seq_len(K), drop = FALSE] %*% diag(sqrt(lam), K, K)
  qr_obj <- qr(t(A_pca))
  A_lower <- t(qr.R(qr_obj))[, seq_len(K), drop = FALSE]
  for (k in seq_len(K)) {
    if (A_lower[k, k] < 0) A_lower[, k] <- -A_lower[, k]
  }
  A_lower[upper.tri(A_lower)] <- 0
  rownames(A_lower) <- colnames(R)
  colnames(A_lower) <- paste0("factor_", seq_len(K))
  A_lower
}

passage_joint_pack_A_for_tmb <- function(A) {
  p <- nrow(A)
  K <- ncol(A)
  log_A_diag <- log(pmax(diag(A[seq_len(K), , drop = FALSE]), 1e-8))
  A_free <- numeric(0)
  for (k in seq_len(K)) {
    if (k < p) {
      A_free <- c(A_free, A[(k + 1L):p, k])
    }
  }
  list(log_A_diag = log_A_diag, A_free = A_free)
}

passage_joint_initialize <- function(Y, X, coords, K) {
  fit0 <- stats::lm.fit(x = X, y = Y)
  fitted0 <- X %*% fit0$coefficients
  R <- Y - fitted0
  A <- passage_joint_lower_triangular_from_pca(R, K)
  S <- crossprod(R) / max(1L, nrow(R))
  resid_cov <- S - A %*% t(A)
  tau2 <- pmax(diag(resid_cov), 0.05 * mean(diag(S)), 1e-4)
  Dinv <- diag(1 / tau2, ncol(Y), ncol(Y))
  M <- passage_joint_near_psd(t(A) %*% Dinv %*% A)
  P <- solve(M, t(A) %*% Dinv)
  U_start <- R %*% t(P)
  d <- as.matrix(stats::dist(coords))
  med_d <- stats::median(d[d > 0])
  list(B = fit0$coefficients, A = A, tau2 = tau2, phi = rep(med_d, K), U = U_start)
}

passage_fit_joint_factor_tmb <- function(Y,
                                         coords,
                                         K_fit = 2L,
                                         X = NULL,
                                         m = 30L,
                                         kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                         site_order = c("maxmin", "input"),
                                         cpp_file = NULL,
                                         tmb_build_dir = file.path(tempdir(), "SpaPath_tmb"),
                                         compile = TRUE,
                                         force_compile = FALSE,
                                         min_tau2 = 1e-3,
                                         nlminb_control = list(iter.max = 200, eval.max = 400),
                                         silent = FALSE) {
  if (!requireNamespace("TMB", quietly = TRUE)) {
    stop("TMB is not installed. Install TMB in the R environment used for joint factor fitting.")
  }
  kernel <- match.arg(kernel)
  site_order <- match.arg(site_order)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  K_fit <- as.integer(K_fit)
  if (K_fit < 1L) stop("K_fit must be at least 1")
  if (K_fit > ncol(Y)) stop("K_fit cannot exceed ncol(Y)")

  ord <- if (site_order == "maxmin") passage_joint_site_order_maxmin(coords) else seq_len(nrow(Y))
  Y_ord <- Y[ord, , drop = FALSE]
  coords_ord <- coords[ord, , drop = FALSE]
  X_ord <- X[ord, , drop = FALSE]
  NN <- passage_joint_vecchia_neighbors(coords_ord, m = m)
  init <- passage_joint_initialize(Y_ord, X_ord, coords_ord, K_fit)
  A_pack <- passage_joint_pack_A_for_tmb(init$A)

  if (compile) {
    dll <- passage_compile_joint_factor_tmb(
      cpp_file = cpp_file,
      build_dir = tmb_build_dir,
      force = force_compile
    )
  } else {
    dll <- tools::file_path_sans_ext(basename(passage_joint_default_cpp(cpp_file)))
  }

  data <- list(
    Y = Y_ord,
    coords = coords_ord,
    X = X_ord,
    NN = NN,
    min_tau2 = min_tau2,
    kernel_code = passage_joint_kernel_code(kernel)
  )
  parameters <- list(
    B = init$B,
    log_A_diag = A_pack$log_A_diag,
    A_free = A_pack$A_free,
    log_phi = log(init$phi),
    log_tau2 = log(pmax(init$tau2 - min_tau2, min_tau2)),
    U = init$U
  )
  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "U",
    DLL = dll,
    silent = silent
  )
  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = nlminb_control
  )
  rep <- obj$report()
  sdrep <- tryCatch(TMB::sdreport(obj), error = function(e) NULL)
  rownames(rep$A) <- colnames(Y)
  colnames(rep$A) <- paste0("factor_", seq_len(K_fit))
  rownames(rep$B) <- colnames(X)
  colnames(rep$B) <- colnames(Y)
  names(rep$phi) <- paste0("factor_", seq_len(K_fit))
  names(rep$tau2) <- colnames(Y)
  colnames(rep$U) <- paste0("factor_", seq_len(K_fit))

  U <- matrix(NA_real_, nrow(Y), K_fit)
  U[ord, ] <- rep$U
  colnames(U) <- colnames(rep$U)
  rownames(U) <- rownames(Y)
  fitted_fixed <- X %*% rep$B
  fitted_spatial <- U %*% t(rep$A)
  fitted_values <- fitted_fixed + fitted_spatial
  colnames(fitted_fixed) <- colnames(Y)
  colnames(fitted_spatial) <- colnames(Y)
  colnames(fitted_values) <- colnames(Y)
  residuals <- Y - fitted_values

  loglik <- -as.numeric(opt$objective)
  npar <- length(obj$par)
  nobs <- length(Y)
  out <- list(
    engine = "joint_factor_tmb_vecchia_v1",
    opt = opt,
    obj = obj,
    sdreport = sdrep,
    B = rep$B,
    A = rep$A,
    phi = rep$phi,
    tau2 = rep$tau2,
    U = U,
    U_ordered = rep$U,
    fitted_fixed = fitted_fixed,
    fitted_spatial = fitted_spatial,
    fitted_values = fitted_values,
    residuals = residuals,
    site_order = ord,
    neighbors = NN,
    init = init,
    data = data,
    X = X,
    coords = coords,
    Y = Y,
    K = K_fit,
    m = m,
    kernel = kernel,
    min_tau2 = min_tau2,
    logLik = loglik,
    npar = npar,
    nobs = nobs,
    AIC = 2 * npar - 2 * loglik,
    BIC = log(nobs) * npar - 2 * loglik,
    call = match.call()
  )
  class(out) <- c("passage_joint_factor_tmb", "list")
  out
}

print.passage_joint_factor_tmb <- function(x, ...) {
  cat("PASSAGE joint factor TMB engine\n")
  cat("  type:", x$engine, "\n")
  cat("  N:", nrow(x$Y), " G:", ncol(x$Y), " K:", x$K, "\n")
  cat("  kernel:", x$kernel, " m:", x$m, "\n")
  cat("  convergence:", x$opt$convergence, " objective:", signif(x$opt$objective, 5), "\n")
  cat("  AIC:", signif(x$AIC, 5), " BIC:", signif(x$BIC, 5), "\n")
  invisible(x)
}

passage_joint_factor_pve <- function(fit,
                                     scale_floor_frac = 0.05,
                                     coherence = c("signed", "absolute", "none"),
                                     ridge = 1e-6) {
  if (!inherits(fit, "passage_joint_factor_tmb")) {
    stop("fit must inherit from passage_joint_factor_tmb")
  }
  coherence <- match.arg(coherence)
  A <- as.matrix(fit$A)
  S_spatial <- A %*% t(A)
  S_residual <- diag(fit$tau2, nrow = nrow(A), ncol = nrow(A))
  rownames(S_residual) <- colnames(S_residual) <- rownames(A)
  tr_sp <- sum(diag(S_spatial))
  tr_res <- sum(diag(S_residual))
  den <- max(tr_sp + tr_res, .Machine$double.eps)
  lambda_sp <- passage_cov_eigen(S_spatial)
  lambda_res <- passage_cov_eigen(S_residual)
  gpve <- passage_generalized_pve_values(S_spatial, S_residual, ridge = ridge)
  prop_sv <- pmax(diag(S_spatial), 0) /
    pmax(pmax(diag(S_spatial), 0) + pmax(diag(S_residual), 0), .Machine$double.eps)

  bbox <- apply(fit$coords, 2L, range)
  tissue_diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  effective_range <- vapply(fit$phi, passage_effective_range, numeric(1), kernel = fit$kernel)
  rel_range <- effective_range / max(tissue_diameter, .Machine$double.eps)
  scale_weight <- pmin(pmax(rel_range / pmax(rel_range + scale_floor_frac, .Machine$double.eps), 0), 1)
  raw_factor <- colSums(A^2)
  coh <- vapply(seq_len(ncol(A)), function(k) {
    passage_factor_coherence(A[, k], method = coherence)
  }, numeric(1))
  epsv_num <- sum(raw_factor * scale_weight, na.rm = TRUE)
  cepsv_num <- sum(raw_factor * scale_weight * coh, na.rm = TRUE)
  C_res <- suppressWarnings(stats::cov2cor(S_residual))
  rv <- C_res[upper.tri(C_res)]
  rv <- rv[is.finite(rv)]
  rho_bar <- if (length(rv)) mean(rv) else NA_real_

  summary <- c(
    pathway_size = nrow(A),
    mean_propSV = mean(prop_sv),
    median_propSV = stats::median(prop_sv),
    cwPVE_trace = tr_sp / den,
    cwPVE_top = if (length(gpve)) max(gpve, na.rm = TRUE) else NA_real_,
    cwPVE_mean = if (length(gpve)) mean(gpve, na.rm = TRUE) else NA_real_,
    ePSV = epsv_num / den,
    cEPSV = cepsv_num / den,
    trace_spatial = tr_sp,
    trace_residual = tr_res,
    spatial_eff_rank = passage_effective_rank_from_values(lambda_sp),
    residual_eff_rank = passage_effective_rank_from_values(lambda_res),
    residual_mean_cor = rho_bar,
    pc1_spatial_fraction = if (sum(lambda_sp) > 0) max(lambda_sp) / sum(lambda_sp) else NA_real_,
    tissue_diameter = tissue_diameter
  )
  list(
    summary = summary,
    propSV = stats::setNames(prop_sv, rownames(A)),
    generalized_pve = gpve,
    spatial_gene_cov = S_spatial,
    residual_gene_cov = S_residual,
    spatial_eigen = lambda_sp,
    residual_eigen = lambda_res,
    factor_table = data.frame(
      factor = seq_len(ncol(A)),
      raw_spatial = raw_factor,
      phi = fit$phi,
      effective_range = effective_range,
      relative_range = rel_range,
      scale_weight = scale_weight,
      coherence = coh,
      ePSV_component = raw_factor * scale_weight,
      cEPSV_component = raw_factor * scale_weight * coh,
      row.names = NULL
    )
  )
}

passage_select_joint_factor_tmb <- function(Y,
                                            coords,
                                            X = NULL,
                                            K_grid = 1:4,
                                            criterion = c("BIC", "AIC"),
                                            ...) {
  criterion <- match.arg(criterion)
  Y <- passage_check_y(Y)
  K_grid <- sort(unique(as.integer(K_grid)))
  K_grid <- K_grid[K_grid >= 1L & K_grid <= ncol(Y)]
  if (!length(K_grid)) stop("K_grid has no valid values for this Y")
  fits <- vector("list", length(K_grid))
  rows <- vector("list", length(K_grid))
  for (ii in seq_along(K_grid)) {
    kk <- K_grid[[ii]]
    fit <- tryCatch(
      passage_fit_joint_factor_tmb(Y = Y, coords = coords, X = X, K_fit = kk, ...),
      error = function(e) e
    )
    fits[[ii]] <- fit
    rows[[ii]] <- if (inherits(fit, "error")) {
      data.frame(K = kk, converged = FALSE, logLik = NA_real_, npar = NA_integer_,
                 AIC = NA_real_, BIC = NA_real_, objective = NA_real_,
                 error_message = conditionMessage(fit))
    } else {
      data.frame(K = kk, converged = fit$opt$convergence == 0L, logLik = fit$logLik,
                 npar = fit$npar, AIC = fit$AIC, BIC = fit$BIC,
                 objective = fit$opt$objective, error_message = NA_character_)
    }
  }
  names(fits) <- paste0("K", K_grid)
  summary <- do.call(rbind, rows)
  ok <- is.finite(summary[[criterion]])
  if (!any(ok)) stop("No joint factor fits completed successfully")
  best_row <- which(ok)[which.min(summary[[criterion]][ok])]
  out <- list(
    best = fits[[best_row]],
    fits = fits,
    summary = summary,
    criterion = criterion,
    best_K = summary$K[[best_row]]
  )
  class(out) <- c("passage_joint_factor_selection", "list")
  out
}

passage_fit_joint_factor_pathway <- function(Y,
                                             coords,
                                             pathway,
                                             X = NULL,
                                             gene_names = NULL,
                                             min_pathway_size = 2L,
                                             K_fit = NULL,
                                             K_grid = NULL,
                                             criterion = c("BIC", "AIC"),
                                             ...) {
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  idx <- passage_resolve_pathway(pathway, gene_names)
  if (length(idx) < min_pathway_size) {
    stop("pathway has fewer than min_pathway_size genes present in Y")
  }
  Yp <- Y[, idx, drop = FALSE]
  if (!is.null(K_grid)) {
    sel <- passage_select_joint_factor_tmb(
      Y = Yp, coords = coords, X = X, K_grid = K_grid,
      criterion = match.arg(criterion), ...
    )
    fit <- sel$best
    selection <- sel$summary
  } else {
    if (is.null(K_fit)) K_fit <- min(2L, ncol(Yp))
    fit <- passage_fit_joint_factor_tmb(Y = Yp, coords = coords, X = X, K_fit = K_fit, ...)
    selection <- NULL
  }
  metrics <- passage_joint_factor_pve(fit)
  list(
    fit = fit,
    metrics = metrics,
    selection = selection,
    genes = colnames(Yp),
    pathway_size = ncol(Yp)
  )
}

passage_fit_joint_factor_pathways <- function(Y,
                                              coords,
                                              pathways,
                                              X = NULL,
                                              min_pathway_size = 2L,
                                              max_pathway_size = 500L,
                                              fdr_method = "BH",
                                              ...) {
  Y <- passage_check_y(Y)
  pathways <- passage_check_pathways(pathways, colnames(Y))
  fits <- vector("list", length(pathways))
  rows <- vector("list", length(pathways))
  names(fits) <- names(pathways)
  for (ii in seq_along(pathways)) {
    pname <- names(pathways)[[ii]]
    genes <- pathways[[ii]]
    q <- length(genes)
    if (q < min_pathway_size || q > max_pathway_size) {
      rows[[ii]] <- data.frame(pathway = pname, pathway_size = q, status = "size_filter",
                               K = NA_integer_, logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                               cEPSV = NA_real_, mean_propSV = NA_real_, error_message = NA_character_)
      next
    }
    ans <- tryCatch(
      passage_fit_joint_factor_pathway(
        Y = Y, coords = coords, pathway = genes, X = X,
        min_pathway_size = min_pathway_size, ...
      ),
      error = function(e) e
    )
    fits[[ii]] <- ans
    if (inherits(ans, "error")) {
      rows[[ii]] <- data.frame(pathway = pname, pathway_size = q, status = "error",
                               K = NA_integer_, logLik = NA_real_, AIC = NA_real_, BIC = NA_real_,
                               cEPSV = NA_real_, mean_propSV = NA_real_,
                               error_message = conditionMessage(ans))
    } else {
      sm <- ans$metrics$summary
      rows[[ii]] <- data.frame(pathway = pname, pathway_size = q, status = "tested",
                               K = ans$fit$K, logLik = ans$fit$logLik,
                               AIC = ans$fit$AIC, BIC = ans$fit$BIC,
                               cEPSV = sm[["cEPSV"]],
                               mean_propSV = sm[["mean_propSV"]],
                               error_message = NA_character_)
    }
  }
  summary <- do.call(rbind, rows)
  summary$cEPSV_rank <- rank(-summary$cEPSV, na.last = "keep", ties.method = "average")
  list(summary = summary, fits = fits, fdr_method = fdr_method)
}

passage_compare_joint_factor_nested <- function(fits) {
  if (!is.list(fits) || length(fits) < 2L) stop("fits must contain at least two fitted models")
  nm <- names(fits)
  if (is.null(nm) || any(!nzchar(nm))) nm <- paste0("M", seq_along(fits))
  rows <- vector("list", length(fits))
  for (ii in seq_along(fits)) {
    fit <- fits[[ii]]
    if (!inherits(fit, "passage_joint_factor_tmb")) stop("all fits must be passage_joint_factor_tmb objects")
    rows[[ii]] <- data.frame(
      model = nm[[ii]],
      K = fit$K,
      q = ncol(fit$X),
      npar = fit$npar,
      logLik = fit$logLik,
      AIC = fit$AIC,
      BIC = fit$BIC,
      cEPSV = passage_joint_factor_pve(fit)$summary[["cEPSV"]],
      mean_propSV = passage_joint_factor_pve(fit)$summary[["mean_propSV"]],
      lr_vs_prev = NA_real_,
      df_vs_prev = NA_integer_,
      p_lrt_vs_prev = NA_real_,
      comparison_note = NA_character_
    )
    if (ii > 1L) {
      prev <- fits[[ii - 1L]]
      lr <- 2 * (fit$logLik - prev$logLik)
      df <- fit$npar - prev$npar
      rows[[ii]]$lr_vs_prev <- lr
      rows[[ii]]$df_vs_prev <- df
      if (fit$K == prev$K && df > 0 && is.finite(lr)) {
        rows[[ii]]$p_lrt_vs_prev <- stats::pchisq(lr, df = df, lower.tail = FALSE)
      } else {
        rows[[ii]]$comparison_note <- "LRT omitted because K changed or df <= 0"
      }
    }
  }
  do.call(rbind, rows)
}

passage_fit_joint_factor_nested <- function(Y,
                                            coords,
                                            X_list,
                                            K_fit = NULL,
                                            K_grid = NULL,
                                            criterion = c("BIC", "AIC"),
                                            ...) {
  if (!is.list(X_list) || length(X_list) == 0L) stop("X_list must be a non-empty list")
  if (is.null(names(X_list)) || any(!nzchar(names(X_list)))) names(X_list) <- paste0("M", seq_along(X_list))
  fits <- vector("list", length(X_list))
  selections <- vector("list", length(X_list))
  names(fits) <- names(selections) <- names(X_list)
  for (nm in names(X_list)) {
    if (!is.null(K_grid)) {
      sel <- passage_select_joint_factor_tmb(
        Y = Y, coords = coords, X = X_list[[nm]], K_grid = K_grid,
        criterion = match.arg(criterion), ...
      )
      fits[[nm]] <- sel$best
      selections[[nm]] <- sel$summary
    } else {
      kk <- if (is.null(K_fit)) min(2L, ncol(as.matrix(Y))) else K_fit
      fits[[nm]] <- passage_fit_joint_factor_tmb(Y = Y, coords = coords, X = X_list[[nm]], K_fit = kk, ...)
      selections[[nm]] <- NULL
    }
  }
  list(fits = fits, selections = selections, comparison = passage_compare_joint_factor_nested(fits))
}

passage_fit_joint_factor_hypotheses <- function(Y,
                                                coords,
                                                pathway = NULL,
                                                X_base = NULL,
                                                Z_cell = NULL,
                                                Z_background = NULL,
                                                gene_names = NULL,
                                                ...) {
  Y <- passage_check_y(Y)
  if (!is.null(pathway)) {
    if (is.null(gene_names)) gene_names <- colnames(Y)
    idx <- passage_resolve_pathway(pathway, gene_names)
    if (!length(idx)) stop("pathway has no genes present in Y")
    Y <- Y[, idx, drop = FALSE]
  }
  X_list <- list(H1 = X_base)
  if (!is.null(Z_cell)) {
    X_list$H2 <- if (is.null(X_base)) Z_cell else cbind(as.matrix(X_base), as.matrix(Z_cell))
  }
  if (!is.null(Z_background)) {
    xb <- X_base
    if (!is.null(Z_cell)) xb <- if (is.null(xb)) Z_cell else cbind(as.matrix(xb), as.matrix(Z_cell))
    X_list$H3 <- if (is.null(xb)) Z_background else cbind(as.matrix(xb), as.matrix(Z_background))
  }
  passage_fit_joint_factor_nested(Y = Y, coords = coords, X_list = X_list, ...)
}

passage_simulate_joint_factor_tmb_data <- function(coords,
                                                   A,
                                                   phi,
                                                   tau2,
                                                   B = NULL,
                                                   X = NULL,
                                                   kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                                   seed = NULL) {
  kernel <- match.arg(kernel)
  if (!is.null(seed)) set.seed(seed)
  coords <- passage_check_coords(coords)
  A <- as.matrix(A)
  n <- nrow(coords)
  p <- nrow(A)
  K <- ncol(A)
  if (length(phi) != K) stop("length(phi) must equal ncol(A)")
  if (length(tau2) != p) stop("length(tau2) must equal nrow(A)")
  if (is.null(X)) X <- matrix(1, nrow = n, ncol = 1L)
  X <- passage_prepare_design(X, n, intercept = TRUE)
  if (is.null(B)) B <- matrix(0, nrow = ncol(X), ncol = p)
  B <- as.matrix(B)
  if (!all(dim(B) == c(ncol(X), p))) stop("B must be ncol(X) by nrow(A)")
  D <- as.matrix(stats::dist(coords))
  U <- matrix(NA_real_, n, K)
  for (k in seq_len(K)) {
    C <- passage_kernel_corr(D, range = phi[[k]], kernel = kernel)
    L <- chol(passage_joint_near_psd(C) + diag(1e-8, n))
    U[, k] <- as.numeric(t(L) %*% stats::rnorm(n))
  }
  E <- matrix(stats::rnorm(n * p, sd = rep(sqrt(tau2), each = n)), nrow = n, ncol = p)
  Y <- X %*% B + U %*% t(A) + E
  colnames(Y) <- rownames(A) %||% paste0("gene_", seq_len(p))
  list(Y = Y, U = U, A = A, phi = phi, tau2 = tau2, B = B, X = X, coords = coords)
}

`%||%` <- function(x, y) if (is.null(x)) y else x
