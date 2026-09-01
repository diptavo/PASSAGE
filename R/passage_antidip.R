# Anti-double-dipping PASSAGE helpers.

passage_make_spot_folds <- function(n, n_folds = 5L, seed = NULL) {
  n <- as.integer(n)
  n_folds <- as.integer(n_folds)
  if (n < 4L) stop("n must be at least 4 for spot cross-fitting")
  n_folds <- max(2L, min(n_folds, n))
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  idx <- sample.int(n)
  split(idx, rep(seq_len(n_folds), length.out = n))
}

passage_align_loadings <- function(A, ref_A) {
  A <- as.matrix(A)
  ref_A <- as.matrix(ref_A)
  K <- min(ncol(A), ncol(ref_A))
  if (K < 1L) stop("A and ref_A must contain at least one factor")
  C <- abs(stats::cor(A[, seq_len(K), drop = FALSE], ref_A[, seq_len(K), drop = FALSE]))
  C[!is.finite(C)] <- 0
  perm <- integer(K)
  used <- logical(K)
  for (kk in seq_len(K)) {
    vals <- C[, kk]
    vals[used] <- -Inf
    jj <- which.max(vals)
    if (!is.finite(vals[jj])) jj <- which(!used)[1L]
    perm[kk] <- jj
    used[jj] <- TRUE
  }
  A2 <- A[, perm, drop = FALSE]
  sign_flip <- vapply(seq_len(K), function(kk) {
    s <- sum(A2[, kk] * ref_A[, kk])
    if (is.finite(s) && s < 0) -1 else 1
  }, numeric(1))
  sweep(A2, 2L, sign_flip, "*")
}

passage_project_scores <- function(R, A, ridge = 1e-8) {
  R <- as.matrix(R)
  A <- as.matrix(A)
  G <- crossprod(A)
  sc <- mean(diag(G))
  if (!is.finite(sc) || sc <= 0) sc <- 1
  as.matrix(R %*% A %*% solve(G + diag(ridge * sc, ncol(G))))
}

passage_normalize_factor_pair <- function(V, A, min_sd = 1e-8) {
  V <- as.matrix(V)
  A <- as.matrix(A)
  mu <- colMeans(V)
  V <- sweep(V, 2L, mu, "-")
  ss <- apply(V, 2L, stats::sd)
  keep <- is.finite(ss) & ss > min_sd
  if (!any(keep)) stop("cross-fitted factors are numerically degenerate")
  V <- V[, keep, drop = FALSE]
  A <- A[, keep, drop = FALSE]
  ss <- ss[keep]
  V <- sweep(V, 2L, ss, "/")
  A <- sweep(A, 2L, ss, "*")
  colnames(V) <- colnames(A) <- paste0("factor_", seq_len(ncol(V)))
  list(V = V, A = A)
}

passage_fit_engine_crossfit_spots <- function(Y,
                                              coords,
                                              X = NULL,
                                              n_folds = 5L,
                                              seed = NULL,
                                              method = c("pca", "spatial_basis", "smoothed_pca", "nmf", "alternating_gp"),
                                              fit_args = list(),
                                              ridge = 1e-8,
                                              verbose = TRUE) {
  method <- match.arg(method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  folds <- passage_make_spot_folds(nrow(Y), n_folds = n_folds, seed = seed)
  R_full <- passage_residualize(Y, X)$resid

  engines <- vector("list", length(folds))
  if (verbose) message("PASSAGE anti-dip: fitting ", length(folds), " spot-crossfit engines")
  for (ff in seq_along(folds)) {
    test_idx <- folds[[ff]]
    train_idx <- setdiff(seq_len(nrow(Y)), test_idx)
    args <- c(
      list(
        Y = Y[train_idx, , drop = FALSE],
        coords = coords[train_idx, , drop = FALSE],
        X = X[train_idx, , drop = FALSE],
        method = method
      ),
      fit_args
    )
    args$verbose <- isTRUE(fit_args$verbose)
    engines[[ff]] <- do.call(passage_fit_factor_engine, args)
  }

  K <- min(vapply(engines, function(e) e$K, integer(1)))
  if (K < 1L) stop("no crossfit engine returned a positive factor rank")
  ref_A <- engines[[1L]]$A[, seq_len(K), drop = FALSE]
  V_cf <- matrix(NA_real_, nrow(Y), K)
  A_sum <- matrix(0, ncol(Y), K, dimnames = list(colnames(Y), paste0("factor_", seq_len(K))))
  for (ff in seq_along(folds)) {
    eng <- engines[[ff]]
    test_idx <- folds[[ff]]
    A_ff <- passage_align_loadings(eng$A[, seq_len(K), drop = FALSE], ref_A)
    R_test <- R_full[test_idx, , drop = FALSE]
    if (!is.null(eng$gene_center)) R_test <- sweep(R_test, 2L, eng$gene_center, "-")
    if (!is.null(eng$gene_scale)) R_test <- sweep(R_test, 2L, pmax(eng$gene_scale, sqrt(.Machine$double.eps)), "/")
    V_cf[test_idx, ] <- passage_project_scores(R_test, A_ff, ridge = ridge)
    A_sum <- A_sum + A_ff
  }
  A_cf <- A_sum / length(folds)
  pair <- passage_normalize_factor_pair(V_cf, A_cf)
  fitted_spatial <- pair$V %*% t(pair$A)
  D <- colMeans((R_full - fitted_spatial)^2)
  D <- pmax(D, .Machine$double.eps)

  base <- engines[[1L]]
  out <- passage_engine_complete(
    Y = Y,
    coords = coords,
    X = X,
    residuals = R_full,
    A = pair$A,
    V = pair$V,
    D = D,
    range_grid = base$range_grid,
    m = base$m,
    kernel = base$kernel,
    ordering = if (!is.null(base$ordering)) base$ordering else "coordinate",
    engine_name = paste0("spot_crossfit_", method, "_vecchia_v1"),
    rank_info = list(K = ncol(pair$V), n_folds = length(folds), fold_sizes = lengths(folds)),
    fitted_spatial = fitted_spatial,
    extra = list(anti_dip = "spot_crossfit", fold_engines = engines),
    verbose = verbose
  )
  out$anti_dip <- "spot_crossfit"
  out
}

passage_fit_engine_pathway_holdout <- function(Y,
                                               coords,
                                               pathway,
                                               X = NULL,
                                               gene_names = NULL,
                                               method = c("pca", "spatial_basis", "smoothed_pca", "nmf", "alternating_gp"),
                                               fit_args = list(),
                                               min_background_genes = 5L,
                                               verbose = TRUE) {
  method <- match.arg(method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  bg <- setdiff(seq_len(ncol(Y)), P)
  if (length(bg) < min_background_genes) {
    stop("not enough non-pathway genes to fit a pathway-holdout engine")
  }
  if (verbose) {
    message("PASSAGE anti-dip: fitting pathway-holdout engine with ", length(bg), " background genes")
  }
  args <- c(
    list(Y = Y[, bg, drop = FALSE], coords = coords, X = X, method = method),
    fit_args
  )
  args$verbose <- isTRUE(fit_args$verbose)
  bg_engine <- do.call(passage_fit_factor_engine, args)
  R_full <- passage_residualize(Y, X)$resid
  fit <- passage_refit_loadings(R_full, bg_engine$V)
  out <- bg_engine
  out$A <- fit$A
  out$V <- fit$V
  out$D <- stats::setNames(fit$D, colnames(Y))
  out$residuals <- R_full
  out$fitted_spatial <- fit$fitted_spatial
  out$gene_names <- colnames(Y)
  out$G <- ncol(Y)
  out$anti_dip <- "pathway_holdout"
  out$holdout_genes <- gene_names[P]
  out$background_genes <- gene_names[bg]
  out$engine <- paste0("pathway_holdout_", out$engine)
  class(out) <- c("passage_engine", "list")
  out
}

passage_score_test_refit_null <- function(Y,
                                          coords,
                                          pathway,
                                          X = NULL,
                                          gene_names = NULL,
                                          method = c("pca", "spatial_basis", "smoothed_pca", "nmf", "alternating_gp"),
                                          fit_args = list(),
                                          weight_schemes = c("equal", "var", "range"),
                                          run_burden = TRUE,
                                          run_spasset = TRUE,
                                          spasset_grid = c(2, 5, 10, 25, 50, 100, 250, 500),
                                          n_perm = 99L,
                                          seed = NULL,
                                          return_engine = FALSE,
                                          verbose = TRUE) {
  method <- match.arg(method)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) stop("pathway has no genes present in Y")
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("n_perm must be at least 1")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  fit_engine <- function(Y_fit) {
    args <- c(list(Y = Y_fit, coords = coords, X = X, method = method), fit_args)
    args$verbose <- isTRUE(fit_args$verbose)
    do.call(passage_fit_factor_engine, args)
  }
  score_engine <- function(engine) {
    pre <- passage_h_precompute(engine, X = X)
    passage_score_test(
      engine = engine,
      Y = Y,
      pathway = P,
      precomp = pre,
      weight_schemes = weight_schemes,
      run_burden = run_burden,
      run_spasset = run_spasset,
      spasset_grid = spasset_grid,
      gene_names = gene_names,
      calibration = "moment"
    )
  }

  if (verbose) message("PASSAGE anti-dip: fitting observed refit-null engine")
  obs_engine <- fit_engine(Y)
  obs <- score_engine(obs_engine)
  obs_stat <- passage_calibration_stat(obs$p_omnibus)
  fit0 <- passage_residualize(Y, X)
  Y_hat <- Y - fit0$resid
  perm_stat <- numeric(n_perm)
  perm_p <- numeric(n_perm)
  for (bb in seq_len(n_perm)) {
    if (verbose && (bb == 1L || bb %% 10L == 0L || bb == n_perm)) {
      message("PASSAGE anti-dip: refit-null permutation ", bb, "/", n_perm)
    }
    Yb <- Y
    ord <- sample.int(nrow(Y))
    Yb[, P] <- Y_hat[, P, drop = FALSE] + fit0$resid[ord, P, drop = FALSE]
    eb <- fit_engine(Yb)
    hb <- score_engine(eb)
    perm_p[bb] <- hb$p_omnibus
    perm_stat[bb] <- passage_calibration_stat(hb$p_omnibus)
  }
  obs$p_omnibus_moment <- obs$p_omnibus
  obs$p_omnibus <- (1 + sum(perm_stat >= obs_stat)) / (1 + n_perm)
  obs$p_omnibus_empirical <- obs$p_omnibus
  obs$calibration <- list(
    method = "refit_null_residual_permutation",
    n_perm = n_perm,
    observed_stat = obs_stat,
    null_stat = perm_stat,
    null_p_moment = perm_p
  )
  obs$anti_dip <- "refit_null"
  if (return_engine) obs$engine <- obs_engine
  class(obs) <- c("passage_score_result", "list")
  obs
}
