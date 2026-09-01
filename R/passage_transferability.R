# Cross-fitted pathway spatial transferability metrics.

passage_make_index_folds <- function(n, n_folds = 5L, seed = NULL) {
  n <- as.integer(n)
  n_folds <- max(2L, min(as.integer(n_folds), n))
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

passage_pst_gene_split <- function(P, seed = NULL) {
  P <- unique(as.integer(P))
  if (length(P) < 4L) stop("PST requires at least four pathway genes")
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  P <- sample(P)
  n_a <- floor(length(P) / 2L)
  list(A = P[seq_len(n_a)], B = P[(n_a + 1L):length(P)])
}

passage_pst_component <- function(R,
                                  train_idx,
                                  eval_idx,
                                  genes,
                                  component = c("pc1", "mean"),
                                  min_sd = 1e-8) {
  component <- match.arg(component)
  genes <- as.integer(genes)
  Z_train <- R[train_idx, genes, drop = FALSE]
  Z_eval <- R[eval_idx, genes, drop = FALSE]
  mu <- colMeans(Z_train)
  ss <- apply(Z_train, 2L, stats::sd)
  keep <- is.finite(ss) & ss > min_sd
  if (sum(keep) < 2L) return(NULL)
  Z_train <- sweep(Z_train[, keep, drop = FALSE], 2L, mu[keep], "-")
  Z_train <- sweep(Z_train, 2L, ss[keep], "/")
  Z_eval <- sweep(Z_eval[, keep, drop = FALSE], 2L, mu[keep], "-")
  Z_eval <- sweep(Z_eval, 2L, ss[keep], "/")
  if (component == "mean") {
    score_train <- rowMeans(Z_train)
    score_eval <- rowMeans(Z_eval)
  } else {
    sv <- tryCatch(svd(Z_train, nu = 0L, nv = 1L), error = function(e) NULL)
    if (is.null(sv) || length(sv$d) < 1L || !is.finite(sv$d[1L]) || sv$d[1L] <= 0) {
      return(NULL)
    }
    w <- sv$v[, 1L]
    score_train <- as.numeric(Z_train %*% w)
    score_eval <- as.numeric(Z_eval %*% w)
  }
  s0 <- stats::sd(score_train)
  if (!is.finite(s0) || s0 <= min_sd) return(NULL)
  m0 <- mean(score_train)
  list(
    train = (score_train - m0) / s0,
    eval = (score_eval - m0) / s0,
    n_genes = sum(keep)
  )
}

passage_pst_neighbor_lag <- function(coords,
                                     train_idx,
                                     test_idx,
                                     train_score,
                                     k_neighbors = 8L,
                                     weight = c("uniform", "inverse_distance"),
                                     eps = 1e-8) {
  weight <- match.arg(weight)
  coords <- passage_check_coords(coords)
  train_idx <- as.integer(train_idx)
  test_idx <- as.integer(test_idx)
  k_neighbors <- max(1L, min(as.integer(k_neighbors), length(train_idx)))
  out <- numeric(length(test_idx))
  for (ii in seq_along(test_idx)) {
    d <- sqrt(rowSums((coords[train_idx, , drop = FALSE] -
      matrix(coords[test_idx[ii], ], nrow = length(train_idx),
             ncol = ncol(coords), byrow = TRUE))^2))
    nb <- order(d)[seq_len(k_neighbors)]
    if (weight == "uniform") {
      ww <- rep(1 / length(nb), length(nb))
    } else {
      ww <- 1 / pmax(d[nb], eps)
      ww <- ww / sum(ww)
    }
    out[ii] <- sum(ww * train_score[nb])
  }
  out
}

passage_pst_r2 <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 4L || stats::sd(x) <= 1e-8 || stats::sd(y) <= 1e-8) {
    return(NA_real_)
  }
  r <- suppressWarnings(stats::cor(x, y))
  if (!is.finite(r)) NA_real_ else r^2
}

passage_pathway_spatial_transferability <- function(Y,
                                                    coords,
                                                    pathway,
                                                    X = NULL,
                                                    gene_names = NULL,
                                                    n_gene_splits = 5L,
                                                    n_folds = 5L,
                                                    k_neighbors = 8L,
                                                    component = c("pc1", "mean"),
                                                    neighbor_weight = c("uniform", "inverse_distance"),
                                                    n_perm = 0L,
                                                    seed = NULL) {
  component <- match.arg(component)
  neighbor_weight <- match.arg(neighbor_weight)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  if (is.null(gene_names)) gene_names <- colnames(Y)
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) < 4L) stop("pathway must contain at least four genes present in Y")
  Xd <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  R <- passage_residualize_with_qr(Y, qr(Xd))
  folds <- passage_make_index_folds(nrow(Y), n_folds = n_folds, seed = seed)
  n_gene_splits <- as.integer(n_gene_splits)
  n_perm <- as.integer(n_perm)
  if (n_gene_splits < 1L) stop("n_gene_splits must be at least 1")
  if (n_perm < 0L) stop("n_perm must be nonnegative")
  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }

  pairs <- list()
  split_rows <- list()
  pp <- 0L
  rr <- 0L
  for (ss in seq_len(n_gene_splits)) {
    spl <- passage_pst_gene_split(P)
    for (ff in seq_along(folds)) {
      test_idx <- folds[[ff]]
      train_idx <- setdiff(seq_len(nrow(Y)), test_idx)
      comp_a <- passage_pst_component(R, train_idx, test_idx, spl$A, component = component)
      comp_b <- passage_pst_component(R, train_idx, test_idx, spl$B, component = component)
      if (is.null(comp_a) || is.null(comp_b)) next
      lag_a <- passage_pst_neighbor_lag(
        coords, train_idx, test_idx, comp_a$train,
        k_neighbors = k_neighbors, weight = neighbor_weight
      )
      lag_b <- passage_pst_neighbor_lag(
        coords, train_idx, test_idx, comp_b$train,
        k_neighbors = k_neighbors, weight = neighbor_weight
      )
      r2_ab <- passage_pst_r2(lag_a, comp_b$eval)
      r2_ba <- passage_pst_r2(lag_b, comp_a$eval)
      if (is.finite(r2_ab)) {
        pp <- pp + 1L
        pairs[[pp]] <- list(x = lag_a, y = comp_b$eval, direction = "A_to_B")
      }
      if (is.finite(r2_ba)) {
        pp <- pp + 1L
        pairs[[pp]] <- list(x = lag_b, y = comp_a$eval, direction = "B_to_A")
      }
      rr <- rr + 1L
      split_rows[[rr]] <- data.frame(
        gene_split = ss,
        fold = ff,
        n_genes_A = comp_a$n_genes,
        n_genes_B = comp_b$n_genes,
        r2_A_to_B = r2_ab,
        r2_B_to_A = r2_ba,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(pairs) == 0L) stop("no valid PST folds were available")
  pair_r2 <- vapply(pairs, function(z) passage_pst_r2(z$x, z$y), numeric(1))
  stat <- mean(pair_r2, na.rm = TRUE)
  null <- numeric(0)
  p_value <- NA_real_
  if (n_perm > 0L) {
    null <- numeric(n_perm)
    for (bb in seq_len(n_perm)) {
      null[bb] <- mean(vapply(pairs, function(z) {
        passage_pst_r2(z$x, sample(z$y))
      }, numeric(1)), na.rm = TRUE)
    }
    p_value <- (1 + sum(null >= stat)) / (1 + length(null))
  }
  out <- list(
    statistic = stat,
    p = p_value,
    pair_r2 = pair_r2,
    null = null,
    split_table = do.call(rbind, split_rows),
    pathway_size = length(P),
    n_pairs = length(pairs),
    n_gene_splits = n_gene_splits,
    n_folds = length(folds),
    k_neighbors = k_neighbors,
    component = component,
    neighbor_weight = neighbor_weight,
    method = "pathway_spatial_transferability"
  )
  class(out) <- c("passage_pst_result", "list")
  out
}
