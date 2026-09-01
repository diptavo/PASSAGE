# PASSAGE Layer 2: pathway score tests.

passage_h_precompute <- function(engine,
                                 X = NULL,
                                 score_kernel = c("precision"),
                                 verbose = FALSE) {
  if (!inherits(engine, "passage_engine")) {
    stop("engine must be from passage_fit_engine_pca()")
  }
  score_kernel <- match.arg(score_kernel)
  X <- passage_prepare_design(X, engine$N, intercept = TRUE)
  moments <- vector("list", engine$K)
  for (k in seq_len(engine$K)) {
    if (verbose) {
      message("precomputing moments for factor ", k, "/", engine$K)
    }
    moments[[k]] <- passage_kernel_moments(engine$K_score[[k]], X[engine$vecchia[[k]]$ord, , drop = FALSE])
  }
  list(
    X = X,
    qr = qr(X),
    moments = moments,
    score_kernel = score_kernel,
    K = engine$K,
    N = engine$N
  )
}

passage_score_test <- function(engine,
                               Y,
                               pathway,
                               precomp,
                               weight_schemes = c("equal", "var", "range", "invvar"),
                               run_burden = TRUE,
                               run_spasset = TRUE,
                               spasset_grid = c(2, 5, 10, 25, 50, 100, 250, 500),
                               gene_names = NULL,
                               calibration = c("permutation", "moment"),
                               n_perm = 199L,
                               seed = NULL) {
  calibration <- match.arg(calibration)
  if (!inherits(engine, "passage_engine")) {
    stop("engine must be from passage_fit_engine_pca()")
  }
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) {
    gene_names <- colnames(Y)
  }
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) < 1L) {
    stop("pathway has no genes present in Y")
  }
  YP <- Y[, P, drop = FALSE]
  A_P <- engine$A[P, , drop = FALSE]
  D_P <- engine$D[P]
  R_P <- passage_residualize_with_qr(YP, precomp$qr)
  out <- passage_score_from_residuals(
    engine = engine,
    R_P = R_P,
    P = P,
    A_P = A_P,
    D_P = D_P,
    precomp = precomp,
    weight_schemes = weight_schemes,
    run_burden = run_burden,
    run_spasset = run_spasset,
    spasset_grid = spasset_grid,
    gene_names = gene_names
  )

  out$p_omnibus_moment <- out$p_omnibus
  if (calibration == "permutation") {
    n_perm <- as.integer(n_perm)
    if (n_perm < 1L) {
      stop("n_perm must be positive when calibration = 'permutation'")
    }
    if (!is.null(seed)) {
      old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
        get(".Random.seed", envir = .GlobalEnv)
      } else {
        NULL
      }
      on.exit({
        if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      }, add = TRUE)
      set.seed(seed)
    }
    obs_stat <- passage_calibration_stat(out$p_omnibus_moment)
    perm_stat <- numeric(n_perm)
    n <- nrow(R_P)
    for (bb in seq_len(n_perm)) {
      R_perm <- R_P[sample.int(n), , drop = FALSE]
      perm_out <- passage_score_from_residuals(
        engine = engine,
        R_P = R_perm,
        P = P,
        A_P = A_P,
        D_P = D_P,
        precomp = precomp,
        weight_schemes = weight_schemes,
        run_burden = run_burden,
        run_spasset = run_spasset,
        spasset_grid = spasset_grid,
        gene_names = gene_names
      )
      perm_stat[bb] <- passage_calibration_stat(perm_out$p_omnibus)
    }
    out$p_omnibus <- (1 + sum(perm_stat >= obs_stat)) / (1 + n_perm)
    out$p_omnibus_empirical <- out$p_omnibus
    out$calibration <- list(
      method = "residual_permutation",
      n_perm = n_perm,
      observed_stat = obs_stat,
      null_stat = perm_stat
    )
  } else {
    out$calibration <- list(method = "moment")
  }
  class(out) <- c("passage_score_result", "list")
  out
}

passage_score_from_residuals <- function(engine,
                                         R_P,
                                         P,
                                         A_P,
                                         D_P,
                                         precomp,
                                         weight_schemes,
                                         run_burden,
                                         run_spasset,
                                         spasset_grid,
                                         gene_names) {
  K <- engine$K
  per_factor <- vector("list", K)
  q_vec <- c_vec <- p_vec <- mean_vec <- var_vec <- numeric(K)
  for (k in seq_len(K)) {
    vc <- engine$vecchia[[k]]
    z <- as.numeric(R_P[vc$ord, , drop = FALSE] %*% A_P[, k])
    qk <- as.numeric(crossprod(z, engine$K_score[[k]] %*% z))
    ck <- sum((A_P[, k]^2) * D_P)
    mean_k <- ck * precomp$moments[[k]]$trace_mk
    var_k <- 2 * ck^2 * precomp$moments[[k]]$trace_mkmk
    pk <- passage_satterthwaite_p(qk, mean_k, var_k)
    q_vec[k] <- qk
    c_vec[k] <- ck
    p_vec[k] <- pk
    mean_vec[k] <- mean_k
    var_vec[k] <- var_k
    per_factor[[k]] <- list(Q = qk, c = ck, p = pk, mean = mean_k, var = var_k)
  }

  joint <- list()
  for (ws in weight_schemes) {
    w <- passage_factor_weights(engine, ws)
    q_joint <- sum(w * q_vec)
    mean_joint <- sum(w * mean_vec)
    var_joint <- sum((w^2) * var_vec)
    joint[[ws]] <- list(
      Q = q_joint,
      p = passage_satterthwaite_p(q_joint, mean_joint, var_joint),
      p_acat = passage_acat(p_vec, weights = w),
      mean = mean_joint,
      var = var_joint,
      weights = w
    )
  }

  burden <- NULL
  if (run_burden) {
    a <- rep(1 / sqrt(length(P)), length(P))
    vc <- engine$vecchia[[1L]]
    z <- as.numeric(R_P[vc$ord, , drop = FALSE] %*% a)
    q_b <- as.numeric(crossprod(z, engine$K_score[[1L]] %*% z))
    c_b <- mean(D_P)
    mean_b <- c_b * precomp$moments[[1L]]$trace_mk
    var_b <- 2 * c_b^2 * precomp$moments[[1L]]$trace_mkmk
    burden <- list(Q = q_b, c = c_b, p = passage_satterthwaite_p(q_b, mean_b, var_b))
  }

  spasset <- NULL
  if (run_spasset && length(P) >= 2L) {
    vc <- engine$vecchia[[1L]]
    R_ord <- R_P[vc$ord, , drop = FALSE]
    gene_q <- colSums(R_ord * as.matrix(engine$K_score[[1L]] %*% R_ord)) /
      pmax(D_P, .Machine$double.eps)
    gene_order <- order(gene_q, decreasing = TRUE)
    grid <- sort(unique(pmin(length(P), spasset_grid)))
    grid <- grid[grid >= 2L]
    p_grid <- numeric(length(grid))
    q_grid <- numeric(length(grid))
    for (ii in seq_along(grid)) {
      keep <- gene_order[seq_len(grid[ii])]
      sub_out <- passage_score_from_residuals(
        engine = engine,
        R_P = R_P[, keep, drop = FALSE],
        P = P[keep],
        A_P = A_P[keep, , drop = FALSE],
        D_P = D_P[keep],
        precomp = precomp,
        weight_schemes = "equal",
        run_burden = FALSE,
        run_spasset = FALSE,
        spasset_grid = spasset_grid,
        gene_names = gene_names
      )
      p_grid[ii] <- sub_out$joint$equal$p_acat
      q_grid[ii] <- sub_out$joint$equal$Q
    }
    best <- which.min(p_grid)
    spasset <- list(
      p = passage_acat(p_grid),
      best_p = p_grid[best],
      best_size = grid[best],
      best_genes = gene_names[P[gene_order[seq_len(grid[best])]]],
      p_by_size = stats::setNames(p_grid, paste0("top", grid)),
      Q_by_size = stats::setNames(q_grid, paste0("top", grid)),
      gene_scores = stats::setNames(gene_q, gene_names[P])
    )
  }

  component_p <- vapply(joint, function(x) x$p_acat, numeric(1))
  if (!is.null(burden)) {
    component_p <- c(component_p, burden = burden$p)
  }
  if (!is.null(spasset)) {
    component_p <- c(component_p, spasset = spasset$p)
  }
  list(
    p_omnibus = passage_acat(component_p),
    p_components = component_p,
    per_factor = per_factor,
    joint = joint,
    burden = burden,
    spasset = spasset,
    pathway_size = length(P),
    pathway_genes = gene_names[P],
    score_kernel = precomp$score_kernel
  )
}

passage_calibration_stat <- function(p) {
  p <- pmin(pmax(as.numeric(p), 1e-300), 1)
  -log(p)
}

passage_factor_weights <- function(engine, scheme = c("equal", "var", "range", "invvar")) {
  scheme <- match.arg(scheme)
  w <- switch(
    scheme,
    equal = rep(1, engine$K),
    var = engine$theta$sigma2,
    range = engine$theta$effective_range,
    invvar = 1 / pmax(engine$theta$sigma2, .Machine$double.eps)
  )
  w <- as.numeric(w)
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) {
    w <- rep(1, engine$K)
  }
  w
}

passage_test_pathway <- function(engine,
                                 Y,
                                 pathway,
                                 Z_CT = NULL,
                                 V_BG = NULL,
                                 hypotheses = c("H1", "H2"),
                                 weight_schemes = c("equal", "var", "range"),
                                 gene_names = NULL,
                                 calibration = c("permutation", "moment"),
                                 n_perm = 199L,
                                 seed = NULL) {
  calibration <- match.arg(calibration)
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) {
    gene_names <- colnames(Y)
  }
  results <- list()
  if ("H1" %in% hypotheses) {
    pre_h1 <- passage_h_precompute(engine, X = matrix(1, nrow(Y), 1L))
    results$H1 <- passage_score_test(
      engine, Y, pathway, pre_h1, weight_schemes,
      gene_names = gene_names, calibration = calibration,
      n_perm = n_perm, seed = seed
    )
  }
  if ("H2" %in% hypotheses) {
    if (is.null(Z_CT)) {
      warning("H2 requested but Z_CT is NULL; skipping H2")
    } else {
      X2 <- passage_prepare_design(Z_CT, nrow(Y), intercept = TRUE)
      pre_h2 <- passage_h_precompute(engine, X = X2)
      results$H2 <- passage_score_test(
        engine, Y, pathway, pre_h2, weight_schemes,
        gene_names = gene_names, calibration = calibration,
        n_perm = n_perm, seed = if (is.null(seed)) NULL else seed + 1L
      )
    }
  }
  if ("H3" %in% hypotheses) {
    if (is.null(V_BG)) {
      warning("H3 requested but V_BG is NULL; skipping H3")
    } else {
      X3 <- if (is.null(Z_CT)) {
        passage_prepare_design(V_BG, nrow(Y), intercept = TRUE)
      } else {
        passage_prepare_design(cbind(Z_CT, V_BG), nrow(Y), intercept = TRUE)
      }
      pre_h3 <- passage_h_precompute(engine, X = X3)
      results$H3 <- passage_score_test(
        engine, Y, pathway, pre_h3, weight_schemes,
        gene_names = gene_names, calibration = calibration,
        n_perm = n_perm, seed = if (is.null(seed)) NULL else seed + 2L
      )
    }
  }

  pve <- passage_pve(engine, Y, pathway, gene_names = gene_names)
  out <- list(
    hypotheses = results,
    decomposition = passage_decomposition(results$H1, results$H2, results$H3),
    pve = pve,
    pathway_size = length(passage_resolve_pathway(pathway, gene_names))
  )
  class(out) <- c("passage_pathway_result", "list")
  out
}

passage_decomposition <- function(h1 = NULL, h2 = NULL, h3 = NULL) {
  q1 <- if (!is.null(h1)) h1$joint$equal$Q else NA_real_
  q2 <- if (!is.null(h2)) h2$joint$equal$Q else NA_real_
  q3 <- if (!is.null(h3)) h3$joint$equal$Q else NA_real_
  q1s <- max(q1, 1e-12, na.rm = TRUE)
  q2s <- max(q2, 1e-12, na.rm = TRUE)
  list(
    Q_H1 = q1,
    Q_H2 = q2,
    Q_H3 = q3,
    p_H1 = if (!is.null(h1)) h1$p_omnibus else NA_real_,
    p_H2 = if (!is.null(h2)) h2$p_omnibus else NA_real_,
    p_H3 = if (!is.null(h3)) h3$p_omnibus else NA_real_,
    cell_type_share = if (is.finite(q1) && is.finite(q2)) pmax(pmin((q1 - q2) / q1s, 1), 0) else NA_real_,
    background_share = if (is.finite(q2) && is.finite(q3)) pmax(pmin((q2 - q3) / q2s, 1), 0) else NA_real_,
    pathway_specific_share = if (is.finite(q1) && is.finite(q3)) pmax(pmin(q3 / q1s, 1), 0) else NA_real_
  )
}
