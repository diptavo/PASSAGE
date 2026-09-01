# Shared utilities for PASSAGE-style spatial pathway testing.

passage_acat <- function(p, weights = NULL) {
  p <- as.numeric(p)
  ok <- is.finite(p) & !is.na(p)
  p <- p[ok]
  if (length(p) == 0L) {
    return(NA_real_)
  }
  eps <- 1e-15
  p <- pmin(pmax(p, eps), 1 - eps)
  if (is.null(weights)) {
    weights <- rep(1 / length(p), length(p))
  } else {
    weights <- as.numeric(weights[ok])
    if (length(weights) != length(p) || any(!is.finite(weights)) || sum(weights) <= 0) {
      stop("weights must be finite, nonnegative, and aligned with p")
    }
    weights <- weights / sum(weights)
  }
  stat <- sum(weights * tan((0.5 - p) * pi))
  pmin(pmax(0.5 - atan(stat) / pi, 0), 1)
}

passage_default_range_grid <- function(coords,
                                       n_grid = 7L,
                                       min_frac = 0.03,
                                       max_frac = 0.60) {
  coords <- passage_check_coords(coords)
  bbox <- apply(coords, 2L, range)
  diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  if (!is.finite(diameter) || diameter <= 0) {
    stop("coords must span a positive spatial domain")
  }
  diameter * exp(seq(log(min_frac), log(max_frac), length.out = n_grid))
}

passage_check_y <- function(Y) {
  Y <- as.matrix(Y)
  storage.mode(Y) <- "double"
  if (!all(is.finite(Y))) {
    stop("Y contains non-finite values")
  }
  if (is.null(colnames(Y))) {
    colnames(Y) <- paste0("gene_", seq_len(ncol(Y)))
  }
  Y
}

passage_check_coords <- function(coords, n = NULL) {
  coords <- as.matrix(coords)
  storage.mode(coords) <- "double"
  if (!is.null(n) && nrow(coords) != n) {
    stop("nrow(coords) must match nrow(Y)")
  }
  if (ncol(coords) < 2L) {
    stop("coords must have at least two columns")
  }
  if (!all(is.finite(coords))) {
    stop("coords contains non-finite values")
  }
  coords
}

passage_prepare_design <- function(X, n, intercept = TRUE) {
  if (is.null(X)) {
    if (!intercept) {
      return(matrix(numeric(0), nrow = n, ncol = 0L))
    }
    X <- matrix(1, nrow = n, ncol = 1L)
    colnames(X) <- "intercept"
    return(X)
  }
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  if (nrow(X) != n) {
    stop("nrow(X) must match nrow(Y)")
  }
  if (!all(is.finite(X))) {
    stop("X contains non-finite values")
  }
  has_intercept <- any(apply(X, 2L, function(x) all(abs(x - 1) < 1e-12)))
  if (intercept && !has_intercept) {
    X <- cbind(intercept = 1, X)
  }
  X
}

passage_residualize <- function(Y, X) {
  Y <- as.matrix(Y)
  X <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  qrx <- qr(X)
  resid <- qr.resid(qrx, Y)
  df <- nrow(Y) - qrx$rank
  if (df <= 0L) {
    stop("design matrix has no residual degrees of freedom")
  }
  sigma2 <- colSums(resid^2) / df
  sigma2 <- pmax(sigma2, .Machine$double.eps)
  list(resid = resid, sigma2 = sigma2, rank = qrx$rank, df = df, qr = qrx, X = X)
}

passage_residualize_with_qr <- function(Y, qrx) {
  qr.resid(qrx, as.matrix(Y))
}

passage_resolve_pathway <- function(pathway, gene_names) {
  if (is.numeric(pathway)) {
    idx <- unique(as.integer(pathway))
    idx <- idx[idx >= 1L & idx <= length(gene_names)]
    return(idx)
  }
  pathway <- unique(as.character(pathway))
  idx <- match(pathway, gene_names)
  idx <- idx[!is.na(idx)]
  unique(idx)
}

passage_check_pathways <- function(pathways, gene_names) {
  if (!is.list(pathways) || length(pathways) == 0L) {
    stop("pathways must be a non-empty list")
  }
  if (is.null(names(pathways)) || any(names(pathways) == "")) {
    names(pathways) <- paste0("pathway_", seq_along(pathways))
  }
  lapply(pathways, function(x) gene_names[passage_resolve_pathway(x, gene_names)])
}

passage_kernel_corr <- function(d, range, kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian")) {
  kernel <- match.arg(kernel)
  x <- pmax(d / range, 0)
  if (kernel == "matern12") {
    return(exp(-x))
  }
  if (kernel == "matern32") {
    z <- sqrt(3) * x
    return((1 + z) * exp(-z))
  }
  if (kernel == "matern52") {
    z <- sqrt(5) * x
    return((1 + z + z^2 / 3) * exp(-z))
  }
  if (kernel == "exponential") {
    return(exp(-x))
  }
  if (kernel == "gaussian") {
    return(exp(-0.5 * x^2))
  }
  stop("unknown kernel")
}

passage_effective_range <- function(range, kernel = "matern32") {
  if (kernel == "matern12") {
    return(3 * range)
  }
  if (kernel == "matern32") {
    return(sqrt(3) * range)
  }
  if (kernel == "matern52") {
    return(sqrt(5) * range)
  }
  if (kernel == "exponential") {
    return(3 * range)
  }
  if (kernel == "gaussian") {
    return(sqrt(6) * range)
  }
  range
}

passage_satterthwaite_p <- function(q, mean, var) {
  if (!is.finite(q) || !is.finite(mean) || !is.finite(var) || mean <= 0 || var <= 0) {
    return(NA_real_)
  }
  scale <- var / (2 * mean)
  df <- 2 * mean^2 / var
  stats::pchisq(q / scale, df = df, lower.tail = FALSE)
}

passage_qform_pvalue <- function(q, lambda, method = c("auto", "davies", "satterthwaite")) {
  method <- match.arg(method)
  lambda <- as.numeric(lambda)
  lambda <- lambda[is.finite(lambda) & lambda > 0]
  if (length(lambda) == 0L) {
    return(NA_real_)
  }
  if (method %in% c("auto", "davies") && requireNamespace("CompQuadForm", quietly = TRUE)) {
    ans <- tryCatch(
      CompQuadForm::davies(q, lambda = lambda),
      error = function(e) NULL
    )
    if (!is.null(ans) && is.finite(ans$Qq)) {
      return(pmin(pmax(ans$Qq, 0), 1))
    }
    if (method == "davies") {
      warning("Davies calculation failed; falling back to Satterthwaite")
    }
  }
  passage_satterthwaite_p(q, mean = sum(lambda), var = 2 * sum(lambda^2))
}

passage_kernel_moments <- function(K, X) {
  X <- passage_prepare_design(X, nrow(K), intercept = TRUE)
  qrx <- qr(X)
  Qx <- qr.Q(qrx)[, seq_len(qrx$rank), drop = FALSE]
  KQ <- as.matrix(K %*% Qx)
  QKQ <- crossprod(Qx, KQ)

  trace_k <- sum(Matrix::diag(K))
  trace_hk <- sum(diag(QKQ))
  trace_mk <- trace_k - trace_hk

  trace_k2 <- sum(K * K)
  trace_hk2 <- sum(KQ * KQ)
  trace_hkhk <- sum(QKQ * QKQ)
  trace_mkmk <- trace_k2 - 2 * trace_hk2 + trace_hkhk

  list(
    trace_mk = as.numeric(max(trace_mk, .Machine$double.eps)),
    trace_mkmk = as.numeric(max(trace_mkmk, .Machine$double.eps)),
    rank_x = qrx$rank
  )
}
