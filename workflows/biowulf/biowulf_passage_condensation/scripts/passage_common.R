gene_symbols_from_obj <- function(obj) {
  gd <- obj$gene_data
  candidates <- c("gene_name", "gene_name_short", "gene_symbol", "symbol", "Symbol", "gene", "gene_id", "ensembl")
  for (cc in intersect(candidates, colnames(gd))) {
    x <- as.character(gd[[cc]])
    if (length(x) == ncol(obj$Y) && sum(!is.na(x) & nzchar(x)) > 0.5 * length(x)) return(toupper(x))
  }
  toupper(colnames(obj$Y))
}

collapse_to_symbols <- function(Y, symbols) {
  symbols <- toupper(symbols)
  ok <- !is.na(symbols) & nzchar(symbols) & is.finite(colSums(Y))
  Y <- Y[, ok, drop = FALSE]
  symbols <- symbols[ok]
  if (anyDuplicated(symbols)) Y <- t(rowsum(t(Y), group = symbols, reorder = FALSE)) else colnames(Y) <- symbols
  keep <- apply(Y, 2L, stats::var) > 1e-8
  Y[, keep, drop = FALSE]
}

zscore_cols <- function(M) {
  M <- as.matrix(M)
  mu <- colMeans(M, na.rm = TRUE)
  S <- sweep(M, 2L, mu, "-")
  sdv <- sqrt(pmax(colMeans(S * S, na.rm = TRUE), 1e-8))
  sweep(S, 2L, sdv, "/")
}

scale_coords01 <- function(coords) {
  coords <- as.matrix(coords)[, seq_len(2L), drop = FALSE]
  coords <- sweep(coords, 2L, apply(coords, 2L, min, na.rm = TRUE), "-")
  denom <- pmax(apply(coords, 2L, max, na.rm = TRUE), .Machine$double.eps)
  coords <- sweep(coords, 2L, denom, "/")
  coords[!is.finite(coords)] <- 0
  coords
}

knn_index <- function(coords, k = 8L) {
  k <- min(k, nrow(coords) - 1L)
  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  t(apply(D, 1L, function(x) order(x)[seq_len(k)]))
}

neighbor_lag <- function(Z, nn) {
  out <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  for (kk in seq_len(ncol(nn))) out <- out + Z[nn[, kk], , drop = FALSE]
  out / ncol(nn)
}

pc1_fraction <- function(M) {
  if (ncol(M) <= 1L) return(1)
  ev <- suppressWarnings(eigen(crossprod(M), symmetric = TRUE, only.values = TRUE)$values)
  ev <- pmax(ev, 0)
  if (!sum(ev) > 0) return(0)
  max(ev) / sum(ev)
}

effective_breadth <- function(w, p) {
  w <- pmax(w, 0)
  if (!sum(w) > 0 || p <= 1L) return(0)
  (sum(w)^2 / sum(w^2)) / p
}

hotspot_score <- function(H, nn) {
  h <- rowMeans(H)
  if (!is.finite(mean(h)) || mean(h) <= 0) return(0)
  hlag <- rep(0, length(h))
  for (kk in seq_len(ncol(nn))) hlag <- hlag + h[nn[, kk]]
  hlag <- hlag / ncol(nn)
  mean(h * hlag) / (mean(h)^2)
}

transport_alignment <- function(absM) {
  denom <- colSums(absM)
  denom[denom <= 0 | !is.finite(denom)] <- 1
  P <- sweep(absM, 2L, denom, "/")
  pbar <- rowMeans(P)
  affinity <- colSums(sqrt(sweep(P, 1L, pbar, "*")))
  concentration <- sqrt(nrow(P) * sum(pbar * pbar))
  mean(affinity) * concentration
}

make_bins <- function(x, n = 4L) {
  x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE))
  if (length(qs) <= 2L) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE))
}

sample_matched <- function(target_bins, bins, size, exclude = integer(0)) {
  out <- integer(0)
  universe <- setdiff(seq_along(bins), exclude)
  tab <- table(target_bins)
  for (bb in names(tab)) {
    pool <- intersect(which(bins == bb), universe)
    n <- as.integer(tab[[bb]])
    if (length(pool) >= n) out <- c(out, sample(pool, n))
  }
  if (length(out) < size) {
    pool <- setdiff(universe, out)
    need <- size - length(out)
    if (!length(pool)) pool <- universe
    out <- c(out, sample(pool, need, replace = length(pool) < need))
  }
  sample(out, size)
}

build_context <- function(Y, X, coords) {
  Z <- qr.resid(qr(as.matrix(X)), Y)
  Z <- zscore_cols(Z)
  coords <- scale_coords01(coords)
  nn8 <- knn_index(coords, 8L)
  nn24 <- knn_index(coords, 24L)
  Lag8 <- neighbor_lag(Z, nn8)
  Lag24 <- neighbor_lag(Z, nn24)
  moran8 <- colSums(Z * Lag8) / pmax(colSums(Z * Z), .Machine$double.eps)
  moran24 <- colSums(Z * Lag24) / pmax(colSums(Z * Z), .Machine$double.eps)
  med <- stats::median(moran8, na.rm = TRUE)
  madv <- stats::mad(moran8, constant = 1.4826, na.rm = TRUE)
  gene_score_z <- if (is.finite(madv) && madv > 1e-8) (moran8 - med) / madv else as.numeric(scale(moran8))
  Hot <- apply(abs(Z), 2L, function(x) x >= stats::quantile(x, 0.90, na.rm = TRUE))
  storage.mode(Hot) <- "double"
  list(Z = Z, Lag8 = Lag8, Lag24 = Lag24, moran8 = moran8, moran24 = moran24,
       gene_score_z = gene_score_z, Hot = Hot, nn8 = nn8, genes = colnames(Y))
}

stat_all <- function(idx, ctx) {
  p <- length(idx)
  M <- ctx$Z[, idx, drop = FALSE]
  L8 <- ctx$Lag8[, idx, drop = FALSE]
  m8 <- pmax(ctx$moran8[idx], 0)
  m24 <- pmax(ctx$moran24[idx], 0)
  pc1 <- pc1_fraction(L8)
  breadth <- effective_breadth(m8, p)
  lowpass <- sum(L8 * L8) / pmax(sum(M * M), .Machine$double.eps)
  range_contrast <- pmax(m8 - m24, 0)
  c(
    CSPS = sum(m8) * pc1 * breadth,
    GSPS = lowpass * pc1 * breadth,
    MMP = sum(M * L8) / pmax(sum(M * M), .Machine$double.eps),
    CSV = mean(range_contrast) * (1 / (1 + stats::mad(range_contrast, constant = 1, na.rm = TRUE))) * breadth,
    HCPS = hotspot_score(ctx$Hot[, idx, drop = FALSE], ctx$nn8),
    OTSAS = transport_alignment(abs(M)) * breadth,
    score_z = mean(ctx$gene_score_z[idx], na.rm = TRUE) * sqrt(p),
    score_z_robust = {
      qq <- stats::quantile(ctx$gene_score_z, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
      mean(pmin(pmax(ctx$gene_score_z[idx], qq[[1L]]), qq[[2L]]), na.rm = TRUE) * sqrt(p)
    }
  )
}

gpd_tail_p <- function(null_values, observed, tail_fraction = 0.20, min_tail = 25L, max_exact_exceedances = 10L) {
  x <- sort(null_values[is.finite(null_values)], decreasing = TRUE)
  n <- length(x)
  n_exceed <- sum(x >= observed)
  emp <- (1 + n_exceed) / (n + 1)
  if (n_exceed > max_exact_exceedances) return(c(p_gpd = emp, p_empirical = emp, gpd_used = 0, gpd_shape = NA, gpd_scale = NA, gpd_threshold = NA, gpd_tail_n = NA))
  if (n < 100L || !is.finite(observed)) return(c(p_gpd = emp, p_empirical = emp, gpd_used = 0, gpd_shape = NA, gpd_scale = NA, gpd_threshold = NA, gpd_tail_n = NA))
  m <- max(min_tail, ceiling(n * tail_fraction))
  m <- min(m, n - 5L)
  threshold <- x[[m + 1L]]
  if (!is.finite(threshold) || observed <= threshold) return(c(p_gpd = emp, p_empirical = emp, gpd_used = 0, gpd_shape = NA, gpd_scale = NA, gpd_threshold = threshold, gpd_tail_n = m))
  y <- x[seq_len(m)] - threshold
  y <- y[y > 0 & is.finite(y)]
  if (length(y) < min_tail || stats::var(y) <= 0) return(c(p_gpd = emp, p_empirical = emp, gpd_used = 0, gpd_shape = NA, gpd_scale = NA, gpd_threshold = threshold, gpd_tail_n = length(y)))
  nll <- function(par) {
    xi <- par[[1L]]
    beta <- exp(par[[2L]])
    z <- 1 + xi * y / beta
    if (beta <= 0 || any(z <= 0) || !is.finite(xi)) return(Inf)
    length(y) * log(beta) + (1 / xi + 1) * sum(log(z))
  }
  nll_exp <- function(log_beta) length(y) * log_beta + sum(y / exp(log_beta))
  fit <- try(stats::optim(c(0.1, log(mean(y))), nll, method = "Nelder-Mead", control = list(maxit = 1000)), silent = TRUE)
  if (inherits(fit, "try-error") || !is.finite(fit$value)) {
    fb <- try(stats::optim(log(mean(y)), nll_exp, method = "Brent", lower = log(mean(y) / 100), upper = log(mean(y) * 100)), silent = TRUE)
    beta <- if (inherits(fb, "try-error")) mean(y) else exp(fb$par)
    surv <- exp(-(observed - threshold) / beta)
    return(c(p_gpd = max(.Machine$double.xmin, min(1, length(y) / n * surv)), p_empirical = emp, gpd_used = 1, gpd_shape = 0, gpd_scale = beta, gpd_threshold = threshold, gpd_tail_n = length(y)))
  }
  xi <- fit$par[[1L]]
  beta <- exp(fit$par[[2L]])
  yy <- observed - threshold
  surv <- if (abs(xi) < 1e-6) exp(-yy / beta) else {
    z <- 1 + xi * yy / beta
    if (z <= 0) 0 else z^(-1 / xi)
  }
  p <- max(.Machine$double.xmin, min(1, length(y) / n * surv))
  c(p_gpd = min(emp, p), p_empirical = emp, gpd_used = 1, gpd_shape = xi, gpd_scale = beta, gpd_threshold = threshold, gpd_tail_n = length(y))
}
