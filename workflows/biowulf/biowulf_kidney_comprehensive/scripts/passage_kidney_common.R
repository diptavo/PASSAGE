gene_symbols_from_obj <- function(obj) {
  gd <- obj$gene_data
  candidates <- c("gene_name", "gene_name_short", "gene_symbol", "symbol", "Symbol",
                  "gene", "gene_id", "ensembl")
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
  if (anyDuplicated(symbols)) {
    Y <- t(rowsum(t(Y), group = symbols, reorder = FALSE))
  } else {
    colnames(Y) <- symbols
  }
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

pc1_loading <- function(M) {
  if (ncol(M) <= 1L) return(1)
  eg <- suppressWarnings(eigen(crossprod(M), symmetric = TRUE))
  v <- abs(eg$vectors[, 1L])
  v / pmax(sum(v), .Machine$double.eps)
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

hotspot_gene_contrib <- function(H, idx, nn) {
  Hp <- H[, idx, drop = FALSE]
  hpath <- rowMeans(Hp)
  lagp <- rep(0, length(hpath))
  for (kk in seq_len(ncol(nn))) lagp <- lagp + hpath[nn[, kk]]
  lagp <- lagp / ncol(nn)
  colMeans(sweep(Hp, 1L, lagp, "*"))
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

transport_gene_contrib <- function(absM) {
  denom <- colSums(absM)
  denom[denom <= 0 | !is.finite(denom)] <- 1
  P <- sweep(absM, 2L, denom, "/")
  pbar <- rowMeans(P)
  colSums(sqrt(sweep(P, 1L, pbar, "*")))
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

build_context <- function(Y, X, coords, keep = NULL) {
  if (is.null(keep)) keep <- seq_len(nrow(Y))
  Yb <- Y[keep, , drop = FALSE]
  Xb <- as.matrix(X[keep, , drop = FALSE])
  coordsb <- scale_coords01(coords[keep, , drop = FALSE])
  Z <- qr.resid(qr(Xb), Yb)
  Z <- zscore_cols(Z)
  nn8 <- knn_index(coordsb, 8L)
  nn24 <- knn_index(coordsb, 24L)
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
  signal <- sum(m8)
  lowpass <- sum(L8 * L8) / pmax(sum(M * M), .Machine$double.eps)
  mult_moran <- sum(M * L8) / pmax(sum(M * M), .Machine$double.eps)
  range_contrast <- pmax(m8 - m24, 0)
  range_agreement <- 1 / (1 + stats::mad(range_contrast, constant = 1, na.rm = TRUE))
  hcps <- hotspot_score(ctx$Hot[, idx, drop = FALSE], ctx$nn8)
  otsas <- transport_alignment(abs(M)) * breadth
  c(
    CSPS = signal * pc1 * breadth,
    GSPS = lowpass * pc1 * breadth,
    MMP = mult_moran,
    CSV = mean(range_contrast) * range_agreement * breadth,
    HCPS = hcps,
    OTSAS = otsas,
    score_z = mean(ctx$gene_score_z[idx], na.rm = TRUE) * sqrt(p),
    score_z_robust = {
      qq <- stats::quantile(ctx$gene_score_z, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
      mean(pmin(pmax(ctx$gene_score_z[idx], qq[[1L]]), qq[[2L]]), na.rm = TRUE) * sqrt(p)
    }
  )
}

driver_scores_for_idx <- function(idx, ctx, statistic) {
  genes <- ctx$genes[idx]
  M <- ctx$Z[, idx, drop = FALSE]
  L8 <- ctx$Lag8[, idx, drop = FALSE]
  m8 <- pmax(ctx$moran8[idx], 0)
  m24 <- pmax(ctx$moran24[idx], 0)
  if (statistic == "CSPS") {
    v <- pc1_loading(L8)
    score <- m8 * v * length(v)
  } else if (statistic == "GSPS") {
    v <- pc1_loading(L8)
    low <- colSums(L8 * L8) / pmax(colSums(M * M), .Machine$double.eps)
    score <- pmax(low, 0) * v * length(v)
  } else if (statistic == "MMP") {
    score <- pmax(colSums(M * L8) / pmax(colSums(M * M), .Machine$double.eps), 0)
  } else if (statistic == "CSV") {
    score <- pmax(m8 - m24, 0)
  } else if (statistic == "HCPS") {
    score <- hotspot_gene_contrib(ctx$Hot, idx, ctx$nn8)
  } else if (statistic == "OTSAS") {
    score <- transport_gene_contrib(abs(M))
  } else if (statistic == "score_z") {
    score <- pmax(ctx$gene_score_z[idx], 0)
  } else if (statistic == "score_z_robust") {
    z <- ctx$gene_score_z
    qq <- stats::quantile(z, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    score <- pmax(pmin(pmax(z[idx], qq[[1L]]), qq[[2L]]), 0)
  } else {
    stop("Unknown statistic: ", statistic)
  }
  stats::setNames(as.numeric(score), genes)
}

asset_select <- function(scores, genes, top_k = 10L, alpha = 2) {
  z <- as.numeric(scores)
  z[!is.finite(z)] <- min(z[is.finite(z)], 0)
  if (length(z) >= 3L) {
    med <- stats::median(z, na.rm = TRUE)
    madv <- stats::mad(z, constant = 1.4826, na.rm = TRUE)
    if (is.finite(madv) && madv > 1e-8) z <- (z - med) / madv
  }
  positive <- z
  positive[positive < 0] <- 0
  w <- exp(alpha * (positive - max(positive, na.rm = TRUE)))
  w[positive <= 0] <- 0
  if (sum(w, na.rm = TRUE) <= 0) w <- rep(1, length(z))
  w <- w / sum(w, na.rm = TRUE)
  ord <- order(w, decreasing = TRUE)
  keep <- ord[seq_len(min(top_k, length(ord)))]
  data.frame(gene = genes[keep], rank = seq_along(keep), driver_weight = w[keep],
             raw_score = scores[keep], stringsAsFactors = FALSE)
}

spatial_block_keep <- function(coords, keep_frac = 0.80, grid_n = 5L) {
  c01 <- scale_coords01(coords)
  gx <- pmin(grid_n, pmax(1L, floor(c01[, 1L] * grid_n) + 1L))
  gy <- pmin(grid_n, pmax(1L, floor(c01[, 2L] * grid_n) + 1L))
  block <- paste(gx, gy, sep = "_")
  blocks <- unique(block)
  n_keep <- max(1L, ceiling(length(blocks) * keep_frac))
  which(block %in% sample(blocks, n_keep))
}

jaccard_mean <- function(sets) {
  if (length(sets) < 2L) return(NA_real_)
  vals <- numeric(0)
  for (i in seq_len(length(sets) - 1L)) {
    for (j in seq.int(i + 1L, length(sets))) {
      u <- union(sets[[i]], sets[[j]])
      vals <- c(vals, if (length(u)) length(intersect(sets[[i]], sets[[j]])) / length(u) else NA_real_)
    }
  }
  mean(vals, na.rm = TRUE)
}
