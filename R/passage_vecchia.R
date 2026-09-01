# Vecchia precision construction for Matérn-like spatial factors.

passage_vecchia_precision <- function(coords,
                                      range,
                                      m = 20L,
                                      kernel = c("matern32", "matern12", "matern52", "exponential", "gaussian"),
                                      ordering = c("coordinate", "none"),
                                      jitter = 1e-8) {
  kernel <- match.arg(kernel)
  ordering <- match.arg(ordering)
  coords <- passage_check_coords(coords)
  n <- nrow(coords)
  m <- min(as.integer(m), n - 1L)
  if (m < 1L) {
    stop("m must be at least 1 and smaller than nrow(coords)")
  }
  if (!is.finite(range) || range <= 0) {
    stop("range must be positive and finite")
  }

  ord <- passage_order_locations(coords, method = ordering)
  inv_ord <- integer(n)
  inv_ord[ord] <- seq_len(n)
  coords_ord <- coords[ord, , drop = FALSE]

  rows <- integer(n + n * m)
  cols <- integer(n + n * m)
  vals <- numeric(n + n * m)
  nn_index <- matrix(NA_integer_, nrow = n, ncol = m)
  cond_var <- rep(1, n)
  ptr <- 0L

  for (i in seq_len(n)) {
    ptr <- ptr + 1L
    rows[ptr] <- i
    cols[ptr] <- i
    vals[ptr] <- 1

    if (i == 1L) {
      next
    }
    prev <- seq_len(i - 1L)
    if (length(prev) <= m) {
      nb <- prev
    } else {
      d2 <- rowSums((coords_ord[prev, , drop = FALSE] -
        matrix(coords_ord[i, ], nrow = length(prev), ncol = ncol(coords_ord), byrow = TRUE))^2)
      nb <- prev[order(d2)[seq_len(m)]]
    }

    d_nn <- as.matrix(stats::dist(coords_ord[nb, , drop = FALSE]))
    c_nn <- as.matrix(passage_kernel_corr(d_nn, range = range, kernel = kernel))
    diag(c_nn) <- 1 + jitter
    d_i <- sqrt(rowSums((coords_ord[nb, , drop = FALSE] -
      matrix(coords_ord[i, ], nrow = length(nb), ncol = ncol(coords_ord), byrow = TRUE))^2))
    c_i <- passage_kernel_corr(d_i, range = range, kernel = kernel)
    b <- tryCatch(
      as.numeric(solve(c_nn, c_i)),
      error = function(e) as.numeric(solve(c_nn + diag(jitter, nrow(c_nn)), c_i))
    )
    d_cond <- 1 - sum(c_i * b)
    d_cond <- max(as.numeric(d_cond), jitter)
    cond_var[i] <- d_cond
    nn_index[i, seq_along(nb)] <- nb

    take <- seq_along(nb)
    j <- ptr + take
    rows[j] <- i
    cols[j] <- nb
    vals[j] <- -b
    ptr <- ptr + length(nb)
  }

  rows <- rows[seq_len(ptr)]
  cols <- cols[seq_len(ptr)]
  vals <- vals[seq_len(ptr)]
  L <- Matrix::sparseMatrix(i = rows, j = cols, x = vals, dims = c(n, n))
  Dinv <- Matrix::Diagonal(n, 1 / cond_var)
  Q <- Matrix::crossprod(L, Dinv %*% L)
  Q <- Matrix::forceSymmetric(Q, uplo = "U")

  out <- list(
    Q = Q,
    L = L,
    cond_var = cond_var,
    log_det_Q = -sum(log(cond_var)),
    nn_index = nn_index,
    ord = ord,
    inv_ord = inv_ord,
    coords_ord = coords_ord,
    m = m,
    range = range,
    kernel = kernel,
    ordering = ordering
  )
  class(out) <- c("passage_vecchia", "list")
  out
}

passage_order_locations <- function(coords, method = c("coordinate", "none")) {
  method <- match.arg(method)
  n <- nrow(coords)
  if (method == "none") {
    return(seq_len(n))
  }
  pc <- tryCatch(
    stats::prcomp(coords, center = TRUE, scale. = FALSE)$x[, 1L],
    error = function(e) rowSums(coords)
  )
  order(pc, coords[, 1L], coords[, 2L], seq_len(n))
}

passage_vecchia_fit_range_grid <- function(v,
                                           coords,
                                           range_grid,
                                           m = 20L,
                                           kernel = "matern32",
                                           ordering = "coordinate") {
  v <- as.numeric(v)
  if (!all(is.finite(v))) {
    stop("factor score vector contains non-finite values")
  }
  fits <- vector("list", length(range_grid))
  obj <- rep(Inf, length(range_grid))
  for (rr in seq_along(range_grid)) {
    vc <- passage_vecchia_precision(
      coords = coords,
      range = range_grid[rr],
      m = m,
      kernel = kernel,
      ordering = ordering
    )
    vo <- v[vc$ord]
    quad <- as.numeric(crossprod(vo, vc$Q %*% vo))
    sigma2 <- max(quad / length(v), .Machine$double.eps)
    obj[rr] <- length(v) * log(sigma2) - vc$log_det_Q
    fits[[rr]] <- list(vecchia = vc, sigma2 = sigma2, objective = obj[rr])
  }
  best <- which.min(obj)
  list(
    range = range_grid[best],
    sigma2 = fits[[best]]$sigma2,
    vecchia = fits[[best]]$vecchia,
    objective = obj[best],
    grid = data.frame(range = range_grid, objective = obj)
  )
}

passage_sparse_covariance_from_vecchia <- function(vecchia) {
  n <- nrow(vecchia$coords_ord)
  nn <- vecchia$nn_index
  from <- rep(seq_len(n), ncol(nn))
  to <- as.vector(nn)
  keep <- !is.na(to)
  from <- from[keep]
  to <- to[keep]
  d <- sqrt(rowSums((vecchia$coords_ord[from, , drop = FALSE] -
    vecchia$coords_ord[to, , drop = FALSE])^2))
  w <- passage_kernel_corr(d, range = vecchia$range, kernel = vecchia$kernel)
  Matrix::sparseMatrix(
    i = c(from, to, seq_len(n)),
    j = c(to, from, seq_len(n)),
    x = c(w, w, rep(1, n)),
    dims = c(n, n)
  )
}
