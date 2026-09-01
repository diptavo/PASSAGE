set.seed(20260526)

root <- getwd()
if (dir.exists(file.path(root, "R"))) {
  for (f in sort(list.files(file.path(root, "R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

n <- 34L
p <- 14L
coords <- cbind(runif(n), runif(n))
D <- as.matrix(dist(coords))
z <- sqrt(3) * D / 0.25
Kc <- (1 + z) * exp(-z)
L <- chol(Kc + diag(1e-6, n))
u1 <- as.numeric(t(L) %*% rnorm(n))
u2 <- sin(2 * pi * coords[, 1]) + 0.5 * cos(2 * pi * coords[, 2])
A <- matrix(rnorm(p * 2, sd = 0.8), p, 2)
Y <- cbind(u1, u2) %*% t(A) + matrix(rnorm(n * p, sd = 0.35), n, p)
colnames(Y) <- paste0("g", seq_len(p))
X <- cbind(batch = rep(c(0, 1), length.out = n))

common <- list(
  K = 2L,
  range_grid = c(0.12, 0.25, 0.5),
  m = 5L,
  kernel = "matern32",
  ordering = "none",
  verbose = FALSE
)

engines <- list(
  spatial_basis = do.call(passage_fit_engine_spatial_basis, c(
    list(Y = Y, coords = coords, X = X, n_basis = 8L, lambda_grid = c(0.1, 1)),
    common
  )),
  smoothed_pca = do.call(passage_fit_engine_smoothed_pca, c(
    list(Y = Y, coords = coords, X = X, n_basis = 8L, lambda_grid = c(0.1, 1)),
    common
  )),
  nmf = do.call(passage_fit_engine_nmf, c(
    list(Y = Y, coords = coords, X = X, n_iter = 25L, n_basis = 8L, lambda_grid = c(0.1, 1)),
    common
  )),
  alternating_gp = do.call(passage_fit_engine_alternating_gp, c(
    list(Y = Y, coords = coords, X = X, n_iter = 1L, smooth_penalty = 0.5),
    common
  ))
)

for (nm in names(engines)) {
  fit <- engines[[nm]]
  stopifnot(inherits(fit, "passage_engine"))
  stopifnot(fit$K >= 1L)
  stopifnot(nrow(fit$A) == p)
  stopifnot(nrow(fit$V) == n)
  stopifnot(all(dim(fit$fitted_spatial) == dim(Y)))
  stopifnot(all(is.finite(fit$D)))
  metric <- passage_pathway_covariance_metrics(fit, Y, pathway = colnames(Y)[1:6], X = X)
  stopifnot(is.finite(metric$summary[["cEPSV"]]))
}

cmp <- passage_compare_factor_engines(
  Y = Y,
  coords = coords,
  X = X,
  methods = c("spatial_basis", "smoothed_pca"),
  fit_args = list(K = 2L, n_basis = 8L, lambda_grid = c(0.1, 1),
                  range_grid = c(0.12, 0.25), m = 5L, verbose = FALSE),
  verbose = FALSE
)
stopifnot(nrow(cmp$summary) == 2L)
stopifnot(all(cmp$summary$status == "fit"))

cat("PASSAGE additional factor engines test passed\n")
