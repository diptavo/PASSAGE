if (dir.exists("R")) {
  source(file.path("R", "spapath.R"))
} else {
  library(PASSAGE)
}

set.seed(42)

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 10),
  y = seq(0, 1, length.out = 10)
))
n <- nrow(coords)
p <- 30

Y <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(Y) <- paste0("g", seq_len(p))

D <- as.matrix(dist(coords))
K <- exp(-D / 0.2) + diag(1e-6, n)
spatial_signal <- as.numeric(t(chol(K)) %*% rnorm(n))
Y[, "g1"] <- Y[, "g1"] + 1.5 * spatial_signal
Y[, "g2"] <- Y[, "g2"] + 1.0 * spatial_signal

pathways <- list(
  signal = paste0("g", 1:5),
  null = paste0("g", 6:15)
)

fit <- spapath_test(
  Y = Y,
  coords = coords,
  pathways = pathways,
  m = 8,
  ranges = c(0.15, 0.30),
  n_sim = 200,
  seed = 7,
  verbose = FALSE
)

stopifnot(inherits(fit, "spapath_result"))
stopifnot(all(c("signal", "null") %in% fit$results$pathway))
stopifnot(all(c("p_value", "fdr", "eSPVE_any", "driver_genes") %in% names(fit$results)))
stopifnot(is.list(fit$results$driver_genes))
stopifnot(fit$results$p_value[fit$results$pathway == "signal"] <= 0.1)

cat("SpaPath smoke test passed\n")
