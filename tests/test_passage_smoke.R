if (dir.exists("R")) {
  for (f in sort(list.files(file.path("R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

set.seed(20260523)

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 9),
  y = seq(0, 1, length.out = 9)
))
n <- nrow(coords)
g <- 50
Y <- matrix(rnorm(n * g), nrow = n, ncol = g)
colnames(Y) <- paste0("g", seq_len(g))

d <- as.matrix(dist(coords))
z <- sqrt(3) * d / 0.25
K <- (1 + z) * exp(-z) + diag(1e-6, n)
v <- as.numeric(t(chol(K)) %*% rnorm(n))
for (jj in 1:6) {
  Y[, jj] <- Y[, jj] + 1.0 * v + rnorm(n, sd = 0.15)
}

pathways <- list(
  signal = paste0("g", 1:10),
  null = paste0("g", 21:35)
)

fit <- passage_run(
  Y = Y,
  coords = coords,
  pathways = pathways,
  K = 3,
  m = 8,
  range_grid = c(0.10, 0.20, 0.35),
  hypotheses = "H1",
  verbose = FALSE
)

stopifnot(inherits(fit, "passage_run"))
stopifnot(inherits(fit$engine, "passage_engine"))
stopifnot(all(c("signal", "null") %in% fit$summary$pathway))
stopifnot(all(c("p_H1", "R2_cca", "PSVS_range", "spasset_genes_H1") %in% names(fit$summary)))
stopifnot(is.list(fit$summary$spasset_genes_H1))
stopifnot(is.finite(fit$summary$p_H1[fit$summary$pathway == "signal"]))
stopifnot(fit$summary$p_H1[fit$summary$pathway == "signal"] <= 0.2)
stopifnot(fit$summary$PSVS_range[fit$summary$pathway == "signal"] >
  fit$summary$PSVS_range[fit$summary$pathway == "null"])

cat("PASSAGE smoke test passed\n")
