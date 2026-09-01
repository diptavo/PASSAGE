# Minimal PASSAGE simulation.
#
# Run from the PASSAGE repository root:
#   Rscript inst/examples/run_passage_simulation.R

if (!"package:PASSAGE" %in% search()) {
  if (dir.exists("R")) {
    for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
      source(f)
    }
  } else {
    library(PASSAGE)
  }
}

set.seed(20260523)

n_side <- 9
coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = n_side),
  y = seq(0, 1, length.out = n_side)
))
n <- nrow(coords)
g <- 50
genes <- paste0("g", seq_len(g))

matern32 <- function(coords, range) {
  d <- as.matrix(stats::dist(coords))
  z <- sqrt(3) * d / range
  (1 + z) * exp(-z)
}

Ksp <- matern32(coords, range = 0.25) + diag(1e-6, n)
v <- as.numeric(t(chol(Ksp)) %*% rnorm(n))

Y <- matrix(rnorm(n * g), nrow = n, ncol = g)
colnames(Y) <- genes
for (jj in 1:6) {
  Y[, jj] <- Y[, jj] + 1.0 * v + rnorm(n, sd = 0.15)
}

pathways <- list(
  spatial_pathway = paste0("g", 1:10),
  null_pathway = paste0("g", 21:35)
)

fit <- passage_run(
  Y = Y,
  coords = coords,
  pathways = pathways,
  K = 3,
  m = 8,
  range_grid = c(0.10, 0.20, 0.35),
  hypotheses = "H1",
  verbose = TRUE
)

print(fit)
