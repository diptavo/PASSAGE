# Minimal PASSAGE simulation.
#
# Run from the PASSAGE repository root:
#   Rscript inst/examples/run_minimal_simulation.R

if (!"package:PASSAGE" %in% search()) {
  if (file.exists(file.path("R", "spapath.R"))) {
    source(file.path("R", "spapath.R"))
  } else {
    library(PASSAGE)
  }
}

set.seed(1)

n_side <- 18
coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = n_side),
  y = seq(0, 1, length.out = n_side)
))
n <- nrow(coords)
p <- 120

feature_names <- paste0("gene", seq_len(p))

matern32_dense <- function(coords, range) {
  D <- as.matrix(stats::dist(coords))
  x <- sqrt(3) * D / range
  (1 + x) * exp(-x)
}

simulate_gp <- function(coords, range, sd = 1) {
  K <- matern32_dense(coords, range)
  K <- K + diag(1e-6, nrow(K))
  as.numeric(t(chol(K)) %*% stats::rnorm(nrow(K))) * sd
}

Y <- matrix(stats::rnorm(n * p, sd = 1), nrow = n, ncol = p)
colnames(Y) <- feature_names

# Sparse pathway signal: three genes drive spatial variability.
sparse_drivers <- c("gene3", "gene7", "gene11")
for (g in sparse_drivers) {
  Y[, g] <- Y[, g] + simulate_gp(coords, range = 0.18, sd = 1.2)
}

# Diffuse pathway signal: many weak spatial genes.
diffuse_drivers <- paste0("gene", 41:55)
for (g in diffuse_drivers) {
  Y[, g] <- Y[, g] + simulate_gp(coords, range = 0.30, sd = 0.45)
}

pathways <- list(
  sparse_pathway = paste0("gene", 1:20),
  null_pathway = paste0("gene", 21:40),
  diffuse_pathway = paste0("gene", 41:70),
  mixed_null_large = paste0("gene", 71:120)
)

fit <- spapath_test(
  Y = Y,
  coords = coords,
  pathways = pathways,
  m = 12,
  ranges = c(0.10, 0.20, 0.35),
  n_sim = 500,
  seed = 99
)

print(fit)

cat("\nDriver genes in sparse_pathway:\n")
print(fit$results$driver_genes[[match("sparse_pathway", fit$results$pathway)]])
