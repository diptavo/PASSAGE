if (dir.exists("R")) {
  for (f in sort(list.files(file.path("R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

set.seed(20260730)

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 12),
  y = seq(0, 1, length.out = 12)
))
n <- nrow(coords)
g <- 80
Y <- matrix(rnorm(n * g, sd = 0.8), nrow = n, ncol = g)
colnames(Y) <- paste0("g", seq_len(g))

d <- as.matrix(dist(coords))
K <- exp(-d / 0.22) + diag(1e-6, n)
v1 <- as.numeric(t(chol(K)) %*% rnorm(n))
v2 <- 0.75 * v1 + as.numeric(t(chol(K)) %*% rnorm(n, sd = 0.45))
for (jj in 1:8) {
  Y[, jj] <- Y[, jj] + 1.4 * v1 + rnorm(n, sd = 0.25)
}
for (jj in 9:16) {
  Y[, jj] <- Y[, jj] + 1.4 * v2 + rnorm(n, sd = 0.25)
}

signal <- paste0("g", 1:16)
null <- paste0("g", 41:60)

pst_signal <- passage_pathway_spatial_transferability(
  Y, coords, signal, n_gene_splits = 4, n_folds = 4,
  k_neighbors = 8, n_perm = 25, seed = 1
)
pst_null <- passage_pathway_spatial_transferability(
  Y, coords, null, n_gene_splits = 4, n_folds = 4,
  k_neighbors = 8, n_perm = 25, seed = 2
)

Y_perm <- Y
Y_perm[, signal] <- Y[sample.int(n), signal]
pst_perm <- passage_pathway_spatial_transferability(
  Y_perm, coords, signal, n_gene_splits = 4, n_folds = 4,
  k_neighbors = 8, n_perm = 25, seed = 3
)

stopifnot(is.finite(pst_signal$statistic))
stopifnot(is.finite(pst_signal$p))
stopifnot(pst_signal$statistic > pst_null$statistic)
stopifnot(pst_signal$statistic > pst_perm$statistic)
stopifnot(nrow(pst_signal$split_table) > 0L)

cat("PASSAGE pathway transferability test passed\n")
