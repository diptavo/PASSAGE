if (dir.exists("R")) {
  for (f in sort(list.files(file.path("R"), pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

set.seed(20260526)

coords <- as.matrix(expand.grid(
  x = seq(0, 1, length.out = 5),
  y = seq(0, 1, length.out = 5)
))
A <- matrix(
  c(
    0.9, 0.0,
    0.5, 0.4,
    0.2, 0.7,
    0.1, 0.3
  ),
  nrow = 4,
  byrow = TRUE
)
rownames(A) <- paste0("g", seq_len(nrow(A)))
sim <- passage_simulate_joint_factor_tmb_data(
  coords = coords,
  A = A,
  phi = c(0.30, 0.18),
  tau2 = rep(0.25, 4),
  kernel = "matern52"
)

stopifnot(all(dim(sim$Y) == c(nrow(coords), 4L)))

fake <- list(
  engine = "joint_factor_tmb_vecchia_v1",
  A = A,
  phi = c(factor_1 = 0.30, factor_2 = 0.18),
  tau2 = stats::setNames(rep(0.25, 4), rownames(A)),
  coords = coords,
  kernel = "matern52"
)
class(fake) <- c("passage_joint_factor_tmb", "list")
pve <- passage_joint_factor_pve(fake)
stopifnot(is.finite(pve$summary[["cEPSV"]]))
stopifnot(is.finite(pve$summary[["mean_propSV"]]))
stopifnot(nrow(pve$factor_table) == 2L)

if (requireNamespace("TMB", quietly = TRUE)) {
  fit <- passage_fit_joint_factor_tmb(
    Y = sim$Y,
    coords = coords,
    K_fit = 2,
    m = 6,
    kernel = "matern52",
    nlminb_control = list(iter.max = 20, eval.max = 40),
    silent = TRUE
  )
  stopifnot(inherits(fit, "passage_joint_factor_tmb"))
  stopifnot(all(dim(fit$fitted_spatial) == dim(sim$Y)))
  stopifnot(is.finite(fit$AIC))

  nested <- passage_fit_joint_factor_hypotheses(
    Y = sim$Y,
    coords = coords,
    X_base = NULL,
    Z_cell = matrix(rnorm(nrow(coords)), ncol = 1),
    K_fit = 1,
    m = 6,
    kernel = "matern32",
    nlminb_control = list(iter.max = 10, eval.max = 20),
    silent = TRUE
  )
  stopifnot(all(c("H1", "H2") %in% nested$comparison$model))
} else {
  cat("TMB not installed; compile/fitting portion skipped\n")
}

cat("PASSAGE joint factor TMB wrapper test passed\n")
