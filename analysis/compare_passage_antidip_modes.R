# Compare PASSAGE anti-double-dipping modes on small local simulations.
#
# Usage:
#   Rscript scripts/compare_passage_antidip_modes.R \
#     --n-reps=8 --n-perm=29 --n-perm-refit=9 \
#     --out=results/passage_antidip_comparison.csv

parse_args <- function(args) {
  cfg <- list(
    n_reps = 8L,
    n_perm = 29L,
    n_perm_refit = 9L,
    seed = 20260730L,
    out = file.path("results", "passage_antidip_comparison.csv")
  )
  for (arg in args) {
    if (grepl("^--n-reps=", arg)) cfg$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--n-perm-refit=", arg)) cfg$n_perm_refit <- as.integer(sub("^--n-perm-refit=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out=", arg)) cfg$out <- sub("^--out=", "", arg)
  }
  cfg
}

if (dir.exists("R")) {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
} else {
  library(PASSAGE)
}

simulate_dataset <- function(rep_id, scenario, seed) {
  set.seed(seed + rep_id * 1000L + match(scenario, c("null", "signal")))
  coords <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = 8),
    y = seq(0, 1, length.out = 8)
  ))
  n <- nrow(coords)
  g <- 80L
  Y <- matrix(rnorm(n * g, sd = 0.8), n, g)
  colnames(Y) <- paste0("g", seq_len(g))
  X <- cbind(
    depth = as.numeric(scale(rnorm(n) + 0.5 * coords[, 1])),
    detected = as.numeric(scale(rnorm(n) + 0.5 * coords[, 2]))
  )
  d <- as.matrix(dist(coords))
  z <- sqrt(3) * d / 0.22
  Kc <- (1 + z) * exp(-z) + diag(1e-6, n)
  u1 <- as.numeric(t(chol(Kc)) %*% rnorm(n))
  u2 <- sin(2 * pi * coords[, 1]) - 0.5 * cos(2 * pi * coords[, 2])

  # Background tissue structure is present in genes outside the tested pathways.
  for (jj in 41:70) {
    Y[, jj] <- Y[, jj] + 0.8 * u1 + rnorm(n, sd = 0.25)
  }
  if (scenario == "signal") {
    for (jj in 1:8) {
      Y[, jj] <- Y[, jj] + 0.9 * u2 + rnorm(n, sd = 0.20)
    }
  }
  list(
    Y = Y,
    coords = coords,
    X = X,
    pathways = list(signal = paste0("g", 1:10), null = paste0("g", 21:35))
  )
}

run_one <- function(dat, mode, rep_id, scenario, cfg) {
  fit_args <- list(
    K = 3L,
    m = 6L,
    range_grid = c(0.10, 0.22, 0.40),
    kernel = "matern32",
    ordering = "none",
    verbose = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- if (mode == "none") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = fit_args$K, m = fit_args$m, range_grid = fit_args$range_grid,
      kernel = fit_args$kernel, hypotheses = "H1", calibration = "permutation",
      n_perm = cfg$n_perm, seed = cfg$seed + rep_id, verbose = FALSE
    )
  } else if (mode == "spot_crossfit") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = fit_args$K, m = fit_args$m, range_grid = fit_args$range_grid,
      kernel = fit_args$kernel, hypotheses = "H1", calibration = "permutation",
      n_perm = cfg$n_perm, seed = cfg$seed + rep_id, anti_dip = "spot_crossfit",
      anti_dip_args = list(n_folds = 4L, seed = cfg$seed + rep_id, fit_args = fit_args),
      verbose = FALSE
    )
  } else if (mode == "pathway_holdout") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = fit_args$K, m = fit_args$m, range_grid = fit_args$range_grid,
      kernel = fit_args$kernel, hypotheses = "H1", calibration = "permutation",
      n_perm = cfg$n_perm, seed = cfg$seed + rep_id, anti_dip = "pathway_holdout",
      anti_dip_args = list(fit_args = fit_args), verbose = FALSE
    )
  } else {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = fit_args$K, m = fit_args$m, range_grid = fit_args$range_grid,
      kernel = fit_args$kernel, hypotheses = "H1", calibration = "moment",
      n_perm = cfg$n_perm_refit, seed = cfg$seed + rep_id, anti_dip = "refit_null",
      anti_dip_args = list(fit_args = fit_args, n_perm = cfg$n_perm_refit),
      verbose = FALSE
    )
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  data.frame(
    rep = rep_id,
    scenario = scenario,
    mode = mode,
    elapsed_sec = elapsed,
    fit$summary[, c("pathway", "pathway_size", "p_H1", "p_H1_moment", "R2_cca", "PSVS_range", "mean_propSV")],
    row.names = NULL
  )
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(dirname(cfg$out), recursive = TRUE, showWarnings = FALSE)
modes <- c("none", "spot_crossfit", "pathway_holdout", "refit_null")
rows <- list()
ii <- 0L
for (scenario in c("null", "signal")) {
  for (rr in seq_len(cfg$n_reps)) {
    dat <- simulate_dataset(rr, scenario, cfg$seed)
    for (mode in modes) {
      message("scenario=", scenario, " rep=", rr, " mode=", mode)
      ii <- ii + 1L
      rows[[ii]] <- run_one(dat, mode, rr, scenario, cfg)
    }
  }
}
out <- do.call(rbind, rows)
write.csv(out, cfg$out, row.names = FALSE)

summary <- aggregate(
  cbind(p_H1, elapsed_sec) ~ scenario + mode + pathway,
  data = out,
  FUN = function(x) c(median = stats::median(x), mean = mean(x))
)
print(summary)
cat("Wrote ", cfg$out, "\n", sep = "")
