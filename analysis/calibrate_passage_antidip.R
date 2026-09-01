# PASSAGE anti-double-dipping calibration checks.
#
# The script runs two bounded calibration panels:
#   1. Simulated null type-I checks with known null pathways.
#   2. Real breast Visium matched-random gene-set checks using prepared data.
#
# Usage:
#   Rscript scripts/calibrate_passage_antidip.R \
#     --sim-reps=20 --sim-refit-reps=8 --sim-n-perm=49 --sim-refit-perm=9 \
#     --real-random=10 --real-n-perm=19 \
#     --out-dir=results/passage_antidip_calibration_quick

parse_args <- function(args) {
  cfg <- list(
    seed = 20260730L,
    out_dir = file.path("results", "passage_antidip_calibration_quick"),
    sim_reps = 20L,
    sim_refit_reps = 8L,
    sim_n_perm = 49L,
    sim_refit_perm = 9L,
    sim_n_side = 8L,
    sim_g = 90L,
    sim_pathway_size = 12L,
    real_random = 10L,
    real_n_perm = 19L,
    real_modes = c("none", "spot_crossfit", "pathway_holdout"),
    K = 3L,
    m = 8L,
    n_folds = 4L
  )
  for (arg in args) {
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--sim-reps=", arg)) cfg$sim_reps <- as.integer(sub("^--sim-reps=", "", arg))
    if (grepl("^--sim-refit-reps=", arg)) cfg$sim_refit_reps <- as.integer(sub("^--sim-refit-reps=", "", arg))
    if (grepl("^--sim-n-perm=", arg)) cfg$sim_n_perm <- as.integer(sub("^--sim-n-perm=", "", arg))
    if (grepl("^--sim-refit-perm=", arg)) cfg$sim_refit_perm <- as.integer(sub("^--sim-refit-perm=", "", arg))
    if (grepl("^--sim-n-side=", arg)) cfg$sim_n_side <- as.integer(sub("^--sim-n-side=", "", arg))
    if (grepl("^--sim-g=", arg)) cfg$sim_g <- as.integer(sub("^--sim-g=", "", arg))
    if (grepl("^--sim-pathway-size=", arg)) cfg$sim_pathway_size <- as.integer(sub("^--sim-pathway-size=", "", arg))
    if (grepl("^--real-random=", arg)) cfg$real_random <- as.integer(sub("^--real-random=", "", arg))
    if (grepl("^--real-n-perm=", arg)) cfg$real_n_perm <- as.integer(sub("^--real-n-perm=", "", arg))
    if (grepl("^--real-modes=", arg)) cfg$real_modes <- strsplit(sub("^--real-modes=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--K=", arg)) cfg$K <- as.integer(sub("^--K=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
  }
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

matern32_dense <- function(coords, range) {
  d <- as.matrix(stats::dist(coords))
  z <- sqrt(3) * d / range
  (1 + z) * exp(-z)
}

simulate_gp <- function(coords, range, sd = 1) {
  K <- matern32_dense(coords, range) + diag(1e-6, nrow(coords))
  as.numeric(t(chol(K)) %*% stats::rnorm(nrow(coords))) * sd
}

make_coords <- function(n_side) {
  as.matrix(expand.grid(x = seq(0, 1, length.out = n_side),
                        y = seq(0, 1, length.out = n_side)))
}

simulate_calibration_null <- function(scenario, rep_id, cfg) {
  set.seed(cfg$seed + 100000L * match(scenario, calibration_scenarios()) + rep_id)
  coords <- make_coords(cfg$sim_n_side)
  n <- nrow(coords)
  g <- cfg$sim_g
  q <- cfg$sim_pathway_size
  Y <- matrix(stats::rnorm(n * g), n, g)
  colnames(Y) <- paste0("g", seq_len(g))
  X <- cbind(tech = as.numeric(scale(stats::rnorm(n))))

  if (scenario == "null_correlated_nonspatial") {
    latent <- matrix(stats::rnorm(n * 4), n, 4)
    loadings <- matrix(stats::rnorm(g * 4, sd = 0.35), g, 4)
    Y <- Y + latent %*% t(loadings)
  } else if (scenario == "null_spatial_background") {
    u1 <- simulate_gp(coords, range = 0.30, sd = 1)
    u2 <- simulate_gp(coords, range = 0.16, sd = 1)
    for (jj in (q + 1L):g) {
      Y[, jj] <- Y[, jj] + stats::rnorm(1, 0.6, 0.15) * u1 +
        stats::rnorm(1, 0.3, 0.10) * u2
    }
  } else if (scenario == "null_covariate_spatial") {
    z <- simulate_gp(coords, range = 0.35, sd = 1)
    X <- cbind(X, spatial_covariate = as.numeric(scale(z)))
    for (jj in seq_len(q)) {
      Y[, jj] <- Y[, jj] + stats::rnorm(1, 0.9, 0.1) * X[, "spatial_covariate"]
    }
  } else if (scenario != "null_independent") {
    stop("unknown scenario: ", scenario)
  }
  list(
    Y = Y,
    coords = coords,
    X = X,
    pathways = list(null_pathway = paste0("g", seq_len(q))),
    scenario = scenario,
    rep = rep_id
  )
}

calibration_scenarios <- function() {
  c("null_independent", "null_correlated_nonspatial",
    "null_spatial_background", "null_covariate_spatial")
}

run_passage_mode <- function(dat, mode, n_perm, cfg) {
  range_grid <- passage_default_range_grid(dat$coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
  fit_args <- list(
    K = cfg$K,
    m = cfg$m,
    range_grid = range_grid,
    kernel = "matern32",
    ordering = "none",
    verbose = FALSE
  )
  t0 <- proc.time()[["elapsed"]]
  fit <- if (mode == "none") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = n_perm,
      seed = cfg$seed + dat$rep, verbose = FALSE
    )
  } else if (mode == "spot_crossfit") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = n_perm,
      seed = cfg$seed + dat$rep, verbose = FALSE, anti_dip = "spot_crossfit",
      anti_dip_args = list(n_folds = cfg$n_folds, seed = cfg$seed + dat$rep, fit_args = fit_args)
    )
  } else if (mode == "pathway_holdout") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "permutation", n_perm = n_perm,
      seed = cfg$seed + dat$rep, verbose = FALSE, anti_dip = "pathway_holdout",
      anti_dip_args = list(fit_args = fit_args)
    )
  } else if (mode == "refit_null") {
    passage_run(
      Y = dat$Y, coords = dat$coords, X = dat$X, pathways = dat$pathways,
      K = cfg$K, m = cfg$m, range_grid = range_grid, kernel = "matern32",
      hypotheses = "H1", calibration = "moment", n_perm = n_perm,
      seed = cfg$seed + dat$rep, verbose = FALSE, anti_dip = "refit_null",
      anti_dip_args = list(fit_args = fit_args, n_perm = n_perm)
    )
  } else {
    stop("unknown mode: ", mode)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  row <- fit$summary[1L, , drop = FALSE]
  data.frame(
    scenario = dat$scenario,
    replicate = dat$rep,
    mode = mode,
    pathway = row$pathway,
    p = row$p_H1,
    p_moment = row$p_H1_moment,
    R2_cca = row$R2_cca,
    PSVS_range = row$PSVS_range,
    mean_propSV = row$mean_propSV,
    elapsed_sec = elapsed,
    n_perm = n_perm,
    stringsAsFactors = FALSE
  )
}

summarize_p <- function(x, group_cols) {
  alphas <- c(0.10, 0.05, 0.01)
  groups <- unique(x[group_cols])
  rows <- list()
  ii <- 0L
  for (gg in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(x))
    for (cc in group_cols) keep <- keep & x[[cc]] == groups[[cc]][gg]
    z <- x[keep, , drop = FALSE]
    for (alpha in alphas) {
      p <- z$p[is.finite(z$p)]
      rate <- mean(p <= alpha)
      n <- length(p)
      se <- sqrt(rate * (1 - rate) / max(1L, n))
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        groups[gg, , drop = FALSE],
        alpha = alpha,
        n = n,
        reject_rate = rate,
        mc_se = se,
        ci95_low = max(0, rate - 1.96 * se),
        ci95_high = min(1, rate + 1.96 * se),
        median_p = stats::median(p),
        min_p = min(p),
        median_elapsed_sec = stats::median(z$elapsed_sec),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

run_simulation_panel <- function(cfg) {
  modes <- c("none", "spot_crossfit", "pathway_holdout")
  rows <- list()
  ii <- 0L
  for (scenario in calibration_scenarios()) {
    for (rr in seq_len(cfg$sim_reps)) {
      dat <- simulate_calibration_null(scenario, rr, cfg)
      for (mode in modes) {
        message("sim scenario=", scenario, " rep=", rr, " mode=", mode)
        ii <- ii + 1L
        rows[[ii]] <- run_passage_mode(dat, mode, cfg$sim_n_perm, cfg)
      }
    }
    for (rr in seq_len(cfg$sim_refit_reps)) {
      dat <- simulate_calibration_null(scenario, rr + 10000L, cfg)
      message("sim scenario=", scenario, " rep=", rr, " mode=refit_null")
      ii <- ii + 1L
      rows[[ii]] <- run_passage_mode(dat, "refit_null", cfg$sim_refit_perm, cfg)
    }
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(cfg$out_dir, "simulation_null_pvalues.csv"), row.names = FALSE)
  summary <- summarize_p(out, c("scenario", "mode"))
  write.csv(summary, file.path(cfg$out_dir, "simulation_null_summary.csv"), row.names = FALSE)
  summary
}

prepared_paths <- function() {
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_perm999_fastcal/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_prepared_data.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_perm999_fastcal/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_prepared_data.rds"
  )
}

make_real_random_pathways <- function(dat, n_random, seed) {
  set.seed(seed)
  sizes <- as.integer(stats::quantile(lengths(dat$pathways), probs = seq(0.1, 0.9, length.out = n_random), names = FALSE))
  sizes <- pmax(5L, pmin(sizes, ncol(dat$Y)))
  genes <- colnames(dat$Y)
  out <- vector("list", length(sizes))
  for (ii in seq_along(sizes)) {
    out[[ii]] <- sample(genes, sizes[[ii]], replace = FALSE)
  }
  names(out) <- paste0("random_", seq_along(out), "_size", sizes)
  out
}

run_real_random_panel <- function(cfg) {
  rows <- list()
  ii <- 0L
  for (dataset in names(prepared_paths())) {
    message("real dataset=", dataset, ": loading prepared data")
    dat <- readRDS(prepared_paths()[[dataset]])
    dat$pathways <- make_real_random_pathways(dat, cfg$real_random, cfg$seed + match(dataset, names(prepared_paths())))
    dat$scenario <- "real_matched_random"
    dat$rep <- 1L
    for (mode in cfg$real_modes) {
      message("real dataset=", dataset, " mode=", mode)
      range_grid <- passage_default_range_grid(dat$coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
      fit_args <- list(K = cfg$K, m = cfg$m, range_grid = range_grid,
                       kernel = "matern32", ordering = "none", verbose = FALSE)
      t0 <- proc.time()[["elapsed"]]
      fit <- if (mode == "none") {
        passage_run(dat$Y, dat$coords, dat$pathways, X = dat$X,
                    K = cfg$K, m = cfg$m, range_grid = range_grid,
                    hypotheses = "H1", calibration = "permutation",
                    n_perm = cfg$real_n_perm, seed = cfg$seed, verbose = FALSE)
      } else if (mode == "spot_crossfit") {
        passage_run(dat$Y, dat$coords, dat$pathways, X = dat$X,
                    K = cfg$K, m = cfg$m, range_grid = range_grid,
                    hypotheses = "H1", calibration = "permutation",
                    n_perm = cfg$real_n_perm, seed = cfg$seed, verbose = FALSE,
                    anti_dip = "spot_crossfit",
                    anti_dip_args = list(n_folds = cfg$n_folds, seed = cfg$seed, fit_args = fit_args))
      } else if (mode == "pathway_holdout") {
        passage_run(dat$Y, dat$coords, dat$pathways, X = dat$X,
                    K = cfg$K, m = cfg$m, range_grid = range_grid,
                    hypotheses = "H1", calibration = "permutation",
                    n_perm = cfg$real_n_perm, seed = cfg$seed, verbose = FALSE,
                    anti_dip = "pathway_holdout",
                    anti_dip_args = list(fit_args = fit_args))
      } else {
        stop("unsupported real-data mode: ", mode)
      }
      elapsed <- proc.time()[["elapsed"]] - t0
      tbl <- fit$summary
      tbl$dataset <- dataset
      tbl$mode <- mode
      tbl$elapsed_sec <- elapsed
      tbl$n_perm <- cfg$real_n_perm
      ii <- ii + 1L
      rows[[ii]] <- tbl[, c("dataset", "mode", "pathway", "pathway_size", "p_H1",
                            "p_H1_moment", "R2_cca", "PSVS_range", "mean_propSV",
                            "elapsed_sec", "n_perm")]
    }
  }
  out <- do.call(rbind, rows)
  names(out)[names(out) == "p_H1"] <- "p"
  write.csv(out, file.path(cfg$out_dir, "real_random_gene_set_pvalues.csv"), row.names = FALSE)
  summary <- summarize_p(out, c("dataset", "mode"))
  write.csv(summary, file.path(cfg$out_dir, "real_random_gene_set_summary.csv"), row.names = FALSE)
  summary
}

write_report <- function(cfg, sim_summary, real_summary) {
  md <- c(
    "# PASSAGE Anti-Dip Calibration",
    "",
    "## Simulation Null Summary",
    "",
    paste(colnames(sim_summary), collapse = " | "),
    paste(rep("---", ncol(sim_summary)), collapse = " | ")
  )
  sim_rows <- apply(sim_summary, 1L, function(x) paste(x, collapse = " | "))
  md <- c(md, sim_rows, "", "## Real Breast Random Gene-Set Summary", "",
          paste(colnames(real_summary), collapse = " | "),
          paste(rep("---", ncol(real_summary)), collapse = " | "))
  real_rows <- apply(real_summary, 1L, function(x) paste(x, collapse = " | "))
  writeLines(c(md, real_rows), file.path(cfg$out_dir, "summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
sim_summary <- run_simulation_panel(cfg)
real_summary <- run_real_random_panel(cfg)
write_report(cfg, sim_summary, real_summary)
message("Wrote calibration outputs to ", cfg$out_dir)
print(sim_summary, row.names = FALSE)
print(real_summary, row.names = FALSE)
