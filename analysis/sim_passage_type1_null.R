# Parallel PASSAGE type-I error simulation for null pathways.
#
# Usage:
#   cd /path/to/PASSAGE
#   Rscript scripts/sim_passage_type1_null.R \
#     --n-reps=1000 \
#     --n-perm=999 \
#     --scenarios=null_independent,null_correlated \
#     --cores=6 \
#     --out-dir=results/passage_type1_null1000_perm999

load_sources <- function() {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
}

parse_args <- function(args) {
  cfg <- list(
    n_reps = 1000L,
    n_perm = 999L,
    seed = 20260523L,
    n_side = 9L,
    g = 60L,
    pathway_size = 12L,
    k_engine = 3L,
    m = 8L,
    cores = max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L)),
    scenarios = c("null_independent", "null_correlated"),
    out_dir = file.path("results", "passage_type1_null1000_perm999")
  )
  for (arg in args) {
    if (grepl("^--n-reps=", arg)) cfg$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--n-side=", arg)) cfg$n_side <- as.integer(sub("^--n-side=", "", arg))
    if (grepl("^--g=", arg)) cfg$g <- as.integer(sub("^--g=", "", arg))
    if (grepl("^--pathway-size=", arg)) cfg$pathway_size <- as.integer(sub("^--pathway-size=", "", arg))
    if (grepl("^--k-engine=", arg)) cfg$k_engine <- as.integer(sub("^--k-engine=", "", arg))
    if (grepl("^--m=", arg)) cfg$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--cores=", arg)) cfg$cores <- as.integer(sub("^--cores=", "", arg))
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--scenarios=", arg)) {
      cfg$scenarios <- strsplit(sub("^--scenarios=", "", arg), ",", fixed = TRUE)[[1]]
    }
  }
  cfg
}

make_coords <- function(n_side) {
  as.matrix(expand.grid(
    x = seq(0, 1, length.out = n_side),
    y = seq(0, 1, length.out = n_side)
  ))
}

simulate_null <- function(scenario, cfg, rep_id) {
  set.seed(cfg$seed + 100000L * match(scenario, cfg$scenarios) + rep_id)
  coords <- make_coords(cfg$n_side)
  n <- nrow(coords)
  g <- cfg$g
  Y <- matrix(stats::rnorm(n * g), nrow = n, ncol = g)
  colnames(Y) <- paste0("g", seq_len(g))

  if (scenario == "null_correlated") {
    latent <- matrix(stats::rnorm(n * 3), nrow = n)
    loadings <- matrix(stats::rnorm(g * 3, sd = 0.35), nrow = g)
    Y <- Y + latent %*% t(loadings)
    colnames(Y) <- paste0("g", seq_len(g))
  } else if (scenario != "null_independent") {
    stop("unsupported null scenario: ", scenario)
  }

  pathways <- list(null_pathway = paste0("g", seq_len(cfg$pathway_size)))
  t0 <- proc.time()[["elapsed"]]
  fit <- passage_run(
    Y = Y,
    coords = coords,
    pathways = pathways,
    K = cfg$k_engine,
    m = cfg$m,
    range_grid = c(0.10, 0.20, 0.35),
    hypotheses = "H1",
    calibration = "permutation",
    n_perm = cfg$n_perm,
    seed = cfg$seed + 500000L * match(scenario, cfg$scenarios) + rep_id,
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  row <- fit$summary[1L, , drop = FALSE]
  data.frame(
    scenario = scenario,
    replicate = rep_id,
    N = nrow(Y),
    G = ncol(Y),
    pathway_size = cfg$pathway_size,
    n_perm = cfg$n_perm,
    p_H1 = row$p_H1,
    p_H1_moment = row$p_H1_moment,
    R2_cca = row$R2_cca,
    PSVS_range = row$PSVS_range,
    elapsed_sec = elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_type1 <- function(results) {
  alphas <- c(0.10, 0.05, 0.01, 0.001)
  rows <- list()
  ii <- 0L
  for (scenario in unique(results$scenario)) {
    x <- results[results$scenario == scenario, , drop = FALSE]
    for (alpha in alphas) {
      ii <- ii + 1L
      p <- x$p_H1
      reject <- p <= alpha
      rate <- mean(reject, na.rm = TRUE)
      n <- sum(is.finite(p))
      se <- sqrt(rate * (1 - rate) / max(n, 1L))
      rows[[ii]] <- data.frame(
        scenario = scenario,
        alpha = alpha,
        n = n,
        type1 = rate,
        mc_se = se,
        ci95_low = max(0, rate - 1.96 * se),
        ci95_high = min(1, rate + 1.96 * se),
        median_p = stats::median(p, na.rm = TRUE),
        ks_p_uniform = stats::ks.test(p, "punif")$p.value,
        median_elapsed_sec = stats::median(x$elapsed_sec, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  load_sources()
  allowed <- c("null_independent", "null_correlated")
  bad <- setdiff(cfg$scenarios, allowed)
  if (length(bad) > 0L) {
    stop("Only null scenarios are supported here: ", paste(bad, collapse = ", "))
  }

  message("PASSAGE type-I null simulation")
  message("  scenarios: ", paste(cfg$scenarios, collapse = ", "))
  message("  reps per scenario: ", cfg$n_reps)
  message("  n_perm: ", cfg$n_perm)
  message("  cores: ", cfg$cores)

  all_rows <- list()
  for (scenario in cfg$scenarios) {
    message("  running ", scenario)
    reps <- seq_len(cfg$n_reps)
    rows <- parallel::mclapply(
      reps,
      function(rep_id) simulate_null(scenario, cfg, rep_id),
      mc.cores = cfg$cores,
      mc.preschedule = FALSE
    )
    scenario_tbl <- do.call(rbind, rows)
    all_rows[[scenario]] <- scenario_tbl
    write.csv(
      scenario_tbl,
      file.path(cfg$out_dir, paste0("type1_", scenario, "_results.csv")),
      row.names = FALSE
    )
  }

  results <- do.call(rbind, all_rows)
  summary <- summarize_type1(results)
  write.csv(results, file.path(cfg$out_dir, "type1_results.csv"), row.names = FALSE)
  write.csv(summary, file.path(cfg$out_dir, "type1_summary.csv"), row.names = FALSE)
  print(summary, row.names = FALSE)
}

main()
