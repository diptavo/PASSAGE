# PASSAGE MVP validation suite.
#
# This script checks whether the current efficient MVP behaves sensibly before
# adding heavier methods-note components. It covers:
#   1. H1 null calibration under independent and correlated non-spatial genes.
#   2. Sparse and diffuse pathway-signal power.
#   3. H1/H2 behavior under a spatial cell-type-like confounder.
#   4. Simple baselines: per-gene ACAT and the original spapath_test prototype.
#   5. Runtime scaling over small N grids.
#
# Usage:
#   cd /path/to/PASSAGE
#   Rscript scripts/sim_passage_validation.R --n-reps=10
#   Rscript scripts/sim_passage_validation.R --n-reps=50 --out-dir=results/passage_validation_50

load_spapath_sources <- function() {
  for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    source(f)
  }
}

parse_args <- function(args) {
  out <- list(
    n_reps = 10L,
    seed = 20260523L,
    out_dir = file.path("results", "passage_validation"),
    n_side = 9L,
    g = 60L,
    k_engine = 3L,
    m = 8L,
    n_perm = 99L,
    run_baseline = TRUE,
    baseline_n_sim = 200L,
    runtime = TRUE,
    scenarios = c(
      "null_independent",
      "null_correlated",
      "sparse_signal",
      "diffuse_signal",
      "celltype_confounded",
      "celltype_plus_specific"
    )
  )
  for (arg in args) {
    if (grepl("^--n-reps=", arg)) out$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--seed=", arg)) out$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--out-dir=", arg)) out$out_dir <- sub("^--out-dir=", "", arg)
    if (grepl("^--n-side=", arg)) out$n_side <- as.integer(sub("^--n-side=", "", arg))
    if (grepl("^--g=", arg)) out$g <- as.integer(sub("^--g=", "", arg))
    if (grepl("^--k-engine=", arg)) out$k_engine <- as.integer(sub("^--k-engine=", "", arg))
    if (grepl("^--m=", arg)) out$m <- as.integer(sub("^--m=", "", arg))
    if (grepl("^--n-perm=", arg)) out$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--baseline-n-sim=", arg)) out$baseline_n_sim <- as.integer(sub("^--baseline-n-sim=", "", arg))
    if (grepl("^--scenarios=", arg)) {
      out$scenarios <- strsplit(sub("^--scenarios=", "", arg), ",", fixed = TRUE)[[1]]
    }
    if (arg == "--no-baseline") out$run_baseline <- FALSE
    if (arg == "--no-runtime") out$runtime <- FALSE
  }
  out
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
  as.matrix(expand.grid(
    x = seq(0, 1, length.out = n_side),
    y = seq(0, 1, length.out = n_side)
  ))
}

scale01 <- function(x) {
  x <- as.numeric(x)
  as.numeric(scale(x))
}

simulate_dataset <- function(scenario,
                             rep_id,
                             n_side = 9L,
                             g = 60L,
                             pathway_size = 12L) {
  coords <- make_coords(n_side)
  n <- nrow(coords)
  genes <- paste0("g", seq_len(g))
  Y <- matrix(stats::rnorm(n * g), nrow = n, ncol = g)
  colnames(Y) <- genes

  signal_genes <- paste0("g", seq_len(pathway_size))
  null_genes <- paste0("g", 31:(30 + pathway_size))
  pathways <- list(signal = signal_genes, null = null_genes)
  Z_CT <- NULL
  true_driver_genes <- character()
  signal_type <- "none"

  if (scenario == "null_independent") {
    signal_type <- "null"
  } else if (scenario == "null_correlated") {
    signal_type <- "null_correlated"
    latent <- matrix(stats::rnorm(n * 3), nrow = n)
    loadings <- matrix(stats::rnorm(g * 3, sd = 0.35), nrow = g)
    Y <- Y + latent %*% t(loadings)
    colnames(Y) <- genes
  } else if (scenario == "sparse_signal") {
    signal_type <- "sparse"
    v <- simulate_gp(coords, range = 0.25, sd = 1)
    true_driver_genes <- paste0("g", 1:4)
    for (gene in true_driver_genes) {
      Y[, gene] <- Y[, gene] + 1.25 * v + stats::rnorm(n, sd = 0.10)
    }
  } else if (scenario == "diffuse_signal") {
    signal_type <- "diffuse"
    true_driver_genes <- paste0("g", 1:pathway_size)
    for (gene in true_driver_genes) {
      Y[, gene] <- Y[, gene] + 0.45 * simulate_gp(coords, range = 0.30, sd = 1)
    }
  } else if (scenario == "celltype_confounded") {
    signal_type <- "celltype_confounded"
    z <- simulate_gp(coords, range = 0.35, sd = 1)
    Z_CT <- cbind(celltype_gradient = scale01(z))
    true_driver_genes <- paste0("g", 1:8)
    for (gene in true_driver_genes) {
      Y[, gene] <- Y[, gene] + 1.0 * Z_CT[, 1]
    }
  } else if (scenario == "celltype_plus_specific") {
    signal_type <- "celltype_plus_specific"
    z <- simulate_gp(coords, range = 0.35, sd = 1)
    v <- simulate_gp(coords, range = 0.18, sd = 1)
    Z_CT <- cbind(celltype_gradient = scale01(z))
    true_driver_genes <- paste0("g", 1:8)
    for (gene in true_driver_genes) {
      Y[, gene] <- Y[, gene] + 0.8 * Z_CT[, 1] + 0.65 * v
    }
  } else {
    stop("unknown scenario: ", scenario)
  }

  list(
    scenario = scenario,
    signal_type = signal_type,
    rep_id = rep_id,
    Y = Y,
    coords = coords,
    pathways = pathways,
    Z_CT = Z_CT,
    true_driver_genes = true_driver_genes
  )
}

flatten_genes <- function(x) {
  if (length(x) == 0L) {
    return("")
  }
  paste(as.character(x), collapse = ";")
}

driver_recall <- function(selected, truth) {
  selected <- unique(as.character(selected))
  truth <- unique(as.character(truth))
  if (length(truth) == 0L) {
    return(NA_real_)
  }
  length(intersect(selected, truth)) / length(truth)
}

extract_passage_row <- function(fit, dat, pathway = "signal", elapsed_sec = NA_real_) {
  row <- fit$summary[fit$summary$pathway == pathway, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("expected exactly one row for pathway ", pathway)
  }
  selected <- row$spasset_genes_H1[[1]]
  data.frame(
    method = "passage",
    scenario = dat$scenario,
    signal_type = dat$signal_type,
    replicate = dat$rep_id,
    pathway = pathway,
    N = nrow(dat$Y),
    G = ncol(dat$Y),
    pathway_size = row$pathway_size,
    p_H1 = row$p_H1,
    p_H2 = if ("p_H2" %in% names(row)) row$p_H2 else NA_real_,
    p_H1_moment = if ("p_H1_moment" %in% names(row)) row$p_H1_moment else NA_real_,
    p_H2_moment = if ("p_H2_moment" %in% names(row)) row$p_H2_moment else NA_real_,
    p_value = row$p_H1,
    R2_cca = row$R2_cca,
    PSVS_range = row$PSVS_range,
    mean_propSV = row$mean_propSV,
    cell_type_share = if ("cell_type_share" %in% names(row)) row$cell_type_share else NA_real_,
    driver_recall = driver_recall(selected, dat$true_driver_genes),
    selected_drivers = flatten_genes(selected),
    elapsed_sec = elapsed_sec,
    stringsAsFactors = FALSE
  )
}

gene_acat_baseline <- function(dat, ranges, m) {
  scores <- spapath_feature_scores(
    Y = dat$Y,
    coords = dat$coords,
    ranges = ranges,
    m = m,
    kernel = "matern32"
  )
  fs <- scores$feature_scores
  p_gene <- stats::pnorm(fs$z, lower.tail = FALSE)
  names(p_gene) <- fs$feature
  rows <- lapply(names(dat$pathways), function(pw) {
    genes <- intersect(dat$pathways[[pw]], names(p_gene))
    p <- passage_acat(p_gene[genes])
    data.frame(
      method = "gene_acat",
      scenario = dat$scenario,
      signal_type = dat$signal_type,
      replicate = dat$rep_id,
      pathway = pw,
      N = nrow(dat$Y),
      G = ncol(dat$Y),
      pathway_size = length(genes),
      p_H1 = NA_real_,
      p_H2 = NA_real_,
      p_H1_moment = NA_real_,
      p_H2_moment = NA_real_,
      p_value = p,
      R2_cca = NA_real_,
      PSVS_range = NA_real_,
      mean_propSV = NA_real_,
      cell_type_share = NA_real_,
      driver_recall = NA_real_,
      selected_drivers = "",
      elapsed_sec = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

spapath_baseline <- function(dat, ranges, m, n_sim) {
  fit <- spapath_test(
    Y = dat$Y,
    coords = dat$coords,
    pathways = dat$pathways,
    ranges = ranges,
    m = m,
    n_sim = n_sim,
    seed = dat$rep_id + 1000L,
    verbose = FALSE
  )
  rows <- lapply(seq_len(nrow(fit$results)), function(i) {
    row <- fit$results[i, , drop = FALSE]
    data.frame(
      method = "spapath",
      scenario = dat$scenario,
      signal_type = dat$signal_type,
      replicate = dat$rep_id,
      pathway = row$pathway,
      N = nrow(dat$Y),
      G = ncol(dat$Y),
      pathway_size = row$q,
      p_H1 = NA_real_,
      p_H2 = NA_real_,
      p_H1_moment = NA_real_,
      p_H2_moment = NA_real_,
      p_value = row$p_value,
      R2_cca = NA_real_,
      PSVS_range = NA_real_,
      mean_propSV = row$eSPVE_any,
      cell_type_share = NA_real_,
      driver_recall = driver_recall(row$driver_genes[[1]], dat$true_driver_genes),
      selected_drivers = flatten_genes(row$driver_genes[[1]]),
      elapsed_sec = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_one <- function(dat, cfg) {
  ranges <- c(0.10, 0.20, 0.35)
  hypotheses <- if (is.null(dat$Z_CT)) "H1" else c("H1", "H2")
  t0 <- proc.time()[["elapsed"]]
  fit <- passage_run(
    Y = dat$Y,
    coords = dat$coords,
    pathways = dat$pathways,
    Z_CT = dat$Z_CT,
    K = cfg$k_engine,
    m = cfg$m,
    range_grid = ranges,
    hypotheses = hypotheses,
    calibration = "permutation",
    n_perm = cfg$n_perm,
    seed = cfg$seed + 10000L * dat$rep_id + match(dat$scenario, c(
      "null_independent",
      "null_correlated",
      "sparse_signal",
      "diffuse_signal",
      "celltype_confounded",
      "celltype_plus_specific"
    )),
    verbose = FALSE
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  passage_rows <- rbind(
    extract_passage_row(fit, dat, pathway = "signal", elapsed_sec = elapsed),
    extract_passage_row(fit, dat, pathway = "null", elapsed_sec = elapsed)
  )

  if (!cfg$run_baseline) {
    return(passage_rows)
  }
  base1 <- gene_acat_baseline(dat, ranges = ranges, m = cfg$m)
  base2 <- spapath_baseline(dat, ranges = ranges, m = cfg$m, n_sim = cfg$baseline_n_sim)
  rbind(passage_rows, base1, base2)
}

summarize_validation <- function(tbl) {
  alpha <- c(0.05, 0.01)
  rows <- list()
  ii <- 0L
  for (method in unique(tbl$method)) {
    for (scenario in unique(tbl$scenario)) {
      for (pathway in unique(tbl$pathway)) {
        sub <- tbl[tbl$method == method & tbl$scenario == scenario & tbl$pathway == pathway, , drop = FALSE]
        if (nrow(sub) == 0L) next
        for (a in alpha) {
          ii <- ii + 1L
          rows[[ii]] <- data.frame(
            method = method,
            scenario = scenario,
            pathway = pathway,
            alpha = a,
            n = sum(is.finite(sub$p_value)),
            reject_rate = mean(sub$p_value <= a, na.rm = TRUE),
            median_p = stats::median(sub$p_value, na.rm = TRUE),
            median_R2_cca = stats::median(sub$R2_cca, na.rm = TRUE),
            median_PSVS_range = stats::median(sub$PSVS_range, na.rm = TRUE),
            mean_driver_recall = mean(sub$driver_recall, na.rm = TRUE),
            median_elapsed_sec = stats::median(sub$elapsed_sec, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  ans <- do.call(rbind, rows)
  ans[order(ans$method, ans$scenario, ans$pathway, ans$alpha), , drop = FALSE]
}

run_runtime_profile <- function(cfg) {
  n_sides <- c(7L, 9L, 11L)
  rows <- vector("list", length(n_sides))
  for (i in seq_along(n_sides)) {
    dat <- simulate_dataset("sparse_signal", rep_id = i, n_side = n_sides[i], g = cfg$g)
    t0 <- proc.time()[["elapsed"]]
    fit <- passage_run(
      Y = dat$Y,
      coords = dat$coords,
      pathways = dat$pathways,
      K = cfg$k_engine,
      m = min(cfg$m, nrow(dat$Y) - 1L),
      range_grid = c(0.10, 0.20, 0.35),
      hypotheses = "H1",
      calibration = "permutation",
      n_perm = cfg$n_perm,
      seed = cfg$seed + 50000L + i,
      verbose = FALSE
    )
    elapsed <- proc.time()[["elapsed"]] - t0
    rows[[i]] <- data.frame(
      n_side = n_sides[i],
      N = nrow(dat$Y),
      G = ncol(dat$Y),
      K = cfg$k_engine,
      m = min(cfg$m, nrow(dat$Y) - 1L),
      n_pathways = length(dat$pathways),
      elapsed_sec = elapsed,
      signal_p_H1 = fit$summary$p_H1[fit$summary$pathway == "signal"],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  load_spapath_sources()
  set.seed(cfg$seed)

  scenarios <- cfg$scenarios
  allowed <- c(
    "null_independent",
    "null_correlated",
    "sparse_signal",
    "diffuse_signal",
    "celltype_confounded",
    "celltype_plus_specific"
  )
  bad <- setdiff(scenarios, allowed)
  if (length(bad) > 0L) {
    stop("unknown scenario(s): ", paste(bad, collapse = ", "))
  }

  message("Running PASSAGE validation: ", cfg$n_reps, " replicate(s) per scenario")
  rows <- list()
  rr <- 0L
  for (scenario in scenarios) {
    message("  scenario: ", scenario)
    for (rep_id in seq_len(cfg$n_reps)) {
      dat <- simulate_dataset(
        scenario = scenario,
        rep_id = rep_id,
        n_side = cfg$n_side,
        g = cfg$g
      )
      rr <- rr + 1L
      rows[[rr]] <- run_one(dat, cfg)
    }
  }
  results <- do.call(rbind, rows)
  summary <- summarize_validation(results)

  result_file <- file.path(cfg$out_dir, "passage_validation_results.csv")
  summary_file <- file.path(cfg$out_dir, "passage_validation_summary.csv")
  write.csv(results, result_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)

  if (cfg$runtime) {
    message("  runtime profile")
    runtime <- run_runtime_profile(cfg)
    write.csv(runtime, file.path(cfg$out_dir, "passage_runtime_profile.csv"), row.names = FALSE)
  }

  message("Wrote: ", result_file)
  message("Wrote: ", summary_file)
  print(summary[summary$method == "passage" & summary$pathway == "signal", ], row.names = FALSE)
}

main()
