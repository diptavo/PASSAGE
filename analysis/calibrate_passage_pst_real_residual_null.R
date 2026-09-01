# Real-data-derived residual-null calibration for Pathway Spatial Transferability.
#
# Usage:
#   Rscript scripts/calibrate_passage_pst_real_residual_null.R \
#     --n-reps=5 --n-perm=99 --max-pathways=4 --pathway-source=hallmark \
#     --out-dir=results/passage_pst_real_residual_null_quick

parse_args <- function(args) {
  cfg <- list(
    seed = 20260730L,
    n_reps = 5L,
    n_perm = 99L,
    max_pathways = 4L,
    n_gene_splits = 3L,
    n_folds = 4L,
    k_neighbors = 8L,
    components = c("pc1"),
    pathway_source = "hallmark",
    out_dir = file.path("results", "passage_pst_real_residual_null_quick")
  )
  for (arg in args) {
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--n-reps=", arg)) cfg$n_reps <- as.integer(sub("^--n-reps=", "", arg))
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--n-gene-splits=", arg)) cfg$n_gene_splits <- as.integer(sub("^--n-gene-splits=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
    if (grepl("^--k-neighbors=", arg)) cfg$k_neighbors <- as.integer(sub("^--k-neighbors=", "", arg))
    if (grepl("^--components=", arg)) cfg$components <- strsplit(sub("^--components=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--pathway-source=", arg)) cfg$pathway_source <- sub("^--pathway-source=", "", arg)
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
  }
  if (!cfg$pathway_source %in% c("hallmark", "random")) {
    stop("--pathway-source must be hallmark or random")
  }
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

prepared_paths <- function() {
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_perm999_fastcal/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_prepared_data.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_perm999_fastcal/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_prepared_data.rds"
  )
}

select_pathways <- function(pathways, max_pathways) {
  preferred <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA"
  )
  keep <- c(intersect(preferred, names(pathways)), setdiff(names(pathways), preferred))
  pathways[keep[seq_len(min(max_pathways, length(keep)))]]
}

make_random_pathways_like <- function(template_pathways, gene_names, seed) {
  set.seed(seed)
  sizes <- pmax(4L, pmin(lengths(template_pathways), length(gene_names)))
  out <- vector("list", length(sizes))
  for (ii in seq_along(sizes)) {
    out[[ii]] <- sample(gene_names, sizes[[ii]], replace = FALSE)
  }
  names(out) <- paste0("random_like_", seq_along(out), "_size", sizes)
  out
}

make_residual_null_y <- function(Y, X, seed) {
  set.seed(seed)
  Xd <- passage_prepare_design(X, nrow(Y), intercept = TRUE)
  qrx <- qr(Xd)
  R <- passage_residualize_with_qr(Y, qrx)
  Yhat <- Y - R
  Yhat + R[sample.int(nrow(R)), , drop = FALSE]
}

summarize_p <- function(x) {
  groups <- unique(x[c("dataset", "pathway_source", "component")])
  rows <- list()
  ii <- 0L
  for (gg in seq_len(nrow(groups))) {
    keep <- rep(TRUE, nrow(x))
    for (cc in colnames(groups)) keep <- keep & x[[cc]] == groups[[cc]][gg]
    z <- x[keep, , drop = FALSE]
    for (alpha in c(0.10, 0.05, 0.01)) {
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
        median_statistic = stats::median(z$statistic),
        median_elapsed_sec = stats::median(z$elapsed_sec),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

write_report <- function(cfg, summary) {
  md <- c(
    "# PASSAGE PST Real Residual Null Calibration",
    "",
    paste0("- Reps per dataset: ", cfg$n_reps),
    paste0("- Pathways per replicate: ", cfg$max_pathways),
    paste0("- PST permutations per pathway: ", cfg$n_perm),
    paste0("- Gene splits: ", cfg$n_gene_splits),
    paste0("- Spot folds: ", cfg$n_folds),
    paste0("- Neighbors: ", cfg$k_neighbors),
    paste0("- Components: ", paste(cfg$components, collapse = ", ")),
    paste0("- Pathway source: ", cfg$pathway_source),
    "",
    paste(colnames(summary), collapse = " | "),
    paste(rep("---", ncol(summary)), collapse = " | "),
    apply(summary, 1L, function(x) paste(x, collapse = " | "))
  )
  writeLines(md, file.path(cfg$out_dir, "summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
ii <- 0L
for (dataset in names(prepared_paths())) {
  message("loading ", dataset)
  dat <- readRDS(prepared_paths()[[dataset]])
  gene_names <- colnames(dat$Y)
  template <- select_pathways(dat$pathways, cfg$max_pathways)
  pathways <- if (cfg$pathway_source == "random") {
    make_random_pathways_like(
      template, gene_names,
      cfg$seed + 10000L * match(dataset, names(prepared_paths()))
    )
  } else {
    template
  }
  for (rr in seq_len(cfg$n_reps)) {
    message("dataset=", dataset, " PST residual-null rep=", rr)
    Y0 <- make_residual_null_y(dat$Y, dat$X, cfg$seed + 1000L * match(dataset, names(prepared_paths())) + rr)
    for (component in cfg$components) {
      for (pp in seq_along(pathways)) {
        message("  component=", component, " pathway=", names(pathways)[pp])
        t0 <- proc.time()[["elapsed"]]
        ans <- passage_pathway_spatial_transferability(
          Y0, dat$coords, pathways[[pp]], X = dat$X, gene_names = gene_names,
          n_gene_splits = cfg$n_gene_splits,
          n_folds = cfg$n_folds,
          k_neighbors = cfg$k_neighbors,
          component = component,
          n_perm = cfg$n_perm,
          seed = cfg$seed + 100000L * rr + 1000L * match(component, cfg$components) + pp
        )
        elapsed <- proc.time()[["elapsed"]] - t0
        ii <- ii + 1L
        rows[[ii]] <- data.frame(
          dataset = dataset,
          pathway_source = cfg$pathway_source,
          replicate = rr,
          component = component,
          pathway = names(pathways)[pp],
          pathway_size = ans$pathway_size,
          statistic = ans$statistic,
          p = ans$p,
          n_pairs = ans$n_pairs,
          n_perm = cfg$n_perm,
          elapsed_sec = elapsed,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

out <- do.call(rbind, rows)
out$fdr <- stats::p.adjust(out$p, method = "BH")
write.csv(out, file.path(cfg$out_dir, "pst_real_residual_null_pvalues.csv"), row.names = FALSE)
summary <- summarize_p(out)
write.csv(summary, file.path(cfg$out_dir, "pst_real_residual_null_summary.csv"), row.names = FALSE)
write_report(cfg, summary)
message("Wrote outputs to ", cfg$out_dir)
print(summary, row.names = FALSE)
