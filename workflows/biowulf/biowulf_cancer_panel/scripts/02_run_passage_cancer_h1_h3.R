#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_perm <- if (length(args) >= 3L) as.integer(args[[3L]]) else 199L
cores <- if (length(args) >= 4L) as.integer(args[[4L]]) else 4L

set.seed(20260803L + task_id)
setwd(root)

source(file.path(root, "scripts", "load_passage.R"))

manifest_path <- file.path(root, "data", "passage_inputs", "passage_input_manifest.csv")
manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  message("task_id ", task_id, " outside manifest rows=", nrow(manifest), "; exiting.")
  quit(save = "no", status = 0)
}
sample_row <- manifest[task_id, , drop = FALSE]
sample_id <- sample_row$sample[[1L]]
message("PASSAGE cancer H1/H3 sample=", sample_id, " task=", task_id)

obj <- readRDS(sample_row$file[[1L]])
gene_symbols_from_obj <- function(obj) {
  gd <- obj$gene_data
  candidates <- c("gene_name", "gene_name_short", "gene_symbol", "symbol", "Symbol",
                  "gene", "gene_id", "ensembl")
  for (cc in intersect(candidates, colnames(gd))) {
    x <- as.character(gd[[cc]])
    if (length(x) == ncol(obj$Y) && sum(!is.na(x) & nzchar(x)) > 0.5 * length(x)) {
      return(toupper(x))
    }
  }
  toupper(colnames(obj$Y))
}

collapse_to_symbols <- function(Y, symbols) {
  symbols <- toupper(symbols)
  ok <- !is.na(symbols) & nzchar(symbols) & is.finite(colSums(Y))
  Y <- Y[, ok, drop = FALSE]
  symbols <- symbols[ok]
  if (anyDuplicated(symbols)) {
    Y <- t(rowsum(t(Y), group = symbols, reorder = FALSE))
  } else {
    colnames(Y) <- symbols
  }
  keep <- apply(Y, 2L, stats::var) > 1e-8
  Y[, keep, drop = FALSE]
}

hallmark_pathways <- local({
  pathway_file <- file.path(root, "refs", "msigdb", "hallmark_human_pathways.rds")
  if (!file.exists(pathway_file)) {
    stop("Missing cached Hallmark pathway file: ", pathway_file,
         ". Create it once before launching the array to avoid parallel msigdbr cache races.")
  }
  x <- readRDS(pathway_file)
  lapply(x, function(g) unique(toupper(g)))
})

Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
hallmark_genes <- sort(unique(unlist(hallmark_pathways, use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), hallmark_genes), drop = FALSE]
pathways <- lapply(hallmark_pathways, intersect, y = colnames(Y))
pathways <- pathways[lengths(pathways) >= 5L & lengths(pathways) <= 500L]
if (length(pathways) == 0L) stop("No Hallmark pathways overlap sample ", sample_id)

coords <- as.matrix(obj$coords)
coords <- scale(coords)
coords <- sweep(coords, 2L, apply(coords, 2L, min), "-")
coords <- sweep(coords, 2L, pmax(apply(coords, 2L, max), .Machine$double.eps), "/")
coords[!is.finite(coords)] <- 0

X_all <- as.matrix(obj$X)
X_tech <- X_all[, seq_len(min(3L, ncol(X_all))), drop = FALSE]
X_h3 <- X_all

range_grid <- passage_default_range_grid(coords, n_grid = 4L, min_frac = 0.05, max_frac = 0.45)
message("Fitting engine: spots=", nrow(Y), " hallmark_genes=", ncol(Y),
        " pathways=", length(pathways), " covariates_H3=", ncol(X_h3))
engine <- passage_fit_engine_pca(
  Y = Y,
  coords = coords,
  X = X_tech,
  K = 6L,
  m = 20L,
  range_grid = range_grid,
  kernel = "matern32",
  verbose = TRUE
)
pre_h1 <- passage_h_precompute(engine, X = X_tech)
pre_h3 <- passage_h_precompute(engine, X = X_h3)

score_one <- function(ii) {
  pname <- names(pathways)[[ii]]
  genes <- pathways[[ii]]
  t0 <- proc.time()[["elapsed"]]
  common <- list(
    engine = engine,
    Y = Y,
    pathway = genes,
    weight_schemes = c("equal", "var", "range"),
    gene_names = colnames(Y),
    calibration = "permutation",
    n_perm = n_perm,
    run_burden = TRUE,
    run_spasset = FALSE
  )
  h1 <- do.call(passage_score_test, c(common, list(precomp = pre_h1, seed = 100000L + 1000L * task_id + ii)))
  h3 <- do.call(passage_score_test, c(common, list(precomp = pre_h3, seed = 200000L + 1000L * task_id + ii)))
  driver <- passage_score_test(
    engine = engine, Y = Y, pathway = genes, precomp = pre_h3,
    weight_schemes = c("equal", "var", "range"), gene_names = colnames(Y),
    calibration = "moment", run_burden = FALSE, run_spasset = TRUE
  )
  spasset <- driver$spasset
  pve <- passage_pve(engine, Y, genes, gene_names = colnames(Y))
  q1 <- h1$joint$equal$Q
  q3 <- h3$joint$equal$Q
  data.frame(
    sample = sample_id,
    cancer = obj$cancer,
    spatial_sample = obj$spatial_sample,
    reference = obj$reference,
    pathway = pname,
    pathway_size = length(genes),
    p_H1 = h1$p_omnibus,
    p_H3 = h3$p_omnibus,
    p_H1_moment = h1$p_omnibus_moment,
    p_H3_moment = h3$p_omnibus_moment,
    Q_H1_equal = q1,
    Q_H3_equal = q3,
    H3_over_H1_Q = q3 / pmax(q1, .Machine$double.eps),
    celltype_adjusted_reduction = pmax(pmin((q1 - q3) / pmax(q1, .Machine$double.eps), 1), 0),
    R2_cca = pve$summary[["R2_cca"]],
    PSVS_range = pve$summary[["PSVS_range"]],
    mean_propSV = pve$summary[["mean_propSV"]],
    p_spasset_H3 = if (!is.null(spasset)) spasset$p else NA_real_,
    spasset_size_H3 = if (!is.null(spasset)) spasset$best_size else NA_integer_,
    spasset_genes_H3 = if (!is.null(spasset)) paste(spasset$best_genes, collapse = ";") else "",
    elapsed_sec = proc.time()[["elapsed"]] - t0,
    stringsAsFactors = FALSE
  )
}

rows <- parallel::mclapply(seq_along(pathways), score_one, mc.cores = cores, mc.preschedule = FALSE)
tbl <- do.call(rbind, rows)
tbl$fdr_H1 <- stats::p.adjust(tbl$p_H1, method = "BH")
tbl$fdr_H3 <- stats::p.adjust(tbl$p_H3, method = "BH")
tbl <- tbl[order(tbl$p_H3, tbl$p_H1), , drop = FALSE]
rownames(tbl) <- NULL

out_dir <- file.path(root, "results", "passage_cancer_h1_h3", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
fit <- list(
  sample = sample_id,
  cancer = obj$cancer,
  spatial_sample = obj$spatial_sample,
  reference = obj$reference,
  summary = tbl,
  engine = engine,
  precomp = list(H1 = pre_h1, H3 = pre_h3),
  data = list(Y = Y, coords = coords, X_tech = X_tech, X_h3 = X_h3, pathways = pathways),
  n_perm = n_perm,
  K = 6L,
  m = 20L,
  range_grid = range_grid,
  primary_celltype_covariate_columns = obj$primary_celltype_covariate_columns,
  detected_celltype_covariate_columns = obj$detected_celltype_covariate_columns
)
saveRDS(fit, file.path(out_dir, "passage_cancer_h1_h3_result.rds"))
write.csv(tbl, file.path(out_dir, "passage_cancer_h1_h3_pathways.csv"), row.names = FALSE)
writeLines(c(
  paste0("# PASSAGE cancer H1/H3: ", sample_id),
  "",
  paste0("- n_perm: ", n_perm),
  paste0("- cancer: ", obj$cancer),
  paste0("- spatial_sample: ", obj$spatial_sample),
  paste0("- reference: ", obj$reference),
  paste0("- spots: ", nrow(Y)),
  paste0("- hallmark genes: ", ncol(Y)),
  paste0("- pathways: ", length(pathways)),
  paste0("- H1 covariates: ", paste(colnames(X_tech), collapse = ", ")),
  paste0("- H3 covariates: ", paste(colnames(X_h3), collapse = ", "))
), file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
