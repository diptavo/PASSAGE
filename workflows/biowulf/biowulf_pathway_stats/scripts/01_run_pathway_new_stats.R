#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
bench_root <- if (length(args) >= 1L) args[[1L]] else getwd()
cohort <- if (length(args) >= 2L) args[[2L]] else "kidney"
task_id <- if (length(args) >= 3L) as.integer(args[[3L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
B <- if (length(args) >= 4L) as.integer(args[[4L]]) else 199L
min_genes <- if (length(args) >= 5L) as.integer(args[[5L]]) else 15L
max_genes <- if (length(args) >= 6L) as.integer(args[[6L]]) else 250L

cohort_key <- match(cohort, c("breast", "kidney", "dlpfc"))
if (is.na(cohort_key)) stop("Unknown cohort: ", cohort)
set.seed(20260814L + task_id + cohort_key * 10000L)

cohort_roots <- c(
  breast = "/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803",
  kidney = "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803",
  dlpfc = "/data/Dutta_lab/SPATH/PASSAGE_spatialDLPFC_20260802"
)
data_root <- cohort_roots[[cohort]]
if (is.null(data_root) || !dir.exists(data_root)) stop("Missing cohort root for ", cohort)

out_dir <- file.path(bench_root, "results", "pathway_new_stats", cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(file.path(data_root, "data", "passage_inputs", "passage_input_manifest.csv"),
                     stringsAsFactors = FALSE)
if (cohort == "breast") manifest <- manifest[manifest$cancer == "breast", , drop = FALSE]
if (task_id < 1L || task_id > nrow(manifest)) {
  message("Task ", task_id, " outside ", cohort, " manifest rows=", nrow(manifest), "; exiting.")
  quit(save = "no", status = 0)
}
row <- manifest[task_id, , drop = FALSE]
sample_id <- row$sample[[1L]]

pathway_file <- file.path(data_root, "refs", "msigdb", "hallmark_human_pathways.rds")
if (!file.exists(pathway_file)) stop("Missing pathway file: ", pathway_file)
pathways <- lapply(readRDS(pathway_file), function(g) unique(toupper(as.character(g))))
names(pathways) <- make.names(names(pathways), unique = TRUE)

gene_symbols_from_obj <- function(obj) {
  gd <- obj$gene_data
  candidates <- c("gene_name", "gene_name_short", "gene_symbol", "symbol", "Symbol",
                  "gene", "gene_id", "ensembl")
  for (cc in intersect(candidates, colnames(gd))) {
    x <- as.character(gd[[cc]])
    if (length(x) == ncol(obj$Y) && sum(!is.na(x) & nzchar(x)) > 0.5 * length(x)) return(toupper(x))
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

zscore_cols <- function(M) {
  M <- as.matrix(M)
  mu <- colMeans(M, na.rm = TRUE)
  S <- sweep(M, 2L, mu, "-")
  sdv <- sqrt(pmax(colMeans(S * S, na.rm = TRUE), 1e-8))
  sweep(S, 2L, sdv, "/")
}

scale_coords01 <- function(coords) {
  coords <- as.matrix(coords)
  coords <- coords[, seq_len(2L), drop = FALSE]
  coords <- sweep(coords, 2L, apply(coords, 2L, min, na.rm = TRUE), "-")
  denom <- pmax(apply(coords, 2L, max, na.rm = TRUE), .Machine$double.eps)
  coords <- sweep(coords, 2L, denom, "/")
  coords[!is.finite(coords)] <- 0
  coords
}

knn_index <- function(coords, k = 8L) {
  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  t(apply(D, 1L, function(x) order(x)[seq_len(k)]))
}

neighbor_lag <- function(Z, nn) {
  out <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  for (kk in seq_len(ncol(nn))) out <- out + Z[nn[, kk], , drop = FALSE]
  out / ncol(nn)
}

pc1_fraction <- function(M) {
  if (ncol(M) <= 1L) return(1)
  C <- crossprod(M)
  ev <- suppressWarnings(eigen(C, symmetric = TRUE, only.values = TRUE)$values)
  ev <- pmax(ev, 0)
  if (!sum(ev) > 0) return(0)
  max(ev) / sum(ev)
}

effective_breadth <- function(w, p) {
  w <- pmax(w, 0)
  if (!sum(w) > 0 || p <= 1L) return(0)
  ((sum(w)^2 / sum(w^2)) / p)
}

hotspot_score <- function(H, nn) {
  h <- rowMeans(H)
  if (!is.finite(mean(h)) || mean(h) <= 0) return(0)
  hlag <- rep(0, length(h))
  for (kk in seq_len(ncol(nn))) hlag <- hlag + h[nn[, kk]]
  hlag <- hlag / ncol(nn)
  mean(h * hlag) / (mean(h)^2)
}

transport_alignment <- function(absM) {
  denom <- colSums(absM)
  denom[denom <= 0 | !is.finite(denom)] <- 1
  P <- sweep(absM, 2L, denom, "/")
  pbar <- rowMeans(P)
  affinity <- colSums(sqrt(sweep(P, 1L, pbar, "*")))
  concentration <- sqrt(nrow(P) * sum(pbar * pbar))
  mean(affinity) * concentration
}

stat_all <- function(idx, Z, Lag8, Lag24, moran8, moran24, Hot, nn8) {
  p <- length(idx)
  M <- Z[, idx, drop = FALSE]
  L8 <- Lag8[, idx, drop = FALSE]
  m8 <- pmax(moran8[idx], 0)
  m24 <- pmax(moran24[idx], 0)
  pc1 <- pc1_fraction(L8)
  breadth <- effective_breadth(m8, p)
  signal <- sum(m8)
  lowpass <- sum(L8 * L8) / pmax(sum(M * M), .Machine$double.eps)
  mult_moran <- sum(M * L8) / pmax(sum(M * M), .Machine$double.eps)
  range_contrast <- pmax(m8 - m24, 0)
  range_agreement <- 1 / (1 + stats::mad(range_contrast, constant = 1, na.rm = TRUE))
  hcps <- hotspot_score(Hot[, idx, drop = FALSE], nn8)
  otsas <- transport_alignment(abs(M)) * breadth
  c(
    CSPS = signal * pc1 * breadth,
    GSPS = lowpass * pc1 * breadth,
    MMP = mult_moran,
    CSV = mean(range_contrast) * range_agreement * breadth,
    HCPS = hcps,
    OTSAS = otsas
  )
}

make_bins <- function(x, n = 4L) {
  x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE))
  if (length(qs) <= 2L) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE))
}

sample_matched <- function(target_bins, bins, size, exclude = integer(0)) {
  out <- integer(0)
  tab <- table(target_bins)
  universe <- setdiff(seq_along(bins), exclude)
  for (bb in names(tab)) {
    pool <- intersect(which(bins == bb), universe)
    n <- as.integer(tab[[bb]])
    if (length(pool) >= n) out <- c(out, sample(pool, n))
  }
  if (length(out) < size) {
    pool <- setdiff(universe, out)
    if (length(pool) < size - length(out)) pool <- universe
    out <- c(out, sample(pool, size - length(out), replace = length(pool) < size - length(out)))
  }
  sample(out, size)
}

t0 <- proc.time()[["elapsed"]]
obj <- readRDS(row$file[[1L]])
Y_raw <- as.matrix(obj$Y)
Y <- collapse_to_symbols(Y_raw, gene_symbols_from_obj(obj))
pathway_genes <- sort(unique(unlist(pathways, use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
if (ncol(Y) < 500L) stop("Too few pathway-universe genes after overlap: ", ncol(Y))

X <- as.matrix(obj$X)
Z <- qr.resid(qr(X), Y)
Z <- zscore_cols(Z)
coords <- scale_coords01(obj$coords)
nn8 <- knn_index(coords, 8L)
nn24 <- knn_index(coords, 24L)
Lag8 <- neighbor_lag(Z, nn8)
Lag24 <- neighbor_lag(Z, nn24)
moran8 <- colSums(Z * Lag8) / pmax(colSums(Z * Z), .Machine$double.eps)
moran24 <- colSums(Z * Lag24) / pmax(colSums(Z * Z), .Machine$double.eps)

absZ <- abs(Z)
Hot <- apply(absZ, 2L, function(x) x >= stats::quantile(x, 0.90, na.rm = TRUE))
storage.mode(Hot) <- "double"

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(moran8, 0)), sep = "_")

pathways2 <- lapply(pathways, intersect, y = colnames(Z))
pathways2 <- pathways2[lengths(pathways2) >= min_genes & lengths(pathways2) <= max_genes]
if (!length(pathways2)) stop("No pathways after size/overlap filtering")

gene_names <- colnames(Z)
rows <- vector("list", length(pathways2) * 6L)
rr <- 0L
for (pw in names(pathways2)) {
  target_genes <- pathways2[[pw]]
  target_idx <- match(target_genes, gene_names)
  target_idx <- target_idx[is.finite(target_idx)]
  size <- length(target_idx)
  obs <- stat_all(target_idx, Z, Lag8, Lag24, moran8, moran24, Hot, nn8)
  null_mat <- matrix(NA_real_, nrow = B, ncol = length(obs), dimnames = list(NULL, names(obs)))
  target_bins <- bins[target_idx]
  for (bb in seq_len(B)) {
    idx <- sample_matched(target_bins, bins, size, exclude = target_idx)
    null_mat[bb, ] <- stat_all(idx, Z, Lag8, Lag24, moran8, moran24, Hot, nn8)
  }
  for (st in names(obs)) {
    rr <- rr + 1L
    rows[[rr]] <- data.frame(
      cohort = cohort,
      sample = sample_id,
      task_id = task_id,
      pathway = pw,
      pathway_size = size,
      statistic = st,
      observed = obs[[st]],
      null_mean = mean(null_mat[, st], na.rm = TRUE),
      null_sd = stats::sd(null_mat[, st], na.rm = TRUE),
      p_value = (1 + sum(null_mat[, st] >= obs[[st]], na.rm = TRUE)) / (B + 1),
      B = B,
      elapsed_sec = proc.time()[["elapsed"]] - t0,
      genes = paste(target_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

out <- do.call(rbind, rows[seq_len(rr)])
out$q_value_sample_stat <- ave(out$p_value, out$sample, out$statistic, FUN = function(x) p.adjust(x, "BH"))
out <- out[order(out$statistic, out$q_value_sample_stat, out$p_value, out$pathway), , drop = FALSE]

outfile <- file.path(out_dir, paste0("pathway_new_stats_", cohort, "_task", task_id, ".csv"))
write.csv(out, outfile, row.names = FALSE)
writeLines(c(
  paste0("cohort=", cohort),
  paste0("sample=", sample_id),
  paste0("spots=", nrow(Z)),
  paste0("genes=", ncol(Z)),
  paste0("pathways=", length(pathways2)),
  paste0("B=", B),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), file.path(out_dir, paste0("pathway_new_stats_", cohort, "_task", task_id, ".log")))
message("Wrote ", outfile)
