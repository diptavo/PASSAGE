#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
n_boot <- if (length(args) >= 3L) as.integer(args[[3L]]) else 40L
n_null <- if (length(args) >= 4L) as.integer(args[[4L]]) else 5L
top_k <- if (length(args) >= 5L) as.integer(args[[5L]]) else 10L

stat_names <- c("CSPS", "GSPS", "MMP", "CSV", "HCPS", "OTSAS", "score_z", "score_z_robust")
set.seed(20260814L + task_id * 1009L)

manifest_file <- file.path(root, "metadata", "driver_stability_manifest.csv")
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  message("Task ", task_id, " outside manifest rows=", nrow(manifest), "; exiting.")
  quit(save = "no", status = 0)
}
task <- manifest[task_id, , drop = FALSE]
cohort <- task$cohort[[1L]]
sample_id <- task$sample[[1L]]
pathways_to_run <- strsplit(task$selected_pathways[[1L]], ";", fixed = TRUE)[[1L]]
pathways_to_run <- pathways_to_run[nzchar(pathways_to_run)]

out_dir <- file.path(root, "results", "driver_stability", cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort_roots <- c(
  breast = "/data/Dutta_lab/SPATH/PASSAGE_cancer_panel_20260803",
  kidney = "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803",
  dlpfc = "/data/Dutta_lab/SPATH/PASSAGE_spatialDLPFC_20260802"
)
pathway_file <- file.path(cohort_roots[[cohort]], "refs", "msigdb", "hallmark_human_pathways.rds")
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
  k <- min(k, nrow(coords) - 1L)
  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  t(apply(D, 1L, function(x) order(x)[seq_len(k)]))
}

neighbor_lag <- function(Z, nn) {
  out <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  for (kk in seq_len(ncol(nn))) out <- out + Z[nn[, kk], , drop = FALSE]
  out / ncol(nn)
}

make_bins <- function(x, n = 4L) {
  x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE))
  if (length(qs) <= 2L) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE))
}

sample_matched <- function(target_bins, bins, size, exclude = integer(0)) {
  out <- integer(0)
  universe <- setdiff(seq_along(bins), exclude)
  tab <- table(target_bins)
  for (bb in names(tab)) {
    pool <- intersect(which(bins == bb), universe)
    n <- as.integer(tab[[bb]])
    if (length(pool) >= n) out <- c(out, sample(pool, n))
  }
  if (length(out) < size) {
    pool <- setdiff(universe, out)
    need <- size - length(out)
    out <- c(out, sample(pool, need, replace = length(pool) < need))
  }
  sample(out, size)
}

asset_select <- function(scores, genes, top_k = 10L, alpha = 2) {
  z <- as.numeric(scores)
  z[!is.finite(z)] <- min(z[is.finite(z)], 0)
  if (length(z) >= 3L) {
    med <- stats::median(z, na.rm = TRUE)
    madv <- stats::mad(z, constant = 1.4826, na.rm = TRUE)
    if (is.finite(madv) && madv > 1e-8) z <- (z - med) / madv
  }
  positive <- z
  positive[positive < 0] <- 0
  w <- exp(alpha * (positive - max(positive, na.rm = TRUE)))
  w[positive <= 0] <- 0
  if (sum(w, na.rm = TRUE) <= 0) w <- rep(1, length(z))
  w <- w / sum(w, na.rm = TRUE)
  ord <- order(w, decreasing = TRUE)
  keep <- ord[seq_len(min(top_k, length(ord)))]
  data.frame(gene = genes[keep], rank = seq_along(keep), driver_weight = w[keep],
             raw_score = scores[keep], stringsAsFactors = FALSE)
}

pc1_loading <- function(M) {
  if (ncol(M) <= 1L) return(1)
  C <- crossprod(M)
  eg <- suppressWarnings(eigen(C, symmetric = TRUE))
  v <- abs(eg$vectors[, 1L])
  v / pmax(sum(v), .Machine$double.eps)
}

hotspot_gene_contrib <- function(H, idx, nn) {
  Hp <- H[, idx, drop = FALSE]
  hpath <- rowMeans(Hp)
  lagp <- rep(0, length(hpath))
  for (kk in seq_len(ncol(nn))) lagp <- lagp + hpath[nn[, kk]]
  lagp <- lagp / ncol(nn)
  colMeans(sweep(Hp, 1L, lagp, "*"))
}

transport_gene_contrib <- function(absM) {
  denom <- colSums(absM)
  denom[denom <= 0 | !is.finite(denom)] <- 1
  P <- sweep(absM, 2L, denom, "/")
  pbar <- rowMeans(P)
  colSums(sqrt(sweep(P, 1L, pbar, "*")))
}

driver_scores_for_idx <- function(idx, ctx, statistic) {
  genes <- ctx$genes[idx]
  M <- ctx$Z[, idx, drop = FALSE]
  L8 <- ctx$Lag8[, idx, drop = FALSE]
  m8 <- pmax(ctx$moran8[idx], 0)
  m24 <- pmax(ctx$moran24[idx], 0)
  if (statistic == "CSPS") {
    v <- pc1_loading(L8)
    score <- m8 * v * length(v)
  } else if (statistic == "GSPS") {
    v <- pc1_loading(L8)
    low <- colSums(L8 * L8) / pmax(colSums(M * M), .Machine$double.eps)
    score <- pmax(low, 0) * v * length(v)
  } else if (statistic == "MMP") {
    score <- pmax(colSums(M * L8) / pmax(colSums(M * M), .Machine$double.eps), 0)
  } else if (statistic == "CSV") {
    score <- pmax(m8 - m24, 0)
  } else if (statistic == "HCPS") {
    score <- hotspot_gene_contrib(ctx$Hot, idx, ctx$nn8)
  } else if (statistic == "OTSAS") {
    score <- transport_gene_contrib(abs(M))
  } else if (statistic == "score_z") {
    score <- pmax(ctx$gene_score_z[idx], 0)
  } else if (statistic == "score_z_robust") {
    z <- ctx$gene_score_z
    qq <- stats::quantile(z, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    score <- pmax(pmin(pmax(z[idx], qq[[1L]]), qq[[2L]]), 0)
  } else {
    stop("Unknown statistic: ", statistic)
  }
  stats::setNames(as.numeric(score), genes)
}

build_context <- function(Y, X, coords, keep = NULL) {
  if (is.null(keep)) keep <- seq_len(nrow(Y))
  Yb <- Y[keep, , drop = FALSE]
  Xb <- as.matrix(X[keep, , drop = FALSE])
  coordsb <- scale_coords01(coords[keep, , drop = FALSE])
  Z <- qr.resid(qr(Xb), Yb)
  Z <- zscore_cols(Z)
  nn8 <- knn_index(coordsb, 8L)
  nn24 <- knn_index(coordsb, 24L)
  Lag8 <- neighbor_lag(Z, nn8)
  Lag24 <- neighbor_lag(Z, nn24)
  moran8 <- colSums(Z * Lag8) / pmax(colSums(Z * Z), .Machine$double.eps)
  moran24 <- colSums(Z * Lag24) / pmax(colSums(Z * Z), .Machine$double.eps)
  z <- moran8
  med <- stats::median(z, na.rm = TRUE)
  madv <- stats::mad(z, constant = 1.4826, na.rm = TRUE)
  gene_score_z <- if (is.finite(madv) && madv > 1e-8) (z - med) / madv else as.numeric(scale(z))
  Hot <- apply(abs(Z), 2L, function(x) x >= stats::quantile(x, 0.90, na.rm = TRUE))
  storage.mode(Hot) <- "double"
  list(Z = Z, Lag8 = Lag8, Lag24 = Lag24, moran8 = moran8, moran24 = moran24,
       gene_score_z = gene_score_z, Hot = Hot, nn8 = nn8, genes = colnames(Y))
}

spatial_block_keep <- function(coords, keep_frac = 0.80, grid_n = 5L) {
  c01 <- scale_coords01(coords)
  gx <- pmin(grid_n, pmax(1L, floor(c01[, 1L] * grid_n) + 1L))
  gy <- pmin(grid_n, pmax(1L, floor(c01[, 2L] * grid_n) + 1L))
  block <- paste(gx, gy, sep = "_")
  blocks <- unique(block)
  n_keep <- max(1L, ceiling(length(blocks) * keep_frac))
  keep_blocks <- sample(blocks, n_keep)
  which(block %in% keep_blocks)
}

jaccard_mean <- function(sets) {
  if (length(sets) < 2L) return(NA_real_)
  vals <- numeric(0)
  for (i in seq_len(length(sets) - 1L)) {
    for (j in seq.int(i + 1L, length(sets))) {
      u <- union(sets[[i]], sets[[j]])
      vals <- c(vals, if (length(u)) length(intersect(sets[[i]], sets[[j]])) / length(u) else NA_real_)
    }
  }
  mean(vals, na.rm = TRUE)
}

t0 <- proc.time()[["elapsed"]]
obj <- readRDS(task$file[[1L]])
Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
pathway_genes <- sort(unique(unlist(pathways, use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
X <- as.matrix(obj$X)
coords <- as.matrix(obj$coords)
full_ctx <- build_context(Y, X, coords)

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(full_ctx$moran8, 0)), sep = "_")

boot_keep <- vector("list", n_boot)
for (bb in seq_len(n_boot)) boot_keep[[bb]] <- spatial_block_keep(coords)
boot_ctx <- vector("list", n_boot)
for (bb in seq_len(n_boot)) boot_ctx[[bb]] <- build_context(Y, X, coords, keep = boot_keep[[bb]])

driver_rows <- list()
summary_rows <- list()
null_rows <- list()
dr <- sr <- nr <- 0L

for (pw in pathways_to_run) {
  target_genes <- intersect(pathways[[pw]], colnames(Y))
  if (length(target_genes) < 10L) next
  target_idx <- match(target_genes, colnames(Y))
  null_sets <- vector("list", n_null)
  for (nn in seq_len(n_null)) {
    null_sets[[nn]] <- sample_matched(bins[target_idx], bins, length(target_idx), exclude = target_idx)
  }
  for (stat in stat_names) {
    selected_sets <- vector("list", n_boot)
    ranks <- list()
    for (bb in seq_len(n_boot)) {
      ctx <- boot_ctx[[bb]]
      scores <- driver_scores_for_idx(target_idx, ctx, stat)
      sel <- asset_select(scores, names(scores), top_k = top_k)
      selected_sets[[bb]] <- sel$gene
      sel$bootstrap <- bb
      sel$cohort <- cohort
      sel$sample <- sample_id
      sel$pathway <- pw
      sel$pathway_size <- length(target_idx)
      sel$statistic <- stat
      dr <- dr + 1L
      driver_rows[[dr]] <- sel[, c("cohort", "sample", "pathway", "pathway_size", "statistic", "bootstrap", "gene", "rank", "driver_weight", "raw_score")]
      ranks[[bb]] <- stats::setNames(sel$rank, sel$gene)
    }
    all_sel <- unlist(selected_sets, use.names = FALSE)
    freq <- sort(table(all_sel) / n_boot, decreasing = TRUE)
    sr <- sr + 1L
    summary_rows[[sr]] <- data.frame(
      cohort = cohort,
      sample = sample_id,
      pathway = pw,
      pathway_size = length(target_idx),
      statistic = stat,
      n_boot = n_boot,
      top_k = top_k,
      n_unique_selected = length(freq),
      max_selection_frequency = if (length(freq)) as.numeric(freq[[1L]]) else NA_real_,
      n_genes_freq_ge_050 = sum(freq >= 0.50),
      n_genes_freq_ge_025 = sum(freq >= 0.25),
      mean_topk_jaccard = jaccard_mean(selected_sets),
      top_genes = paste(names(head(freq, 20L)), collapse = ";"),
      top_gene_frequencies = paste(sprintf("%.3f", as.numeric(head(freq, 20L))), collapse = ";"),
      stringsAsFactors = FALSE
    )
    for (nn in seq_len(n_null)) {
      nsets <- vector("list", n_boot)
      for (bb in seq_len(n_boot)) {
        ctx <- boot_ctx[[bb]]
        nscores <- driver_scores_for_idx(null_sets[[nn]], ctx, stat)
        nsets[[bb]] <- asset_select(nscores, names(nscores), top_k = top_k)$gene
      }
      nf <- sort(table(unlist(nsets, use.names = FALSE)) / n_boot, decreasing = TRUE)
      nr <- nr + 1L
      null_rows[[nr]] <- data.frame(
        cohort = cohort,
        sample = sample_id,
        pathway = pw,
        statistic = stat,
        null_id = nn,
        null_size = length(null_sets[[nn]]),
        null_max_selection_frequency = if (length(nf)) as.numeric(nf[[1L]]) else NA_real_,
        null_n_genes_freq_ge_050 = sum(nf >= 0.50),
        null_mean_topk_jaccard = jaccard_mean(nsets),
        stringsAsFactors = FALSE
      )
    }
  }
}

drivers <- if (length(driver_rows)) do.call(rbind, driver_rows) else data.frame()
summ <- if (length(summary_rows)) do.call(rbind, summary_rows) else data.frame()
nulls <- if (length(null_rows)) do.call(rbind, null_rows) else data.frame()
if (nrow(summ) && nrow(nulls)) {
  key <- paste(summ$cohort, summ$sample, summ$pathway, summ$statistic, sep = "\r")
  nkey <- paste(nulls$cohort, nulls$sample, nulls$pathway, nulls$statistic, sep = "\r")
  summ$null_p_max_frequency <- vapply(seq_len(nrow(summ)), function(i) {
    z <- nulls$null_max_selection_frequency[nkey == key[[i]]]
    (1 + sum(z >= summ$max_selection_frequency[[i]], na.rm = TRUE)) / (1 + length(z))
  }, numeric(1))
  summ$null_p_jaccard <- vapply(seq_len(nrow(summ)), function(i) {
    z <- nulls$null_mean_topk_jaccard[nkey == key[[i]]]
    (1 + sum(z >= summ$mean_topk_jaccard[[i]], na.rm = TRUE)) / (1 + length(z))
  }, numeric(1))
}

prefix <- file.path(out_dir, paste0("driver_stability_", cohort, "_task", task_id))
write.csv(drivers, paste0(prefix, "_drivers.csv"), row.names = FALSE)
write.csv(summ, paste0(prefix, "_summary.csv"), row.names = FALSE)
write.csv(nulls, paste0(prefix, "_nulls.csv"), row.names = FALSE)
writeLines(c(
  paste0("cohort=", cohort),
  paste0("sample=", sample_id),
  paste0("pathways=", length(pathways_to_run)),
  paste0("statistics=", paste(stat_names, collapse = ",")),
  paste0("n_boot=", n_boot),
  paste0("n_null=", n_null),
  paste0("top_k=", top_k),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", prefix)
