#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
task_id <- if (length(args) >= 2L) as.integer(args[[2L]]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))
top_k <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L
n_controls <- if (length(args) >= 4L) as.integer(args[[4L]]) else 10L

stat_names <- c("CSPS", "GSPS", "MMP", "CSV", "HCPS", "OTSAS", "score_z", "score_z_robust")
set.seed(20260815L + task_id * 997L)

manifest <- read.csv(file.path(root, "metadata", "driver_stability_manifest.csv"), stringsAsFactors = FALSE)
if (task_id < 1L || task_id > nrow(manifest)) {
  message("Task ", task_id, " outside manifest rows=", nrow(manifest), "; exiting.")
  quit(save = "no", status = 0)
}
task <- manifest[task_id, , drop = FALSE]
cohort <- task$cohort[[1L]]
sample_id <- task$sample[[1L]]
pathways_to_run <- strsplit(task$selected_pathways[[1L]], ";", fixed = TRUE)[[1L]]
pathways_to_run <- pathways_to_run[nzchar(pathways_to_run)]

gene_freq_file <- file.path(root, "results", "driver_stability_summary", "driver_stability_gene_frequencies.csv")
gene_freq <- read.csv(gene_freq_file, stringsAsFactors = FALSE)
gene_freq <- gene_freq[gene_freq$cohort == cohort & gene_freq$sample == sample_id, , drop = FALSE]

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
  candidates <- c("gene_name", "gene_name_short", "gene_symbol", "symbol", "Symbol", "gene", "gene_id", "ensembl")
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
  if (anyDuplicated(symbols)) Y <- t(rowsum(t(Y), group = symbols, reorder = FALSE)) else colnames(Y) <- symbols
  keep <- apply(Y, 2L, stats::var) > 1e-8
  Y[, keep, drop = FALSE]
}

zscore_cols <- function(M) {
  mu <- colMeans(M, na.rm = TRUE)
  S <- sweep(as.matrix(M), 2L, mu, "-")
  sdv <- sqrt(pmax(colMeans(S * S, na.rm = TRUE), 1e-8))
  sweep(S, 2L, sdv, "/")
}

scale_coords01 <- function(coords) {
  coords <- as.matrix(coords)[, seq_len(2L), drop = FALSE]
  coords <- sweep(coords, 2L, apply(coords, 2L, min, na.rm = TRUE), "-")
  denom <- pmax(apply(coords, 2L, max, na.rm = TRUE), .Machine$double.eps)
  coords <- sweep(coords, 2L, denom, "/")
  coords[!is.finite(coords)] <- 0
  coords
}

knn_index <- function(coords, k = 8L) {
  D <- as.matrix(stats::dist(coords))
  diag(D) <- Inf
  t(apply(D, 1L, function(x) order(x)[seq_len(min(k, length(x) - 1L))]))
}

neighbor_lag <- function(Z, nn) {
  out <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  for (kk in seq_len(ncol(nn))) out <- out + Z[nn[, kk], , drop = FALSE]
  out / ncol(nn)
}

pc1_fraction <- function(M) {
  if (ncol(M) <= 1L) return(1)
  ev <- suppressWarnings(eigen(crossprod(M), symmetric = TRUE, only.values = TRUE)$values)
  ev <- pmax(ev, 0)
  if (!sum(ev) > 0) return(0)
  max(ev) / sum(ev)
}

effective_breadth <- function(w, p) {
  w <- pmax(w, 0)
  if (!sum(w) > 0 || p <= 1L) return(0)
  (sum(w)^2 / sum(w^2)) / p
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

build_context <- function(Y, X, coords) {
  Z <- qr.resid(qr(as.matrix(X)), Y)
  Z <- zscore_cols(Z)
  coords <- scale_coords01(coords)
  nn8 <- knn_index(coords, 8L)
  nn24 <- knn_index(coords, 24L)
  Lag8 <- neighbor_lag(Z, nn8)
  Lag24 <- neighbor_lag(Z, nn24)
  moran8 <- colSums(Z * Lag8) / pmax(colSums(Z * Z), .Machine$double.eps)
  moran24 <- colSums(Z * Lag24) / pmax(colSums(Z * Z), .Machine$double.eps)
  med <- stats::median(moran8, na.rm = TRUE)
  madv <- stats::mad(moran8, constant = 1.4826, na.rm = TRUE)
  gene_score_z <- if (is.finite(madv) && madv > 1e-8) (moran8 - med) / madv else as.numeric(scale(moran8))
  Hot <- apply(abs(Z), 2L, function(x) x >= stats::quantile(x, 0.90, na.rm = TRUE))
  storage.mode(Hot) <- "double"
  list(Z = Z, Lag8 = Lag8, Lag24 = Lag24, moran8 = moran8, moran24 = moran24,
       gene_score_z = gene_score_z, Hot = Hot, nn8 = nn8, genes = colnames(Y))
}

pathway_stat <- function(idx, ctx, statistic) {
  p <- length(idx)
  if (p < 2L) return(NA_real_)
  M <- ctx$Z[, idx, drop = FALSE]
  L8 <- ctx$Lag8[, idx, drop = FALSE]
  m8 <- pmax(ctx$moran8[idx], 0)
  m24 <- pmax(ctx$moran24[idx], 0)
  if (statistic == "CSPS") {
    return(sum(m8) * pc1_fraction(L8) * effective_breadth(m8, p))
  }
  if (statistic == "GSPS") {
    lowpass <- sum(L8 * L8) / pmax(sum(M * M), .Machine$double.eps)
    return(lowpass * pc1_fraction(L8) * effective_breadth(m8, p))
  }
  if (statistic == "MMP") {
    return(sum(M * L8) / pmax(sum(M * M), .Machine$double.eps))
  }
  if (statistic == "CSV") {
    rc <- pmax(m8 - m24, 0)
    return(mean(rc) * (1 / (1 + stats::mad(rc, constant = 1, na.rm = TRUE))) * effective_breadth(m8, p))
  }
  if (statistic == "HCPS") return(hotspot_score(ctx$Hot[, idx, drop = FALSE], ctx$nn8))
  if (statistic == "OTSAS") return(transport_alignment(abs(M)) * effective_breadth(m8, p))
  if (statistic == "score_z") return(mean(ctx$gene_score_z[idx], na.rm = TRUE) * sqrt(p))
  if (statistic == "score_z_robust") {
    z <- ctx$gene_score_z
    qq <- stats::quantile(z, probs = c(0.01, 0.99), na.rm = TRUE, names = FALSE)
    zr <- pmin(pmax(z[idx], qq[[1L]]), qq[[2L]])
    return(mean(zr, na.rm = TRUE) * sqrt(p))
  }
  stop("Unknown statistic: ", statistic)
}

make_bins <- function(x, n = 4L) {
  x[!is.finite(x)] <- stats::median(x[is.finite(x)], na.rm = TRUE)
  qs <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1L), na.rm = TRUE))
  if (length(qs) <= 2L) return(rep(1L, length(x)))
  as.integer(cut(x, breaks = qs, include.lowest = TRUE, labels = FALSE))
}

t0 <- proc.time()[["elapsed"]]
obj <- readRDS(task$file[[1L]])
Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
pathway_genes <- sort(unique(unlist(pathways, use.names = FALSE)))
Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
ctx <- build_context(Y, as.matrix(obj$X), as.matrix(obj$coords))

expr_mean <- colMeans(Y, na.rm = TRUE)
detect_rate <- colMeans(Y > 0, na.rm = TRUE)
expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(ctx$moran8, 0)), sep = "_")

rows <- list()
rr <- 0L
for (pw in pathways_to_run) {
  genes <- intersect(pathways[[pw]], colnames(Y))
  if (length(genes) < 10L) next
  idx <- match(genes, colnames(Y))
  for (stat in stat_names) {
    gf <- gene_freq[gene_freq$pathway == pw & gene_freq$statistic == stat, , drop = FALSE]
    gf <- gf[order(-gf$selection_frequency, gf$mean_rank), , drop = FALSE]
    drivers <- head(gf$gene, top_k)
    drivers <- intersect(drivers, genes)
    non_drivers <- setdiff(genes, drivers)
    full <- pathway_stat(idx, ctx, stat)
    for (gene in drivers) {
      gi <- match(gene, colnames(Y))
      drop <- full - pathway_stat(setdiff(idx, gi), ctx, stat)
      rr <- rr + 1L
      rows[[rr]] <- data.frame(cohort = cohort, sample = sample_id, pathway = pw, statistic = stat,
                               gene = gene, type = "driver", selection_frequency = gf$selection_frequency[match(gene, gf$gene)],
                               full_stat = full, loo_stat = full - drop, stat_drop = drop,
                               stringsAsFactors = FALSE)
    }
    control_pool <- non_drivers
    if (length(control_pool)) {
      control_score <- abs(ctx$gene_score_z[match(control_pool, colnames(Y))])
      control_pool <- control_pool[order(control_score, decreasing = TRUE)]
      controls <- head(control_pool, min(n_controls, length(control_pool)))
      for (gene in controls) {
        gi <- match(gene, colnames(Y))
        drop <- full - pathway_stat(setdiff(idx, gi), ctx, stat)
        rr <- rr + 1L
        rows[[rr]] <- data.frame(cohort = cohort, sample = sample_id, pathway = pw, statistic = stat,
                                 gene = gene, type = "control", selection_frequency = NA_real_,
                                 full_stat = full, loo_stat = full - drop, stat_drop = drop,
                                 stringsAsFactors = FALSE)
      }
    }
  }
}
out <- if (length(rows)) do.call(rbind, rows) else data.frame()
out_dir <- file.path(root, "results", "leave_one_out", cohort)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
prefix <- file.path(out_dir, paste0("leave_one_out_", cohort, "_task", task_id))
write.csv(out, paste0(prefix, ".csv"), row.names = FALSE)
writeLines(c(
  paste0("cohort=", cohort),
  paste0("sample=", sample_id),
  paste0("rows=", nrow(out)),
  paste0("elapsed_sec=", round(proc.time()[["elapsed"]] - t0, 3))
), paste0(prefix, ".log"))
message("Wrote ", paste0(prefix, ".csv"))
