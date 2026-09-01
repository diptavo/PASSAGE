#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
B <- if (length(args) >= 2L) as.integer(args[[2L]]) else 999L
top_k <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L
gene_loc <- if (length(args) >= 4L) args[[4L]] else "/data/Dutta_lab/tools/NCBI38.gene.loc"
gene_result_root <- if (length(args) >= 5L) args[[5L]] else "/data/Dutta_lab/SPATH/PASSAGE_kidney_RCC_GWAS_20260803/results/magma_gene"

source(file.path(root, "scripts", "passage_kidney_common.R"))
set.seed(20260816L)
out_dir <- file.path(root, "results", "gwas_driver_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sample_matched_pool <- function(target_bins, bins, pool_idx, size) {
  out <- integer(0)
  tab <- table(target_bins)
  for (bb in names(tab)) {
    pool <- intersect(which(bins == bb), pool_idx)
    n <- as.integer(tab[[bb]])
    if (length(pool) >= n) out <- c(out, sample(pool, n))
  }
  if (length(out) < size) {
    pool <- setdiff(pool_idx, out)
    need <- size - length(out)
    if (!length(pool)) pool <- pool_idx
    out <- c(out, sample(pool, need, replace = length(pool) < need))
  }
  sample(out, size)
}

drivers <- read.csv(file.path(root, "results", "driver_stability_summary", "driver_stability_gene_frequencies.csv"), stringsAsFactors = FALSE)
manifest <- read.csv(file.path(root, "metadata", "driver_manifest.csv"), stringsAsFactors = FALSE)
pathways <- readRDS(file.path(root, "refs", "msigdb", "msigdb_human_pathways_filtered.rds"))

loc <- read.table(gene_loc, stringsAsFactors = FALSE, fill = TRUE, quote = "", comment.char = "")
if (ncol(loc) < 6L) stop("Unexpected gene location file: ", gene_loc)
symbol_map <- unique(data.frame(entrez = as.character(loc[[1L]]), symbol = toupper(as.character(loc[[6L]])), stringsAsFactors = FALSE))

read_magma <- function(pheno) {
  f <- file.path(gene_result_root, paste0(pheno, ".genes.out"))
  z <- read.table(f, header = TRUE, stringsAsFactors = FALSE)
  z$entrez <- as.character(z$GENE)
  z <- merge(z, symbol_map, by = "entrez")
  z <- z[!duplicated(z$symbol), , drop = FALSE]
  rownames(z) <- z$symbol
  z
}
phenos <- c("RCC", "CC2", "PRCC")
magma <- setNames(lapply(phenos, read_magma), phenos)

rows <- list()
rr <- 0L
for (ii in seq_len(nrow(manifest))) {
  task <- manifest[ii, , drop = FALSE]
  obj <- readRDS(task$file[[1L]])
  Y <- collapse_to_symbols(as.matrix(obj$Y), gene_symbols_from_obj(obj))
  pws <- strsplit(task$selected_pathways[[1L]], ";", fixed = TRUE)[[1L]]
  pws <- pws[nzchar(pws) & pws %in% names(pathways)]
  pathway_genes <- sort(unique(unlist(pathways[pws], use.names = FALSE)))
  Y <- Y[, intersect(colnames(Y), pathway_genes), drop = FALSE]
  ctx <- build_context(Y, as.matrix(obj$X), as.matrix(obj$coords))
  expr_mean <- colMeans(Y, na.rm = TRUE)
  detect_rate <- colMeans(Y > 0, na.rm = TRUE)
  expr_var <- apply(Y, 2L, stats::var, na.rm = TRUE)
  bins <- paste(make_bins(expr_mean), make_bins(detect_rate), make_bins(expr_var), make_bins(pmax(ctx$moran8, 0)), sep = "_")

  for (pw in pws) {
    genes <- intersect(pathways[[pw]], colnames(Y))
    if (length(genes) < 15L) next
    idx <- match(genes, colnames(Y))
    for (stat in c("score_z", "score_z_robust", "CSPS", "GSPS", "HCPS")) {
      gf <- drivers[drivers$sample == task$sample[[1L]] & drivers$pathway == pw & drivers$statistic == stat, , drop = FALSE]
      gf <- gf[order(-gf$selection_frequency, gf$mean_rank), , drop = FALSE]
      driver_genes <- head(gf$gene[gf$selection_frequency >= 0.50], top_k)
      driver_genes <- intersect(driver_genes, genes)
      if (length(driver_genes) < 2L) next
      driver_idx <- match(driver_genes, colnames(Y))
      target_bins <- bins[driver_idx]
      pool_idx <- setdiff(idx, driver_idx)
      if (length(pool_idx) < length(driver_idx)) next

      for (pheno in phenos) {
        gw <- magma[[pheno]]
        d <- intersect(driver_genes, rownames(gw))
        if (length(d) < 2L) next
        gwas_score <- -log10(pmax(gw$P[match(d, rownames(gw))], .Machine$double.xmin))
        obs_mean_logp <- mean(gwas_score, na.rm = TRUE)
        obs_min_p <- min(gw$P[match(d, rownames(gw))], na.rm = TRUE)
        null_mean <- rep(NA_real_, B)
        for (bb in seq_len(B)) {
          nidx <- sample_matched_pool(target_bins, bins, pool_idx, length(driver_idx))
          ngenes <- colnames(Y)[nidx]
          ngenes <- intersect(ngenes, rownames(gw))
          null_mean[[bb]] <- if (length(ngenes) >= 2L) {
            mean(-log10(pmax(gw$P[match(ngenes, rownames(gw))], .Machine$double.xmin)), na.rm = TRUE)
          } else NA_real_
        }
        rr <- rr + 1L
        rows[[rr]] <- data.frame(
          sample = task$sample[[1L]],
          reference = task$reference[[1L]],
          pathway = pw,
          statistic = stat,
          phenotype = pheno,
          n_driver_genes = length(driver_genes),
          n_driver_genes_in_magma = length(d),
          observed_mean_neglog10p = obs_mean_logp,
          observed_min_gene_p = obs_min_p,
          null_mean = mean(null_mean, na.rm = TRUE),
          null_sd = stats::sd(null_mean, na.rm = TRUE),
          empirical_p = (1 + sum(null_mean >= obs_mean_logp, na.rm = TRUE)) / (1 + sum(is.finite(null_mean))),
          driver_genes = paste(driver_genes, collapse = ";"),
          driver_genes_in_magma = paste(d, collapse = ";"),
          B = B,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

out <- if (length(rows)) do.call(rbind, rows) else data.frame()
if (nrow(out)) {
  out$q_value_reference_stat_pheno <- ave(out$empirical_p, out$reference, out$statistic, out$phenotype, FUN = function(x) p.adjust(x, "BH"))
  out <- out[order(out$phenotype, out$reference, out$statistic, out$q_value_reference_stat_pheno, out$empirical_p), , drop = FALSE]
}
write.csv(out, file.path(out_dir, "gwas_driver_matched_enrichment_results.csv"), row.names = FALSE)

if (nrow(out)) {
  meta <- do.call(rbind, lapply(split(out, list(out$reference, out$statistic, out$phenotype), drop = TRUE), function(z) {
    p <- z$empirical_p[z$empirical_p > 0 & z$empirical_p <= 1]
    data.frame(
      reference = z$reference[[1L]],
      statistic = z$statistic[[1L]],
      phenotype = z$phenotype[[1L]],
      n_tests = nrow(z),
      fisher_p = if (length(p)) stats::pchisq(-2 * sum(log(p)), df = 2 * length(p), lower.tail = FALSE) else NA_real_,
      min_p = min(z$empirical_p, na.rm = TRUE),
      median_p = stats::median(z$empirical_p, na.rm = TRUE),
      n_q05 = sum(z$q_value_reference_stat_pheno <= 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  meta$q_value_global <- p.adjust(meta$fisher_p, "BH")
} else {
  meta <- data.frame()
}
write.csv(meta, file.path(out_dir, "gwas_driver_matched_enrichment_meta.csv"), row.names = FALSE)

md <- c(
  "# RCC GWAS Validation of PASSAGE Spatial Driver Genes",
  "",
  paste0("- result rows: ", nrow(out)),
  paste0("- phenotypes: ", paste(phenos, collapse = ", ")),
  paste0("- matched null draws per test: ", B),
  "",
  "The empirical test asks whether stable PASSAGE drivers have stronger MAGMA gene-level association than matched non-driver genes from the same spatial dataset/pathway universe."
)
if (nrow(meta)) {
  top <- head(meta[order(meta$q_value_global, meta$fisher_p), , drop = FALSE], 30L)
  md <- c(md, "", "## Top Meta-Signals", "",
          paste(c("reference", "statistic", "phenotype", "fisher_p", "q", "n_tests", "n_q05"), collapse = " | "),
          paste(rep("---", 7), collapse = " | "))
  for (ii in seq_len(nrow(top))) {
    r <- top[ii, ]
    md <- c(md, paste(c(r$reference, r$statistic, r$phenotype,
                        sprintf("%.4g", r$fisher_p), sprintf("%.4g", r$q_value_global),
                        r$n_tests, r$n_q05), collapse = " | "))
  }
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", out_dir)
