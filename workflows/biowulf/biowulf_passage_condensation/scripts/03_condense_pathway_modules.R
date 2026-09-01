#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) args[[1L]] else getwd()
q_cut <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 0.10
jaccard_cut <- if (length(args) >= 3L) as.numeric(args[[3L]]) else 0.25
overlap_cut <- if (length(args) >= 4L) as.numeric(args[[4L]]) else 0.50
max_candidates <- if (length(args) >= 5L) as.integer(args[[5L]]) else 2500L
primary_stats <- if (length(args) >= 6L) strsplit(args[[6L]], ",", fixed = TRUE)[[1L]] else c("score_z", "score_z_robust", "CSPS", "GSPS", "HCPS")

if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package is required")
meta <- read.csv(file.path(root, "results", "pathway_testing_summary", "pathway_reference_meta.csv"), stringsAsFactors = FALSE)
meta <- meta[meta$statistic %in% primary_stats, , drop = FALSE]
out_dir <- file.path(root, "results", "pathway_modules")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  unique(unlist(strsplit(x[nzchar(x)], ";", fixed = TRUE), use.names = FALSE))
}

components_from_edges <- function(n, i, j) {
  parent <- seq_len(n)
  find <- function(x) {
    while (parent[[x]] != x) {
      parent[[x]] <<- parent[[parent[[x]]]]
      x <- parent[[x]]
    }
    x
  }
  union <- function(a, b) {
    ra <- find(a); rb <- find(b)
    if (ra != rb) parent[[rb]] <<- ra
  }
  if (length(i)) for (k in seq_along(i)) union(i[[k]], j[[k]])
  vapply(seq_len(n), find, integer(1))
}

module_rows <- list()
member_rows <- list()
gene_sets <- list()
diag_rows <- list()
mr <- gr <- dr <- 0L

group_key <- interaction(meta$cohort, meta$reference, meta$statistic, drop = TRUE, sep = "\r")
for (key in levels(group_key)) {
  z0 <- meta[group_key == key, , drop = FALSE]
  sig <- z0[z0$q_value_cohort_reference_stat <= q_cut, , drop = FALSE]
  n_sig <- nrow(sig)
  truncated <- FALSE
  if (nrow(sig) > max_candidates) {
    sig$rank_score <- rank(sig$q_value_cohort_reference_stat, ties.method = "average") + rank(-sig$median_effect_z, ties.method = "average")
    sig <- sig[order(sig$rank_score), , drop = FALSE]
    sig <- head(sig, max_candidates)
    truncated <- TRUE
  }
  if (!nrow(sig)) {
    dr <- dr + 1L
    diag_rows[[dr]] <- data.frame(cohort = z0$cohort[[1L]], reference = z0$reference[[1L]], statistic = z0$statistic[[1L]],
                                  n_tested = nrow(z0), n_significant = n_sig, n_condensed_candidates = 0,
                                  n_modules = 0, compression_ratio = NA_real_, truncated = truncated, stringsAsFactors = FALSE)
    next
  }
  genes_list <- lapply(sig$genes, split_genes)
  names(genes_list) <- sig$pathway
  all_genes <- sort(unique(unlist(genes_list, use.names = FALSE)))
  gi <- match(unlist(genes_list, use.names = FALSE), all_genes)
  pj <- rep(seq_along(genes_list), lengths(genes_list))
  X <- Matrix::sparseMatrix(i = gi, j = pj, x = 1L, dims = c(length(all_genes), length(genes_list)))
  inter <- Matrix::summary(Matrix::crossprod(X))
  inter <- inter[inter$i < inter$j & inter$x > 0, , drop = FALSE]
  sizes <- lengths(genes_list)
  keep <- rep(FALSE, nrow(inter))
  if (nrow(inter)) {
    jac <- inter$x / (sizes[inter$i] + sizes[inter$j] - inter$x)
    ov <- inter$x / pmin(sizes[inter$i], sizes[inter$j])
    keep <- jac >= jaccard_cut | ov >= overlap_cut
  }
  edges <- inter[keep, , drop = FALSE]
  comp <- components_from_edges(length(genes_list), edges$i, edges$j)
  comp_ids <- as.integer(factor(comp))
  for (cc in sort(unique(comp_ids))) {
    idx <- which(comp_ids == cc)
    sub <- sig[idx, , drop = FALSE]
    gl <- genes_list[idx]
    freq <- sort(table(unlist(gl, use.names = FALSE)), decreasing = TRUE)
    min_count <- if (length(idx) <= 2L) 1L else max(2L, ceiling(0.20 * length(idx)))
    module_genes <- names(freq)[freq >= min_count]
    if (length(module_genes) < 10L) module_genes <- names(freq)
    module_id <- sprintf("%s__%s__%s__M%04d", sub$cohort[[1L]], sub$reference[[1L]], sub$statistic[[1L]], cc)
    module_id <- make.names(module_id)
    representative <- sub$pathway[[which.min(sub$q_value_cohort_reference_stat)]]
    mr <- mr + 1L
    module_rows[[mr]] <- data.frame(
      module_id = module_id,
      cohort = sub$cohort[[1L]],
      reference = sub$reference[[1L]],
      source_statistic = sub$statistic[[1L]],
      n_member_pathways = length(idx),
      n_module_genes = length(module_genes),
      representative_pathway = representative,
      best_q = min(sub$q_value_cohort_reference_stat, na.rm = TRUE),
      best_p = min(sub$fisher_p, na.rm = TRUE),
      median_effect_z = stats::median(sub$median_effect_z, na.rm = TRUE),
      member_pathways = paste(sub$pathway, collapse = ";"),
      module_genes = paste(module_genes, collapse = ";"),
      stringsAsFactors = FALSE
    )
    gene_sets[[module_id]] <- module_genes
    for (jj in seq_len(nrow(sub))) {
      gr <- gr + 1L
      member_rows[[gr]] <- data.frame(module_id = module_id, pathway = sub$pathway[[jj]], collection = sub$collection[[jj]],
                                      pathway_q = sub$q_value_cohort_reference_stat[[jj]], pathway_effect_z = sub$median_effect_z[[jj]],
                                      stringsAsFactors = FALSE)
    }
  }
  dr <- dr + 1L
  diag_rows[[dr]] <- data.frame(cohort = z0$cohort[[1L]], reference = z0$reference[[1L]], statistic = z0$statistic[[1L]],
                                n_tested = nrow(z0), n_significant = n_sig, n_condensed_candidates = nrow(sig),
                                n_modules = length(unique(comp_ids)),
                                compression_ratio = if (length(unique(comp_ids))) n_sig / length(unique(comp_ids)) else NA_real_,
                                truncated = truncated, stringsAsFactors = FALSE)
}

modules <- if (length(module_rows)) do.call(rbind, module_rows) else data.frame()
members <- if (length(member_rows)) do.call(rbind, member_rows) else data.frame()
diag <- if (length(diag_rows)) do.call(rbind, diag_rows) else data.frame()
write.csv(modules, file.path(out_dir, "pathway_modules.csv"), row.names = FALSE)
write.csv(members, file.path(out_dir, "pathway_module_members.csv"), row.names = FALSE)
write.csv(diag, file.path(out_dir, "pathway_module_diagnostics.csv"), row.names = FALSE)
saveRDS(gene_sets, file.path(out_dir, "pathway_module_gene_sets.rds"))

md <- c(
  "# PASSAGE Redundancy-Condensed Pathway Modules",
  "",
  paste0("- modules: ", nrow(modules)),
  paste0("- q_cut: ", q_cut),
  paste0("- jaccard_cut: ", jaccard_cut),
  paste0("- overlap_cut: ", overlap_cut),
  "",
  "## Compression Diagnostics",
  "",
  paste(c("cohort", "reference", "statistic", "significant", "modules", "compression", "truncated"), collapse = " | "),
  paste(rep("---", 7), collapse = " | ")
)
if (nrow(diag)) {
  for (ii in seq_len(nrow(diag))) {
    r <- diag[ii, ]
    md <- c(md, paste(c(r$cohort, r$reference, r$statistic, r$n_significant, r$n_modules,
                        sprintf("%.2f", r$compression_ratio), r$truncated), collapse = " | "))
  }
}
writeLines(md, file.path(out_dir, "summary.md"))
message("Wrote ", nrow(modules), " modules")
