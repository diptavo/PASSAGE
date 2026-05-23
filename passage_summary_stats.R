# ============================================================================
# passage_summary_stats.R
#
# PASSAGE: Layer 5 / Summary-statistics fallback.
#
# Given only per-gene fitted parameters from a Vecchia-GP fit (e.g. nnSVG
# output: sigma_j^2, phi_j, tau_j^2) plus a gene-correlation matrix from the
# pathway, derives an approximate pathway-level H1 p-value.
#
# This is the analogue of MAGMA / S-PrediXcan-style summary-statistics
# methods in genetics: it lets researchers apply PASSAGE without re-fitting
# the joint multivariate engine, using only pre-computed gene-level outputs
# from existing tools.
#
# Approach:
#   For each gene j in pathway P, nnSVG's score statistic against the no-
#   spatial null is approximately
#       Z_j = (Y_j^T R(phi_j) Y_j - tau_j^2 N) / sqrt(2 * tau_j^4 * sum_eigvals)
#   under the null Y_j ~ N(0, tau_j^2 I).  Pathway aggregation:
#       T = sum_j w_j Z_j   (weighted by gene importance)
#   with a Cauchy / Liu-style p-value combining via the gene correlation
#   matrix.
#
# Limitations: this is an approximation; the joint model captures cross-gene
# spatial coherence which the summary-statistics version cannot.  Use as a
# fast deployable fallback when the joint engine is impractical.
#
# Dependencies: CompQuadForm (Davies), passage_score_h1.R (acat_combine)
# ============================================================================


# ----------------------------------------------------------------------------
# Build a pathway-level summary p-value from per-gene Vecchia-GP outputs.
#
# Inputs
#   gene_stats : data.frame with columns:
#                  gene       (character; gene name)
#                  Z          (numeric; per-gene standardised spatial score
#                              statistic, e.g. nnSVG's log-likelihood ratio)
#                  p_value    (numeric; per-gene p-value)
#                  propSV     (numeric; nnSVG-style per-gene spatial fraction)
#                              -- optional
#   pathway     : character vector of gene names defining the pathway
#   gene_cor    : optional G x G gene-gene correlation matrix used for
#                 dependence-aware combination. If NULL, assumes independence
#                 (anticonservative for highly correlated pathways).
#   combine     : "acat" (default, Cauchy combination, dependence-robust),
#                 "stouffer" (sqrt-N weighted Z), "fisher" (chi^2 combination,
#                 dependence-naive).
#
# Returns a list with:
#   p_pathway   : combined p-value
#   T_stat      : underlying combination statistic
#   n_genes_in  : number of pathway genes found in gene_stats
#   method      : combination method used
# ----------------------------------------------------------------------------
passage_summary_h1 <- function(gene_stats, pathway,
                               gene_cor = NULL,
                               combine = c("acat", "stouffer", "fisher"),
                               weights = NULL) {
  combine <- match.arg(combine)
  stopifnot(all(c("gene", "p_value") %in% names(gene_stats)))

  in_path <- gene_stats$gene %in% pathway
  if (sum(in_path) < 2L) {
    return(list(p_pathway = NA_real_, n_genes_in = sum(in_path),
                method = combine,
                note = "Fewer than 2 pathway genes in summary stats"))
  }

  pstats <- gene_stats[in_path, , drop = FALSE]
  p_vals <- pmin(pmax(pstats$p_value, .Machine$double.eps), 1 - 1e-15)

  if (is.null(weights)) {
    weights <- rep(1, length(p_vals))
    if ("propSV" %in% names(pstats)) {
      weights <- pmax(pstats$propSV, 1e-3)  # weight by spatial fraction
    }
  }

  T_stat <- NA_real_
  p_pathway <- NA_real_

  if (combine == "acat") {
    p_pathway <- acat_combine(p_vals, weights = weights)
    T_stat <- sum(weights * tan((0.5 - p_vals) * pi)) / sum(weights)
  } else if (combine == "stouffer") {
    # Stouffer: T = sum(w_j Z_j) / sqrt(sum(w_j^2))
    # If correlated: variance inflation factor = sum_{j,j'} w_j w_j' cor_jj'
    Z <- qnorm(1 - p_vals)
    if (is.null(gene_cor)) {
      var_T <- sum(weights^2)
    } else {
      pw_genes_in_cor <- pstats$gene[pstats$gene %in% rownames(gene_cor)]
      sub_cor <- gene_cor[pw_genes_in_cor, pw_genes_in_cor, drop = FALSE]
      idx_match <- match(pw_genes_in_cor, pstats$gene)
      w_sub <- weights[idx_match]
      var_T <- as.numeric(t(w_sub) %*% sub_cor %*% w_sub)
    }
    T_stat <- sum(weights * Z) / sqrt(max(var_T, 1e-12))
    p_pathway <- pnorm(T_stat, lower.tail = FALSE)
  } else if (combine == "fisher") {
    # Fisher's combined: -2 sum log(p) ~ chi^2(2k) under independence
    T_stat <- -2 * sum(log(p_vals))
    p_pathway <- pchisq(T_stat, df = 2 * length(p_vals), lower.tail = FALSE)
  }

  list(
    p_pathway = p_pathway,
    T_stat = T_stat,
    n_genes_in = nrow(pstats),
    method = combine,
    weights_used = weights,
    pathway_genes_in_stats = pstats$gene
  )
}


# ----------------------------------------------------------------------------
# Bulk runner: apply passage_summary_h1 to a list of pathways and assemble
# into a tidy data.frame.  Mirrors the run_passage interface but is much
# faster (no engine fit needed).
# ----------------------------------------------------------------------------
run_passage_summary <- function(gene_stats, pathways,
                                gene_cor = NULL,
                                combine = "acat",
                                min_pathway_size = 5L,
                                max_pathway_size = 500L,
                                adjust_method = "BH",
                                verbose = TRUE) {
  if (is.null(names(pathways))) {
    names(pathways) <- paste0("pw_", seq_along(pathways))
  }
  sizes <- vapply(pathways, length, integer(1))
  keep <- sizes >= min_pathway_size & sizes <= max_pathway_size
  pathways <- pathways[keep]
  if (verbose) {
    cat(sprintf("run_passage_summary: scoring %d pathways (filtered from %d)\n",
                length(pathways), length(keep)))
  }

  rows <- list()
  for (nm in names(pathways)) {
    res <- passage_summary_h1(gene_stats, pathways[[nm]],
                              gene_cor = gene_cor, combine = combine)
    rows[[nm]] <- data.frame(
      pathway = nm,
      pathway_size = length(pathways[[nm]]),
      n_genes_in = res$n_genes_in,
      T_stat = res$T_stat,
      p_pathway = res$p_pathway,
      stringsAsFactors = FALSE
    )
  }
  tbl <- do.call(rbind, rows)
  if (!is.null(adjust_method) && nrow(tbl) > 0L) {
    tbl$p_adj <- p.adjust(tbl$p_pathway, method = adjust_method)
  }
  tbl <- tbl[order(tbl$p_pathway), , drop = FALSE]
  rownames(tbl) <- NULL
  tbl
}


# ----------------------------------------------------------------------------
# Convenience: convert nnSVG output to the gene_stats format expected here.
# ----------------------------------------------------------------------------
nnSVG_to_gene_stats <- function(nnSVG_results,
                                gene_col = "gene_name",
                                pval_col = "padj",
                                propSV_col = "prop_sv") {
  out <- data.frame(
    gene = nnSVG_results[[gene_col]],
    p_value = nnSVG_results[[pval_col]],
    stringsAsFactors = FALSE
  )
  if (propSV_col %in% names(nnSVG_results)) {
    out$propSV <- nnSVG_results[[propSV_col]]
  }
  # Compute Z from p-value
  out$Z <- qnorm(1 - pmin(pmax(out$p_value, 1e-300), 1 - 1e-15))
  out
}
