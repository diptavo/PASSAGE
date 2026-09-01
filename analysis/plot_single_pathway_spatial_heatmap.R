# Render spatial heatmaps for one pathway from a saved PASSAGE conditional run.

parse_args <- function(args) {
  cfg <- list(
    result_dir = NULL,
    pathway = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    condition = "H3",
    out = NULL
  )
  for (arg in args) {
    if (grepl("^--result-dir=", arg)) cfg$result_dir <- sub("^--result-dir=", "", arg)
    if (grepl("^--pathway=", arg)) cfg$pathway <- sub("^--pathway=", "", arg)
    if (grepl("^--condition=", arg)) cfg$condition <- sub("^--condition=", "", arg)
    if (grepl("^--out=", arg)) cfg$out <- sub("^--out=", "", arg)
  }
  if (is.null(cfg$result_dir)) stop("--result-dir is required")
  if (is.null(cfg$out)) {
    cfg$out <- file.path(cfg$result_dir, paste0("spatial_heatmap_", cfg$pathway, "_", cfg$condition, ".png"))
  }
  cfg
}

zscale <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x)
  if (!is.finite(s) || s <= 1e-10) return(x - mean(x))
  (x - mean(x)) / s
}

cap_score <- function(x, lo = 0.02, hi = 0.98) {
  q <- stats::quantile(x, c(lo, hi), na.rm = TRUE)
  pmax(pmin(x, q[[2L]]), q[[1L]])
}

score_to_col <- function(x, pal) {
  z <- cap_score(x)
  rng <- range(z, finite = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) {
    return(rep(pal[[ceiling(length(pal) / 2)]], length(x)))
  }
  ii <- round(1 + (z - rng[[1L]]) / diff(rng) * (length(pal) - 1L))
  pal[pmax(1L, pmin(length(pal), ii))]
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
fit <- readRDS(file.path(cfg$result_dir, "passage_msigdb_conditional_result.rds"))
dat <- readRDS(file.path(cfg$result_dir, "passage_msigdb_conditional_prepared_data.rds"))
if (!cfg$pathway %in% names(dat$pathways)) stop("Pathway not found: ", cfg$pathway)
genes <- dat$pathways[[cfg$pathway]]
idx <- passage_resolve_pathway(genes, colnames(dat$Y))
if (length(idx) < 2L) stop("Fewer than 2 pathway genes present")
pre <- fit$precomp[[cfg$condition]]
if (is.null(pre)) stop("Condition not found in fit$precomp: ", cfg$condition)

Yp <- dat$Y[, idx, drop = FALSE]
raw_score <- rowMeans(scale(Yp))
resid <- passage_residualize_with_qr(Yp, pre$qr)
resid_score <- rowMeans(scale(resid))
fitted <- fit$engine$V %*% t(fit$engine$A[idx, , drop = FALSE])
fitted <- passage_residualize_with_qr(fitted, pre$qr)
fitted_score <- rowMeans(scale(fitted))

coords <- as.data.frame(dat$coords)
if (ncol(coords) < 2L) stop("dat$coords must have at least two columns")
x <- coords[[1L]]
y <- coords[[2L]]
y_plot <- max(y, na.rm = TRUE) - y

comp_path <- file.path(cfg$result_dir, paste0("passage_competitive_metrics_", cfg$condition, ".csv"))
comp <- if (file.exists(comp_path)) read.csv(comp_path, stringsAsFactors = FALSE) else NULL
comp_row <- if (!is.null(comp)) comp[match(cfg$pathway, comp$pathway), , drop = FALSE] else NULL
base_path <- file.path(cfg$result_dir, "passage_msigdb_conditional_pathways.csv")
base <- if (file.exists(base_path)) read.csv(base_path, stringsAsFactors = FALSE) else NULL
base_row <- if (!is.null(base)) base[match(cfg$pathway, base$pathway), , drop = FALSE] else NULL

subtitle <- paste0(
  "genes=", length(idx),
  if (!is.null(base_row) && nrow(base_row) == 1L) {
    p_col <- paste0("p_", cfg$condition)
    fdr_col <- paste0("fdr_", cfg$condition)
    if (all(c(p_col, fdr_col) %in% names(base_row))) {
      paste0("; ", cfg$condition, " p=", signif(base_row[[p_col]], 3), "; ", cfg$condition, " FDR=", signif(base_row[[fdr_col]], 3))
    } else {
      ""
    }
  } else "",
  if (!is.null(comp_row) && nrow(comp_row) == 1L) paste0("; cEPSV=", signif(comp_row$cEPSV, 3), "; competitive cEPSV FDR=", signif(comp_row$competitive_cEPSV_fdr_perm, 3)) else ""
)

has_fitted <- any(is.finite(fitted_score)) && stats::sd(fitted_score, na.rm = TRUE) > 1e-10
png(cfg$out, width = if (has_fitted) 2100 else 1500, height = 760, res = 150)
layout(matrix(seq_len(if (has_fitted) 3L else 2L), nrow = 1), widths = rep(1, if (has_fitted) 3L else 2L))
pal <- hcl.colors(101, "Blue-Red 3")
par(mar = c(4, 4, 4, 1), oma = c(0, 0, 3, 0))
plot_panel <- function(score, title) {
  plot(x, y_plot, asp = 1, pch = 16, cex = 0.75, col = score_to_col(score, pal),
       xlab = "spatial x", ylab = "spatial y", main = title)
  box()
}
plot_panel(raw_score, "Mean z-scored pathway expression")
plot_panel(resid_score, paste0(cfg$condition, " residual pathway score"))
if (has_fitted) plot_panel(fitted_score, paste0(cfg$condition, " fitted spatial component"))
mtext(cfg$pathway, outer = TRUE, line = 1.5, cex = 1.1, font = 2)
mtext(subtitle, outer = TRUE, line = 0.2, cex = 0.8)
dev.off()

write.csv(
  data.frame(
    barcode = dat$spot_metadata$barcode,
    x = x,
    y = y,
    raw_pathway_score = raw_score,
    residual_pathway_score = resid_score,
    fitted_spatial_score = fitted_score
  ),
  sub("[.]png$", ".csv", cfg$out),
  row.names = FALSE
)
message("Wrote ", cfg$out)
