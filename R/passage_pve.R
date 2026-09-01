# PASSAGE Layer 3: pathway-level effect sizes.

passage_pve <- function(engine,
                        Y,
                        pathway,
                        gene_names = NULL,
                        compute = c("cca", "range", "meangene")) {
  Y <- passage_check_y(Y)
  if (is.null(gene_names)) {
    gene_names <- colnames(Y)
  }
  P <- passage_resolve_pathway(pathway, gene_names)
  if (length(P) == 0L) {
    stop("pathway has no genes present in Y")
  }
  out <- list(pathway_size = length(P))
  if ("cca" %in% compute) {
    out$cca <- passage_pve_cca(engine, Y, P)
  }
  if ("range" %in% compute) {
    out$range <- passage_pve_range(engine, P)
  }
  if ("meangene" %in% compute) {
    out$meangene <- passage_pve_meangene(engine, P)
  }
  out$summary <- c(
    R2_cca = if (!is.null(out$cca)) out$cca$R2_cca else NA_real_,
    PSVS_range = if (!is.null(out$range)) out$range$PSVS_range else NA_real_,
    mean_propSV = if (!is.null(out$meangene)) out$meangene$mean_propSV else NA_real_
  )
  class(out) <- c("passage_pve_result", "list")
  out
}

passage_pve_cca <- function(engine, Y, P) {
  YP <- engine$residuals[, P, drop = FALSE]
  Yhat <- engine$V %*% t(engine$A[P, , drop = FALSE])
  YP <- sweep(YP, 2L, colMeans(YP), "-")
  Yhat <- sweep(Yhat, 2L, colMeans(Yhat), "-")
  ok_y <- apply(YP, 2L, stats::sd) > 1e-10
  ok_h <- apply(Yhat, 2L, stats::sd) > 1e-10
  if (sum(ok_y) < 2L || sum(ok_h) < 2L) {
    return(list(R2_cca = NA_real_, rho_squared = numeric(0)))
  }
  cc <- tryCatch(stats::cancor(YP[, ok_y, drop = FALSE], Yhat[, ok_h, drop = FALSE]),
                 error = function(e) NULL)
  if (is.null(cc)) {
    return(list(R2_cca = NA_real_, rho_squared = numeric(0)))
  }
  rho2 <- pmin(pmax(cc$cor, 0), 1)^2
  list(R2_cca = mean(rho2), rho_squared = rho2)
}

passage_pve_range <- function(engine, P, tissue_diameter = NULL) {
  A_P <- engine$A[P, , drop = FALSE]
  if (is.null(tissue_diameter)) {
    bbox <- apply(engine$coords, 2L, range)
    tissue_diameter <- sqrt(sum((bbox[2L, ] - bbox[1L, ])^2))
  }
  rel <- engine$theta$effective_range / max(tissue_diameter, .Machine$double.eps)
  g <- pmin(pmax(rel, 0), 1)
  spatial <- sum(g * engine$theta$sigma2 * colSums(A_P^2))
  den <- spatial + sum(engine$D[P])
  list(
    PSVS_range = spatial / max(den, .Machine$double.eps),
    factor_contributions = g * engine$theta$sigma2 * colSums(A_P^2),
    tissue_diameter = tissue_diameter
  )
}

passage_pve_meangene <- function(engine, P) {
  A_P <- engine$A[P, , drop = FALSE]
  spatial_gene <- as.numeric((A_P^2) %*% engine$theta$sigma2)
  prop <- spatial_gene / pmax(spatial_gene + engine$D[P], .Machine$double.eps)
  list(mean_propSV = mean(prop), median_propSV = stats::median(prop), propSV = prop)
}
