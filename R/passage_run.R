# Top-level PASSAGE workflow.

passage_run <- function(Y,
                        coords,
                        pathways,
                        Z_CT = NULL,
                        X = NULL,
                        K = 6L,
                        m = 20L,
                        range_grid = NULL,
                        kernel = c("matern32", "exponential", "gaussian"),
                        hypotheses = c("H1", "H2"),
                        weight_schemes = c("equal", "var", "range"),
                        min_pathway_size = 2L,
                        max_pathway_size = 500L,
                        fdr_method = "BH",
                        calibration = c("permutation", "moment"),
                        n_perm = 199L,
                        seed = NULL,
                        verbose = TRUE,
                        anti_dip = c("none", "spot_crossfit", "pathway_holdout", "refit_null"),
                        anti_dip_args = list()) {
  kernel <- match.arg(kernel)
  calibration <- match.arg(calibration)
  anti_dip <- match.arg(anti_dip)
  Y <- passage_check_y(Y)
  coords <- passage_check_coords(coords, nrow(Y))
  pathways <- passage_check_pathways(pathways, colnames(Y))

  engine_fit_args <- list(
    K = K,
    m = m,
    range_grid = range_grid,
    kernel = kernel,
    verbose = verbose
  )
  engine <- NULL
  if (anti_dip == "none") {
    engine <- do.call(passage_fit_engine_pca, c(
      list(Y = Y, coords = coords, X = X),
      engine_fit_args
    ))
  } else if (anti_dip == "spot_crossfit") {
    cf_args <- anti_dip_args
    if (is.null(cf_args$fit_args)) cf_args$fit_args <- engine_fit_args
    if (is.null(cf_args$method)) cf_args$method <- "pca"
    if (is.null(cf_args$verbose)) cf_args$verbose <- verbose
    engine <- do.call(passage_fit_engine_crossfit_spots, c(
      list(Y = Y, coords = coords, X = X),
      cf_args
    ))
  }

  if (verbose && !is.null(engine)) {
    message("PASSAGE: precomputing H1/H2 score moments")
  }
  precomp <- list()
  if (!is.null(engine) && "H1" %in% hypotheses) {
    precomp$H1 <- passage_h_precompute(engine, X = passage_prepare_design(X, nrow(Y)))
  }
  if (!is.null(engine) && "H2" %in% hypotheses && !is.null(Z_CT)) {
    X2 <- if (is.null(X)) {
      Z_CT
    } else {
      cbind(as.matrix(X), as.matrix(Z_CT))
    }
    precomp$H2 <- passage_h_precompute(engine, X = passage_prepare_design(X2, nrow(Y)))
  }
  if ("H2" %in% hypotheses && is.null(Z_CT)) {
    warning("H2 requested but Z_CT is NULL; skipping H2")
    hypotheses <- setdiff(hypotheses, "H2")
  }

  if (verbose) {
    message("PASSAGE: scoring pathways")
  }
  rows <- vector("list", length(pathways))
  per_pathway <- vector("list", length(pathways))
  names(per_pathway) <- names(pathways)
  for (ii in seq_along(pathways)) {
    pname <- names(pathways)[ii]
    genes <- pathways[[ii]]
    q <- length(genes)
    if (q < min_pathway_size || q > max_pathway_size) {
      rows[[ii]] <- passage_empty_row(pname, q, "size_filter")
      next
    }
    engine_i <- engine
    precomp_i <- precomp
    if (anti_dip == "pathway_holdout") {
      hold_args <- anti_dip_args
      if (is.null(hold_args$fit_args)) hold_args$fit_args <- engine_fit_args
      if (is.null(hold_args$method)) hold_args$method <- "pca"
      if (is.null(hold_args$verbose)) hold_args$verbose <- verbose
      engine_i <- do.call(passage_fit_engine_pathway_holdout, c(
        list(Y = Y, coords = coords, pathway = genes, X = X, gene_names = colnames(Y)),
        hold_args
      ))
      precomp_i <- list()
      if ("H1" %in% hypotheses) {
        precomp_i$H1 <- passage_h_precompute(engine_i, X = passage_prepare_design(X, nrow(Y)))
      }
      if ("H2" %in% hypotheses) {
        X2 <- if (is.null(X)) Z_CT else cbind(as.matrix(X), as.matrix(Z_CT))
        precomp_i$H2 <- passage_h_precompute(engine_i, X = passage_prepare_design(X2, nrow(Y)))
      }
    }
    hres <- list()
    if ("H1" %in% hypotheses) {
      if (anti_dip == "refit_null") {
        refit_args <- anti_dip_args
        if (is.null(refit_args$fit_args)) refit_args$fit_args <- engine_fit_args
        if (is.null(refit_args$method)) refit_args$method <- "pca"
        if (is.null(refit_args$n_perm)) refit_args$n_perm <- n_perm
        if (is.null(refit_args$return_engine)) refit_args$return_engine <- TRUE
        if (is.null(refit_args$verbose)) refit_args$verbose <- verbose
        hres$H1 <- do.call(passage_score_test_refit_null, c(
          list(
            Y = Y,
            coords = coords,
            pathway = genes,
            X = passage_prepare_design(X, nrow(Y)),
            gene_names = colnames(Y),
            weight_schemes = weight_schemes,
            run_burden = TRUE,
            run_spasset = TRUE,
            seed = if (is.null(seed)) NULL else seed + 1000L * ii + 1L
          ),
          refit_args
        ))
        engine_i <- hres$H1$engine
      } else {
        hres$H1 <- passage_score_test(
          engine_i, Y, genes, precomp_i$H1,
          weight_schemes = weight_schemes,
          gene_names = colnames(Y),
          calibration = calibration,
          n_perm = n_perm,
          seed = if (is.null(seed)) NULL else seed + 1000L * ii + 1L
        )
      }
    }
    if ("H2" %in% hypotheses) {
      if (anti_dip == "refit_null") {
        X2 <- if (is.null(X)) Z_CT else cbind(as.matrix(X), as.matrix(Z_CT))
        refit_args <- anti_dip_args
        if (is.null(refit_args$fit_args)) refit_args$fit_args <- engine_fit_args
        if (is.null(refit_args$method)) refit_args$method <- "pca"
        if (is.null(refit_args$n_perm)) refit_args$n_perm <- n_perm
        if (is.null(refit_args$return_engine)) refit_args$return_engine <- FALSE
        if (is.null(refit_args$verbose)) refit_args$verbose <- verbose
        hres$H2 <- do.call(passage_score_test_refit_null, c(
          list(
            Y = Y,
            coords = coords,
            pathway = genes,
            X = passage_prepare_design(X2, nrow(Y)),
            gene_names = colnames(Y),
            weight_schemes = weight_schemes,
            run_burden = TRUE,
            run_spasset = TRUE,
            seed = if (is.null(seed)) NULL else seed + 1000L * ii + 2L
          ),
          refit_args
        ))
      } else {
        hres$H2 <- passage_score_test(
          engine_i, Y, genes, precomp_i$H2,
          weight_schemes = weight_schemes,
          gene_names = colnames(Y),
          calibration = calibration,
          n_perm = n_perm,
          seed = if (is.null(seed)) NULL else seed + 1000L * ii + 2L
        )
      }
    }
    pve <- passage_pve(engine_i, Y, genes, gene_names = colnames(Y))
    decomp <- passage_decomposition(hres$H1, hres$H2, NULL)
    per_pathway[[pname]] <- list(hypotheses = hres, decomposition = decomp, pve = pve)
    rows[[ii]] <- data.frame(
      pathway = pname,
      pathway_size = q,
      p_H1 = decomp$p_H1,
      p_H2 = decomp$p_H2,
      p_H1_moment = passage_moment_p(hres$H1),
      p_H2_moment = passage_moment_p(hres$H2),
      Q_H1 = decomp$Q_H1,
      Q_H2 = decomp$Q_H2,
      cell_type_share = decomp$cell_type_share,
      R2_cca = pve$summary["R2_cca"],
      PSVS_range = pve$summary["PSVS_range"],
      mean_propSV = pve$summary["mean_propSV"],
      spasset_genes_H1 = I(list(passage_spasset_best_genes(hres$H1))),
      status = "tested",
      stringsAsFactors = FALSE
    )
  }
  tbl <- do.call(rbind, rows)
  if ("p_H1" %in% names(tbl)) {
    tbl$p_H1_adj <- stats::p.adjust(tbl$p_H1, method = fdr_method)
  }
  if ("p_H2" %in% names(tbl)) {
    tbl$p_H2_adj <- stats::p.adjust(tbl$p_H2, method = fdr_method)
  }
  sort_col <- if ("p_H2" %in% names(tbl) && any(is.finite(tbl$p_H2))) "p_H2" else "p_H1"
  tbl <- tbl[order(tbl[[sort_col]], na.last = TRUE), , drop = FALSE]
  rownames(tbl) <- NULL

  out <- list(
    summary = tbl,
    per_pathway = per_pathway,
    engine = engine,
    precomp = precomp,
    hypotheses = hypotheses,
    calibration = calibration,
    anti_dip = anti_dip,
    n_perm = if (calibration == "permutation") n_perm else 0L,
    call = match.call()
  )
  class(out) <- c("passage_run", "list")
  out
}

print.passage_run <- function(x, n = 10L, ...) {
  cat("PASSAGE run\n")
  cat("  hypotheses:", paste(x$hypotheses, collapse = ", "), "\n")
  cat("  pathways:", nrow(x$summary), "\n\n")
  cols <- intersect(
    c("pathway", "pathway_size", "p_H1", "p_H1_adj", "p_H2", "p_H2_adj",
      "cell_type_share", "R2_cca", "PSVS_range", "status"),
    names(x$summary)
  )
  print(utils::head(x$summary[, cols, drop = FALSE], n), row.names = FALSE)
  invisible(x)
}

passage_empty_row <- function(pathway, q, status) {
  data.frame(
    pathway = pathway,
    pathway_size = q,
    p_H1 = NA_real_,
    p_H2 = NA_real_,
    p_H1_moment = NA_real_,
    p_H2_moment = NA_real_,
    Q_H1 = NA_real_,
    Q_H2 = NA_real_,
    cell_type_share = NA_real_,
    R2_cca = NA_real_,
    PSVS_range = NA_real_,
    mean_propSV = NA_real_,
    spasset_genes_H1 = I(list(character())),
    status = status,
    stringsAsFactors = FALSE
  )
}

passage_moment_p <- function(h) {
  if (is.null(h) || is.null(h$p_omnibus_moment)) {
    return(NA_real_)
  }
  h$p_omnibus_moment
}

passage_spasset_best_genes <- function(h) {
  if (is.null(h) || is.null(h$spasset) || !is.list(h$spasset)) {
    return(character())
  }
  if (is.null(h$spasset$best_genes)) {
    return(character())
  }
  as.character(h$spasset$best_genes)
}
