# Optional descriptive GP fits. These do not replace the unsmoothed scores
# consumed by the pathway tests and do not supply calibrated test p-values.
passage_fit_pathway_gp <- function(fit, coords,
                                   covfun_name = "matern15_isotropic",
                                   m_seq = c(10L, 30L), max_iter = 40L,
                                   seed = 1L) {
  if (!inherits(fit, "passage_pathway_fit") ||
      !identical(fit$scores_source, "unsmoothed_expression")) {
    stop("fit must be a passage_pathway_fit with unsmoothed expression scores")
  }
  if (length(fit$method) != 1L || !fit$method %in% c("pca", "factor_analysis")) {
    stop("unsupported pathway estimation method")
  }
  scores <- as.matrix(fit$scores)
  X <- as.matrix(fit$X)
  if (!is.numeric(scores) || !is.numeric(X) || !all(is.finite(scores)) ||
      !all(is.finite(X)) || nrow(scores) < 5L || ncol(scores) < 1L ||
      nrow(X) != nrow(scores) || ncol(X) < 1L) {
    stop("fit needs finite numeric scores and an aligned adjustment design")
  }
  if (length(fit$K) != 1L || !is.finite(fit$K) || fit$K != ncol(scores)) {
    stop("fit$K must equal the number of score columns")
  }
  valid_ids <- function(x, expected) {
    !is.null(x) && length(x) == expected && !anyNA(x) &&
      all(nzchar(x)) && !anyDuplicated(x)
  }
  ids <- rownames(scores)
  components <- colnames(scores)
  if (!valid_ids(ids, nrow(scores)) || !valid_ids(components, ncol(scores)) ||
      !identical(rownames(X), ids)) {
    stop("scores need unique spot/component names; X must have exactly aligned spot IDs")
  }
  coords <- passage_check_coords(coords, nrow(scores))
  if (!identical(rownames(coords), ids)) stop("coordinate spot IDs are not aligned with the scores")
  if (anyDuplicated(as.data.frame(coords))) stop("duplicate spatial coordinates")
  qx <- qr(X)
  if (qx$rank != ncol(X) || nrow(scores) - qx$rank < 3L) {
    stop("GP adjustment design must have full rank and at least three residual degrees of freedom")
  }
  if (any(colSums(qr.resid(qx, scores)^2) <= .Machine$double.eps)) {
    stop("each component must have nonzero variation after adjustment")
  }
  if (!is.character(covfun_name) || length(covfun_name) != 1L ||
      is.na(covfun_name) || !nzchar(covfun_name)) stop("invalid covariance-function name")
  positive_integer <- function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x)) &&
      all(x >= 1 & x <= .Machine$integer.max & x == floor(x))
  }
  if (!positive_integer(m_seq) || is.unsorted(m_seq)) {
    stop("m_seq must be a nondecreasing sequence of positive integers")
  }
  if (length(max_iter) != 1L || !positive_integer(max_iter)) {
    stop("max_iter must be a positive integer")
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("seed must be a nonnegative integer")
  }
  m_used <- unique(pmin(as.integer(m_seq), nrow(scores) - 1L))
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    } else assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }, add = TRUE)
  package_warnings <- character()
  available <- withCallingHandlers(requireNamespace("GpGp", quietly = TRUE),
    warning = function(w) {
      package_warnings <<- c(package_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!available) stop("Optional package GpGp is required for descriptive GP fitting")
  set.seed(as.integer(seed))
  gp_fits <- warnings <- vector("list", ncol(scores))
  names(gp_fits) <- names(warnings) <- components
  rows <- parameters <- vector("list", ncol(scores))
  for (k in seq_len(ncol(scores))) {
    component_warnings <- character()
    gp <- tryCatch(withCallingHandlers(
      GpGp::fit_model(y = scores[, k], locs = coords, X = X,
        covfun_name = covfun_name, m_seq = m_used,
        max_iter = as.integer(max_iter), silent = TRUE),
      warning = function(w) {
        component_warnings <<- c(component_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) {
        stop("GpGp fit failed for component ", components[k], ": ", conditionMessage(e), call. = FALSE)
      })
    cp <- as.numeric(gp$covparms)
    if (!length(cp) || any(!is.finite(cp))) {
      stop("GpGp returned nonfinite covariance parameters for ", components[k])
    }
    converged <- if (is.logical(gp$conv) && length(gp$conv) == 1L) gp$conv else NA
    loglik <- if (is.numeric(gp$loglik) && length(gp$loglik) == 1L) gp$loglik else NA_real_
    default_matern <- identical(covfun_name, "matern15_isotropic")
    if (default_matern && (length(cp) != 3L || any(cp[1:2] <= 0) || cp[3] < 0)) {
      stop("invalid Matern 3/2 covariance parameters from GpGp for ", components[k])
    }
    parameter_names <- if (default_matern) c("variance", "alpha", "relative_nugget") else names(gp$covparms)
    if (is.null(parameter_names) || length(parameter_names) != length(cp)) {
      parameter_names <- paste0("parameter_", seq_along(cp))
    }
    parameters[[k]] <- data.frame(component = components[k],
      parameter = parameter_names, estimate = cp, stringsAsFactors = FALSE)
    rows[[k]] <- data.frame(component = components[k], converged = converged,
      loglik = loglik, n_warnings = length(component_warnings),
      spatial_variance = if (default_matern) cp[1L] else NA_real_,
      range_gp_alpha = if (default_matern) cp[2L] else NA_real_,
      range_passage_rho = if (default_matern) sqrt(3) * cp[2L] else NA_real_,
      relative_nugget = if (default_matern) cp[3L] else NA_real_,
      noise_variance = if (default_matern) cp[1L] * cp[3L] else NA_real_,
      stringsAsFactors = FALSE)
    gp_fits[[k]] <- gp
    warnings[[k]] <- component_warnings
  }
  summary <- do.call(rbind, rows)
  rownames(summary) <- NULL
  list(gp_fits = gp_fits, summary = summary,
    covariance_parameters = do.call(rbind, parameters), warnings = warnings,
    package_warnings = package_warnings,
    configuration = list(interface_version = "1.0.0",
      GpGp_version = as.character(utils::packageVersion("GpGp")),
      estimation_method = fit$method, covfun_name = covfun_name,
      m_seq_requested = m_seq, m_seq_used = m_used,
      max_iter = as.integer(max_iter), seed = as.integer(seed),
      coordinates = "supplied coordinate units; no automatic rescaling",
      covariance_parameterization = if (identical(covfun_name, "matern15_isotropic")) {
        "sigma2 * (1 + distance/alpha) * exp(-distance/alpha); diagonal noise sigma2 * relative_nugget; PASSAGE rho = sqrt(3) * alpha"
      } else "GpGp covariance-function parameterization; inspect covariance_parameters and GpGp documentation"),
    scores_source = "unsmoothed_expression",
    role = "Descriptive plug-in per-component GP fits; no calibrated p-values or predictions are returned",
    inference_note = paste(
      "Do not substitute GP-smoothed predictions or learned covariance kernels into the fixed-kernel test without recalibration.",
      "PCA/FA and GP parameter uncertainty is not propagated; returned package standard errors are not full-pipeline uncertainty.",
      "Inspect convergence and warnings for every component before interpreting fits."))
}
