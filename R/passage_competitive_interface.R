# Statistic-independent competitive inference. A statistic is a deterministic
# function of a gene-index set; it may perform its own pathway-local fitting.
passage_matching_diagnostics <- function(pathway, gene_bins, gene_names = NULL,
                                          reference = c("complement", "universe")) {
  reference <- match.arg(reference)
  if (!is.data.frame(gene_bins) || !all(c("gene", "bin") %in% names(gene_bins)))
    stop("gene_bins must contain gene and bin columns")
  if (is.null(gene_names)) gene_names <- as.character(gene_bins$gene)
  if (anyNA(gene_names) || anyDuplicated(gene_names) || any(!nzchar(gene_names)) ||
      !identical(as.character(gene_bins$gene), as.character(gene_names)))
    stop("gene_bins$gene must exactly match unique gene_names")
  if (is.numeric(pathway) && (anyNA(pathway) || any(!is.finite(pathway)) ||
      any(pathway != as.integer(pathway)) || any(pathway < 1 | pathway > length(gene_names))))
    stop("pathway indices must be integers within the gene universe")
  P <- passage_resolve_pathway(pathway, gene_names)
  if (!length(P)) stop("pathway has no genes in the supplied universe")
  strata <- as.character(gene_bins$bin)
  if ("module" %in% names(gene_bins)) {
    if (anyNA(gene_bins$module)) stop("module matching labels contain missing values")
    strata <- paste(strata, gene_bins$module, sep = "|")
  }
  if (anyNA(gene_bins$bin) || anyNA(strata) || any(!nzchar(strata)))
    stop("matching strata must be nonempty and nonmissing")
  required <- table(strata[P])
  pool <- if (reference == "complement") setdiff(seq_along(strata), P) else seq_along(strata)
  available <- table(factor(strata[pool], levels = names(required)))
  support <- data.frame(stratum = names(required), required = as.integer(required),
    available = as.integer(available), stringsAsFactors = FALSE)
  support$supported <- support$available >= support$required
  log_set_count <- sum(lchoose(support$available, support$required))
  list(supported = all(support$supported), strata = support,
    pathway = P, gene_strata = strata, n_genes = length(P),
    universe_size = length(gene_names), reference = reference,
    log_reference_set_count = log_set_count,
    reference_has_variation = is.finite(log_set_count) && log_set_count > 0,
    matching = if ("module" %in% names(gene_bins))
      "unique_bin_and_module" else "unique_bin")
}

# This helper also permits independent scoring engines to supply a score matrix.
# Row 1 is observed. Every later row is one shared background gene set.
passage_calibrate_competitive_scores <- function(scores, combine = "none",
                                                alternative = "greater",
                                                mc_intervals = FALSE,
                                                conf_level = 0.95) {
  scores <- as.matrix(scores)
  if (!is.numeric(scores) || nrow(scores) < 2L || ncol(scores) < 1L ||
      any(!is.finite(scores))) stop("scores must be a finite numeric matrix with at least two rows")
  if (is.null(colnames(scores))) colnames(scores) <- paste0("statistic", seq_len(ncol(scores)))
  nm <- colnames(scores)
  if (anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm)) stop("statistic names must be unique and nonempty")
  combine <- unique(as.character(combine))
  if (!length(combine) || any(!combine %in% c("none", "cauchy", "minp", "fisher")) ||
      ("none" %in% combine && length(combine) != 1L)) stop("choose none, cauchy, minp and/or fisher")
  if (length(alternative) == 1L) alternative <- rep(alternative, ncol(scores)) else {
    if (!is.null(names(alternative))) {
      if (anyDuplicated(names(alternative)) || !setequal(names(alternative), nm))
        stop("named alternatives must match every statistic")
      alternative <- alternative[nm]
    }
  }
  if (length(alternative) != ncol(scores) || anyNA(alternative) ||
      any(!alternative %in% c("greater", "less"))) stop("alternative must be greater/less for each statistic")
  if (length(conf_level) != 1L || !is.finite(conf_level) || conf_level <= 0 || conf_level >= 1)
    stop("conf_level must lie between zero and one")
  oriented <- sweep(scores, 2L, ifelse(alternative == "greater", 1, -1), "*")
  B <- nrow(scores) - 1L
  counts <- colSums(sweep(oriented[-1L, , drop = FALSE], 2L, oriented[1L, ], ">="))
  intervals <- matrix(NA_real_, ncol(scores), 2L)
  if (isTRUE(mc_intervals)) for (j in seq_len(ncol(scores)))
    intervals[j, ] <- stats::binom.test(counts[j], B, conf.level = conf_level)$conf.int
  components <- data.frame(statistic = nm, observed = scores[1L, ],
    alternative = as.character(alternative), exceedances = as.integer(counts), B = B,
    p_competitive = (1 + counts) / (B + 1), mc_lower = intervals[, 1L],
    mc_upper = intervals[, 2L], mc_floor = 1 / (B + 1), at_floor = counts == 0,
    row.names = NULL, stringsAsFactors = FALSE)
  # All rows, including observed, receive exactly the same transformation.
  # Ties use the conservative upper-tail rank. No independent-p approximation.
  rank_p <- matrix(NA_real_, nrow(scores), ncol(scores), dimnames = dimnames(scores))
  for (j in seq_len(ncol(scores)))
    rank_p[, j] <- rank(-oriented[, j], ties.method = "max") / nrow(scores)
  combined_scores <- matrix(numeric(0), nrow(scores), 0L)
  combinations <- data.frame(method = character(), observed = numeric(),
    exceedances = integer(), B = integer(), p_competitive = numeric(),
    mc_lower = numeric(), mc_upper = numeric(), mc_floor = numeric(), at_floor = logical())
  if (!identical(combine, "none")) {
    combined_scores <- vapply(combine, function(method) {
      switch(method,
        cauchy = {
          # Finite endpoint convention is symmetric across observed/null rows.
          eps <- 0.5 / nrow(scores)
          pp <- pmin(pmax(rank_p, eps), 1 - eps)
          rowMeans(tan(pi * (0.5 - pp)))
        },
        minp = -apply(rank_p, 1L, min),
        fisher = -2 * rowSums(log(rank_p)))
    }, numeric(nrow(scores)))
    colnames(combined_scores) <- combine
    cc <- colSums(sweep(combined_scores[-1L, , drop = FALSE], 2L, combined_scores[1L, ], ">="))
    combinations <- data.frame(method = combine, observed = combined_scores[1L, ],
      exceedances = as.integer(cc), B = B, p_competitive = (1 + cc) / (B + 1),
      mc_lower = NA_real_, mc_upper = NA_real_, mc_floor = 1 / (B + 1),
      at_floor = cc == 0, row.names = NULL)
  }
  list(components = components, combinations = combinations, scores = scores,
    rank_p = rank_p, combined_scores = combined_scores,
    calibration = list(method = "shared_background_empirical_tail",
      combination = "symmetric pooled ranks followed by empirical tail",
      assumptions = "target and matched-set score distributions must be comparable; study-specific calibration required",
      mc_intervals = if (isTRUE(mc_intervals)) "exact binomial intervals for individual conditional background-tail probabilities" else "not available for externally supplied sets",
      combination_mc_intervals = "not binomial: pooled ranks share all background rows",
      conf_level = conf_level,
      selection = "choose the primary statistic/combiner before analysis; returning several does not correct selection across combiners"))
}

passage_competitive_test <- function(pathway, statistic, gene_bins = NULL,
                                     gene_names = NULL, B = 9999L, seed = 1L,
                                     matched_sets = NULL, combine = "none",
                                     alternative = "greater", conf_level = 0.95,
                                     keep_matched_sets = TRUE,
                                     reference = c("complement", "universe")) {
  reference <- match.arg(reference)
  if (!is.function(statistic)) stop("statistic must be a function of gene indices")
  if (is.null(gene_names) && !is.null(gene_bins)) gene_names <- as.character(gene_bins$gene)
  if (is.null(gene_names)) stop("supply gene_names or gene_bins to define the background universe")
  if (is.null(gene_bins)) gene_bins <- data.frame(gene = gene_names, bin = "all")
  support <- passage_matching_diagnostics(pathway, gene_bins, gene_names, reference)
  P <- support$pathway
  if (!support$supported) stop("insufficient unique background genes in matching strata; inspect passage_matching_diagnostics() and prespecify coarser bins")
  old <- if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) get(".Random.seed", .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  if (!is.null(seed)) set.seed(seed)
  generated <- is.null(matched_sets)
  if (generated) matched_sets <- passage_sample_unique_sets(P, support$gene_strata, B, seed,
                                                          reference = reference)
  if (!is.list(matched_sets) || !length(matched_sets)) stop("matched_sets must be a nonempty list")
  # External sets must obey the same declared comparison; no silent replacement.
  for (b in seq_along(matched_sets)) {
    S <- matched_sets[[b]]
    if (is.character(S)) S <- match(S, gene_names)
    if (!is.numeric(S) || anyNA(S) || any(!is.finite(S)) || any(S != as.integer(S)) ||
        length(S) != length(P) || anyDuplicated(S) || any(S < 1L | S > length(gene_names)) ||
        (reference == "complement" && any(S %in% P)))
      stop("invalid matched set ", b, ": require unique genes of the same size; complement references must exclude the target")
    S <- as.integer(S)
    if (!identical(sort(support$gene_strata[S]), sort(support$gene_strata[P])))
      stop("matched set ", b, " does not preserve matching strata")
    matched_sets[[b]] <- S
  }
  sets <- c(list(P), matched_sets)
  if (is.function(attr(statistic, "reset_diagnostics"))) attr(statistic, "reset_diagnostics")()
  batch <- attr(statistic, "batch")
  if (is.function(batch)) {
    scores <- batch(sets)
    if (!is.matrix(scores) || nrow(scores) != length(sets)) stop("statistic batch function returned incompatible rows")
  } else {
    first <- statistic(P)
    if (!is.numeric(first) || !length(first) || !is.null(dim(first))) stop("statistic must return a numeric vector")
    if (is.null(names(first)) && length(first) == 1L) names(first) <- "statistic"
    if (is.null(names(first)) || anyNA(names(first)) || any(!nzchar(names(first))) || anyDuplicated(names(first)))
      stop("statistic must return unique named components (or one unnamed number)")
    scores <- matrix(NA_real_, length(sets), length(first), dimnames = list(NULL, names(first)))
    scores[1L, ] <- first
    for (b in seq_along(matched_sets)) {
      value <- statistic(matched_sets[[b]])
      if (is.null(names(value)) && length(value) == 1L) names(value) <- "statistic"
      if (!is.numeric(value) || !is.null(dim(value)) || !identical(names(value), names(first)))
        stop("statistic components changed in background set ", b)
      scores[b + 1L, ] <- value
    }
  }
  rownames(scores) <- c("observed", paste0("background_", seq_along(matched_sets)))
  ans <- passage_calibrate_competitive_scores(scores, combine, alternative,
    mc_intervals = generated, conf_level = conf_level)
  support$gene_strata <- NULL
  support$n_draws <- length(matched_sets)
  support$n_distinct_sets <- length(unique(vapply(matched_sets, function(x) paste(sort(x), collapse = ","), character(1))))
  support$set_source <- if (generated) "independent uniform draws within fixed strata" else "user supplied; independence not assumed"
  support$mean_target_overlap <- mean(vapply(matched_sets, function(S) sum(S %in% P), integer(1)))
  # Matching summaries in actual units, in addition to exact bin composition.
  cols <- intersect(c("mean_expr", "detection", "variance"), names(gene_bins))
  support$balance <- if (length(cols)) do.call(rbind, lapply(cols, function(column) {
    v <- gene_bins[[column]]
    obs <- mean(v[P]); nul <- vapply(matched_sets, function(S) mean(v[S]), numeric(1))
    sd_all <- stats::sd(v)
    data.frame(characteristic = column, observed = obs, background_mean = mean(nul),
      background_sd = stats::sd(nul), difference_in_gene_sd = if (is.finite(sd_all) && sd_all > 0) (obs - mean(nul)) / sd_all else NA_real_)
  })) else data.frame()
  ans$matching <- support
  ans$matched_sets <- if (isTRUE(keep_matched_sets)) matched_sets else NULL
  if (is.function(attr(statistic, "diagnostics"))) ans$estimation <- attr(statistic, "diagnostics")()
  ans$calibration$reference <- reference
  ans$calibration$assumptions <- if (reference == "universe")
    paste("Conditional gene-label randomization: the target must be uniform over all sets with the declared stratum counts,",
          "conditional on expression, design, strata and any shared fit. Every set is scored identically.",
          "This assumption is not guaranteed for curated correlated biological modules.") else
    "Target-excluded matched-set comparison; finite-universe exchangeability is not guaranteed; empirical calibration required."
  ans$configuration <- list(B = length(matched_sets), seed = seed, reference = reference,
    statistic = attr(statistic, "configuration"), combine = combine,
    calibration_status = "study-specific calibration required")
  class(ans) <- c("passage_competitive_result", "list")
  ans
}
