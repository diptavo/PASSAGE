passage_empirical_tail_calibrate <- function(x, p_col,
                                             group_cols,
                                             strata_cols = character(),
                                             leaveout_cols = character(),
                                             out_col = NULL,
                                             n_col = NULL) {
  if (!is.data.frame(x)) stop("x must be a data.frame")
  if (!p_col %in% names(x)) stop("p_col not found in x: ", p_col)
  needed <- unique(c(group_cols, strata_cols, leaveout_cols))
  missing <- setdiff(needed, names(x))
  if (length(missing)) stop("columns not found in x: ", paste(missing, collapse = ", "))
  if (is.null(out_col)) {
    suffix <- if (length(strata_cols)) paste(strata_cols, collapse = "_") else "pooled"
    out_col <- paste0(p_col, "_empirical_", suffix)
    if (length(leaveout_cols)) out_col <- paste0(out_col, "_leaveout")
  }

  p <- as.numeric(x[[p_col]])
  out <- rep(NA_real_, length(p))
  n_ref <- rep(NA_integer_, length(p))
  match_cols <- unique(c(group_cols, strata_cols))

  same_column_value <- function(v, i) {
    out <- v == v[[i]]
    out[is.na(out)] <- is.na(v[is.na(out)]) & is.na(v[[i]])
    out
  }

  for (ii in seq_along(p)) {
    if (!is.finite(p[[ii]])) next
    keep <- rep(TRUE, nrow(x))
    for (cc in match_cols) {
      keep <- keep & same_column_value(x[[cc]], ii)
    }
    if (length(leaveout_cols)) {
      leave_same <- rep(TRUE, nrow(x))
      for (cc in leaveout_cols) {
        leave_same <- leave_same & same_column_value(x[[cc]], ii)
      }
      keep <- keep & !leave_same
    }
    null <- p[keep & is.finite(p)]
    if (!length(null)) next
    n_ref[[ii]] <- length(null)
    out[[ii]] <- (1 + sum(null <= p[[ii]])) / (1 + length(null))
  }
  x[[out_col]] <- out
  if (!is.null(n_col)) x[[n_col]] <- n_ref
  x
}

passage_empirical_competitive_calibration <- function(x,
                                                      p_col = "competitive_score_p",
                                                      group_cols = c("dataset", "mode", "matching", "mc_sampler"),
                                                      size_col = "pathway_size",
                                                      replicate_col = "replicate",
                                                      prefix = "competitive_score") {
  if (!p_col %in% names(x)) stop("p_col not found in x: ", p_col)
  missing_group <- setdiff(group_cols, names(x))
  if (length(missing_group)) {
    stop("group_cols not found in x: ", paste(missing_group, collapse = ", "))
  }

  x <- passage_empirical_tail_calibrate(
    x, p_col = p_col, group_cols = group_cols,
    out_col = paste0(prefix, "_empirical_pooled_p"),
    n_col = paste0(prefix, "_empirical_pooled_n")
  )
  if (replicate_col %in% names(x)) {
    x <- passage_empirical_tail_calibrate(
      x, p_col = p_col, group_cols = group_cols,
      leaveout_cols = replicate_col,
      out_col = paste0(prefix, "_empirical_leave_rep_p"),
      n_col = paste0(prefix, "_empirical_leave_rep_n")
    )
  }
  if (size_col %in% names(x)) {
    x <- passage_empirical_tail_calibrate(
      x, p_col = p_col, group_cols = group_cols,
      strata_cols = size_col,
      out_col = paste0(prefix, "_empirical_size_p"),
      n_col = paste0(prefix, "_empirical_size_n")
    )
    if (replicate_col %in% names(x)) {
      x <- passage_empirical_tail_calibrate(
        x, p_col = p_col, group_cols = group_cols,
        strata_cols = size_col, leaveout_cols = replicate_col,
        out_col = paste0(prefix, "_empirical_leave_rep_size_p"),
        n_col = paste0(prefix, "_empirical_leave_rep_size_n")
      )
    }
  }
  x
}
