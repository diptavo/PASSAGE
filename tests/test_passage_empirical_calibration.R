if (file.exists(file.path("R", "passage_empirical_calibration.R"))) {
  source(file.path("R", "passage_empirical_calibration.R"))
} else {
  library(PASSAGE)
}

x <- data.frame(
  dataset = rep("d1", 6),
  mode = rep("m1", 6),
  matching = rep("match", 6),
  mc_sampler = rep("module", 6),
  replicate = rep(1:3, each = 2),
  pathway_size = rep(c(10, 20), 3),
  competitive_score_p = c(0.01, 0.50, 0.04, 0.60, 0.20, 0.90)
)

y <- passage_empirical_competitive_calibration(x)
stopifnot("competitive_score_empirical_pooled_p" %in% names(y))
stopifnot("competitive_score_empirical_leave_rep_p" %in% names(y))
stopifnot("competitive_score_empirical_size_p" %in% names(y))
stopifnot("competitive_score_empirical_leave_rep_size_p" %in% names(y))
stopifnot("competitive_score_empirical_leave_rep_size_n" %in% names(y))
stopifnot(abs(y$competitive_score_empirical_pooled_p[[1]] - 2 / 7) < 1e-12)
stopifnot(abs(y$competitive_score_empirical_leave_rep_p[[1]] - 1 / 5) < 1e-12)
stopifnot(abs(y$competitive_score_empirical_size_p[[1]] - 2 / 4) < 1e-12)
stopifnot(abs(y$competitive_score_empirical_leave_rep_size_p[[1]] - 1 / 3) < 1e-12)
stopifnot(y$competitive_score_empirical_leave_rep_size_n[[1]] == 2L)

z <- passage_empirical_tail_calibrate(
  x, p_col = "competitive_score_p",
  group_cols = c("dataset", "mode", "matching", "mc_sampler"),
  strata_cols = "pathway_size",
  leaveout_cols = "replicate",
  out_col = "p_emp"
)
stopifnot(all.equal(z$p_emp, y$competitive_score_empirical_leave_rep_size_p))

message("PASSAGE empirical calibration test passed")
