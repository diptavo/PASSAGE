#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) stop("Usage: passage_cli.R BUNDLE.rds PATHWAYS.rds OUTPUT_DIR [B=9999] [SEED=1]")
suppressPackageStartupMessages(library(PASSAGE))
out <- args[[3L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
B <- if (length(args) >= 4L) as.integer(args[[4L]]) else 9999L
seed <- if (length(args) >= 5L) as.integer(args[[5L]]) else 1L
fit <- passage_analyze_bundle(args[[1L]], pathways = readRDS(args[[2L]]), B = B, seed = seed)
saveRDS(fit, file.path(out, "passage_result.rds"))
write.csv(fit$summary, file.path(out, "pathways.csv"), row.names = FALSE)
write.csv(fit$coverage, file.path(out, "pathway_coverage.csv"), row.names = FALSE)
capture.output(fit$session_info, file = file.path(out, "sessionInfo.txt"))
writeLines("PASS", file.path(out, "analysis.complete"))
