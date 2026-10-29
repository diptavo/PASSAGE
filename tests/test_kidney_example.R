library(PASSAGE)

input <- system.file("extdata", "kidney_KC2.rds", package = "PASSAGE")
runner <- system.file("examples", "run_kidney.R", package = "PASSAGE")
stopifnot(nzchar(input), nzchar(runner),
          unname(tools::md5sum(input)) == "755944e9845d917b6bdccadecef3492f")
b <- passage_validate_bundle(input)
stopifnot(identical(dim(b$Y), c(2166L, 17630L)), ncol(b$X) == 13L,
          identical(rownames(b$raw_counts), colnames(b$Y)),
          identical(colnames(b$raw_counts), rownames(b$Y)),
          identical(b$spot_qc$barcode, rownames(b$Y)),
          identical(unname(lengths(b$pathway_sets)), c(200L, 36L, 200L)))
coverage <- vapply(b$pathway_sets, function(g) sum(g %in% colnames(b$Y)), integer(1))
stopifnot(identical(unname(coverage), c(190L, 32L, 197L)))
# A pathway-only export would break competitive inference even if it passed
# structural validation. Require the frozen full expression universe.
stopifnot(length(setdiff(colnames(b$Y), unique(unlist(b$pathway_sets)))) > 17000L)
probe <- unique(c(1L, 100L, 1000L, ncol(b$Y)))
normalized <- log1p(sweep(t(as.matrix(b$raw_counts[probe, , drop = FALSE])),
                         1L, b$spot_qc$library_size / 10000, "/"))
stopifnot(max(abs(normalized - b$Y[, probe, drop = FALSE])) < 1e-12)

# Sourcing the installed runner defines a function without starting an analysis.
env <- new.env()
source(runner, local = env)
stopifnot(is.function(env$run_kidney))
out <- tempfile("kidney-protected-")
dir.create(out)
writeLines("keep", file.path(out, "existing.txt"))
err <- try(env$run_kidney(mode = "gp", B = 1L, output = out, input = input), silent = TRUE)
stopifnot(inherits(err, "try-error"), grepl("new or empty", as.character(err)),
          identical(readLines(file.path(out, "existing.txt")), "keep"))
unlink(out, recursive = TRUE)
cat("Kidney example data and installed-runner contracts passed\n")
