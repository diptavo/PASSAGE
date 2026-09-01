# Run Pathway Spatial Transferability on prepared 10x breast Visium data.
#
# Usage:
#   Rscript scripts/run_10x_breast_pst_smoke.R \
#     --n-perm=199 --max-pathways=6 --components=mean,pc1 \
#     --out-dir=results/passage_10x_breast_pst_smoke

parse_args <- function(args) {
  cfg <- list(
    seed = 20260730L,
    n_perm = 199L,
    max_pathways = 6L,
    n_gene_splits = 5L,
    n_folds = 5L,
    k_neighbors = 8L,
    components = c("mean", "pc1"),
    out_dir = file.path("results", "passage_10x_breast_pst_smoke")
  )
  for (arg in args) {
    if (grepl("^--seed=", arg)) cfg$seed <- as.integer(sub("^--seed=", "", arg))
    if (grepl("^--n-perm=", arg)) cfg$n_perm <- as.integer(sub("^--n-perm=", "", arg))
    if (grepl("^--max-pathways=", arg)) cfg$max_pathways <- as.integer(sub("^--max-pathways=", "", arg))
    if (grepl("^--n-gene-splits=", arg)) cfg$n_gene_splits <- as.integer(sub("^--n-gene-splits=", "", arg))
    if (grepl("^--n-folds=", arg)) cfg$n_folds <- as.integer(sub("^--n-folds=", "", arg))
    if (grepl("^--k-neighbors=", arg)) cfg$k_neighbors <- as.integer(sub("^--k-neighbors=", "", arg))
    if (grepl("^--components=", arg)) cfg$components <- strsplit(sub("^--components=", "", arg), ",", fixed = TRUE)[[1]]
    if (grepl("^--out-dir=", arg)) cfg$out_dir <- sub("^--out-dir=", "", arg)
  }
  cfg
}

for (f in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
  source(f)
}

prepared_paths <- function() {
  c(
    Visium_FFPE_Human_Breast_Cancer =
      "results/passage_10x_hallmark_perm999_fastcal/Visium_FFPE_Human_Breast_Cancer/passage_hallmark_prepared_data.rds",
    V1_Breast_Cancer_Block_A_Section_1 =
      "results/passage_10x_hallmark_perm999_fastcal/V1_Breast_Cancer_Block_A_Section_1/passage_hallmark_prepared_data.rds"
  )
}

select_pathways <- function(pathways, max_pathways) {
  preferred <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_ESTROGEN_RESPONSE_LATE",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_HYPOXIA"
  )
  keep <- c(intersect(preferred, names(pathways)), setdiff(names(pathways), preferred))
  pathways[keep[seq_len(min(max_pathways, length(keep)))]]
}

write_report <- function(cfg, tbl) {
  top <- tbl[order(tbl$p, -tbl$statistic, na.last = TRUE), , drop = FALSE]
  lines <- c(
    "# PASSAGE PST Breast Smoke",
    "",
    paste0("- Pathways per dataset: ", cfg$max_pathways),
    paste0("- PST permutations per pathway: ", cfg$n_perm),
    paste0("- Gene splits: ", cfg$n_gene_splits),
    paste0("- Spot folds: ", cfg$n_folds),
    paste0("- Neighbors: ", cfg$k_neighbors),
    "",
    paste(colnames(top), collapse = " | "),
    paste(rep("---", ncol(top)), collapse = " | "),
    apply(top, 1L, function(x) paste(x, collapse = " | "))
  )
  writeLines(lines, file.path(cfg$out_dir, "summary.md"))
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
ii <- 0L
for (dataset in names(prepared_paths())) {
  message("loading ", dataset)
  dat <- readRDS(prepared_paths()[[dataset]])
  gene_names <- colnames(dat$Y)
  pathways <- select_pathways(dat$pathways, cfg$max_pathways)
  for (component in cfg$components) {
    for (pp in seq_along(pathways)) {
      message("dataset=", dataset, " component=", component, " pathway=", names(pathways)[pp])
      t0 <- proc.time()[["elapsed"]]
      ans <- passage_pathway_spatial_transferability(
        dat$Y, dat$coords, pathways[[pp]], X = dat$X, gene_names = gene_names,
        n_gene_splits = cfg$n_gene_splits,
        n_folds = cfg$n_folds,
        k_neighbors = cfg$k_neighbors,
        component = component,
        n_perm = cfg$n_perm,
        seed = cfg$seed + 10000L * match(dataset, names(prepared_paths())) +
          1000L * match(component, cfg$components) + pp
      )
      ii <- ii + 1L
      rows[[ii]] <- data.frame(
        dataset = dataset,
        component = component,
        pathway = names(pathways)[pp],
        pathway_size = ans$pathway_size,
        statistic = ans$statistic,
        p = ans$p,
        n_pairs = ans$n_pairs,
        elapsed_sec = proc.time()[["elapsed"]] - t0,
        stringsAsFactors = FALSE
      )
    }
  }
}

out <- do.call(rbind, rows)
out$fdr <- ave(out$p, out$dataset, out$component, FUN = function(p) stats::p.adjust(p, method = "BH"))
write.csv(out, file.path(cfg$out_dir, "breast_pst_pathway_results.csv"), row.names = FALSE)
write_report(cfg, out)
message("Wrote outputs to ", cfg$out_dir)
print(out[order(out$p, -out$statistic), ], row.names = FALSE)
