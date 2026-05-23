# ============================================================================
# load_passage.R
#
# Convenience loader.  Sources all PASSAGE files in the correct dependency
# order so users can do `source("load_passage.R")` and have everything
# available.
#
# Required packages: Matrix, GpGp, ucminf, CompQuadForm.
# Optional: parallel (multi-core in run_passage), sparseinv (faster v_k
#           updates in CAVI engine; v1 falls back to dense inverse).
# ============================================================================

# Layer 1 (engines)
source("engine_pca.R")
source("engine_cavi.R")

# Layer 2 (score tests)
source("passage_score_h1.R")
source("passage_score_h2.R")
source("passage_score_h3.R")
source("passage_score_h4.R")

# Layer 3 (PVE estimators)
source("passage_pve.R")

# Layer 4 (omnibus reporting)
source("passage_omnibus.R")

# Layer 5 (summary-statistics fallback)
source("passage_summary_stats.R")

# Top-level driver
source("run_passage.R")

cat("PASSAGE v0.1.0 loaded.\n")
cat("Main entry points:\n")
cat("  fit_engine_pca(Y, locs, K)         - PCA two-stage engine (E1)\n")
cat("  fit_engine_cavi(Y, locs, K)        - sparse CAVI engine (E3)\n")
cat("  passage_h1 / h2 / h3 / h4(...)     - per-hypothesis tests\n")
cat("  passage_pve(engine, Y, pathway)    - PVE estimators (b, c, d)\n")
cat("  passage_test_pathway(...)          - all H + PVE for one pathway\n")
cat("  run_passage(Y, locs, pathways)     - top-level driver\n")
cat("  run_passage_summary(gene_stats,    - summary-stats fallback\n")
cat("                      pathways)\n")
