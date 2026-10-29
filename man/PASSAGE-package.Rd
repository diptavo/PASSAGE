\name{PASSAGE-package}
\alias{PASSAGE-package}
\alias{PASSAGE}
\docType{package}
\alias{passage_acat}
\alias{passage_default_range_grid}
\alias{passage_fit_engine_pca}
\alias{passage_fit_engine_crossfit_spots}
\alias{passage_fit_engine_pathway_holdout}
\alias{passage_fit_engine_alternating_gp}
\alias{passage_fit_engine_nmf}
\alias{passage_fit_engine_smoothed_pca}
\alias{passage_fit_engine_spatial_basis}
\alias{passage_fit_factor_engine}
\alias{passage_compare_factor_engines}
\alias{passage_conditional_competitive_test}
\alias{passage_competitive_fast_context}
\alias{passage_empirical_competitive_calibration}
\alias{passage_empirical_tail_calibrate}
\alias{passage_gene_spatial_stats}
\alias{passage_gene_score_z}
\alias{passage_gene_pERSA}
\alias{passage_gene_lse}
\alias{passage_gene_activity_pathway_score_stat}
\alias{passage_effective_pathway_size}
\alias{passage_pathway_activity_hotspot_stat}
\alias{passage_pathway_factor_hsic_stat}
\alias{passage_pERSA_pathway_score_stat}
\alias{passage_fast_cEPSV}
\alias{passage_fast_ePSV}
\alias{passage_fast_pc1_spatial_fraction}
\alias{passage_h_precompute}
\alias{passage_make_gene_bins}
\alias{passage_make_expression_coherence_balance}
\alias{passage_make_factor_coherence_balance}
\alias{passage_combine_balance_functions}
\alias{passage_pathway_coherence_test}
\alias{passage_pathway_covariance_metrics}
\alias{passage_pathway_score_stat}
\alias{passage_standardize_gene_stat_by_bin}
\alias{passage_robust_pathway_score_stat}
\alias{passage_sparse_topk_pathway_score_stat}
\alias{passage_pathway_spatial_transferability}
\alias{passage_pve}
\alias{passage_region_enrichment_test}
\alias{passage_resolve_pathway}
\alias{passage_run}
\alias{passage_score_test}
\alias{passage_score_test_refit_null}
\alias{passage_competitive_gene_stat_z}
\alias{passage_competitive_permutation_test}
\alias{passage_compile_joint_factor_tmb}
\alias{passage_compare_joint_factor_nested}
\alias{passage_fit_joint_factor_hypotheses}
\alias{passage_fit_joint_factor_nested}
\alias{passage_fit_joint_factor_pathway}
\alias{passage_fit_joint_factor_pathways}
\alias{passage_fit_joint_factor_tmb}
\alias{passage_joint_factor_pve}
\alias{passage_select_joint_factor_tmb}
\alias{passage_simulate_joint_factor_tmb_data}
\alias{passage_test_pathway}
\alias{passage_vecchia_precision}
\alias{spapath_acat}
\alias{spapath_default_ranges}
\alias{spapath_feature_scores}
\alias{spapath_knn}
\alias{spapath_pathway_tests}
\alias{spapath_sparse_kernel_from_knn}
\alias{spapath_test}
\title{PASSAGE: pathway-level spatial association testing}
\description{
PASSAGE tests predefined gene sets for residual spatial association. It
provides pathway PCA and factor analysis, fixed-kernel subspace tests,
GP covariance comparisons, and a separate global competitive analysis with
matched gene sets.
}
\details{
For pathway estimation use \code{\link{passage_fit_pathway}}. Its unsmoothed
scores can be tested with \code{\link{passage_test_pathway_gp}} or
\code{\link{passage_test_pathway_fit}}. For competitive analysis use
\code{\link{passage_analyze_bundle}}, which compares pathway scores with
unique matched gene sets using a global factor engine.

The pathway tests assume independent Gaussian row errors with common gene
covariance conditional on a fixed adjustment design. Heteroskedasticity and
expression-derived covariates require further calibration. The competitive
reference depends on exchangeability within the matching strata.

The remaining functions in this consolidated topic support older engines and
research analyses. The \code{spapath_*} functions retain compatibility with
the earlier subset-adaptive prototype. See the repository methods and software
guides for the current workflows.
}
\keyword{package}
