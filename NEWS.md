# PASSAGE 0.6.1

- Added `reference = "universe"` to generic competitive testing, matching
  diagnostics, set sampling and bundle analysis. It samples the complete
  eligible gene universe within strata, including possible target overlap.
  The legacy `"complement"` default remains available unchanged. Full-universe
  inference has a conditional uniform-gene-label null; it does not establish
  validity for arbitrary correlated biological modules.
- Added `passage_prepare_nuisance()` for measured covariates and optional
  prespecified negative-control PCs, with design/leverage diagnostics and an
  analysis universe excluding controls.
- Added `passage_test_pathway_gp_wild()`: an HC2 residual wild bootstrap that
  preserves spot-specific residual magnitudes and refits PCA/FA, rank selection
  and GP searches. This is an approximate self-contained calibration option;
  the observed GP statistic is unchanged. Gaussian rotation remains separate.
- Added exhaustive finite-universe rank checks and regression tests for the
  bootstrap, RNG restoration, adjustment alignment and negative-control isolation.

# PASSAGE 0.6.0

- Separated statistic construction from competitive inference.
  `passage_make_statistic()` supports all 17 existing global choices plus
  fitted-profile absolute correlation; `passage_competitive_test()` accepts
  any deterministic numeric statistic function of gene indices. Legacy
  conditional tests now use the shared statistic factory.
- Exposed directional burden, sparse, and sign-insensitive choices through
  the main bundle workflow, with cached global quantities and batched default
  scoring.
- Added pathway-local competitive GP statistics with sum, maximum, and top-k
  component aggregation. Local estimation is repeated for every matched set.
  Self-contained GP tests also repeat the selected aggregation in reference
  draws.
- Added Cauchy, minimum-p, and Fisher combinations calibrated against shared
  background score ranks. These do not use independent-p Fisher or analytic
  ACAT reference distributions and do not correct post hoc combiner selection.
- Added matching preflight, exact stratum support and balance diagnostics,
  distinct-set counts, Monte Carlo resolution flags, and component-level
  binomial intervals when background sets were independently generated.
- Exposed global rank rules, prespecified rank sensitivity panels, scaling,
  and estimation diagnostics in the bundle workflow. Known nuisance
  covariates remain part of the estimation design.
- Corrected scoring to use the global engine's gene scales and default
  covariate design. Corrected pathway-holdout preprocessing metadata after
  its final unscaled loading refit. Rerun previously scaled analyses.
- Added a complete interactive R guide for statistics, custom callbacks,
  combined tests, local/global choices, and interpretation.

- Added a complete prepared KC2 example, including the competitive background,
  three Hallmark pathways and a runner for estimation and both testing questions.
- Declared optional-engine dependencies needed for GpGp fitting and TMB
  compilation in clean installations.

These are research capabilities. Matching assumptions, local-GP calibration,
and spatial heteroskedasticity still require study-specific assessment; no
general heteroskedasticity correction was added. Bundle BH/BY summaries adjust
across pathways within each named test, not selection across statistics or
combiners.

# PASSAGE 0.5.0

- Added a pathway GP covariance comparison. The test fits component spatial
  fractions and ranges and calibrates their likelihood improvement using joint
  Gaussian rotations, repeating the parameter search for every null draw.
- Added reusable kernel caches, Nyström approximations and a dense reference.
  Batched spectral calculations avoid repeated GP optimization.
- Added component parameter estimates, boundary flags, settings and timing to
  the test output. The pathway CLI now supports both GP and fixed-kernel tests.
- Reorganized documentation into methods and software guides.

The GP reference assumes a fixed adjustment design and independent Gaussian
rows with common gene covariance. A small pilot showed substantial inflation
under spatial heteroskedasticity; correction for that setting remains open.

# PASSAGE 0.4.0

- Added pathway-only covariance PCA, with 95% variance retention or a fixed
  rank, and fixed-rank factor analysis with regression scores.
- Added testing of the fitted score span using scale-free trace and largest-root
  statistics. These functions do not require a global factor model.
- Added optional descriptive GpGp fits and a single-pathway CLI.

# PASSAGE 0.3.0

- Corrected the PCA/varimax reconstruction so singular-value scale is applied
  once. Results from earlier versions require rerunning after this correction.
- Fixed sampling from singleton matched pools and changed the default legacy
  competitive statistic to `score_z`.
- Added bundle validation, competitive testing with unique matched sets,
  empirical p-values, BH/BY summaries, coverage and candidate driver outputs.
- Added batched set-score calculations and spatial subspace reference tests.
