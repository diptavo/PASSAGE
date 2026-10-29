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
