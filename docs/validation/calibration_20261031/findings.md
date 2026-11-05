# Unequal noise and strong nuisance variation: PASSAGE 0.6.1

The new calibration options substantially improve the two targeted scenarios.
The main experiment completed **6,800 datasets**, followed by **600 independent
weak-signal datasets**. There were no computational failures or runtime warnings
in either experiment. All remote numerical work used Biowulf batch jobs.
These findings support explicit options under stated assumptions; they do not
establish universal calibration for curated pathways or spatially dependent
spot errors.

## What changed

1. **Competitive reference:** `reference = "universe"` samples the full eligible
   gene universe within matching strata, including possible overlap with the
   target. The same statistic, estimation and selection rules apply to every
   set. This corrects the asymmetry of removing the target from every reference
   draw under a conditional uniform-gene-label null. The legacy complement
   comparison remains the default for reproducibility; users must select the
   new reference explicitly.
2. **Nuisance adjustment:** `passage_prepare_nuisance()` supports measured
   covariates and optional PCs from prespecified negative-control genes. It
   reports leverage and residual degrees of freedom, and returns an analysis
   gene universe excluding the controls. Known covariates are preferred.
3. **Self-contained GP calibration:** `passage_test_pathway_gp_wild()` uses
   HC2-adjusted residuals and one Rademacher multiplier per spot, applied jointly
   across pathway genes. Each draw refits covariate adjustment, PCA/FA, selected
   rank and the GP parameter search. Spot-specific residual magnitudes stay at
   their original locations. The observed GP statistic is unchanged.

The wild bootstrap is approximate. It does not fit a heteroskedastic GP variance
surface; spatial fractions/ranges from the constant-nugget profile remain
descriptive, potentially misspecified parameter summaries.

## Competitive calibration

The main count simulation used 64 spots, 240 genes, 12-gene pathways, global
rank 4 and 399 reference draws. Pathway membership was uniformly random before
gene parameters and expression were generated. The noisy scenario varied NB
dispersion across space; the nuisance scenario included eight strong latent
factors. Expression was log1p CPM10k. Four existing statistics and three
empirical combinations were evaluated without tuning to the new results.

The following comparisons use **the same matching-supported datasets** for
both references. The minimum-p combination contains loading, directional,
robust and sparse statistics. Rates are at a nominal 5% threshold.

| Null scenario | Paired N | Complement minimum-p | Full-universe minimum-p |
|---|---:|---:|---:|
| Unequal noise | 961 | 7.08% | 4.06% |
| Strong nuisance | 986 | 8.01% | 4.87% |
| Both | 491 | 7.13% | 4.89% |

On all full-universe datasets, the minimum-p rates were 4.1% (41/1,000;
95% CI 3.0–5.5%), 4.9% (49/1,000; 3.6–6.4%) and 4.8% (24/500; 3.1–7.1%).
All seven unadjusted full-universe outputs ranged from 4.0–4.9% under unequal
noise, 4.8–6.1% under strong nuisance and 4.2–6.0% with both. Every one of those
pointwise intervals included 5%.

The complement reference had insufficient matching support in **108 of 3,700
competitive datasets**. These cases were retained as unsupported, not dropped
silently or repaired by post hoc coarsening. Full-universe randomization was
available in all 3,700. Its all-case rates and the paired comparisons are
reported separately; changing denominators is not presented as a paired gain.

The validity argument is conditional: hold expression, design, strata and any
target-independent shared fit fixed, and assume uniform target membership
within the declared stratum counts. Target and reference scores are then
exchangeable. This can tolerate heterogeneous spot noise because it randomizes
gene labels, not spot errors. It does **not** establish that a fixed curated
co-regulated pathway has uniform labels relative to unrelated genes. A
target-specific shared fit must instead be refitted for each candidate set.

## Nuisance adjustment and detection

With the full-universe reference, the directional statistic gave:

| Scenario | No adjustment | Known nuisance factors | Negative-control PCs |
|---|---:|---:|---:|
| Strong-nuisance null, N=1,000 | 4.9% | 4.8% | 4.7% |
| Strong-nuisance alternative, N=200 | 54.5% | 72.5% | 71.0% |
| Both-noise-and-nuisance null, N=500 | 4.2% | 4.6% | 4.2% |
| Both-noise-and-nuisance alternative, N=200 | 58.5% | 79.5% | 77.0% |

These are paired method comparisons on the same generated datasets. The controls
were idealized, independent auxiliary measurements sharing the nuisance factors
but no target spatial signal. Their quality is part of the simulation design,
not evidence that arbitrary housekeeping genes are valid real-data controls.
Known-factor adjustment is a favorable measured-covariate benchmark. In the count
simulation, the linear adjustment is applied after normalization; it is not an
oracle fit to the entire nonlinear count mean.

**Residual warnings are retained.** Under the combined null, the adjusted
robust statistic rejected 7.2% with measured factors (95% CI 5.1–9.8%). Sparse
top-k rejected 7.2% with control PCs and 7.4% with measured factors (5.3–10.1%
for the latter). These are pointwise, unadjusted flags among many correlated
checks; they do not establish systematic inflation, but neither are they
discarded. No general adjusted robust/sparse default is recommended from this
experiment. There was no tuning or repeated simulation until these flags vanished.

## Self-contained GP calibration

This experiment used 500 datasets per null scenario, 64 spots, 12 correlated
genes, fixed pathway PCA rank 2, the sum statistic, dense kernels at ranges 0.1
and 0.3, and 399 draws. The mean design was known and correctly specified.
Heteroskedasticity produced an approximately 20-fold spot-variance gradient.

| Null scenario | Gaussian rotation | Wild bootstrap | Wild 95% CI |
|---|---:|---:|---:|
| Gaussian, common variance | 5.0% | 5.0% | 3.3–7.3% |
| Unequal noise | 14.4% | 4.2% | 2.6–6.3% |
| Strong measured nuisance | 4.0% | 4.2% | 2.6–6.3% |
| Unequal noise plus nuisance | 12.4% | 6.2% | 4.3–8.7% |
| Heavy tails plus unequal noise | 17.8% | 5.6% | 3.8–8.0% |

The wild reference substantially reduces the clear heteroskedastic failures of
Gaussian rotations. The combined case remains an uncertainty: 6.2% is compatible
with 5% here, but its interval also allows meaningful inflation. This is not a
certificate of exact finite-sample control.

The primary alternatives were too strong: both methods detected all 200 datasets
in each scenario. Therefore a **separately labeled exploratory follow-up** reduced
the added coefficient from 0.45 to 0.15, used new independent seeds, and left all
methods/settings unchanged. With 200 datasets per scenario:

| Weak alternative | Gaussian rotation | Wild bootstrap |
|---|---:|---:|
| Common noise variance | 32.0% | 31.5% |
| Unequal noise | 25.5% | 14.5% |
| Unequal noise plus nuisance | 19.5% | 12.0% |

The lower heteroskedastic detection rate must be interpreted alongside the
rotation's inflated null rate. It is not a fair comparison of power at equal
type-I error. The common-variance comparison shows essentially unchanged
detection in this small experiment. The weak-signal follow-up is exploratory,
not a prespecified part of the original 6,800-dataset design.

## Practical recommendation and remaining limits

- For competitive inference under the stated conditional gene-label null,
  explicitly use `reference = "universe"`. Preserve the prespecified statistic
  or panel and combiner. Report matching/reference diagnostics.
- Adjust measured nuisance covariates during estimation. Use control PCs only
  with defensible negative controls, and keep those genes out of the tested
  universe. Neither a larger global rank nor empirical p-value combinations
  automatically removes nuisance.
- For self-contained GP testing with unequal noise, prefer the wild-refit
  option over the Gaussian rotation reference within the independently
  distributed spot-error setting studied here. Keep its approximate status
  and the combined-case uncertainty visible.
- Further work is needed for correlated curated pathway nulls, unknown spatial
  dependence between residual errors, skewed normalized-count self-contained
  tests, estimated-covariate uncertainty, high-leverage designs, rank/FA choices,
  smaller p-value tails and study-wide FDR. The calibration experiments do not
  establish external-method superiority or real biological discoveries.

## Software verification and evidence

PASSAGE 0.6.1 passed a full `R CMD check --no-manual --no-vignettes` with
`Status: OK` and all 18 test files/examples. Optional suggests were not forced.
Final documentation/build checks and both installed interactive examples passed.
The package is installed in isolated project libraries locally and on Biowulf.
Frozen simulation R sources, final release R sources and the local repository
match. At the time of this October 31 audit, no public release or Git push had
been performed.

- [Calibration figure](main/calibration.pdf), [PNG](main/calibration.png)
- [All main rates and exact intervals](main/rejection_summary.csv)
- [Paired reference comparisons](main/paired_reference_comparison.csv)
- [Main audit](main/audit.txt), [availability](main/availability.csv)
- [Weak-signal rates](weak_power/rejection_summary.csv), [follow-up audit](weak_power/audit.txt)
- [Original full package check](package_check.txt)
- [Source verification](source_release_verification.txt)

Original October 31 source package SHA256:
`a2ef955e4053f41f0021ce4a19ba5c63862e9bf100ea840f54c179cd9d8c0141`.

The tables, audits, source hashes and package-check record linked above retain
the public evidence for this experiment. The package checksum refers to the
archived audit build, before subsequent publication metadata and documentation
updates.
