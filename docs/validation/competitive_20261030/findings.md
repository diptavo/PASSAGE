# Interpretation of the competitive interface validation

The new interface successfully exposes complementary signal patterns, but the pilot does **not** establish a universally superior statistic or a fully calibrated combined default. The clearest benefit is the directional pathway-average (`score_activity_hotspot`) option for same-direction effects. Sparse and opposing effects show no improvement over the existing `score_z` reference in this particular pilot. Calibration remains an important limitation, especially with heteroskedasticity and strong nuisance variation.

The completed run contains 1,200 datasets. Of these, 1,147 are jointly evaluable across all seven outputs and 53 have unsupported matching. The batch audit reports no computational failures or runtime warnings. Unsupported cases are retained; bins were not relaxed. All rates below are conditional on matching support and use identical replicates and background sets across methods within a scenario.

## Null calibration

At a nominal 5% threshold, the observed rejection percentages were:

| Statistic | IID null, n=190 | Spatial null, n=185 | Heteroskedastic null, n=195 | Strong-nuisance null, n=193 |
|---|---:|---:|---:|---:|
| Existing `score_z` | 4.7 | 8.1 | **8.7** | 5.7 |
| Directional activity hotspot | 4.2 | 7.0 | 6.7 | **9.3** |
| Robust gene-score aggregate | 4.7 | 5.9 | 7.2 | 7.3 |
| Sparse top-k | 8.4 | 6.5 | 8.2 | 6.7 |
| Cauchy combination | 7.4 | 8.1 | 8.2 | 8.3 |
| Minimum-p combination | 5.3 | 7.0 | 8.2 | 7.3 |
| Fisher combination | 7.4 | 7.6 | 7.7 | 7.3 |

Two **pointwise, unadjusted exact 95% binomial intervals** exclude 5%:

- `score_z` under heteroskedasticity: 17/195 = **8.7%**, interval **5.2–13.6%**.
- Activity hotspot under strong nuisance variation: 18/193 = **9.3%**, interval **5.6–14.3%**.

These are calibration warnings, not a claim of familywise significance across the many checks. Other intervals include 5%, but they are broad enough to allow meaningful inflation. For example, sparse top-k under the IID null is 16/190 = 8.4%, interval 4.9–13.3%. A name such as “robust” does not establish robustness to spatially varying noise.

The combination null rates are mostly 7–8%. Their intervals all include 5%, so this pilot does not distinguish modest inflation from simulation uncertainty. It also does not certify control at 5%. The minimum-p combination is closest to nominal under the IID null, but reaches 8.2% under heteroskedasticity, as does Cauchy. Combining statistics does not repair a mismatch between the target and matched-background score distributions.

Do **not** pool these four null scenarios as independent replicates: they share random variates at matching replicate numbers. Exact intervals for each separate scenario are available in [rejection_summary.csv](rejection_summary.csv).

## Detection under alternatives

The table reports rejection percentages at the nominal 5% cutoff. These are **unadjusted detection rates**, not power after equalizing actual type-I error.

| Statistic | Dense, n=96 | Sparse, n=93 | Opposing, n=95 | Dense + strong nuisance, n=100 |
|---|---:|---:|---:|---:|
| Existing `score_z` | 18.8 | 20.4 | 24.2 | 5.0 |
| Directional activity hotspot | **89.6** | 11.8 | 3.2 | **52.0** |
| Robust gene-score aggregate | 21.9 | 17.2 | 20.0 | 5.0 |
| Sparse top-k | 15.6 | 19.4 | 21.1 | 5.0 |
| Cauchy combination | 69.8 | 20.4 | 21.1 | 28.0 |
| Minimum-p combination | 71.9 | 17.2 | 20.0 | 34.0 |
| Fisher combination | 40.6 | 22.6 | 20.0 | 9.0 |

For dense same-direction effects, activity hotspot detects 86/96 datasets versus 18/96 for `score_z`: a **70.8 percentage-point** increase. On the paired data, hotspot gains 69 detections and loses one. Its exact detection interval is 81.7–94.9%, compared with 11.5–28.0% for `score_z`. This is strong evidence of complementary sensitivity in this setting, although the corresponding spatial-null rates are 7.0% and 8.1%, so it is not a comparison at demonstrated equal size.

With strong nuisance variation, hotspot detects 52/100 versus 5/100 for `score_z`, gaining 49 detections and losing two. However, hotspot's corresponding null rate is 9.3% with the warning above. The observed increase is worth pursuing, but it must **not** be described as a validated calibrated-power advantage or as proof that nuisance sensitivity has been solved.

The directional score is weak under opposing effects, consistent with cancellation when pathway genes move in opposite directions. The sparse and robust options do not outperform the reference in the sparse scenario here. These results support a small complementary statistic panel, rather than expecting every new option to improve every scenario.

## Which combination looks most useful?

Minimum-p and Cauchy retain much of the dense-signal sensitivity that Fisher loses when only one component of the panel provides strong evidence. Minimum-p has somewhat higher observed detection here: 71.9% versus 69.8% for dense signal and 34% versus 28% with nuisance. Those differences are only two and six datasets, respectively; the pilot does not establish that minimum-p is generally better. Fisher has the highest sparse detection point estimate, 22.6%, but this differs from the reference by only two datasets.

All three should remain optional, with the statistic panel and primary combiner selected **before** examining pathway results. Reporting the smallest of the three combined p-values would introduce another uncorrected selection step. These implementations use pooled component ranks followed by empirical recalibration on shared background sets; the findings do not directly compare the usual analytic Cauchy or independent-p Fisher reference distributions. With B=99, component p-values have 0.01 resolution, and pooled ranks and ties can materially influence combination behavior.

## Scope and next validation

This is a small software-and-method pilot: 64 spots, 240 genes, pathways of size 12, fixed global rank 4, and 99 backgrounds. It differs from the earlier pilot with 144 spots, 600 genes and 30-gene pathways. Rejection rates across the two studies must not be interpreted as a before/after software comparison. Pathway-local GP sum/max testing was checked in smoke runs only; no local-GP power or calibration result is established by this main run.

The next scientific priority is independent, larger null validation under heteroskedasticity and nuisance variation, followed by comparisons at matched empirical false-positive rates. Predeclare the panel and combination, use more background draws, retain unsupported-matching counts, and include correlated module nulls before claiming broad competitive validity. No automatic default change is justified by this pilot alone.

Evidence: [full rate table and exact intervals](rejection_summary.csv), [paired detection counts](paired_differences.csv), and the completed-run audit and availability reports.

## Software verification and provenance

PASSAGE 0.6.0 passed `R CMD check --no-manual --no-vignettes` with optional suggests not forced: all 17 test files and package examples passed. Optional `fastICA` and `TMB` were unavailable locally, so their optional branches were not fully covered. The installed package passed six focused regression files and the interactive example on Biowulf. The 24-case smoke run included successful global and pathway-local analyses in all eight scenarios; one case had unsupported matching.

All remote numerical work used batch jobs. The main array, audit, installed-package smoke and evidence export completed successfully. Frozen simulation source hashes matched the checked local R source at this audit. The tables and audit reports alongside this document retain the public summary evidence.

Package SHA256: `a5f23524e2bcb4d98744e2407c5c43d0579ad3cfbd3a30a8e314e2a57f80af38`. Downloaded evidence archive SHA256: `baf26d96e1b4584da7b00abf7b5dc2e6f2bb2bef3dfdb3bc49d542b3edc55c12`.
