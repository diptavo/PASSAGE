# PASSAGE Analysis Plan: spatialDLPFC

## Main Question

Identify pathways with spatial variation in human DLPFC Visium sections after
controlling for technical covariates and local cell-type composition.

## Data

Use LIBD spatialDLPFC because it provides:

- 30 Visium DLPFC samples from neurotypical controls.
- 113,927 Visium spots in the full processed SpatialExperiment object.
- Same-project snRNA-seq reference with 77,604 nuclei.
- Known spatial structure: cortical layers and white matter.

## Cell-Type Adjustment

The H3 covariate matrix should include:

- Intercept.
- Spot-level technical covariates: log library size, detected genes, and any
  available QC covariates.
- Estimated broad cell-type proportions from the DLPFC snRNA-seq reference.

Because cell-type proportions are compositional, use either:

- Drop-one encoding of broad cell-type proportions; or
- CLR-transformed proportions with a small pseudocount, followed by QR cleanup.

Initial PASSAGE H3 should use broad cell classes, not fine subclusters, to avoid
overfitting and collinearity. Fine subclusters can be used in sensitivity runs.

## Pathway Sets

Start with:

- MSigDB Hallmark.
- Reactome immune/synaptic/metabolic subsets.
- Curated brain layer/cell-type marker pathways only as positive controls, not
  as discovery evidence.

## PASSAGE Tests

Primary inference:

- Competitive H3 test.
- Test statistic: `score_z` initially, then `score_robust_z` sensitivity.
- P-value calibration: module matched competitive null plus empirical
  leave-replicate calibration.

Secondary properties:

- `score_sparse_topk_z`: driver concentration.
- `score_activity_hotspot`: coherent pathway activity localization.
- `score_coherence`: pathway spatial coherence.
- cEPSV: descriptive only, not primary p-value statistic.

## Validation

Expected positive controls:

- Synaptic/neuron projection pathways vary across cortical layers.
- Myelination/oligodendrocyte pathways enrich in white matter.
- Immune/inflammatory pathways should be weaker in neurotypical DLPFC.

Expected H3 behavior:

- H1 should find many layer/cell-type driven pathways.
- H3 should reduce pathways explained by cell-type composition.
- Residual H3 hits should reflect within-cell-type spatial programs,
  laminar gradients not fully explained by broad cell types, or local tissue
  microenvironment effects.

## Deliverables

1. Per-section PASSAGE input files.
2. Cell-type covariate matrix per section.
3. H1 vs H3 pathway result tables.
4. Competitive H3 calibration summary.
5. Driver gene tables for significant H3 pathways.
