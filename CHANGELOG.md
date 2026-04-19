# Changelog

## [3.1.0] — 2026-04-18

### Updates
- Expanded manifest.json with full CI, tier counts, and sensitivity results
- Added scikit-learn and tqdm to requirements.txt with pinned versions
- Documentation alignment across README, AUDIT.md, and manifest.json

## [3.0.0] — 2026-04-13

### Predictor Benchmark and Concordance (Notebook 10)
- Four-predictor benchmark: AlphaMissense, REVEL, CADD, PrimateAI
- Within-germline concordance analysis (AM × REVEL)
- Leave-one-gene-out sensitivity
- Exploratory cancer cell fraction analysis
- Six-tier biological characterization (TRUE_HRD, PROBABLE_HRD, BIALLELIC_NO_SCAR, MONOALLELIC, NO_LOH_EVALUABLE, BENIGN)

## [2.0.0] — 2026-03-17

### Extended Robustness Suite
- **Notebook7_Extended_Robustness.py**: 5 sensitivity analyses
  - Germline vs somatic gene-identity proxy
  - Stage-adjusted Cox
  - RMST at three truncation points
  - Progressive purity filtering
  - E-value for unmeasured confounding
- **Notebook8_BRCA1_Methylation.py**: Promoter methylation independence analysis

### New Data Files
- `results/robustness/` — Summary CSVs and full details JSON
- `results/methylation/` — BRCA1 probe-level data
- `results/immune/` — Thorsson feature comparison

### Documentation
- Added AUDIT.md with data provenance and decision justifications
- Updated README.md, METHODS.md, manifest.json

## [1.2.0] — 2026-03-16

### Robustness Analyses
- **Notebook5_Synthetic_Validation.ipynb**: Pipeline unit test, power analysis, adversarial stress test
- **Notebook6_Robustness_Analyses.ipynb**: Age-adjusted Cox, Schoenfeld PH test, bootstrap, permutation, ridge-penalized Cox

## [1.1.0] — 2026-03-16

### Repository Organization
- Renamed notebooks, removed debug cells, added METHODS.md and CHANGELOG.md

## [1.0.0] — 2026-02-16

### Initial Release
- 4 notebooks (NB1-NB4), 25 HRR genes, 31 TCGA tumor types
- ClinVar kappa = 0.733, VUS reclassification 90.1%
- Stratified Cox, meta-analysis, LOH/HRDsum characterization
