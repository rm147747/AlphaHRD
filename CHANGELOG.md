# Changelog

## [2.0.0] — 2026-03-17

### New Analyses (8 Robustness Suite)
- **Notebook7_Extended_Robustness.py**: 5 new sensitivity analyses
  - 7.1 Germline vs Somatic: HR=0.62 (p<0.001) for germline BRCA1/2/ATM/PALB2/CHEK2
  - 7.2 Stage-adjusted Cox: HR attenuates to 0.87 (p=0.16), delta=3.1%
  - 7.3 RMST at 3 tau: +3.1 to +7.3 months survival difference, all p<0.01
  - 7.4 Purity filter: HR strengthens monotonically 0.80→0.70 with increasing purity
  - 7.5 E-value: 1.80 (CI bound 1.31) for baseline HR
- **Notebook8_BRCA1_Methylation.py**: BRCA1 promoter methylation (Layer 2B)
  - 3 CpG probes (cg13601799, cg08047457, cg19531713) from 450K array via Xena API
  - 85.4% overlap (1,608/1,883 patients)
  - BRCA1 methylation independent of AM status (OR=0.91-1.27)
### New Data Files
- `results/robustness/` — 5-analysis and 8-analysis summary CSVs + full details JSON
- `results/methylation/` — BRCA1 probe-level data + merged analysis
- `results/immune/` — Thorsson feature comparison
- `figures/Fig_robustness_forest_plot.png` — Forest plot for analyses 1-5

### Documentation
- Added AUDIT.md: complete audit trail with data provenance, decision justifications, negative findings, AI tool disclosure
- Updated README.md: reflects 9 notebooks, 8 robustness analyses, full file tree
- Updated METHODS.md: sections 8.6-8.13 for all new analyses
- Updated manifest.json: includes extended robustness results

## [1.2.0] — 2026-03-16

### New Notebooks
- **Notebook5_Synthetic_Validation.ipynb**: Pipeline unit test (87.3%), power analysis, adversarial stress test (12/12 survive)
- **Notebook6_Robustness_Analyses.ipynb**: Age-adjusted Cox, Schoenfeld (p=0.870), bootstrap HR (0.802, CI 0.669-0.942), permutation (p<0.001), ridge-penalized Cox

## [1.1.0] — 2026-03-16

### Repository Cleanup (Audit-Ready)
- Renamed notebooks, removed debug cells, added METHODS.md and CHANGELOG.md
- HR corrected from 0.87 (sensitivity) to 0.85 (primary) in manuscript

## [1.0.0] — 2026-02-16

### Initial Release
- 4 notebooks (NB1-NB4), 25 HRR genes, 31 TCGA tumor types
- ClinVar kappa=0.733, VUS reclassification 90.1%
- Stratified Cox HR=0.85, meta-analytic HR=0.86
