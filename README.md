# Concordant Pathogenicity Predictors Identify Survival-Associated HRR Variants of Uncertain Significance in Germline-Enriched Genes Across 31 Cancer Types

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

## Overview

**AlphaHRD** is a three-layer computational framework for interpreting HRR missense variants of uncertain significance, combining pathogenicity prediction with allele-specific LOH and genomic scar characterization. We benchmarked four predictors (AlphaMissense, REVEL, CADD, PrimateAI) against overall survival within germline-enriched HRR genes across 31 TCGA pan-cancer types.

> Raphael B. Moreira, Indianara V. Brandão, Felipe Batalini, Vamsidhar Velcheti
>
> Target: JCO Precision Oncology

## Key Findings

| Finding | Value |
|---------|-------|
| ClinVar concordance (κ) | 0.733 (CI 0.712–0.754) |
| VUS reclassified | 90.1% (66,912 / 74,246) |
| **Concordant pathogenic (AM + REVEL) HR** | **0.563 (0.390–0.814, p = 0.0022)** |
| AM within-germline HR | 0.646 (0.467–0.893, p = 0.0083) |
| REVEL within-germline HR | 0.679 (0.523–0.882, p = 0.0037) |
| Gene-stratified HR (within-gene effect) | 0.644 (p = 0.013) |
| E-value (concordance) | 2.95 (CI bound 1.76) |
| Stratified Cox HR (pan-cancer) | 0.868 (0.737–1.023, p = 0.091) |
| Biallelic vs monoallelic HRDsum | 27.0 vs 11.0 (p = 4.63 × 10⁻⁹) |
| Purity-filtered HR (≥ 0.5) | 0.699 (0.573–0.853, p = 0.0004) |
| Tumors analyzed | 31 TCGA PanCancer Atlas |
| Patients | 1,939 (1,883 with OS) |
| Variants | 4,301 HRR missense |

## Central Result

Within germline-enriched HRR genes (BRCA1/2, ATM, PALB2, CHEK2), patients classified as pathogenic by **both** AlphaMissense and REVEL had the strongest survival association (HR = 0.563, p = 0.0022). When the two predictors disagreed, neither group reached significance. Cross-predictor concordance defines the clinically informative VUS subset.

CADD trended without significance (p = 0.091). PrimateAI did not discriminate (p = 0.53).

## Pipeline

```
Layer 1: Pathogenicity classification
         AlphaMissense (≥0.564) + benchmark vs REVEL (≥0.5), CADD (≥25), PrimateAI (≥0.803)
         → Concordance analysis: AM × REVEL agreement

Layer 2: Allele-specific LOH (ASCAT3)
         → Biallelic vs monoallelic classification

Layer 3: Genomic scar characterization
         → HRDsum (LOH + TAI + LST)
         → Six-tier classification (biological characterization):
           TRUE_HRD, PROBABLE_HRD, BIALLELIC_NO_SCAR,
           MONOALLELIC, NO_LOH_EVALUABLE, BENIGN
```

## Notebooks

| # | Script | Purpose | Key Output |
|---|--------|---------|------------|
| NB1 | `Notebook1_Data_Acquisition.ipynb` | TCGA-PRAD download, HRR filter, AM annotation (prostate subset) | `annotated_hrr_variants.csv` (52 variants; pan-cancer data in `pancancer_hrr_variants.csv` from NB4) |
| NB2 | `Notebook2_Concordance_VUS.ipynb` | ClinVar κ, VUS reclassification | `concordance_results.csv` |
| NB3 | `Notebook3_mCRPC_Validation.py` | MSK-IMPACT external validation (negative) | `mcrpc_validation_results.csv` |
| NB4 | `Notebook4_PanCancer_Survival.ipynb` | 31-tumor Cox, LOH, HRDsum, meta-analysis | `pancancer_hrr_variants.csv` |
| NB5 | `Notebook5_Synthetic_Validation.ipynb` | Simulation accuracy 87%, power analysis | `synthetic_module1-3.csv` |
| NB6 | `Notebook6_Robustness_Analyses.ipynb` | Schoenfeld, bootstrap, permutation | `reviewer_robustness.csv` |
| NB7 | `Notebook7_Extended_Robustness.py` | Germline/somatic proxy, RMST, purity, E-value | `robustness_5analyses_FINAL.csv` |
| NB8 | `Notebook8_BRCA1_Methylation.py` | Promoter methylation independence | `analysis6_results.csv` |
| **NB10** | **`Notebook10_Benchmark_Concordance.py`** | **Predictor benchmark, AM×REVEL concordance, LOGO, CCF** | **`benchmark_predictors.csv`, `concordance_survival.csv`** |

## Sensitivity Analyses

| # | Analysis | Result |
|---|----------|--------|
| 1a | Germline-enriched vs somatic-enriched | Germline HR = 0.617 (p < .0001); Somatic HR = 0.914 (p = .333) |
| 1b | Within-germline predictor concordance | Both-pathogenic HR = 0.563 (p = 0.0022) |
| 2 | RMST | +3.1 to +7.3 months across τ = 60, 81, 120 (all p < 0.01) |
| 3 | Purity filter | HR strengthens: 0.806 → 0.699 |
| 4 | E-value | Concordance: 2.95 (CI 1.76); Pan-cancer unadjusted: 1.80 (CI 1.31); Stage-adjusted: 1.56 (CI 1.00) |
| 5 | BRCA1 methylation | Independent (OR ≈ 1.0, p = 0.66) |
| 6 | TMB + immune | TMB 3× higher; immune infiltrate identical |
| 7 | Stage-adjusted Cox | HR = 0.87 (p = 0.156) |

## Important Limitations

- ATM concentrates the within-germline signal (excluding ATM: p = 0.24)
- Concordance analysis was added during pre-submission audit, not pre-specified
- AlphaMissense does not predict HRD status independently
- TCGA is treatment-naive; PARPi-treated validation needed
- Six-tier classification is biological characterization, not a validated prognostic tool
- TRUE_HRD subgroup (n = 36 with OS data) is underpowered for survival analysis

## Data Sources

- Somatic mutations: cBioPortal REST API, TCGA PanCan Atlas 2018 (MC3)
- AlphaMissense: Zenodo (doi:10.5281/zenodo.8208688)
- REVEL, CADD, PrimateAI: myvariant.info batch API (GRCh37)
- ASCAT3 segments: GDC allele-specific copy number
- HRDsum: Knijnenburg et al. Cell Reports 2018
- ClinVar: NCBI FTP (February 2026)

## Documentation

- [AUDIT.md](AUDIT.md) — Pre-submission audit trail
- [CHANGELOG.md](CHANGELOG.md) — Version history
- `manifest.json` — Key results for reproducibility verification

## Reproducibility

Python 3.12.3, seed 42. All dependencies pinned in `requirements.txt`.
Run NB1 → NB2 → NB3 → NB4 → NB5–NB8 → NB10 in order. Codespaces-ready.


## License

[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Contact

- Raphael B. Moreira, MD — rb@firstsaude.com
