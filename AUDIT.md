# Audit Trail — AlphaHRD Repository

This document provides a complete audit trail for the AlphaHRD project.
Every analytical result is traceable to a specific script/notebook cell,
with explicit documentation of data provenance, decisions, and limitations.

## Data Provenance

| Data Source | Access Method | Date | Identifier |
|-------------|--------------|------|------------|
| TCGA PanCan Atlas mutations | cBioPortal REST API | 2026-02-16 | 31 `_tcga_pan_can_atlas_2018` studies |
| TCGA PanCan Atlas clinical | cBioPortal REST API | 2026-02-16 | OS, age, stage, purity (ABSOLUTE) |
| AlphaMissense scores | hegelab.org per-protein TSVs | 2026-02-16 | Cheng et al., Science 2023 |
| ClinVar classifications | NCBI E-utilities API | 2026-02-16 | Pathogenic/LP + Benign/LB only |
| MSK-IMPACT clinical sequencing | cBioPortal (`msk_impact_2017`) | 2026-02-16 | Zehir et al., Nat Med 2017 |
| TCGA LOH (ASCAT) | Google BigQuery ISB-CGC | 2026-03-16 | `isb-cgc-bq.pancancer_atlas.Filtered_all_ASCAT_data_by_Sample` |
| HRDsum (Knijnenburg) | GDC PanCan Atlas | 2026-03-16 | Knijnenburg et al., Cell Rep 2018 |
| 450K methylation probes | UCSC Xena API (xenaPython) | 2026-03-17 | `jhu-usc.edu_PANCAN_HumanMethylation450...xena` |
| Thorsson immune features | Cell.com supplementary mmc2 | 2026-03-17 | Thorsson et al., Immunity 2018, Table S1 |

## Analysis-to-Code Mapping

| Analysis | Script/Notebook | Output File(s) |
|----------|----------------|-----------------|
| Data acquisition | NB1 | `results/annotated_hrr_variants.csv` |
| ClinVar concordance | NB2 | Figures 1-5, kappa statistics |
| MSK-IMPACT external validation | NB3 | `results/mcrpc_validation_results.csv`, Figure 6 |
| Pan-cancer survival | NB4 | `results/pancancer_cox_results.csv` |
| Synthetic validation | NB5 | `results/synthetic_module1-3.csv` |
| Original robustness (5 tests) | NB6 | `results/reviewer_robustness_analyses.csv` |
| **Extended robustness (5 new)** | **NB7** | **`results/robustness/robustness_5analyses_FINAL.csv`** |
| **BRCA1 methylation** | **NB8** | **`results/methylation/analysis6_results.csv`** |
| **Predictor benchmark + concordance** | **NB10** | **`results/benchmark_predictors.csv`, `results/concordance_survival.csv`** |

## Key Statistical Decisions (with justification)

### 1. Germline-gene proxy (NB7, Analysis 1)
**Decision:** Use gene identity (BRCA1/2, ATM, PALB2, CHEK2) as proxy for germline origin.
**Justification:** TCGA PanCan Atlas does not reliably separate germline from somatic calls in all tumor types. These 5 genes have >80% germline origin in the Huang et al. (Cell 2018) TCGA germline analysis. This is a conservative proxy; some somatic variants in these genes will be misclassified as germline.
**Alternative tested:** None available without controlled germline data.
**Impact if wrong:** Would dilute the germline signal toward null, making our estimate conservative.

### 2. Stage adjustment with complete-case analysis (NB7, Analysis 2)
**Decision:** Complete-case analysis (patients with available stage data only).
**Justification:** 30% of patients lack AJCC stage data. Multiple imputation was considered but not pursued because (a) the primary analysis is associative, not causal, (b) the missing mechanism may be MNAR (some tumor types do not use AJCC staging), and (c) the complete-case sensitivity analysis is more conservative.
**Impact:** Reduced sample size (1,883 → 1,311) widened confidence intervals and contributed to loss of significance. The HR direction was preserved (0.87 < 1.0).

### 3. RMST truncation points (NB7, Analysis 3)
**Decision:** Three tau values: 60 months (5 years, clinical standard), 81 months (90th percentile of observed follow-up), 120 months (10 years).
**Justification:** Multiple tau values demonstrate robustness across follow-up horizons. The 90th percentile avoids extrapolation beyond observed data. Pre-specified, not selected post-hoc.

### 4. Purity thresholds (NB7, Analysis 4)
**Decision:** Four thresholds: 0.2, 0.3, 0.4, 0.5.
**Justification:** Progressive filtering from liberal (0.2, loses 25 patients) to restrictive (0.5, loses 547 patients). Demonstrates dose-response relationship between purity and signal strength.

### 5. BRCA1 methylation probes (NB8, Analysis 6)
**Decision:** Three probes (cg13601799, cg08047457, cg19531713) averaged for composite promoter methylation score.
**Justification:** cg13601799 is the most cited BRCA1 promoter probe in the literature. The other two are adjacent CpGs in the same promoter region. Averaging reduces noise from individual probe outliers.
**Threshold:** Beta > 0.2 (moderate), > 0.3 (specific). Both reported.

## Negative and Null Findings (explicitly documented)

| Finding | Implication |
|---------|------------|
| Stage adjustment attenuates HR to 0.87 (p=0.156) | Confounding by stage is real |
| E-value = 1.80 (CI bound 1.31) | Moderate robustness to unmeasured confounding |
| BRCA1 methylation does not differ by AM group | Pipeline captures different HRD pathway |
| Somatic-only carriers: HR=0.914 (p=0.333) | No survival benefit for somatic HRR variants |
| mCRPC validation: opposite direction (HR=1.15) | Referral bias in advanced cohort |

## Reproducibility Checklist

- [x] Random seed: 42 (set in every notebook/script)
- [x] Python 3.12.3 with pinned packages (requirements.txt)
- [x] All data from public sources (no controlled/restricted access needed for core analyses)
- [x] Every figure traceable to specific code cell
- [x] Every number in manuscript cross-referenceable to results CSV
- [x] Notebooks run in order: NB1 → NB2 → NB3 → NB4 → NB5 → NB6 → NB7 → NB8 → NB10
- [x] NB3 depends on NB1 outputs (msk_impact_v72_scored.csv)
- [x] NB7-NB8 depend on NB4 outputs (analysis_dataset_robustness.csv)
- [x] NB10 depends on NB4 outputs (variants_with_scores.csv, pancancer_hrr_variants.csv)
- [x] manifest.json contains key results for automated verification


### What AI (LLM) assisted with
- Code support for data processing, statistical analysis,
- Pipeline automation (API queries, data harmonization)
- Figure generation (forest plots, KM curves)
- Repository organization and documentation

### What was human-led 
- Code generation for data processing, statistical analysis, and visualization (Python and R)
- Study design and research question formulation
- Gene panel selection (25 HRR genes)
- Literature search and reference verification
- Statistical method selection and justification
- Clinical interpretation of all results
- All data source decisions
- Review, validation, and approval of every output
- Senior author communications and collaboration decisions

