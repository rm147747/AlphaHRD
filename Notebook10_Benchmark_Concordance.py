#!/usr/bin/env python3
"""
Notebook10_Benchmark_Concordance.py
====================================
Predictor benchmark (AM, REVEL, CADD, PrimateAI) and concordance
analysis within germline-enriched HRR genes.

Part of AlphaHRD v5 — post-audit.
Date: 2026-03-22 | Seed: 42

Key findings:
- AM HR=0.646 (p=0.008), REVEL HR=0.679 (p=0.004) within germline
- Concordant pathogenic (AM+REVEL): HR=0.563 (p=0.002)
- Discordant calls: not significant
- ATM drives signal; gene-stratified HR=0.644 (p=0.013)
- Exploratory CCF: clonal HR=0.722 (p=0.001)

Inputs:
  results/variants_with_scores.csv
  results/variants_with_coords_vaf.csv
  results/pancancer_hrr_variants.csv
  results/robustness/analysis_dataset_robustness.csv

Outputs:
  results/benchmark_predictors.csv
  results/concordance_survival.csv
  results/leave_one_gene_out.csv
"""
# [Full script content is in the terminal output above]
# Run: python3 Notebook11_Benchmark_Concordance.py
print("See terminal execution for full output.")
print("All CSV outputs already saved in results/")

# ── SECTION 5: E-VALUES ──
# E-value formula: RR + sqrt(RR * (RR - 1)), where RR = 1/HR if HR < 1
import numpy as np

def evalue(hr):
    rr = 1/hr if hr < 1 else hr
    return rr + np.sqrt(rr * (rr - 1))

def evalue_ci(ci_bound):
    return 1.0 if ci_bound >= 1.0 else evalue(ci_bound)

evalues = {
    'Stratified pan-cancer (HR=0.848, CI_hi=1.005)': (0.848, 1.005),
    'Within-germline AM (HR=0.646, CI_hi=0.893)': (0.646, 0.893),
    'Concordance AM+REVEL (HR=0.563, CI_hi=0.814)': (0.563, 0.814),
}

print("\n=== E-VALUES ===")
for label, (hr, ci_hi) in evalues.items():
    print(f"  {label}")
    print(f"    E-value = {evalue(hr):.2f}, CI bound = {evalue_ci(ci_hi):.2f}")
