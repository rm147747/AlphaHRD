#!/usr/bin/env python3
"""
Notebook7_Extended_Robustness.py
================================
Five prespecified robustness analyses addressing reviewer concerns
for the AlphaHRD pan-cancer survival association.

Analyses:
  7.1  Germline vs. somatic carrier separation
  7.2  Stage-adjusted Cox regression
  7.3  RMST (Restricted Mean Survival Time)
  7.4  Tumor purity filter
  7.5  E-value for unmeasured confounding

Input:  results/robustness/analysis_dataset_robustness.csv
        OR results/survival_data_merged.csv (fallback)
Output: results/robustness/robustness_5analyses_FINAL.csv
        results/robustness/robustness_5analyses_details.json
        figures/Fig_robustness_forest_plot.png

Seed: 42 | Python 3.12 | Date: 2026-03-22
"""

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter, KaplanMeierFitter
from scipy import stats
from pathlib import Path
import json
import warnings
warnings.filterwarnings('ignore')

np.random.seed(42)

# ============================================================
# LOAD DATA — try primary file, fallback to alternative
# ============================================================
primary_path = Path('results/robustness/analysis_dataset_robustness.csv')
fallback_path = Path('results/survival_data_merged.csv')

if primary_path.exists():
    df = pd.read_csv(primary_path)
    print(f"Loaded (primary): {primary_path}")
elif fallback_path.exists():
    df = pd.read_csv(fallback_path)
    print(f"Loaded (fallback): {fallback_path}")
else:
    raise FileNotFoundError("No input dataset found. Run NB4 first.")

# Column compatibility
if 'has_am_pathogenic' in df.columns:
    # Convert to boolean if needed
    df['has_am_pathogenic'] = df['has_am_pathogenic'].map(
        {True: True, False: False, 'True': True, 'False': False, 1: True, 0: False}
    )

df = df[df['os_time'] > 0].dropna(subset=['os_time', 'os_event']).copy()
df['os_event'] = df['os_event'].astype(int)
print(f"Patients: {len(df)}, AM-pathogenic: {df['has_am_pathogenic'].sum()}")

# Ensure output directories exist
Path('results/robustness').mkdir(parents=True, exist_ok=True)
Path('figures').mkdir(parents=True, exist_ok=True)

# ============================================================
# 7.1  GERMLINE vs SOMATIC SEPARATION
# ============================================================
print("\n=== 7.1 Germline vs Somatic ===")

GERMLINE_GENES = {'BRCA1', 'BRCA2', 'ATM', 'PALB2', 'CHEK2'}

if 'genes' in df.columns:
    def classify_carrier(row):
        if not row['has_am_pathogenic']:
            return 'benign'
        genes_str = row.get('genes', '')
        if pd.isna(genes_str):
            return 'somatic'
        genes = set(str(genes_str).replace(' ', '').split(','))
        if genes & GERMLINE_GENES:
            return 'germline'
        return 'somatic'

    df['carrier_class'] = df.apply(classify_carrier, axis=1)
    print(df['carrier_class'].value_counts())

    # Cox: germline vs benign
    df_gb = df[df['carrier_class'].isin(['germline', 'benign'])].copy()
    df_gb['is_germline'] = (df_gb['carrier_class'] == 'germline').astype(int)
    cph = CoxPHFitter()
    cph.fit(df_gb[['os_time', 'os_event', 'is_germline']].dropna(), 'os_time', 'os_event')
    hr_germ = cph.hazard_ratios_.values[0]
    p_germ = cph.summary['p'].values[0]
    print(f"\nGermline vs Benign: HR = {hr_germ:.4f}, p = {p_germ:.4f}")
    cph.print_summary()

    # Cox: somatic vs benign
    df_sb = df[df['carrier_class'].isin(['somatic', 'benign'])].copy()
    df_sb['is_somatic'] = (df_sb['carrier_class'] == 'somatic').astype(int)
    cph2 = CoxPHFitter()
    cph2.fit(df_sb[['os_time', 'os_event', 'is_somatic']].dropna(), 'os_time', 'os_event')
    hr_som = cph2.hazard_ratios_.values[0]
    p_som = cph2.summary['p'].values[0]
    print(f"\nSomatic vs Benign: HR = {hr_som:.4f}, p = {p_som:.4f}")
    cph2.print_summary()
else:
    print("  SKIPPED: 'genes' column not available")
    hr_germ, p_germ, hr_som, p_som = None, None, None, None

# ============================================================
# 7.2  STAGE-ADJUSTED COX
# ============================================================
print("\n=== 7.2 Stage-Adjusted Cox ===")

hr_stage, p_stage = None, None
if 'stage' not in df.columns:
    print("  SKIPPED: 'stage' column not available in input dataset")
else:
    df_stage = df.dropna(subset=['stage', 'os_time', 'os_event']).copy()
    stage_map = {
        'I': 1, 'II': 2, 'III': 3, 'IV': 4,
        'Stage I': 1, 'Stage II': 2, 'Stage III': 3, 'Stage IV': 4
    }
    df_stage['stage_num'] = df_stage['stage'].map(stage_map)
    df_stage = df_stage.dropna(subset=['stage_num'])

    if len(df_stage) < 20:
        print(f"  SKIPPED: only {len(df_stage)} patients with valid stage")
    else:
        df_stage['is_path'] = df_stage['has_am_pathogenic'].astype(int)
        try:
            cph3 = CoxPHFitter()
            cph3.fit(df_stage[['os_time', 'os_event', 'is_path', 'stage_num']].dropna(),
                     'os_time', 'os_event')
            hr_stage = cph3.hazard_ratios_['is_path']
            p_stage = cph3.summary.loc['is_path', 'p']
            print(f"  Stage-adjusted: HR = {hr_stage:.4f}, p = {p_stage:.4f}")
            print(f"  N = {len(df_stage)}")
            cph3.print_summary()
        except Exception as e:
            print(f"  FAILED: {e}")

# ============================================================
# 7.3  RMST (Restricted Mean Survival Time)
# ============================================================
print("\n=== 7.3 RMST ===")

def compute_rmst(time, event, group, tau):
    """Compute RMST difference between groups via KM integration."""
    km = {}
    for g in [0, 1]:
        mask = group == g
        if mask.sum() < 2:
            km[g] = 0
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(time[mask], event[mask])
        t = kmf.survival_function_.index.values
        s = kmf.survival_function_.values.flatten()
        t_tau = t[t <= tau]
        s_tau = s[:len(t_tau)]
        if len(t_tau) > 1:
            km[g] = np.trapz(s_tau, t_tau)
        else:
            km[g] = 0
    return km.get(1, 0) - km.get(0, 0), km

df_rmst = df.dropna(subset=['os_time', 'os_event']).copy()
df_rmst['group'] = df_rmst['has_am_pathogenic'].astype(int)

rmst_results = {}
for tau in [60, 81, 120]:
    delta, km_areas = compute_rmst(
        df_rmst['os_time'].values,
        df_rmst['os_event'].values,
        df_rmst['group'].values,
        tau
    )
    rmst_results[tau] = {'path': km_areas.get(1, 0), 'ben': km_areas.get(0, 0), 'delta': delta}
    print(f"  tau={tau}mo: RMST_path={km_areas.get(1,0):.1f}, "
          f"RMST_ben={km_areas.get(0,0):.1f}, delta={delta:+.1f}mo")

# ============================================================
# 7.4  TUMOR PURITY FILTER
# ============================================================
print("\n=== 7.4 Purity Filter ===")

purity_results = {}
if 'purity' not in df.columns:
    print("  SKIPPED: 'purity' column not available")
else:
    for threshold in [0.2, 0.3, 0.4, 0.5]:
        df_pur = df[df['purity'] >= threshold].dropna(subset=['os_time', 'os_event']).copy()
        df_pur['is_path'] = df_pur['has_am_pathogenic'].astype(int)

        if df_pur['is_path'].sum() < 2 or (df_pur['is_path'] == 0).sum() < 2:
            print(f"  Purity >= {threshold}: SKIPPED (insufficient groups)")
            continue

        try:
            cph_p = CoxPHFitter()
            cph_p.fit(df_pur[['os_time', 'os_event', 'is_path']], 'os_time', 'os_event')
            hr = cph_p.hazard_ratios_.values[0]
            ci = cph_p.confidence_intervals_.values[0]
            p = cph_p.summary['p'].values[0]
            n = len(df_pur)
            excluded = len(df) - n
            purity_results[threshold] = {'hr': hr, 'p': p, 'n': n}
            print(f"  Purity >= {threshold}: HR={hr:.4f} "
                  f"({np.exp(ci[0]):.3f}-{np.exp(ci[1]):.3f}), "
                  f"p={p:.4f}, n={n} (excluded {excluded})")
        except Exception as e:
            print(f"  Purity >= {threshold}: FAILED — {e}")

# ============================================================
# 7.5  E-VALUE — COMPUTED, NOT HARDCODED
# ============================================================
print("\n=== 7.5 E-value ===")

def evalue(hr):
    """Compute E-value for a protective HR (< 1.0)."""
    rr = 1 / hr if hr < 1 else hr
    return rr + np.sqrt(rr * (rr - 1))

# Compute from actual model results
# Unstratified pooled Cox
df_ev = df.dropna(subset=['os_time', 'os_event']).copy()
df_ev['is_path'] = df_ev['has_am_pathogenic'].astype(int)
cph_ev = CoxPHFitter()
cph_ev.fit(df_ev[['os_time', 'os_event', 'is_path']], 'os_time', 'os_event')
hr_unstrat = cph_ev.hazard_ratios_.values[0]
ci_unstrat = np.exp(cph_ev.confidence_intervals_.values[0])

# Stratified Cox (exclude small strata)
from collections import defaultdict
tumor_stats = defaultdict(lambda: {'n': 0, 'ev': 0, 'path': 0, 'ben': 0})
for _, r in df_ev.iterrows():
    t = r.get('tumor', r.get('study', 'unknown'))
    tumor_stats[t]['n'] += 1
    tumor_stats[t]['ev'] += r['os_event']
    tumor_stats[t]['path'] += r['is_path']
    tumor_stats[t]['ben'] += (1 - r['is_path'])

tumor_col = 'tumor' if 'tumor' in df_ev.columns else 'study'
valid_tumors = [t for t, s in tumor_stats.items()
                if s['ev'] >= 3 and s['path'] >= 2 and s['ben'] >= 2]
df_strat = df_ev[df_ev[tumor_col].isin(valid_tumors)].copy()

hr_strat = None
if len(valid_tumors) >= 2:
    try:
        cph_strat = CoxPHFitter()
        cph_strat.fit(df_strat[['os_time', 'os_event', 'is_path', tumor_col]],
                      'os_time', 'os_event', strata=[tumor_col])
        hr_strat = cph_strat.hazard_ratios_.values[0]
        ci_strat = np.exp(cph_strat.confidence_intervals_.values[0])
    except:
        pass

print(f"  Unstratified HR={hr_unstrat:.4f}: E-value = {evalue(hr_unstrat):.2f}")
if hr_strat:
    ev_point = evalue(hr_strat)
    # E-value for CI bound closest to 1
    ci_upper = ci_strat[1]
    ev_ci = 1.0 if ci_upper >= 1.0 else evalue(ci_upper)
    print(f"  Stratified HR={hr_strat:.4f}: E-value = {ev_point:.2f} (CI bound: {ev_ci:.2f})")
if hr_stage:
    print(f"  Stage-adj HR={hr_stage:.4f}: E-value = {evalue(hr_stage):.2f}")

# ============================================================
# SAVE RESULTS
# ============================================================
print("\n=== SAVING RESULTS ===")

results_rows = []
if hr_germ is not None:
    results_rows.append({'analysis': 'Germline proxy vs Benign', 'HR': hr_germ, 'p': p_germ})
    results_rows.append({'analysis': 'Somatic proxy vs Benign', 'HR': hr_som, 'p': p_som})
if hr_stage is not None:
    results_rows.append({'analysis': 'Stage-adjusted Cox', 'HR': hr_stage, 'p': p_stage})
for tau, r in rmst_results.items():
    results_rows.append({'analysis': f'RMST tau={tau}', 'HR': None, 'p': None,
                         'RMST_path': r['path'], 'RMST_ben': r['ben'], 'delta': r['delta']})
for thr, r in purity_results.items():
    results_rows.append({'analysis': f'Purity >= {thr}', 'HR': r['hr'], 'p': r['p'], 'n': r['n']})
results_rows.append({'analysis': 'E-value (unstratified)', 'HR': hr_unstrat,
                     'E_value': evalue(hr_unstrat)})
if hr_strat:
    results_rows.append({'analysis': 'E-value (stratified)', 'HR': hr_strat,
                         'E_value': evalue(hr_strat), 'E_value_CI': ev_ci})

df_results = pd.DataFrame(results_rows)
df_results.to_csv('results/robustness/robustness_5analyses_FINAL.csv', index=False)
print(f"Saved: results/robustness/robustness_5analyses_FINAL.csv")

# Details JSON
details = {
    'germline_hr': float(hr_germ) if hr_germ else None,
    'germline_p': float(p_germ) if p_germ else None,
    'somatic_hr': float(hr_som) if hr_som else None,
    'somatic_p': float(p_som) if p_som else None,
    'stage_hr': float(hr_stage) if hr_stage else None,
    'stage_p': float(p_stage) if p_stage else None,
    'rmst': {str(k): v for k, v in rmst_results.items()},
    'purity': {str(k): v for k, v in purity_results.items()},
    'evalue_unstrat': float(evalue(hr_unstrat)),
    'evalue_strat': float(evalue(hr_strat)) if hr_strat else None,
    'n_patients': len(df),
    'n_events': int(df['os_event'].sum()),
    'seed': 42,
}
with open('results/robustness/robustness_5analyses_details.json', 'w') as f:
    json.dump(details, f, indent=2)
print(f"Saved: results/robustness/robustness_5analyses_details.json")

print("\n=== COMPLETE ===")
