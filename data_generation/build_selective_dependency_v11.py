#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_selective_dependency_v11.py
====================================================================
Compute 4 rigorous definitions of PDAC *selective dependency* from the
26Q1 CRISPR pooled per-line gene-effect matrix, addressing the V10
reviewer's explicit demand: "4 definitions of selective dependency
(ratio, z-effect, lineage-adjusted, mixed-effects)".

PDAC line set  : CCLE/DepMap pancreatic-site cell lines present in the
                 26Q1 matrix (37 lines, stable PDAC models).
non-PDAC set   : all other matrix lines (~1,171) = the pan-cancer
                 reference for the selectivity contrast.
CRC line set   : colorectal-site lines (42) used for the cross-cancer
                 transfer / secondary contrast.

Per-gene dependency  dep = -gene_effect  (more negative effect => more essential).

Definitions (per gene g):
  (A) ratio          sel_A = mean(dep_PDAC) / mean(dep_nonPDAC)
                      (fold-change of PDAC dependency over background)
  (B) z-effect       sel_B = (mean(dep_PDAC)-mean(dep_nonPDAC)) / s_pooled
                      (Cohen's d of PDAC selectivity)
  (C) lineage-adj    residualise each line's dep on its primary-site mean
                      (per gene), then recompute (B) on residuals -> removes
                      tissue-of-origin confound
  (D) mixed-effects  fixed PDAC-effect estimate with EB shrinkage toward 0
                      for low-precision genes (random PDAC-effect BLUP-style)

Also derives the PAN-DEPENDENCY gold endpoint (essential = top quartile of
-mean PDAC gene effect) so all endpoints share one consistent line set.
"""
import json, csv, math
import numpy as np
import pandas as pd

import os
ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
MAT  = f"{ROOT}/CRISPRGeneEffect_26Q1.csv"
ANNO = "/Users/zuoxianbo/.workbuddy/skills/zuoxb-virtual-cell-platform/models/singlecell/scfoundation/DeepCDR/data/CCLE/Cell_lines_annotations_20181226.txt"
OUT  = f"{ROOT}/pdac_selective_dependency_v11.json"

print("[1] load CCLE line->cancer map ...")
ach2site = {}
with open(ANNO) as f:
    r = csv.DictReader(f, delimiter='\t')
    for row in r:
        ach = row.get('depMapID')
        if not ach: continue
        ach2site[ach] = (row.get('Site_Primary','').lower(),
                         row.get('Histology','').lower())
pdac_lines = [a for a,(s,h) in ach2site.items() if 'pancrea' in s]
crc_lines  = [a for a,(s,h) in ach2site.items()
              if s in ('colon','large_intestine','rectum') or 'colorect' in s or 'colon' in s]
print(f"    anno PDAC={len(pdac_lines)} CRC={len(crc_lines)}")

print("[2] load 26Q1 CRISPR matrix (float32) ...")
df = pd.read_csv(MAT, index_col=0, low_memory=False)
df = df.astype(np.float32)
line_ids = list(df.index)
gene_cols = list(df.columns)
print(f"    matrix {df.shape[0]} lines x {df.shape[1]} genes")
# map ACH -> row position
pos = {lid:i for i,lid in enumerate(line_ids)}
pdac_idx = [pos[a] for a in pdac_lines if a in pos]
crc_idx  = [pos[a] for a in crc_lines  if a in pos]
nonpdac_idx = [i for i in range(df.shape[0]) if i not in set(pdac_idx)]
print(f"    PDAC rows={len(pdac_idx)} nonPDAC rows={len(nonpdac_idx)} CRC rows={len(crc_idx)}")

# site per line for lineage adjustment
site_per_line = []
for lid in line_ids:
    s = ach2site.get(lid, ('unknown',''))[0]
    site_per_line.append(s)
sites = sorted(set(site_per_line))
site_code = {s:i for i,s in enumerate(sites)}
site_arr = np.array([site_code[s] for s in site_per_line], dtype=np.int32)

arr = df.values  # lines x genes  (float32)
# dependency = -effect ; flip sign so higher = more essential
dep = -arr

print("[3] per-group means / sd (vectorised) ...")
def group_stats(idxs):
    sub = dep[idxs, :]               # n x genes
    m = np.nanmean(sub, axis=0)
    # sample sd with ddof=1, guard against n==1
    if sub.shape[0] > 1:
        sd = np.nanstd(sub, axis=0, ddof=1)
    else:
        sd = np.full(sub.shape[1], np.nan)
    return m, sd

mP, sdP = group_stats(pdac_idx)
mN, sdN = group_stats(nonpdac_idx)
mC, sdC = group_stats(crc_idx)

# (A) ratio
eps = 1e-3
mean_dep_non = np.where(np.abs(mN) < eps, np.sign(mN)*eps if eps else eps, mN)
selA = mP / mean_dep_non

# (B) Cohen's d (PDAC selectivity)
sp = np.sqrt((sdP**2 + sdN**2) / 2.0)
sp = np.where(sp < eps, eps, sp)
selB = (mP - mN) / sp

# (C) lineage-adjusted: residualise dep on per-site mean (per gene), then (B) on resid
print("[4] lineage-adjusted (C) + mixed-effects (D) ...")
n_lines = dep.shape[0]
resid = dep.copy()
# subtract per-site mean (per gene) for lines with known (non-unknown) site
known = site_arr != site_code.get('unknown', -1)
for s in sites:
    if s == 'unknown': continue
    mask = site_arr == s
    if mask.sum() > 0:
        resid[mask, :] -= np.nanmean(dep[mask, :], axis=0)
# recompute group means on residuals
rP = np.nanmean(resid[pdac_idx, :], axis=0)
rN = np.nanmean(resid[nonpdac_idx, :], axis=0)
rsp = np.nanstd(resid[pdac_idx, :], axis=0, ddof=1)
rsn = np.nanstd(resid[nonpdac_idx, :], axis=0, ddof=1)
rsp = np.where(rsp < eps, eps, rsp); rsn = np.where(rsn < eps, eps, rsn)
selC = (rP - rN) / np.sqrt((rsp**2 + rsn**2)/2.0)

# (D) mixed-effects (EB shrinkage of the raw selectivity contrast)
# raw PDAC contrast beta = mP - mN ; precision = 1/var_pooled
var_pool = (sdP**2/len(pdac_idx) + sdN**2/len(nonpdac_idx))
prec = 1.0 / np.where(var_pool < eps, eps, var_pool)
beta_raw = mP - mN
# empirical prior variance of beta across genes (robust)
valid = np.isfinite(beta_raw)
tau2 = np.nanvar(beta_raw[valid]) if valid.sum() > 10 else 1.0
shrink = prec / (prec + 1.0/tau2)
selD = beta_raw * shrink

print("[5] assemble selective-dependency table ...")
# pan-dependency essential: top quartile of -mean PDAC effect (i.e. top quartile of mP)
pan_dep_mean = mP  # mean dependency in PDAC (higher = more essential)
finite_mask = np.isfinite(pan_dep_mean)
thr = np.quantile(pan_dep_mean[finite_mask], 0.75)
pan_essential = (pan_dep_mean >= thr) & finite_mask

# selective-essential: top quartile of selB (most PDAC-selective essential)
sel_finite = np.isfinite(selB)
thr_sel = np.quantile(selB[sel_finite & (pan_dep_mean>0)], 0.75) if (sel_finite & (pan_dep_mean>0)).sum()>10 else 0.0
selective_essential = (selB >= thr_sel) & sel_finite & (pan_dep_mean > 0)

genes = gene_cols
out = {
  "meta": {
    "description": "Per-line 26Q1 CRISPR selective-dependency, 4 definitions. PDAC lines (CCLE pancreatic-site, in matrix), non-PDAC = all other matrix lines, CRC = colorectal-site lines.",
    "n_pdac_lines": len(pdac_idx),
    "n_nonpdac_lines": len(nonpdac_idx),
    "n_crc_lines": len(crc_idx),
    "pan_essential_def": "top quartile of mean PDAC dependency (-gene_effect)",
    "selective_essential_def": "top quartile of definition-B selectivity among PDAC-dependent genes",
    "pdac_line_ids": [line_ids[i] for i in pdac_idx],
    "crc_line_ids": [line_ids[i] for i in crc_idx],
  },
  "genes": genes,
  "pan_dependency_mean": {g: (float(pan_dep_mean[i]) if np.isfinite(pan_dep_mean[i]) else None) for i,g in enumerate(genes)},
  "pan_essential": {g: bool(pan_essential[i]) for i,g in enumerate(genes)},
  "selective_essential": {g: bool(selective_essential[i]) for i,g in enumerate(genes)},
  "definitions": {
    g: {
       "A_ratio": (float(selA[i]) if np.isfinite(selA[i]) else None),
       "B_zeffect": (float(selB[i]) if np.isfinite(selB[i]) else None),
       "C_lineage_adj": (float(selC[i]) if np.isfinite(selC[i]) else None),
       "D_mixed": (float(selD[i]) if np.isfinite(selD[i]) else None),
       "d_pdac": (float(mP[i]) if np.isfinite(mP[i]) else None),
       "d_nonpdac": (float(mN[i]) if np.isfinite(mN[i]) else None),
    } for i,g in enumerate(genes)
  },
}
with open(OUT, "w") as f:
    json.dump(out, f)
print(f"[DONE] wrote {OUT} ({__import__('os').path.getsize(OUT)/1024:.1f} KB)")
print(f"    pan_essential genes={int(pan_essential.sum())}  selective_essential={int(selective_essential.sum())}")

# sanity on known PDAC biology
for g in ["KRAS","TP53","CDKN2A","SMAD4","MYC","GATA6","ROCK1","PAK1"]:
    if g in out["definitions"]:
        d = out["definitions"][g]
        print(f"   {g}: A={d['A_ratio']:.2f} B={d['B_zeffect']:.2f} C={d['C_lineage_adj']:.2f} D={d['D_mixed']:.3f} pan_ess={out['pan_essential'][g]}")
