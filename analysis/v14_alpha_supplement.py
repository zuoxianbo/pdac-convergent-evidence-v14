# -*- coding: utf-8 -*-
"""
V14 lightweight supplement: compute ECS alpha-sensitivity (the analysis script did not
emit this key, but make_v14_figures.py requires B["alpha_sensitivity"]).

Replicates the EXACT ECS formula from build_v14_analysis.py (D weighted 0.80/0.10/0.10;
PHI = mean of 5 support layers) for a grid of alpha, computes AUROC on E1 & E3, validates
against the known benchmark (alpha=0.6 -> E3 = 0.815), then patches benchmark_v14.json.
No bootstrap / permutation -> fast.
"""
import json
import numpy as np

import os
ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
RES = os.environ.get("V14_RES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))

EV = json.load(open(f"{ROOT}/evidence_layers_v11.json"))
LAYERS = EV["layers"]
ALL_GENES = EV["genes"]

DEP_W = {"string_centrality": 0.80, "mutation_freq": 0.10, "impc_animal_ko": 0.10}
SUPPORT = ["cancer_driver", "ot_genetics_pdac", "druggability", "hpa_pdac_prognostic", "hpa_rna_tissue_spec"]

def lv(g, L):
    rec = LAYERS.get(L, {}).get(g)
    return rec["norm"] if (rec and rec.get("present")) else None

def ecs_alpha(g, alpha):
    num = 0.0; den = 0.0
    for L, w in DEP_W.items():
        v = lv(g, L)
        if v is not None:
            num += w * v; den += w
    if den == 0:
        return None
    D = num / den
    vs = [lv(g, L) for L in SUPPORT if lv(g, L) is not None]
    PHI = float(np.mean(vs)) if vs else 0.0
    return D * (1.0 + alpha * PHI)

# endpoints
dpd = json.load(open(f"{ROOT}/depmap_pdac_dependency.json"))
E1 = {g: (bool(d.get("essential")) if isinstance(d, dict) else False) for g, d in dpd.items()}
DRUG = {g: lv(g, "druggability") for g in ALL_GENES}
E3 = {g: (E1.get(g, False) and DRUG.get(g) is not None) for g in ALL_GENES}
# Paired universe for E3 (same as the main analysis): ECS and STRING both present.
PAIRED_E3 = set(g for g in ALL_GENES if ecs_alpha(g, 0.6) is not None and lv(g, "string_centrality") is not None)

def auroc(scores, labels):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    if len(scores) < 2: return None
    order = np.argsort(scores, kind="mergesort")
    s = scores[order]; lab = labels[order]
    ranks = np.empty(len(s)); i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[i:j + 1] = (i + j) / 2.0 + 1.0
        i = j + 1
    npos = int(lab.sum()); nneg = len(lab) - npos
    if npos == 0 or nneg == 0: return None
    sumr = ranks[lab == 1].sum()
    return (sumr - npos * (npos + 1) / 2.0) / (npos * nneg)

def auc_of_scorer_dict(sc, ep):
    genes = PAIRED_E3 if ep == "E3_actionable_target" else ALL_GENES
    S = []; L = []
    for g in genes:
        s = sc.get(g); l = E3.get(g) if ep == "E3_actionable_target" else E1.get(g)
        if s is None or l is None: continue
        S.append(s); L.append(1 if l else 0)
    if len(S) < 10: return None
    return auroc(np.asarray(S, float), np.asarray(L, int))

alpha_grid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
alpha_sens = {}
for a in alpha_grid:
    sc = {g: ecs_alpha(g, a) for g in ALL_GENES}
    alpha_sens[f"a{a}"] = {
        "E1_pan_dependency": auc_of_scorer_dict(sc, "E1_pan_dependency"),
        "E3_actionable_target": auc_of_scorer_dict(sc, "E3_actionable_target"),
    }
    print(f"  alpha={a:.1f}  E1_AUROC={alpha_sens[f'a{a}']['E1_pan_dependency']:.3f}  E3_AUROC={alpha_sens[f'a{a}']['E3_actionable_target']:.3f}")

# VALIDATION: alpha=0.6 E3 must match the paired benchmark ECS E3 (~0.812)
val = alpha_sens["a0.6"]["E3_actionable_target"]
assert abs(val - 0.812) < 0.005, f"alpha=0.6 E3={val:.4f} inconsistent with paired benchmark 0.812!"
print(f"  [ok] alpha=0.6 E3_AUROC={val:.4f} matches paired benchmark 0.812")

# patch benchmark_v14.json
bp = json.load(open(f"{RES}/benchmark_v14.json"))
bp["alpha_sensitivity"] = alpha_sens
json.dump(bp, open(f"{RES}/benchmark_v14.json", "w"), indent=1, default=str)
print("Patched benchmark_v14.json with alpha_sensitivity.")
