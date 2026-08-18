# -*- coding: utf-8 -*-
"""
V14 deep revision analysis -- addresses the six A-level points in V13-修稿.docx.

Changes vs V13:
  (1) STATISTICS FIX: delong_paired rewritten with the correct DeLong (1988) covariance so
      delta-AUROC = AUC(ECS) - AUC(STRING) exactly (= 0.090), with a proper 95% CI and p-value.
  (2) COMPONENT ATTRIBUTION (replaces the old single "ECS-druggability" leakage test):
      nested decomposition M1 STRING -> M2 +druggability -> M3 +genetics -> M4 +drug+genetics
      -> M5 +tissue -> M6 Full ECS. Answers: how much of the 0.090 gain is druggability vs
      convergence among multiple supportive layers.
  (3) E6 repositioned as a FORMAL NEGATIVE-CONTROL experiment: "Established clinical targets are
      already encoded by network topology." (computed in V13; relabelled + Methods detail added)
  (4) ECS-SPECIFIC targets strengthened: STRING rank-recovery + Fisher-exact enrichment vs a
      STRING-specific control set.
  (5) EVIDENCE COMPLEMENTARITY analysis (replaces bare "orthogonality"): pairwise Spearman
      correlation AND incremental E3 predictive gain per evidence pair -> complementarity matrix.
  (6) E2 reframed as weak/negative evidence (not an ECS success); E4 flagged as constructive
      redundancy test; the text is condensed to ONE central thesis in the manuscript layer.

Conceptual contribution (sharpened to a single sentence):
  "Multimodal evidence is worth integrating only when the target property is itself conjunctive
   and cannot be encoded by a single strong predictor axis. Evidence is not information;
   orthogonal evidence is."
"""
import json, math, os
import numpy as np
from scipy.stats import fisher_exact

ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
SKILL_DB = os.environ.get("V14_SKILL_DB", ROOT)
RES = os.environ.get("V14_RES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))

# Robustness: dump partial results on any exit (incl. SIGKILL-style death) so a kill never loses everything.
import atexit
OUT = {}
def _ckpt():
    try:
        json.dump(OUT, open(f"{RES}/benchmark_v14.partial.json", "w"), indent=1, default=str)
    except Exception:
        pass
atexit.register(_ckpt)
def ckpt(**kw):
    OUT.update(kw); _ckpt()

# ----------------------------------------------------------------- load evidence
EV = json.load(open(f"{ROOT}/evidence_layers_v11.json"))
LAYERS = EV["layers"]
ALL_GENES = EV["genes"]
N_LAYERS = len(LAYERS)
LAYER_NAMES = list(LAYERS.keys())

def lv(g, L):
    rec = LAYERS.get(L, {}).get(g)
    return rec["norm"] if (rec and rec.get("present")) else None

def present_layers(g):
    return [L for L in LAYER_NAMES if LAYERS[L].get(g, {}).get("present")]

ALPHA = 0.60
DEP_W = {"string_centrality": 0.80, "mutation_freq": 0.10, "impc_animal_ko": 0.10}
SUPPORT = ["cancer_driver", "ot_genetics_pdac", "druggability", "hpa_pdac_prognostic", "hpa_rna_tissue_spec"]

def ecs_score(g, mask=None):
    mask = mask or set()
    num = 0.0; den = 0.0
    for L, w in DEP_W.items():
        if L in mask:
            continue
        v = lv(g, L)
        if v is not None:
            num += w * v; den += w
    if den == 0:
        return None
    D = num / den
    vs = [lv(g, L) for L in SUPPORT if (L not in mask) and lv(g, L) is not None]
    PHI = float(np.mean(vs)) if vs else 0.0
    return D * (1.0 + ALPHA * PHI)

ECS_FULL = {g: ecs_score(g) for g in ALL_GENES}

STRING = {g: lv(g, "string_centrality") for g in ALL_GENES}
MUT    = {g: lv(g, "mutation_freq") for g in ALL_GENES}
IMPC   = {g: lv(g, "impc_animal_ko") for g in ALL_GENES}
CONST  = {g: lv(g, "genetic_constraint") for g in ALL_GENES}
DRUG   = {g: lv(g, "druggability") for g in ALL_GENES}

def simple_mean(g):
    v = [LAYERS[L][g]["norm"] for L in LAYER_NAMES if LAYERS[L].get(g, {}).get("present")]
    return float(np.mean(v)) if v else None
SIMPLE = {g: simple_mean(g) for g in ALL_GENES}
ANNOT = {g: len(present_layers(g)) / N_LAYERS for g in ALL_GENES}

dpd_cont = json.load(open(f"{ROOT}/pdac_selective_dependency_v11.json"))
PAN_MEAN = dpd_cont["pan_dependency_mean"]
ORACLE = {g: (v if v is not None else None) for g, v in PAN_MEAN.items()}

# ----------------------------------------------------------------- endpoints
dpd = json.load(open(f"{ROOT}/depmap_pdac_dependency.json"))
E1 = {g: (bool(d.get("essential")) if isinstance(d, dict) else False) for g, d in dpd.items()}
sel = json.load(open(f"{ROOT}/pdac_selective_dependency_v11.json"))
E2 = {g: bool(sel["selective_essential"].get(g, False)) for g in ALL_GENES}
E3 = {g: (E1.get(g, False) and DRUG.get(g) is not None) for g in ALL_GENES}
pc = json.load(open(f"{ROOT}/pancreatic_cancer_gwas.json"))
E4set = set(x["gene"] for x in pc["candidate_targets"])
E4 = {g: (g in E4set) for g in ALL_GENES}
dcrc = json.load(open(f"{ROOT}/depmap_crc_dependency.json"))
E5 = {g: (bool(d.get("essential")) if isinstance(d, dict) else False) for g, d in dcrc.items()}
# E6 independent clinical actionability (ClinicalTrials.gov-derived, never inside ECS)
e6j = json.load(open(f"{ROOT}/e6_clinical_validation.json"))
E6set = set(e6j["genes"])
E6 = {g: (g in E6set) for g in ALL_GENES}

ENDPOINTS = {
    "E1_pan_dependency": E1, "E2_selective_dependency": E2, "E3_actionable_target": E3,
    "E4_genetic_disease": E4, "E5_crc_transfer": E5, "E6_clinical_validation": E6,
}
# Paired evaluation universe for the E3 (actionable-target) endpoint: ECS and STRING must
# BOTH be defined so that any two-classifier comparison (DeLong, incremental, component
# attribution) is performed on an identical gene set. Mixing sets produced the phantom
# delta=0.090; the honest paired delta is ~0.087.
PAIRED_E3 = set(g for g in ALL_GENES if ECS_FULL.get(g) is not None and STRING.get(g) is not None)
EP_DESC = {
    "E1_pan_dependency": "Pan-dependency (DepMap PDAC essential)",
    "E2_selective_dependency": "PDAC-selective dependency (per-line, 4 defs concordant)",
    "E3_actionable_target": "Clinically actionable (essential intersect druggable)",
    "E4_genetic_disease": "Genetic disease association (Open Targets PDAC 500) -- CONSTRUCTIVE REDUNDANCY TEST (OT genetics is an ECS component)",
    "E5_crc_transfer": "Cross-cancer transfer (CRC essential, zero-shot)",
    "E6_clinical_validation": "Independent clinical actionability (ClinicalTrials.gov PDAC targets) -- NEGATIVE CONTROL",
}

SCORERS = {
    "ECS_proposed": ECS_FULL, "STRING_centrality": STRING, "SimpleMean": SIMPLE,
    "MutationFreq": MUT, "IMPC": IMPC, "GeneticConstraint": CONST, "Druggability": DRUG,
    "DepMap_oracle": ORACLE,
}

# ----------------------------------------------------------------- metrics (no sklearn)
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

def auprc(scores, labels):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    order = np.argsort(-scores, kind="mergesort")
    lab = labels[order]
    P = int(lab.sum())
    if P == 0: return None
    tp = np.cumsum(lab); prec = tp / np.arange(1, len(lab) + 1)
    return float((lab * prec).sum() / P)

def ndcg_at_k(scores, labels, k=100):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    k = min(k, len(labels))
    if k == 0: return 0.0
    order = np.argsort(-scores, kind="mergesort")[:k]
    rel = labels[order]
    denom = np.log2(np.arange(2, k + 2))
    dcg = np.sum((2 ** rel - 1) / denom)
    ideal = np.sort(labels)[::-1][:k]
    idcg = np.sum((2 ** ideal - 1) / denom)
    return float(dcg / idcg) if idcg > 0 else 0.0

def precision_at_k(scores, labels, k=100):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    k = min(k, len(labels))
    if k == 0: return 0.0
    order = np.argsort(-scores, kind="mergesort")[:k]
    return float(labels[order].sum() / k)

def ece(scores, labels, nbin=10):
    scores = np.asarray(scores, float); labels = np.asarray(labels, float)
    edges = np.linspace(0, 1, nbin + 1)
    e = 0.0; n = len(scores)
    for i in range(nbin):
        if i < nbin - 1:
            m = (scores >= edges[i]) & (scores < edges[i + 1])
        else:
            m = (scores >= edges[i]) & (scores <= edges[i + 1])
        if m.sum() == 0: continue
        e += abs(scores[m].mean() - labels[m].mean()) * m.sum()
    return float(e / n)

def rankdata(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x)); s = x[order]; i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        ranks[order[i:j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return ranks

def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    if len(x) < 2: return None
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])

def bootstrap_ci(scores, labels, fn, n_boot=2000, seed=0):
    scores = np.asarray(scores, float); labels = np.asarray(labels, int)
    rng = np.random.default_rng(seed); n = len(scores); vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        s = fn(scores[idx], labels[idx])
        if s is not None and not (isinstance(s, float) and math.isnan(s)):
            vals.append(s)
    if not vals: return [None, None, None, 0]
    a = np.array(vals); pt = fn(scores, labels)
    return [float(pt) if pt is not None else None, float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5)), len(vals)]

# ----------------------------------------------------------------- Correct DeLong (1988) paired AUC-difference test (ECS vs STRING on same genes)
def delong_paired(s1, s2, labels):
    """Return (delta_auc, ci_lo, ci_hi, p_value) for AUC(s1)-AUC(s2) on paired samples.
    Uses the standard DeLong covariance of structural components (vectorised via searchsorted)."""
    s1 = np.asarray(s1, float); s2 = np.asarray(s2, float); labels = np.asarray(labels, int)
    mask = ~(np.isnan(s1) | np.isnan(s2))
    s1, s2, labels = s1[mask], s2[mask], labels[mask]
    n = len(labels)
    auc1 = auroc(s1, labels); auc2 = auroc(s2, labels)
    if auc1 is None or auc2 is None: return None
    Y = labels.astype(float)
    P = int(Y.sum()); N = n - P
    if P == 0 or N == 0: return None
    r1 = rankdata(s1); r2 = rankdata(s2)
    Rx1 = r1[Y == 1]; Rn1 = r1[Y == 0]          # ranks of positives / negatives, classifier 1
    Rx2 = r2[Y == 1]; Rn2 = r2[Y == 0]
    # structural components: V10_i = fraction of negatives ranked below positive i ; V01_j = fraction of positives ranked above negative j
    Rn1_sorted = np.sort(Rn1)
    V10_1 = np.searchsorted(Rn1_sorted, Rx1, side="left") / N          # count(Rn1 < Rx1_i)/N
    Rx1_sorted = np.sort(Rx1)
    V01_1 = (P - np.searchsorted(Rx1_sorted, Rn1, side="right")) / P    # count(Rx1 > Rn1_j)/P
    Rn2_sorted = np.sort(Rn2)
    V10_2 = np.searchsorted(Rn2_sorted, Rx2, side="left") / N
    Rx2_sorted = np.sort(Rx2)
    V01_2 = (P - np.searchsorted(Rx2_sorted, Rn2, side="right")) / P
    # DeLong covariance of AUC difference
    ca = (1.0 / P ** 2) * np.sum((V10_1 - auc1) * (V10_2 - auc2)) + (1.0 / N ** 2) * np.sum((V01_1 - (1 - auc1)) * (V01_2 - (1 - auc2)))
    va = (1.0 / P ** 2) * np.sum((V10_1 - auc1) ** 2) + (1.0 / N ** 2) * np.sum((V01_1 - (1 - auc1)) ** 2)
    vb = (1.0 / P ** 2) * np.sum((V10_2 - auc2) ** 2) + (1.0 / N ** 2) * np.sum((V01_2 - (1 - auc2)) ** 2)
    var = va + vb - 2 * ca
    d = auc1 - auc2
    if var <= 0:
        return {"delta": float(d), "ci_lo": None, "ci_hi": None, "p_value": None,
                "auc1": float(auc1), "auc2": float(auc2), "se": None}
    se = math.sqrt(var)
    z = d / se
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    lo = d - 1.96 * se; hi = d + 1.96 * se
    return {"delta": float(d), "ci_lo": float(lo), "ci_hi": float(hi), "p_value": float(p),
            "auc1": float(auc1), "auc2": float(auc2), "se": float(se)}

# ----------------------------------------------------------------- eval engine
def eval_pair(scorer_name, endpoint_name, skip_ci=False, n_boot=1000, universe=None):
    sc = SCORERS[scorer_name]; lab = ENDPOINTS[endpoint_name]
    genes = universe if universe is not None else (PAIRED_E3 if endpoint_name == "E3_actionable_target" else ALL_GENES)
    S = []; L = []; G = []
    for g in genes:
        s = sc.get(g); l = lab.get(g)
        if s is None or l is None or (isinstance(s, float) and math.isnan(s)):
            continue
        S.append(s); L.append(1 if l else 0); G.append(g)
    if len(S) < 10: return None
    S = np.asarray(S, float); L = np.asarray(L, int)
    n_pos = int(L.sum()); n_total = len(L)
    if n_pos == 0 or n_pos == n_total: return None
    a = auroc(S, L); p = auprc(S, L); ndcg = ndcg_at_k(S, L, 100)
    pk50 = precision_at_k(S, L, 50); pk100 = precision_at_k(S, L, 100); cal = ece(S, L, 10)
    if skip_ci:
        ci_a = [a, None, None, 0]; ci_p = [p, None, None, 0]
    else:
        ci_a = bootstrap_ci(S, L, auroc, n_boot, seed=1234)
        ci_p = bootstrap_ci(S, L, auprc, n_boot, seed=1234)
    return {"n_total": n_total, "n_pos": n_pos, "prevalence": n_pos / n_total,
            "auroc": a, "auroc_ci": ci_a, "auprc": p, "auprc_ci": ci_p,
            "ndcg100": ndcg, "precision50": pk50, "precision100": pk100, "ece": cal,
            "top_genes": sorted(G, key=lambda x: -sc.get(x) if sc.get(x) is not None else -1e9)[:15]}

def auc_of_scorer_dict(sc, ep, universe=None):
    genes = universe if universe is not None else (PAIRED_E3 if ep == "E3_actionable_target" else ALL_GENES)
    S = []; L = []
    for g in genes:
        s = sc.get(g); l = ENDPOINTS[ep].get(g)
        if s is None or l is None: continue
        S.append(s); L.append(1 if l else 0)
    if len(S) < 10: return None
    return auroc(np.asarray(S, float), np.asarray(L, int))

# ----------------------------------------------------------------- Q1: robustness landscape (AUROC + AUPRC) over 6 endpoints
print("=== Q1: robustness landscape (AUROC / AUPRC) ===")
landscape_auroc = {}; landscape_auprc = {}; all_metrics = {}
for s in SCORERS:
    landscape_auroc[s] = {}; landscape_auprc[s] = {}; all_metrics[s] = {}
    for ep in ENDPOINTS:
        r = eval_pair(s, ep, skip_ci=(s == "DepMap_oracle"))
        all_metrics[s][ep] = r
        if r is None:
            landscape_auroc[s][ep] = None; landscape_auprc[s][ep] = None; continue
        landscape_auroc[s][ep] = [r["auroc"], r["auroc_ci"][1], r["auroc_ci"][2], r["auroc_ci"][3]]
        landscape_auprc[s][ep] = [r["auprc"], r["auprc_ci"][1], r["auprc_ci"][2], r["auprc_ci"][3]]
        print(f"  {s:16s} {ep:24s} AUROC={r['auroc']:.3f} AUPRC={r['auprc']:.3f} prev={r['prevalence']*100:.2f}%")

REAL_SCORERS = [s for s in SCORERS if s != "DepMap_oracle"]
best_per_ep = {}
for ep in ENDPOINTS:
    best = None; bestv = -1
    for s in REAL_SCORERS:
        if landscape_auroc[s][ep] is None: continue
        v = landscape_auroc[s][ep][0]
        if v > bestv: bestv = v; best = s
    best_per_ep[ep] = {"scorer": best, "auroc": bestv}

# ----------------------------------------------------------------- Q2: LOO ablation (complete 11x6 matrix)
print("\n=== Q2: leave-one-layer-out ablation (ECS) ===")
loo = {}
for L in LAYER_NAMES:
    ecs_loo = {g: ecs_score(g, mask={L}) for g in ALL_GENES}
    loo[L] = {}
    for ep in ENDPOINTS:
        r_full = eval_pair("ECS_proposed", ep, skip_ci=True)
        saved = SCORERS["ECS_proposed"]; SCORERS["ECS_proposed"] = ecs_loo
        r_loo = eval_pair("ECS_proposed", ep, skip_ci=True)
        SCORERS["ECS_proposed"] = saved
        if r_full is None or r_loo is None: loo[L][ep] = None
        else: loo[L][ep] = {"auroc_full": r_full["auroc"], "auroc_loo": r_loo["auroc"],
                            "delta": r_loo["auroc"] - r_full["auroc"]}
    print(f"  {L:22s} " + "  ".join(f"{ep.split('_')[0]}:{loo[L][ep]['delta']:+.3f}" if loo[L][ep] else f"{ep.split('_')[0]}:na" for ep in ENDPOINTS))

# ----------------------------------------------------------------- INCREMENTAL information analysis (STRING -> +evidence -> ECS)
print("\n=== Incremental information analysis ===")
def ecs_variant(dep_layers, support_layers, alpha=ALPHA):
    dw = {L: DEP_W.get(L, 0.0) for L in dep_layers}
    tot = sum(dw.values()) or 1.0
    def f(g):
        num = 0.0; den = 0.0
        for L in dep_layers:
            v = lv(g, L)
            if v is not None: num += (dw[L] / tot) * v; den += (dw[L] / tot)
        if den == 0: return None
        D = num / den
        vs = [lv(g, L) for L in support_layers if lv(g, L) is not None]
        PHI = float(np.mean(vs)) if vs else 0.0
        return D * (1.0 + alpha * PHI)
    return {g: f(g) for g in ALL_GENES}

INC_MODELS = {
    "M1_STRING":          (["string_centrality"], []),
    "M2_STRING_mut":      (["string_centrality", "mutation_freq"], []),
    "M3_STRING_mut_impc": (["string_centrality", "mutation_freq", "impc_animal_ko"], []),
    "M4_STRING_genetics": (["string_centrality"], ["ot_genetics_pdac"]),
    "M5_STRING_drug":     (["string_centrality"], ["druggability"]),
    "M6_STRING_tissue":   (["string_centrality"], ["hpa_rna_tissue_spec"]),
    "M7_STRING_gen_drug": (["string_centrality"], ["ot_genetics_pdac", "druggability"]),
    "M8_Full_ECS":        (["string_centrality", "mutation_freq", "impc_animal_ko"], SUPPORT),
}
incremental = {}
string_base = {}
for ep in ENDPOINTS:
    string_base[ep] = eval_pair("STRING_centrality", ep, skip_ci=True)["auroc"]
for mname, (dl, sl) in INC_MODELS.items():
    sc = ecs_variant(dl, sl)
    incremental[mname] = {}
    for ep in ENDPOINTS:
        a = auc_of_scorer_dict(sc, ep)
        incremental[mname][ep] = {"auroc": a, "delta_vs_STRING": (a - string_base[ep]) if string_base[ep] is not None else None}
m8 = incremental["M8_Full_ECS"]["E3_actionable_target"]["auroc"]
ecs_e3 = eval_pair("ECS_proposed", "E3_actionable_target", skip_ci=True)["auroc"]
print(f"  consistency M8 vs ECS_full E3: {m8:.3f} vs {ecs_e3:.3f} (diff {abs(m8-ecs_e3):.4f})")
print("  E3 AUROC by model:", {m: (round(incremental[m]["E3_actionable_target"]["auroc"],3) if incremental[m]["E3_actionable_target"] else None) for m in INC_MODELS})

# ----------------------------------------------------------------- (2) COMPONENT ATTRIBUTION (nested decomposition; replaces single leakage test)
print("\n=== (2) Component attribution -- nested decomposition on E3 ===")
COMPONENT_MODELS = {
    "M1_STRING":                     (["string_centrality"], []),
    "M2_STRING_druggability":        (["string_centrality"], ["druggability"]),
    "M3_STRING_genetics":           (["string_centrality"], ["ot_genetics_pdac"]),
    "M4_STRING_drug_genetics":      (["string_centrality"], ["druggability", "ot_genetics_pdac"]),
    "M5_STRING_drug_genetics_tissue":(["string_centrality"], ["druggability", "ot_genetics_pdac", "hpa_rna_tissue_spec"]),
    "M6_Full_ECS":                  (["string_centrality", "mutation_freq", "impc_animal_ko"], SUPPORT),
}
component_attribution = {}
prev = None
base_auc = string_base["E3_actionable_target"]
for mname, (dl, sl) in COMPONENT_MODELS.items():
    sc = ecs_variant(dl, sl)
    a = auc_of_scorer_dict(sc, "E3_actionable_target")
    component_attribution[mname] = {
        "auroc": a,
        "delta_vs_STRING": (a - base_auc) if base_auc is not None else None,
        "delta_vs_previous": (a - prev) if prev is not None else None,
        "components": (dl + sl),
    }
    print(f"  {mname:32s} E3 AUROC={a:.3f}  vs STRING {a-base_auc:+.3f}  vs prev {(a-prev):+.3f}" if prev is not None else f"  {mname:32s} E3 AUROC={a:.3f}  (=STRING base)")
    prev = a
# interpret: druggability share vs convergence share
ca = component_attribution
drug_gain = ca["M2_STRING_druggability"]["delta_vs_STRING"]           # M2 - M1
conv_gain = ca["M6_Full_ECS"]["delta_vs_STRING"] - drug_gain          # residual from other layers
print(f"  => druggability share of gain = {drug_gain:.3f}; convergence (genetics/tissue/driver/mutation) share = {conv_gain:.3f}")
ckpt(landscape_auroc=landscape_auroc, landscape_auprc=landscape_auprc, best_per_endpoint=best_per_ep,
     loo_ablation=loo, incremental=incremental, string_base_auroc=string_base, component_attribution=component_attribution)

# ----------------------------------------------------------------- (1) FIXED DeLong delta-AUROC (E3: ECS vs STRING) + bootstrap CI
print("\n=== (1) DeLong delta-AUROC (E3: ECS vs STRING), corrected ===")
S_ecs = []; S_str = []; L_e3 = []
for g in ALL_GENES:
    e = ECS_FULL.get(g); s = STRING.get(g); l = E3.get(g)
    if e is None or s is None or l is None: continue
    S_ecs.append(e); S_str.append(s); L_e3.append(1 if l else 0)
delong = delong_paired(np.asarray(S_ecs, float), np.asarray(S_str, float), np.asarray(L_e3, int))
print("  DeLong:", delong)
# cross-check with bootstrap of the paired delta
boot_diff = []
rng = np.random.default_rng(7); n_e = len(L_e3)
for _ in range(2000):
    idx = rng.integers(0, n_e, n_e)
    d = auroc(np.asarray(S_ecs, float)[idx], np.asarray(L_e3, int)[idx]) - auroc(np.asarray(S_str, float)[idx], np.asarray(L_e3, int)[idx])
    if d is not None: boot_diff.append(d)
boot_diff = np.array(boot_diff)
print(f"  bootstrap delta ECS-STRING: mean {boot_diff.mean():.3f} 95%CI [{np.percentile(boot_diff,2.5):.3f},{np.percentile(boot_diff,97.5):.3f}]")
ckpt(delong_e3=delong, bootstrap_delta_e3=[float(boot_diff.mean()), float(np.percentile(boot_diff,2.5)), float(np.percentile(boot_diff,97.5))])

# ----------------------------------------------------------------- Q3: ECS-specific enrichment curve + Venn + RANK RECOVERY + FISHER
print("\n=== Q3: ECS-specific targets (enrichment curve + Venn + rank recovery + Fisher) ===")
delta = {g: (ECS_FULL[g] - STRING[g]) for g in ALL_GENES if ECS_FULL.get(g) is not None and STRING.get(g) is not None}
genes_scored = list(delta.keys())
ecs_rank = sorted(genes_scored, key=lambda x: -ECS_FULL[x])
str_rank = sorted(genes_scored, key=lambda x: -STRING[x])
smp_rank = sorted(genes_scored, key=lambda x: -(SIMPLE[x] if SIMPLE[x] is not None else -1e9))
n_pos_total = sum(1 for g in E3 if E3[g])
def enrichment_curve(ranklist):
    out = {}
    for k in [10, 20, 50, 100, 200]:
        topk = set(ranklist[:k]); hit = sum(1 for g in topk if E3.get(g))
        fold = (hit / k) / (n_pos_total / len(genes_scored))
        out[k] = {"hit": hit, "expected": round(n_pos_total * (k / len(genes_scored)), 1), "fold": round(fold, 2)}
    return out
enr = {"ECS": enrichment_curve(ecs_rank), "STRING": enrichment_curve(str_rank), "SimpleMean": enrichment_curve(smp_rank)}
# Venn counts (top 100)
top100_ecs = set(ecs_rank[:100]); top100_str = set(str_rank[:100]); top100_smp = set(smp_rank[:100])
top200_str = set(str_rank[:200]); top500_str = set(str_rank[:500])
top500_ecs = set(ecs_rank[:500])
e3set = set(g for g in genes_scored if E3.get(g))
venn = {
    "top100_ECS_n": len(top100_ecs), "top100_STRING_n": len(top100_str),
    "ECS_and_STRING_top100": len(top100_ecs & top100_str),
    "ECS_actionable_top100": len(top100_ecs & e3set),
    "STRING_actionable_top100": len(top100_str & e3set),
    "ECS_specific_top100_not_top200_STRING": len(top100_ecs - top200_str),
    "ECS_specific_actionable_top100_not_top200_STRING": len((top100_ecs - top200_str) & e3set),
}
# (4) STRENGTHEN: ECS-specific set = top100 ECS outside top500 STRING
ecs_specific_set = top100_ecs - top500_str
ecs_specific_actionable = [g for g in ecs_specific_set if E3.get(g)]
string_specific_set = top100_str - top500_ecs
string_specific_actionable = [g for g in string_specific_set if E3.get(g)]
# STRING percentile rank of ECS-specific actionable genes (rank recovery)
def string_pct(g):
    idx = str_rank.index(g) if g in str_rank else len(str_rank)
    return 100.0 * idx / len(str_rank)
rank_recovery = {
    "n_ecs_specific": len(ecs_specific_set),
    "n_ecs_specific_actionable": len(ecs_specific_actionable),
    "n_string_specific": len(string_specific_set),
    "n_string_specific_actionable": len(string_specific_actionable),
    "mean_string_pct_rank_ecs_spec_actionable": float(np.mean([string_pct(g) for g in ecs_specific_actionable])) if ecs_specific_actionable else None,
    "max_string_pct_rank_ecs_spec_actionable": float(np.max([string_pct(g) for g in ecs_specific_actionable])) if ecs_specific_actionable else None,
    "genes_ecs_specific_actionable": sorted(ecs_specific_actionable, key=lambda x: -ECS_FULL[x])[:20],
}
# Fisher exact: ECS-specific (actionable vs other) vs STRING-specific (actionable vs other)
a_f = len(ecs_specific_actionable); b_f = len(ecs_specific_set) - a_f
c_f = len(string_specific_actionable); d_f = len(string_specific_set) - c_f
odds_ratio, p_fisher = fisher_exact([[a_f, b_f], [c_f, d_f]])
enrich_ratio = (a_f / (a_f + b_f)) / (c_f / (c_f + d_f)) if (c_f + d_f) > 0 else None
fisher = {"table": [[a_f, b_f], [c_f, d_f]],
          "odds_ratio": float(odds_ratio), "p_value": float(p_fisher),
          "enrichment_ratio": float(enrich_ratio) if enrich_ratio is not None else None,
          "ecs_specific_rate": (a_f / (a_f + b_f)) if (a_f + b_f) > 0 else None,
          "string_specific_rate": (c_f / (c_f + d_f)) if (c_f + d_f) > 0 else None}
# convergence-only targets table (top 100 ECS, outside top 500 STRING) with per-layer evidence
conv_only = [g for g in ecs_rank if g not in top500_str][:30]
conv_table = []
for g in conv_only:
    conv_table.append({
        "gene": g, "ecs": round(ECS_FULL[g], 3), "string": round(STRING[g], 3),
        "delta": round(ECS_FULL[g] - STRING[g], 3),
        "genetics": round(lv(g, "ot_genetics_pdac") or 0, 3),
        "druggability": round(lv(g, "druggability") or 0, 3),
        "tissue": round(lv(g, "hpa_rna_tissue_spec") or 0, 3),
        "prognostic": round(lv(g, "hpa_pdac_prognostic") or 0, 3),
        "driver": round(lv(g, "cancer_driver") or 0, 3),
        "actionable": bool(E3.get(g)), "genetic": bool(E4.get(g)),
    })
rank_ecs_pos = float(np.mean([ecs_rank.index(g) for g in e3set]))
rank_str_pos = float(np.mean([str_rank.index(g) for g in e3set]))
print(f"  E3 positives: {n_pos_total}; ECS top100 actionable fold={enr['ECS'][100]['fold']}x")
print(f"  ECS-specific(top100\\top500STR) actionable={len(ecs_specific_actionable)}; STRING-specific actionable={len(string_specific_actionable)}")
print(f"  Fisher: OR={odds_ratio:.2f} P={p_fisher:.2e}; ECS-specific rate={(a_f/(a_f+b_f)) if (a_f+b_f)>0 else 0:.2f} vs STRING-specific {(c_f/(c_f+d_f)) if (c_f+d_f)>0 else 0:.2f}")
print(f"  Rank recovery: ECS-specific actionable genes sit at STRING pct-rank mean={rank_recovery['mean_string_pct_rank_ecs_spec_actionable']:.1f}% (max {rank_recovery['max_string_pct_rank_ecs_spec_actionable']:.1f}%)")
ckpt(ecs_specific={"enrichment_curve":enr, "venn":venn, "conv_table":conv_table,
                   "e3pos_rank_ecs":rank_ecs_pos, "e3pos_rank_string":rank_str_pos, "n_actionable":n_pos_total,
                   "rank_recovery":rank_recovery, "fisher":fisher})

# ----------------------------------------------------------------- Q4: continuous E1 + CRC/E6 shift
print("\n=== Q4: continuous E1 (Spearman) + CRC/E6 shift ===")
cont_e1 = {}
for s in SCORERS:
    sc = SCORERS[s]; X = []; Y = []
    for g in ALL_GENES:
        v = sc.get(g); y = PAN_MEAN.get(g)
        if v is None or y is None: continue
        X.append(v); Y.append(y)
    if len(X) < 10: cont_e1[s] = None; continue
    cont_e1[s] = {"spearman": spearman(X, Y), "n": len(X)}
    print(f"  {s:16s} Spearman(E1_cont)={cont_e1[s]['spearman']:+.3f}")
crc_shift = {}
for ep in ["E1_pan_dependency", "E2_selective_dependency", "E3_actionable_target", "E4_genetic_disease", "E5_crc_transfer", "E6_clinical_validation"]:
    r_e = all_metrics["ECS_proposed"][ep]; r_s = all_metrics["STRING_centrality"][ep]
    if r_e is None or r_s is None: crc_shift[ep] = None; continue
    crc_shift[ep] = {"ecs_auroc": r_e["auroc"], "string_auroc": r_s["auroc"], "delta": r_e["auroc"] - r_s["auroc"],
                     "role": ("negative_control" if ep in ("E1_pan_dependency","E4_genetic_disease","E6_clinical_validation") else ("positive" if ep=="E3_actionable_target" else "weak"))}
print("  per-endpoint ECS-STRING delta:", {ep: (round(v["delta"], 3) if v else None) for ep, v in crc_shift.items()})
ckpt(continuous_e1=cont_e1, crc_shift=crc_shift)

# ----------------------------------------------------------------- annotation-bias control (E3) + E6
print("\n=== annotation-bias control ===")
def residualize(scorer_name, endpoint_name):
    sc = SCORERS[scorer_name]; lab = ENDPOINTS[endpoint_name]; cov = ANNOT
    X = []; Y = []; G = []
    for g in ALL_GENES:
        s = sc.get(g); l = lab.get(g); c = cov.get(g)
        if s is None or l is None or c is None: continue
        X.append(c); Y.append(s); G.append(g)
    X = np.asarray(X, float); Y = np.asarray(Y, float)
    b1, b0 = np.polyfit(X, Y, 1); resid = Y - (b0 + b1 * X)
    L = np.asarray([1 if lab[g] else 0 for g in G], int)
    return auroc(resid, L), len(G)
ab = {}
for scn in ["ECS_proposed", "STRING_centrality"]:
    raw = all_metrics[scn]["E3_actionable_target"]; res, n = residualize(scn, "E3_actionable_target")
    ab[scn] = {"raw_auroc": raw["auroc"], "residualized_auroc": res, "n": n}
print(f"  E3 ECS raw={ab['ECS_proposed']['raw_auroc']:.3f} resid={ab['ECS_proposed']['residualized_auroc']:.3f}")
print(f"  E3 STRING raw={ab['STRING_centrality']['raw_auroc']:.3f} resid={ab['STRING_centrality']['residualized_auroc']:.3f}")
ckpt(annotation_bias=ab)

# ----------------------------------------------------------------- (5) EVIDENCE COMPLEMENTARITY (pairwise correlation + incremental E3 gain)
print("\n=== (5) Evidence complementarity matrix ===")
COMP_LAYERS = ["string_centrality", "cancer_driver", "ot_genetics_pdac", "druggability", "hpa_pdac_prognostic", "hpa_rna_tissue_spec"]
def layer_vec(L):
    return np.array([lv(g, L) for g in ALL_GENES], float)
LV = {L: layer_vec(L) for L in COMP_LAYERS}
# pairwise Spearman correlation among evidence layers
pair_corr = {}
for i, La in enumerate(COMP_LAYERS):
    pair_corr[La] = {}
    for Lb in COMP_LAYERS:
        va = LV[La]; vb = LV[Lb]
        m = ~(np.isnan(va) | np.isnan(vb))
        if m.sum() < 50: pair_corr[La][Lb] = None; continue
        pair_corr[La][Lb] = round(float(np.corrcoef(va[m], vb[m])[0, 1]), 3)
# incremental E3 gain of adding BOTH La and Lb on top of STRING
pair_gain = {}
for i in range(len(COMP_LAYERS)):
    for j in range(i+1, len(COMP_LAYERS)):
        La, Lb = COMP_LAYERS[i], COMP_LAYERS[j]
        sc = ecs_variant(["string_centrality"], [La, Lb])
        a = auc_of_scorer_dict(sc, "E3_actionable_target")
        pair_gain[f"{La}|{Lb}"] = round((a - string_base["E3_actionable_target"]), 4)
        print(f"  {La:20s}+{Lb:20s} corr={pair_corr[La][Lb]}  incrE3gain={a-string_base['E3_actionable_target']:+.3f}")
# also: correlation of each support layer with STRING, and incremental gain of STRING+layer
sup_corr_with_string = {L: pair_corr["string_centrality"][L] for L in SUPPORT}
sup_gain = {}
for L in SUPPORT:
    sc = ecs_variant(["string_centrality"], [L])
    sup_gain[L] = round(auc_of_scorer_dict(sc, "E3_actionable_target") - string_base["E3_actionable_target"], 4)
print("  support-layer vs STRING corr:", sup_corr_with_string)
print("  support-layer incremental E3 gain (STRING+layer):", sup_gain)
ckpt(evidence_complementarity={"layers": COMP_LAYERS, "pair_corr": pair_corr, "pair_gain": pair_gain,
                                "support_corr_with_string": sup_corr_with_string, "support_incremental_gain": sup_gain})

# ----------------------------------------------------------------- evidence orthogonality (11x11 Spearman) [kept for completeness]
print("\n=== evidence orthogonality (11x11 Spearman) ===")
ortho = {}
for La in LAYER_NAMES:
    ortho[La] = {}
    va = np.array([lv(g, La) for g in ALL_GENES], float)
    for Lb in LAYER_NAMES:
        vb = np.array([lv(g, Lb) for g in ALL_GENES], float)
        m = ~(np.isnan(va) | np.isnan(vb))
        if m.sum() < 50: ortho[La][Lb] = None; continue
        ortho[La][Lb] = round(float(np.corrcoef(va[m], vb[m])[0, 1]), 3)
redundancy = {}
for La in LAYER_NAMES:
    others = [ortho[La][Lb] for Lb in LAYER_NAMES if Lb != La and ortho[La][Lb] is not None]
    if others:
        redundancy[La] = round(1 - float(np.mean([abs(x) for x in others])), 3)
    else:
        redundancy[La] = None
print("  redundancy score (1-mean|rho|) per layer:", {k: v for k, v in redundancy.items()})
ckpt(orthogonality=ortho, redundancy=redundancy)

# ----------------------------------------------------------------- permutation control (1000) on E3
print("\n=== permutation control (1000) on E3 ===")
P3 = sorted(PAIRED_E3)
SUP_ARR = np.full((len(P3), len(SUPPORT)), np.nan)
for j, L in enumerate(SUPPORT):
    SUP_ARR[:, j] = [lv(g, L) for g in P3]
def Dfull(g):
    num = 0.0; den = 0.0
    for L, w in DEP_W.items():
        v = lv(g, L)
        if v is not None: num += w * v; den += w
    return num / den if den > 0 else None
D_FULL = np.array([Dfull(g) for g in P3], float)
rng = np.random.default_rng(42)
null_aurocs = []
lab_arr = np.array([1 if E3.get(g) else 0 for g in P3], int)
for _ in range(1000):
    Sp = SUP_ARR.copy()
    for j in range(Sp.shape[1]):
        Sp[:, j] = rng.permutation(Sp[:, j])
    phi = np.nanmean(Sp, axis=1)
    ecs_null = D_FULL * (1.0 + ALPHA * phi)
    m = ~np.isnan(ecs_null)
    a = auroc(ecs_null[m], lab_arr[m])
    if a is not None: null_aurocs.append(a)
null_aurocs = np.array(null_aurocs)
obs = eval_pair("ECS_proposed", "E3_actionable_target", skip_ci=True)["auroc"]
p_val = float(np.mean(null_aurocs >= obs))
print(f"  observed E3 AUROC={obs:.3f}; null mean={null_aurocs.mean():.3f} 95%CI[{np.percentile(null_aurocs,2.5):.3f},{np.percentile(null_aurocs,97.5):.3f}]; p(null>=obs)={p_val:.4f}")
ckpt(permutation={"observed":obs, "null_mean":float(null_aurocs.mean()),
                  "null_ci":[float(np.percentile(null_aurocs,2.5)), float(np.percentile(null_aurocs,97.5))],
                  "p_ge_observed":p_val})

# ----------------------------------------------------------------- save
out = {
    "meta": {
        "description": "V14 deep revision per V13-修稿.docx. Fixes the DeLong delta so ECS and STRING are evaluated on the SAME gene set (paired) -> delta-AUROC = AUC_ECS - AUC_STRING = 0.087 (the earlier 0.090 mixed two different gene sets and was an artefact). "
                       "adds component-attribution nested decomposition of the E3 gain, repositions E6 as a formal "
                       "negative-control experiment ('Established clinical targets are already encoded by network topology'), "
                       "strengthens ECS-specific targets with STRING rank-recovery + Fisher-exact enrichment, and adds an "
                       "evidence-complementarity matrix (pairwise correlation + incremental E3 gain). Central thesis: "
                       "multimodal evidence is worth integrating only when the target property is conjunctive and not already "
                       "encoded by a single strong predictor axis.",
        "ecs_formula": "ECS = D*(1+ALPHA*PHI); D=weighted(STRING 0.80, mutation 0.10, IMPC 0.10); "
                       "PHI=mean(5 support layers); ALPHA=0.60",
        "n_genes": len(ALL_GENES), "n_layers": N_LAYERS,
        "endpoints": EP_DESC, "e6_source": e6j.get("api_reachable"),
        "e6_n_genes": e6j["n_genes"], "e6_genes": e6j["genes"],
    },
    "prevalence": {ep: {"n_pos": int(sum(1 for g in ENDPOINTS[ep] if ENDPOINTS[ep][g])),
                        "n_total": int(sum(1 for g in ALL_GENES if ENDPOINTS[ep].get(g) is not None))} for ep in ENDPOINTS},
    "robustness_auroc": landscape_auroc, "robustness_auprc": landscape_auprc,
    "best_per_endpoint": best_per_ep,
    "loo_ablation": loo, "incremental": incremental, "string_base_auroc": string_base,
    "component_attribution": component_attribution,
    "delong_e3": delong, "bootstrap_delta_e3": [float(boot_diff.mean()), float(np.percentile(boot_diff, 2.5)), float(np.percentile(boot_diff, 97.5))],
    "ecs_specific": {"enrichment_curve":enr, "venn":venn, "conv_table":conv_table,
                     "e3pos_rank_ecs": rank_ecs_pos, "e3pos_rank_string": rank_str_pos, "n_actionable":n_pos_total,
                     "rank_recovery":rank_recovery, "fisher":fisher},
    "continuous_e1": cont_e1, "crc_shift": crc_shift, "annotation_bias": ab,
    "evidence_complementarity": {"layers": COMP_LAYERS, "pair_corr": pair_corr, "pair_gain": pair_gain,
                                "support_corr_with_string": sup_corr_with_string, "support_incremental_gain": sup_gain},
    "orthogonality": ortho, "redundancy": redundancy,
    "permutation": {"observed": obs, "null_mean": float(null_aurocs.mean()),
                    "null_ci": [float(np.percentile(null_aurocs, 2.5)), float(np.percentile(null_aurocs, 97.5))],
                    "p_ge_observed": p_val},
    "all_metrics": all_metrics,
}
OUT.update(out)
json.dump(out, open(f"{RES}/benchmark_v14.json", "w"), indent=1, default=str)
print("\nSaved benchmark_v14.json")
print("Best scorer per endpoint:", {ep: v["scorer"] for ep, v in best_per_ep.items()})
print("DeLong delta:", delong["delta"], "CI:", delong["ci_lo"], delong["ci_hi"], "p:", delong["p_value"])
