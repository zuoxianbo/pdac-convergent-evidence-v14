#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_evidence_v11.py
====================================================================
Unified extraction of ~10 normalized evidence layers for the PDAC
(Evidence Convergence Score) V11 revision, directly addressing the
V10 reviewer's core demand: "exhaust ALL possible data for analysis".

Each layer is stored with a CONSISTENT schema:
    {gene_symbol: {"raw": <raw value>, "norm": <0..1 or signed>, "present": bool}}
and a per-layer meta block recording how normalization was done
(min/max/transform), so the downstream analysis + manuscript never
hard-codes a number.

Sources integrated (all already downloaded & on disk):
  1. STRING PPI degree-centrality            (network)      string_centrality.json
  2. PDAC somatic mutation frequency          (genetic)      paad_mutfreq_genomewide.json
  3. IMPC animal KO viability/lethality        (animal)       impc_layer.json
  4. Open Targets gnomAD genetic constraint   (constraint)   opentargets parquet -> ensg_symbol_map
  5. Open Targets mouse KO score              (animal)       opentargets parquet
  6. Open Targets cancer-driver flag          (driver)       opentargets parquet
  7. Open Targets tissue specificity          (tissue)       opentargets parquet
  8. HPA RNA tissue-specificity score         (tissue)       proteinatlas.tsv
  9. HPA PDAC cancer-prognostic signal         (cancer)       proteinatlas.tsv
 10. Open Targets genetics PDAC 500 targets    (genetic)      pancreatic_cancer.json
 11. Druggability (approved/small-molecule)    (drug)         drug_targets.json + OT parquet

Output: outputs/data_v9/evidence_layers_v11.json
"""
import json, os, csv, math, collections
import numpy as np

ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
SKILL_DB = os.environ.get("V14_SKILL_DB", ROOT)
PY3 = "/Users/zuoxianbo/.workbuddy/binaries/python/versions/3.13.12/bin/python3"

def load_json(p):
    with open(p) as f:
        return json.load(f)

# ---------------------------------------------------------------- helpers
def to_float(x):
    try:
        if x is None: return None
        if isinstance(x, str):
            x = x.strip()
            if x == "" or x.lower() in ("na","nan","none","null"): return None
        return float(x)
    except Exception:
        return None

def minmax_norm(vals):
    """vals: list of floats (may contain None). Returns (lo, hi, n)."""
    vs = [v for v in vals if v is not None]
    if not vs: return (0.0, 1.0, 0)
    return (min(vs), max(vs), len(vs))

def norm01(v, lo, hi):
    if hi == lo: return 0.5
    return (v - lo) / (hi - lo)

def build_layer_from_dict(raw_dict, key_for_value, transform="minmax"):
    """raw_dict: gene->payload. Build normalized layer dict[str,int]."""
    out = {}
    vals = []
    for g, payload in raw_dict.items():
        v = to_float(payload[key_for_value]) if isinstance(payload, dict) else to_float(payload)
        if v is None: 
            out[g] = {"raw": None, "norm": None, "present": False}
        else:
            vals.append(v)
            out[g] = {"raw": v, "norm": None, "present": True}
    lo, hi, n = minmax_norm(vals)
    for g, rec in out.items():
        if rec["present"]:
            rec["norm"] = norm01(rec["raw"], lo, hi)
    return out, {"transform": transform, "lo": lo, "hi": hi, "n": n}

# ---------------------------------------------------------------- 1. STRING centrality
print("[1/11] STRING PPI centrality ...")
sc = load_json(os.path.join(ROOT, "string_centrality.json"))
string_layer = {}
for g, p in sc.items():
    c = to_float(p.get("centrality"))
    if c is None: continue
    string_layer[g] = {"raw": c, "norm": None, "present": True}
lo, hi, n = minmax_norm([r["raw"] for r in string_layer.values()])
for g, r in string_layer.items():
    r["norm"] = norm01(r["raw"], lo, hi)

# ---------------------------------------------------------------- 2. PDAC mutation frequency
print("[2/11] PDAC somatic mutation frequency ...")
mf = load_json(os.path.join(ROOT, "paad_mutfreq_genomewide.json"))
mut_layer = {}
for g, p in mf.items():
    f = to_float(p.get("freq_pct"))
    if f is None: continue
    mut_layer[g] = {"raw": f, "norm": None, "present": True}
lo, hi, n = minmax_norm([r["raw"] for r in mut_layer.values()])
for g, r in mut_layer.items():
    r["norm"] = norm01(r["raw"], lo, hi)

# ---------------------------------------------------------------- 3. IMPC animal KO
print("[3/11] IMPC animal KO ...")
impc = load_json(os.path.join(ROOT, "impc_layer.json"))
impc_layer = {}
for g, p in impc.items():
    s = to_float(p.get("score"))
    if s is None: continue
    impc_layer[g] = {"raw": s, "norm": None, "present": True}
lo, hi, n = minmax_norm([r["raw"] for r in impc_layer.values()])
for g, r in impc_layer.items():
    r["norm"] = norm01(r["raw"], lo, hi)

# ---------------------------------------------------------------- 4-7. Open Targets parquet via ENSG->symbol
print("[4/11] Open Targets parquet (constraint/animal/tissue/driver) ...")
import pyarrow.parquet as pq
pf = pq.ParquetFile(os.path.join(ROOT, "opentargets_target_prioritisation_26.06.parquet"))
ensg_map = load_json(os.path.join(ROOT, "downloads_v11", "ensg_symbol_map.json"))
# reverse: one symbol per ensg (map is ensg->symbol already)
cols = ['targetId','geneticConstraint','mouseKOScore','isCancerDriverGene','tissueSpecificity']
t = pf.read(columns=cols)
d = t.to_pydict()

constraint_layer, animalko_layer, driver_layer, tissue_layer = {}, {}, {}, {}
constraint_vals, animalko_vals, tissue_vals = [], [], []
for ensg, gc, mk, drv, ts in zip(d['targetId'], d['geneticConstraint'], d['mouseKOScore'],
                                  d['isCancerDriverGene'], d['tissueSpecificity']):
    sym = ensg_map.get(ensg)
    if not sym: continue
    # genetic constraint (gnomAD o/e z, already signed; higher = more constrained = more intolerant)
    if gc is not None:
        constraint_layer[sym] = {"raw": float(gc), "norm": None, "present": True}
        constraint_vals.append(float(gc))
    # mouse KO score (signed; more negative = more lethal/essential)
    if mk is not None:
        animalko_layer[sym] = {"raw": float(mk), "norm": None, "present": True}
        animalko_vals.append(float(mk))
    # cancer driver flag (-1 sentinel = driver)
    if drv is not None:
        driver_layer[sym] = {"raw": 1.0, "norm": 1.0, "present": True}
    # tissue specificity (signed; higher = more specific)
    if ts is not None:
        tissue_layer[sym] = {"raw": float(ts), "norm": None, "present": True}
        tissue_vals.append(float(ts))
# normalize constraint / animalko / tissue (signed -> 0..1 by minmax of observed)
for L, V in [(constraint_layer, constraint_vals), (animalko_layer, animalko_vals), (tissue_layer, tissue_vals)]:
    lo, hi, _ = minmax_norm(V)
    for g, r in L.items():
        r["norm"] = norm01(r["raw"], lo, hi)
print(f"      constraint={len(constraint_layer)} animalko={len(animalko_layer)} driver={len(driver_layer)} tissue={len(tissue_layer)}")

# ---------------------------------------------------------------- 8-9. HPA proteinatlas
print("[5/11] HPA protein atlas (tissue specificity + PDAC prognostic) ...")
hpa_path = os.path.join(ROOT, "downloads_v11", "proteinatlas.tsv")
hpa_tiss = {}   # RNA tissue specificity score (only present for specific genes)
hpa_prog = {}   # PDAC prognostic: {raw: 'unfavorable (p)' / 'favorable (p)' / 'unprognostic (p)', norm, present}
prog_dir = {"unfavorable": 1.0, "favorable": 1.0, "unprognostic": 0.0}
with open(hpa_path, newline='') as f:
    r = csv.reader(f, delimiter='\t')
    header = next(r)
    # index columns by substring (header cells are NOT quoted in the raw tsv)
    def find_col(sub, required=True):
        for i, h in enumerate(header):
            if sub in h:
                return i
        if required:
            raise RuntimeError(f"HPA column not found: {sub}")
        return None
    i_sym = 0
    i_rna_spec = find_col('RNA tissue specificity score')
    i_rna_spec_ntpm = find_col('RNA tissue specific nTPM')
    i_panc_tcga = find_col('Pancreatic Adenocarcinoma (TCGA)', required=False)
    i_panc_val = find_col('Pancreatic Adenocarcinoma (validation)', required=False)
    for row in r:
        max_i = max(i_rna_spec or 0, i_panc_tcga or 0, i_panc_val or 0)
        if len(row) <= max_i: continue
        sym = row[i_sym]
        # RNA tissue specificity score (empty => not tissue-specific)
        sv = row[i_rna_spec].strip().strip('"')
        if sv and sv.lower() not in ("na","",):
            v = to_float(sv)
            if v is not None:
                hpa_tiss[sym] = {"raw": v, "norm": None, "present": True}
        # PDAC prognostic (TCGA primary; fall back to validation)
        prog_txt = row[i_panc_tcga].strip().strip('"') if i_panc_tcga is not None else ""
        if (not prog_txt) and i_panc_val is not None:
            prog_txt = row[i_panc_val].strip().strip('"')
        if prog_txt:
            low = prog_txt.lower()
            # pattern: "<un/favorable/prognostic> (<p>)"
            if "potential prognostic unfavorable" in low or "unfavorable" in low:
                direction = "unfavorable"
            elif "favorable" in low:
                direction = "favorable"
            else:
                direction = "unprognostic"
            # extract p-value
            pv = None
            import re
            m = re.search(r"\(([0-9.eE+-]+)\)", prog_txt)
            if m: pv = to_float(m.group(1))
            significant = (pv is not None and pv < 0.05)
            # norm: prognostic-unfavorable in PDAC = strong cancer-relevance signal (1.0)
            #       prognostic-favorable = relevance too (1.0); unprognostic = 0.0
            norm = 1.0 if direction in ("unfavorable","favorable") else 0.0
            hpa_prog[sym] = {"raw": prog_txt, "norm": norm, "present": True,
                             "direction": direction, "pvalue": pv, "significant": significant}
lo, hi, n = minmax_norm([r["raw"] for r in hpa_tiss.values()]) if hpa_tiss else (0,1,0)
for g, r in hpa_tiss.items():
    r["norm"] = norm01(r["raw"], lo, hi)
print(f"      hpa_tissue_specificity={len(hpa_tiss)} hpa_pdac_prognostic={len(hpa_prog)}")

# ---------------------------------------------------------------- 10. OT genetics PDAC 500
print("[6/11] Open Targets genetics PDAC 500 targets ...")
pc = load_json(os.path.join(SKILL_DB, "disease_gwas_major", "pancreatic_cancer.json"))
pc_layer = {}
for item in pc.get("candidate_targets", []):
    g = item.get("gene")
    s = to_float(item.get("score"))
    if not g or s is None: continue
    pc_layer[g] = {"raw": s, "norm": None, "present": True}
lo, hi, n = minmax_norm([r["raw"] for r in pc_layer.values()])
for g, r in pc_layer.items():
    r["norm"] = norm01(r["raw"], lo, hi)

# ---------------------------------------------------------------- 11. Druggability
print("[7/11] Druggability (approved + small-molecule binders) ...")
dt = load_json(os.path.join(SKILL_DB, "drug_targets.json"))
drug_approved = set()
for item in dt.get("drugs", []):
    tg = item.get("target")
    if tg: drug_approved.add(tg.upper())
# also OT parquet small-molecule binder flags (already loaded? reload needed)
drug_layer = {}
for g in drug_approved:
    drug_layer[g] = {"raw": 1.0, "norm": 1.0, "present": True, "approved": True}
# add OT small-molecule binder as continuous druggability (if not already present)
t2 = pf.read(columns=['targetId','hasSmallMoleculeBinder','hasPocket','isInMembrane','maxClinicalStage'])
d2 = t2.to_pydict()
for ensg, smb, pocket, mem, stage in zip(d2['targetId'], d2['hasSmallMoleculeBinder'],
                                          d2['hasPocket'], d2['isInMembrane'], d2['maxClinicalStage']):
    sym = ensg_map.get(ensg)
    if not sym: continue
    # druggability score: small-molecule binder=1, pocket=0.6, membrane=0.3, clinical stage present adds
    score = 0.0
    if smb: score = max(score, 1.0)
    elif pocket: score = max(score, 0.6)
    if mem: score = max(score, 0.3)
    if stage is not None:
        score = max(score, float(stage))  # 0.01..1.0 clinical advancement
    if score > 0.0 and sym not in drug_layer:
        drug_layer[sym] = {"raw": score, "norm": score, "present": True, "approved": False}
print(f"      druggable genes={len(drug_layer)} (approved-known={len(drug_approved)})")

# ---------------------------------------------------------------- assemble
all_genes = set()
for L in [string_layer, mut_layer, impc_layer, constraint_layer, animalko_layer,
          driver_layer, tissue_layer, hpa_tiss, hpa_prog, pc_layer, drug_layer]:
    all_genes.update(L.keys())
all_genes = sorted(all_genes)
print(f"[*] TOTAL unique genes across all layers: {len(all_genes)}")

evidence = {
    "meta": {
        "description": "Normalized multimodal evidence layers for PDAC ECS V11. Each layer: gene->{raw,norm,present}. norm in [0,1] unless signed-raw preserved.",
        "n_genes_total": len(all_genes),
        "layers": {
            "string_centrality": {"type":"network","n":len(string_layer),"direction":"higher=better connected"},
            "mutation_freq": {"type":"genetic","n":len(mut_layer),"direction":"higher=more recurrently mutated"},
            "impc_animal_ko": {"type":"animal","n":len(impc_layer),"direction":"higher=more lethal/essential in mouse"},
            "genetic_constraint": {"type":"constraint","n":len(constraint_layer),"direction":"higher=more intolerant of LoF"},
            "mouse_ko_score": {"type":"animal","n":len(animalko_layer),"direction":"higher=more viable (less lethal)"},
            "cancer_driver": {"type":"driver","n":len(driver_layer),"direction":"1=cancer driver gene"},
            "tissue_specificity": {"type":"tissue","n":len(tissue_layer),"direction":"higher=more tissue-specific"},
            "hpa_rna_tissue_spec": {"type":"tissue","n":len(hpa_tiss),"direction":"higher=more RNA tissue-specific"},
            "hpa_pdac_prognostic": {"type":"cancer","n":len(hpa_prog),"direction":"1=prognostic in PDAC"},
            "ot_genetics_pdac": {"type":"genetic","n":len(pc_layer),"direction":"higher=genetic association score"},
            "druggability": {"type":"drug","n":len(drug_layer),"direction":"1=approved/druggable target"},
        }
    },
    "genes": all_genes,
    "layers": {
        "string_centrality": string_layer,
        "mutation_freq": mut_layer,
        "impc_animal_ko": impc_layer,
        "genetic_constraint": constraint_layer,
        "mouse_ko_score": animalko_layer,
        "cancer_driver": driver_layer,
        "tissue_specificity": tissue_layer,
        "hpa_rna_tissue_spec": hpa_tiss,
        "hpa_pdac_prognostic": hpa_prog,
        "ot_genetics_pdac": pc_layer,
        "druggability": drug_layer,
    }
}

out_path = os.path.join(ROOT, "evidence_layers_v11.json")
with open(out_path, "w") as f:
    json.dump(evidence, f, indent=1)
print(f"[DONE] wrote {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")

# quick sanity: coverage of key PDAC genes
key = ["KRAS","TP53","CDKN2A","SMAD4","STK11","PALB2","BRCA2","ATM","TGFBR2","GATA6"]
print("[sanity] key-gene layer coverage:")
for g in key:
    cov = [L for L in evidence["layers"] if g in evidence["layers"][L] and evidence["layers"][L][g]["present"]]
    print(f"   {g}: {len(cov)}/{len(evidence['layers'])} layers -> {cov}")
