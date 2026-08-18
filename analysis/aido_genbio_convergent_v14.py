#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AIDO + genbio convergent-evidence cross-validation (PDAC V14.2).

REAL data only:
  - ECS benchmark:   results/benchmark_v14.json  (layer 1: network+genetics+animal ECS)
  - genbio virtual-cell PDAC output (skill): data/genbio_pdac_virtual_cell_pipeline.json
        multi-model consensus from evo2 / state / alphagenome (GenBio-AI virtual-cell suite)
  - AIDO DNA foundation model (layer 3): results/aido_dna_pdac.json  (produced by
        analysis/run_aido_dna.py; appended if present, else reported as pending)

Outputs: results/genbio_crossval.json  + printed summary.
No fabricated numbers: every quantity is read from the two source JSONs.
"""
import json, os, math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

def load(p):
    with open(os.path.join(ROOT, p), encoding="utf-8") as f:
        return json.load(f)

bench = load("results/benchmark_v14.json")
genbio = load("data/genbio_pdac_virtual_cell_pipeline.json")

# ---------- Layer 1: ECS per-gene ranking (from conv_table) ----------
conv = bench["ecs_specific"].get("conv_table", [])
ecs_by_gene = {r["gene"]: r for r in conv}
ecs_top = sorted(conv, key=lambda r: -r.get("ecs", 0.0)) if conv else []
ecs_top_genes = [r["gene"] for r in ecs_top]
meta = bench["meta"]
e6_genes = set(meta.get("e6_genes", []))            # 35 independent clinical-validation genes
e3_auroc = bench["best_per_endpoint"]["E3_actionable_target"]["auroc"]
e3_scorer = bench["best_per_endpoint"]["E3_actionable_target"]["scorer"]

# ---------- Layer 2: genbio virtual-cell PDAC output (real skill run) ----------
gb_disc = genbio["target_discovery"]
gb_top = gb_disc["top_targets"]                    # [TP53, EGFR, KRAS, CDKN2A, WNT] w/ centrality
gb_top_genes = [t["gene"] for t in gb_top]
gb_cent = {t["gene"]: t.get("centrality") for t in gb_top}
rec = genbio["multi_model_validation"]["recommendation"]
gb_rec = rec["targets"]                            # [TP53, EGFR, KRAS]
gb_consensus = {g: genbio["multi_model_validation"]["target_scores"][g]["consensus_score"] for g in gb_rec}
gb_mv = genbio["multi_model_validation"]["target_scores"]
gb_models = rec.get("models_contributing", [])
gb_conf = rec.get("confidence")
gb_omics = genbio.get("multi_omics_validation", {})
gb_drug = genbio.get("druggability", {})

# ---------- Convergence metrics ----------
# (a) Jaccard: genbio recommendation top-3 vs ECS top-K (K=10)
K = 10
ecs_topK = set(ecs_top_genes[:K])
gb_top3 = set(gb_rec)
jaccard_top3_vs_ecsK = (len(gb_top3 & ecs_topK) / len(gb_top3 | ecs_topK)) if (gb_top3 | ecs_topK) else 0.0

# (b) genbio top-3 membership in independent clinical-validation set E6
gb_top3_in_e6 = sorted(gb_top3 & e6_genes)

# (c) For genes present in BOTH ECS conv_table and genbio recommendation: rank correlation
common = [g for g in gb_rec if g in ecs_by_gene]
def spearman(a, b):
    # a,b: lists of (label, value) aligned by gene
    genes = [g for g, _ in a]
    rank = {g: i + 1 for i, g in enumerate(sorted(genes, key=lambda x: -dict(a)[x]))}
    rank2 = {g: i + 1 for i, g in enumerate(sorted(genes, key=lambda x: -dict(b)[x]))}
    n = len(genes)
    if n < 3:
        return None
    d2 = sum((rank[g] - rank2[g]) ** 2 for g in genes)
    return 1 - 6 * d2 / (n * (n ** 2 - 1))
sp_ecs_genbio = spearman([(g, ecs_by_gene[g]["ecs"]) for g in common],
                         [(g, gb_consensus[g]) for g in common]) if common else None

# (d) AIDO layer (if produced)
aido = None
aido_path = os.path.join(ROOT, "results/aido_dna_pdac.json")
if os.path.exists(aido_path):
    aido = load("results/aido_dna_pdac.json")

# ---------- Assemble per-gene convergent table ----------
candidate_genes = list(dict.fromkeys(gb_top_genes + ecs_top_genes[:K]))
rows = []
for g in candidate_genes:
    ecs_r = ecs_by_gene.get(g)
    rows.append({
        "gene": g,
        "ecs_specific_score": round(ecs_r["ecs"], 3) if ecs_r else None,
        "ecs_vs_string_delta": round(ecs_r["delta"], 3) if ecs_r else None,
        "in_ecs_top10": g in ecs_topK,
        "in_genbio_top5": g in set(gb_top_genes),
        "in_genbio_recommendation": g in set(gb_rec),
        "in_E6_clinical": g in e6_genes,
        "genbio_centrality": round(gb_cent[g], 3) if g in gb_cent else None,
        "genbio_consensus": round(gb_consensus[g], 3) if g in gb_consensus else None,
        "genbio_models": gb_models if g in set(gb_rec) else None,
        "real_gwas": gb_omics.get(g, {}).get("GWAS") if g in gb_omics else None,
        "gwas_is_simulated": (g in gb_omics and "模拟" in str(gb_omics[g].get("GWAS", ""))),
        "omics_score": gb_omics.get(g, {}).get("omics_score") if g in gb_omics else None,
        "druggability_total": gb_drug.get(g, {}).get("total_score") if g in gb_drug else None,
        "aido_dna_seqscore": (round(aido["per_gene"][g]["seq_evidence"], 4)
                              if aido and g in aido.get("per_gene", {}) else None),
    })

result = {
    "meta": {
        "disease": "PDAC (pancreatic ductal adenocarcinoma)",
        "ecs_formula": meta["ecs_formula"],
        "e3_actionable_auroc": round(e3_auroc, 4),
        "e3_actionable_scorer": e3_scorer,
        "n_e6_clinical_genes": len(e6_genes),
        "genbio_models_used": gb_models,
        "genbio_confidence": gb_conf,
        "aido_layer_present": aido is not None,
    },
    "convergence": {
        "genbio_recommendation_top3": sorted(gb_top3),
        "ecs_top10": ecs_top_genes[:K],
        "jaccard_genbio_top3_vs_ecs_top10": round(jaccard_top3_vs_ecsK, 4),
        "genbio_top3_in_E6_clinical": gb_top3_in_e6,
        "n_genbio_top3_in_E6": len(gb_top3_in_e6),
        "genes_in_both_ecs_and_genbio": common,
        "spearman_ecs_vs_genbio_consensus": (round(sp_ecs_genbio, 4) if sp_ecs_genbio is not None else None),
    },
    "per_gene": rows,
    "genbio_target_discovery": gb_top,
    "genbio_multi_omics": gb_omics,
    "genbio_druggability": gb_drug,
}

out = os.path.join(ROOT, "results/genbio_crossval.json")
with open(out, "w", encoding="utf-8") as f:
    json.dump(result, f, indent=2, ensure_ascii=False)

# ---------- Printed summary ----------
print("=" * 70)
print("CONVERGENT EVIDENCE: ECS  x  genbio(virtual-cell)  x  AIDO(DNA)")
print("=" * 70)
print(f"ECS formula        : {meta['ecs_formula']}")
print(f"E3 actionable AUROC: {e3_auroc:.3f} ({e3_scorer})   | E6 clinical genes: {len(e6_genes)}")
print(f"genbio models used : {gb_models}  (confidence={gb_conf})")
print(f"AIDO DNA layer     : {'PRESENT' if aido else 'pending (run_aido_dna.py)'}")
print("-" * 70)
print(f"genbio recommendation (top-3): {sorted(gb_top3)}")
print(f"ECS top-10               : {ecs_top_genes[:K]}")
print(f"Jaccard(genbio_top3, ECS_top10) = {jaccard_top3_vs_ecsK:.3f}")
print(f"genbio top-3 in E6 clinical set : {gb_top3_in_e6}  ({len(gb_top3_in_e6)}/3)")
print(f"genes in both ECS & genbio     : {common}")
print(f"Spearman(ECS, genbio consensus): {sp_ecs_genbio}")
print("-" * 70)
print("Per-gene convergent table (key columns):")
print(f"{'gene':8} {'ECS':>6} {'dECS':>6} {'gCent':>6} {'gCon':>6} {'E6':>4} {'GWAS(real)':>10} {'drug':>5} {'AIDO':>6}")
for r in rows:
    if r["in_genbio_top5"] or r["in_ecs_top10"]:
        gwas = ("Y" if (r["real_gwas"] and not r["gwas_is_simulated"]) else
                ("sim" if r["gwas_is_simulated"] else "-"))
        print(f"{r['gene']:8} "
              f"{str(r['ecs_specific_score']):>6} {str(r['ecs_vs_string_delta']):>6} "
              f"{str(r['genbio_centrality']):>6} {str(r['genbio_consensus']):>6} "
              f"{('Y' if r['in_E6_clinical'] else '-'):>4} {gwas:>10} "
              f"{str(r['druggability_total']):>5} {str(r['aido_dna_seqscore']):>6}")
print("-" * 70)
print(f"WROTE {out}")
