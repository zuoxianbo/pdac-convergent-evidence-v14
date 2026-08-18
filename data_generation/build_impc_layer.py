#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build genome-wide IMPC in-vivo functional layer (D3) mapped to human symbols.
Reads IMPC genotype-phenotype (mouse genes) + MGI HOM file (mouse->human orthologs).
Output: impc_layer.json {human_symbol: {score, n_sig, lethal, source}}
Score (transparent, in [0,1]):
  score = clip(0.5*min(n_sig/5,1) + 0.5*lethal_flag, 0, 1)
  where n_sig = # phenotype associations with p < 1e-4; lethal_flag = any lethality MP term.
"""
import os, json, gzip, re
OUT = os.path.dirname(os.path.abspath(__file__))
IMP = "/Users/zuoxianbo/Desktop/zuoxb-虚拟科研系统-MAC/Zuoxb-Data-Medicine-platform/models/animal_models_db/impc_genotype_phenotype.json"
HOM = os.path.join(OUT, "HOM_MouseHumanSequence.rpt")

# lethality-related MP terms
LETHAL_MP = {
    "MP:0008765",  # embryonic lethal
    "MP:0011110",  # preweaning lethality, incomplete penetrance
    "MP:0006021",  # postnatal lethality
    "MP:0009265",  # perinatal lethality
    "MP:0005172",  # embryonic growth retardation (proxy)
    "MP:0002082",  # decreased viability
}

def mouse_to_human():
    """Group HOM file by DB Class Key; collect mouse & human symbols."""
    m2h = {}
    grp = {}
    with open(HOM) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c:i for i,c in enumerate(header)}
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(idx.values()): continue
            key = p[idx.get("DB Class Key",0)]
            taxon = p[idx.get("NCBI Taxon ID",2)]
            sym = p[idx.get("Symbol",3)]
            if not key or not sym: continue
            grp.setdefault(key, {})[taxon] = sym
    for key, d in grp.items():
        ms = d.get("10090"); hs = d.get("9606")
        if ms and hs:
            m2h[ms] = hs
    return m2h

m2h = mouse_to_human()
print("mouse->human orthologs:", len(m2h), flush=True)

d = json.load(open(IMP))
genes = d["genes"]  # dict mouse_symbol -> {gene, phenotypes:[{mp_term,mp_id,top_level,p_value,sex}]}
print("IMPC genes (mouse):", len(genes), flush=True)

impc = {}
for mg, rec in genes.items():
    pheno = rec.get("phenotypes", [])
    if not pheno: continue
    n_sig = sum(1 for x in pheno if isinstance(x.get("p_value"),(int,float)) and x["p_value"] < 1e-4)
    lethal = any(x.get("mp_id") in LETHAL_MP for x in pheno)
    if n_sig == 0 and not lethal: continue
    score = min(1.0, 0.5*min(n_sig/5.0,1.0) + 0.5*(1.0 if lethal else 0.0))
    hs = m2h.get(mg)
    if not hs: continue
    impc[hs] = {"score": round(score,4), "n_sig": n_sig, "lethal": int(lethal), "mouse_gene": mg}

with open(os.path.join(OUT,"impc_layer.json"),"w") as fh:
    json.dump(impc, fh, indent=1)
print("WROTE impc_layer.json :", len(impc), "human genes with IMPC in-vivo score", flush=True)
top = sorted(impc.items(), key=lambda kv: kv[1]["score"], reverse=True)[:10]
for s,v in top:
    print(f"  {s}: score={v['score']} n_sig={v['n_sig']} lethal={v['lethal']}")
