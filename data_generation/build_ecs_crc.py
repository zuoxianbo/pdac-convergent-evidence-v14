#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_ecs_crc.py  -  Colorectal (CRC) ZERO-SHOT transfer of the PDAC-trained ECS (Major 10).

Frozen, NO re-tuning:
  * weights W        loaded from paad_ecs_v9.json (the PDAC-trained dependence weights)
  * lambda, C_THR    loaded from paad_ecs_v9.json
  * evidence definition, missingness policy, ranking rule  identical to PAAD

Only D1 (mutation frequency) is swapped to TCGA-COAD+READ; D2 (STRING) and D3 (IMPC)
are cancer-agnostic and therefore reused unchanged.

Input : crc_mutfreq_genomewide.json, string_centrality.json, impc_layer.json,
        paad_ecs_v9.json (for frozen params)
Output: crc_ecs_v9.json  (same schema as paad_ecs_v9.json)
"""
import json, os, math, random, itertools
import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))
random.seed(20260816); np.random.seed(20260816)

def load(p):
    with open(os.path.join(OUT,p)) as f: return json.load(f)

src = load("paad_ecs_v9.json")
W   = src["effective_weights_W"]
LAMBDA = src["meta"]["lambda"]
C_THR  = src["meta"]["consistency_thr"]
DEP_THR= src["meta"]["dep_thr"]
DISCOVERY = src["meta"]["discovery_layers"]
KNOWN = src["meta"]["known_positives"]
N_BOOT = src["meta"]["n_boot"]

mut = load("crc_mutfreq_genomewide.json")
net = load("string_centrality.json")
imp = load("impc_layer.json")

def norm_mut():
    vals=[v["freq_pct"] for g,v in mut.items() if isinstance(v,dict) and v.get("freq_pct") is not None]
    mx=max(vals) if vals else 1.0
    return {g:(v["freq_pct"]/mx) for g,v in mut.items() if isinstance(v,dict) and v.get("freq_pct") is not None}
MUT=norm_mut(); NET={g:v["centrality"] for g,v in net.items()}; IMP={g:v["score"] for g,v in imp.items()}
# identical per-layer [0,1] normalization as PDAC engine (STRING/IMPC sets are
# cancer-agnostic so normalization is identical -> genuinely frozen evidence def.)
def to01(d):
    vals=[v for v in d.values() if v is not None]
    if not vals: return d
    lo,hi=min(vals),max(vals); rng=(hi-lo) if hi>lo else 1.0
    return {g:((v-lo)/rng if v is not None else None) for g,v in d.items()}
MUT=to01(MUT); NET=to01(NET); IMP=to01(IMP)
universe=set(MUT)|set(NET)|set(IMP)
X={g:{"D1_mut":MUT.get(g),"D2_net":NET.get(g),"D3_impc":IMP.get(g)} for g in universe}

def boot_cv(g):
    ms=np.array([X[g][a] for a in DISCOVERY if X[g][a] is not None])
    if len(ms)<2: return 0.0
    means=[np.random.choice(ms,len(ms),replace=True).mean() for _ in range(N_BOOT)]
    means=np.array(means); return min(1.0,float(means.std()/(means.mean()+1e-9)))

def ecs_of(g):
    ms=[(a,X[g][a]) for a in DISCOVERY if X[g][a] is not None]
    if not ms: return None
    wsum=sum(W[a] for a,_ in ms); S=sum(W[a]*v for a,v in ms)/wsum
    vals=[v for _,v in ms]; sd=float(np.std(vals)) if len(vals)>1 else 0.0
    C=math.exp(-sd/C_THR); B=len(ms)/len(DISCOVERY); cv=boot_cv(g)
    U=min(1.0,max(0.0,0.5*(1.0-B)+0.5*cv))
    return {"ECS":S*C*(1.0-LAMBDA*U),"S":S,"C":C,"B":B,"U":U,"cv":cv,
            "n_measured":len(ms),"measured":[a for a,_ in ms]}

scored={g:ecs_of(g) for g in universe}; scored={g:v for g,v in scored.items() if v is not None}
ranked=sorted(scored,key=lambda g:scored[g]["ECS"],reverse=True)
rank={g:i+1 for i,g in enumerate(ranked)}

CANCER_GENES=set(KNOWN)|set(mut.keys())
novel=[g for g in ranked if g not in CANCER_GENES][:50]
out={"meta":{
        "method":"Evidence Convergence Score (ECS) v9 — colorectal ZERO-SHOT transfer (frozen PDAC params)",
        "transfer":"D1 swapped to TCGA-COAD+READ; D2/D3 reused; W,lambda,C_THR FROZEN from PAAD (no re-tuning)",
        "discovery_layers":DISCOVERY,"lambda":LAMBDA,"consistency_thr":C_THR,
        "universe_size":len(universe),"scored_genes":len(scored),
        "known_positives":KNOWN,"frozen_from":"paad_ecs_v9.json"},
     "effective_weights_W":W,
     "ranking":[{"gene":g,"rank":rank[g],"ECS":round(scored[g]["ECS"],5),
                "S":round(scored[g]["S"],4),"C":round(scored[g]["C"],4),
                "B":round(scored[g]["B"],3),"U":round(scored[g]["U"],4),
                "n_measured":scored[g]["n_measured"]} for g in ranked[:200]],
     "novel_candidates_top50":novel,"top10":ranked[:10]}
with open(os.path.join(OUT,"crc_ecs_v9.json"),"w") as f: json.dump(out,f,indent=2,ensure_ascii=False)
with open(os.path.join(OUT,"crc_ecs_scores.json"),"w") as f:
    json.dump({g:round(scored[g]["ECS"],6) for g in scored}, f)
print("WROTE crc_ecs_v9.json + crc_ecs_scores.json  top10:",ranked[:10],flush=True)
