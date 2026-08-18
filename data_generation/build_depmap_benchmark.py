#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_depmap_benchmark.py  -  v9 held-out validation & benchmark suite.

Addresses:
  Major 3/4/11/13 : DepMap/CRISPR external validation (Spearman, AUROC, AUPRC,
                   P@K, NDCG@K, enrichment@K) comparing ECS vs baselines.
  Major 6         : data-driven lambda/tau selection on a DEV split, locked &
                   evaluated on TEST (honest generalization, not post-hoc tuning).
  Major 12        : low-mutation conditional incremental value (Delta AUPRC
                   ECS - MutFreq across mutation bins).
  Major 20/21     : novel discovery endpoint (top genes outside driver/cancer set,
                   supported by DepMap).
  Major 10        : colorectal zero-shot transfer (frozen PDAC params) + failure analysis.

Scorers compared:
  ECS (proposed), Mutation frequency, STRING centrality, IMPC, SimpleMean,
  NetworkPropagation (k-NN diffusion over STRING).
  NOTE: OpenTargets & legacy 'v6' baselines require the Open Targets platform /
  prior-version scores not re-run in this analysis; this is stated honestly.
"""
import json, os, math, random, gzip
import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))
random.seed(20260816); np.random.seed(20260816)

def load(p):
    with open(os.path.join(OUT,p)) as f: return json.load(f)

# ---------- metric helpers ----------
def rankdata(a):
    a=np.asarray(a,float); n=len(a); order=np.argsort(a,kind="mergesort")
    ranks=np.empty(n); r=1
    while r<=n:
        j=r
        while j<n and a[order[j-1]]==a[order[j]]: j+=1
        avg=(r+j-1)/2.0
        for k in range(r-1,j): ranks[order[k]]=avg
        r=j+1
    return ranks

def spearman(a,b):
    a=np.asarray(a,float); b=np.asarray(b,float); n=len(a)
    if n<3: return float("nan")
    ra,rb=rankdata(a),rankdata(b); d=ra-rb
    return 1.0-6.0*float(np.sum(d*d))/(n*(n*n-1))

def auroc(y,s):
    y=np.asarray(y,int); s=np.asarray(s,float); n=len(y); r=rankdata(s)
    P=int(y.sum()); 
    if P==0 or P==n: return float("nan")
    sp=float(r[y==1].sum())
    return (sp-P*(P+1)/2.0)/(P*(n-P))

def auprc(y,s):
    y=np.asarray(y,int); s=np.asarray(s,float); n=len(y); P=int(y.sum())
    if P==0: return float("nan")
    order=np.argsort(-s,kind="mergesort"); tp=0; prev=0.0; ap=0.0
    for i,idx in enumerate(order):
        if y[idx]==1:
            tp+=1; prec=tp/(i+1.0); rec=tp/P
            ap+=prec*(rec-prev); prev=rec
    return ap

def p_at_k(y,s,k):
    order=np.argsort(-np.asarray(s,float),kind="mergesort")[:k]
    return float(np.asarray(y,int)[order].sum())/k

def ndcg_at_k(y,s,k):
    y=np.asarray(y,int); s=np.asarray(s,float); n=len(y)
    order=np.argsort(-s,kind="mergesort")[:k]
    dcg=sum(y[order[i]]/math.log2(2+i) for i in range(len(order)))
    ideal=sorted(y.tolist(),reverse=True)[:k]
    idcg=sum(ideal[i]/math.log2(2+i) for i in range(len(ideal)))
    return float(dcg/idcg) if idcg>0 else 0.0

def enrich_at_k(y,s,k):
    base=float(np.asarray(y,int).sum())/len(y)
    return p_at_k(y,s,k)/base if base>0 else float("nan")

def full_metrics(y,s,ks=(10,50,100)):
    return {"spearman":spearman(s,y),"auroc":auroc(y,s),"auprc":auprc(y,s),
            **{f"P@{k}":p_at_k(y,s,k) for k in ks},
            **{f"NDCG@{k}":ndcg_at_k(y,s,k) for k in ks},
            **{f"enrich@{k}":enrich_at_k(y,s,k) for k in ks},
            "base_rate":float(np.asarray(y,int).sum())/len(y)}

# ---------- load layers & ground truth ----------
gt = load("depmap_pdac_dependency.json")           # gene -> {gene_effect_mean, dependency_score, essential}
mut = load("paad_mutfreq_genomewide.json")
net = load("string_centrality.json")
imp = load("impc_layer.json")
ecs_scores = load("paad_ecs_scores.json")
ecs_comp   = load("paad_ecs_components.json")
meta = load("paad_ecs_v9.json")["meta"]
LAMBDA0, C_THR0 = meta["lambda"], meta["consistency_thr"]
KNOWN = set(meta["known_positives"])

# common universe = genes with DepMap PDAC label AND in STRING universe
string_syms = set(net.keys())
U = [g for g in gt if g in string_syms]
print("benchmark universe (DepMap PDAC & STRING):", len(U), flush=True)

# normalized mutfreq (identical to build_ecs_v9)
mvals=[v["freq_pct"] for g,v in mut.items() if isinstance(v,dict) and v.get("freq_pct") is not None]
mmax=max(mvals) if mvals else 1.0
def mut_n(g):
    v=mut.get(g)
    return (v["freq_pct"]/mmax) if isinstance(v,dict) and v.get("freq_pct") is not None else 0.0
def net_n(g): return net[g]["centrality"] if g in net else 0.0
def imp_n(g): return imp[g]["score"] if g in imp else 0.0

# build scorer vectors over U (missing -> 0, consistent with ECS missingness)
y_labels = np.array([1 if gt[g]["essential"] else 0 for g in U], float)
y_dep    = np.array([gt[g]["dependency_score"] for g in U], float)
S_mut  = np.array([mut_n(g) for g in U], float)
S_net  = np.array([net_n(g) for g in U], float)
S_imp  = np.array([imp_n(g) for g in U], float)
S_ecs  = np.array([ecs_scores.get(g,0.0) for g in U], float)
S_mean = np.array([np.mean([x for x in (mut_n(g),net_n(g),imp_n(g)) if x>0]) if any(x>0 for x in (mut_n(g),net_n(g),imp_n(g))) else 0.0 for g in U], float)

# ---------- Network propagation baseline (k-NN diffusion over STRING) ----------
def build_prop(seed, K=15, iters=12, alpha=0.6):
    # adjacency restricted to U, keep top-K neighbours per node by combined score
    pid2sym={}
    with gzip.open(os.path.join(OUT,"string_info.txt.gz"),"rt") as fh:
        fh.readline()
        for line in fh:
            p=line.rstrip("\n").split("\t")
            if len(p)>=2 and p[1]: pid2sym[p[0]]=p[1]
    adj={g:{} for g in U}
    with gzip.open(os.path.join(OUT,"string_links.txt.gz"),"rt") as fh:
        for line in fh:
            if line.startswith("#"): continue
            p=line.rstrip("\n").split()
            if len(p)<3: continue
            sa,sb=pid2sym.get(p[0]),pid2sym.get(p[1])
            if sa in adj and sb in adj and sa!=sb:
                try: s=float(p[2])
                except: continue
                if s>0:
                    if sa not in adj[sb] or s>adj[sa][sb]: adj[sa][sb]=s; adj[sb][sa]=s
    # keep top-K per node
    for g in adj:
        if len(adj[g])>K:
            top=sorted(adj[g].items(),key=lambda kv:-kv[1])[:K]
            adj[g]={k:v for k,v in top}
    # row-normalize
    for g in adj:
        tot=sum(adj[g].values())
        if tot>0: adj[g]={k:v/tot for k,v in adj[g].items()}
    f0=np.array([seed.get(g,0.0) for g in U],float)
    idx={g:i for i,g in enumerate(U)}
    f=f0.copy()
    for _ in range(iters):
        nf=np.zeros_like(f)
        for g in U:
            i=idx[g]
            for nb,w in adj[g].items():
                nf[i]+=w*f[idx[nb]]
        f=alpha*nf+(1-alpha)*f0
    return f

print("building network propagation baseline (k-NN diffusion)...", flush=True)
S_prop = build_prop({g:mut_n(g) for g in U})
# also a 'raw mutfreq-seeded' baseline is just S_mut

scorers = {
    "ECS_proposed": S_ecs,
    "MutationFreq": S_mut,
    "STRING_centrality": S_net,
    "IMPC": S_imp,
    "SimpleMean": S_mean,
    "NetworkPropagation": S_prop,
}

# ---------- primary DepMap benchmark table ----------
bench = {}
for name,s in scorers.items():
    m = full_metrics(y_labels, s)
    bench[name] = {"vs_essentiality": m,
                   "vs_continuous_dependency_spearman": spearman(s, y_dep)}
print("\n=== DepMap PDAC held-out validation (essentiality AUROC/AUPRC) ===", flush=True)
for name in scorers:
    b=bench[name]["vs_essentiality"]
    print(f"  {name:22s} AUROC={b['auroc']:.3f} AUPRC={b['auprc']:.3f} "
          f"P@10={b['P@10']:.3f} NDCG@10={b['NDCG@10']:.3f} "
          f"rho(dep)={bench[name]['vs_continuous_dependency_spearman']:.3f}", flush=True)

# ---------- Major 6: data-driven lambda/tau on DEV, locked on TEST ----------
dev_idx=[i for i,g in enumerate(U) if g not in KNOWN]
random.shuffle(dev_idx)
cut=len(dev_idx)//2
dev,test=dev_idx[:cut],dev_idx[cut:]
dev_set=set(dev); test_set=set(test)
print(f"\nMajor6 dev/test split: dev={len(dev)} test={len(test)} (genes outside KNOWN)", flush=True)

# ECS recompute for arbitrary (lam,tau) using stored components
def ecs_lamtau(comp, lam, tau):
    S=comp["S"]; sd=comp["sd"]; Uu=comp["U"]
    C=math.exp(-sd/tau)
    return S*C*(1.0-lam*Uu)

best=None
for lam in (0.3,0.6,0.9):
    for tau in (0.2,0.5,1.0):
        ds=np.array([ecs_lamtau(ecs_comp[U[i]],lam,tau) for i in dev],float)
        dy=np.array([y_dep[i] for i in dev],float)
        r=spearman(ds,dy)
        if best is None or (not math.isnan(r) and r>best[0]):
            best=(r,lam,tau)
print(f"  data-driven best (DEV Spearman): lam={best[1]} tau={best[2]} rho={best[0]:.3f}  (heuristic lam={LAMBDA0} tau={C_THR0})", flush=True)
# evaluate on TEST
def ecs_vec(lam,tau):
    return np.array([ecs_lamtau(ecs_comp[U[i]],lam,tau) for i in range(len(U))],float)
S_ecs_dd = ecs_vec(best[1],best[2])
bench["ECS_dataDriven"]={"params":{"lambda":best[1],"tau":best[2]},
    "TEST_vs_essentiality":full_metrics(np.array([y_labels[i] for i in test],float), S_ecs_dd[test]),
    "TEST_spearman_dependency":spearman(S_ecs_dd[test], np.array([y_dep[i] for i in test],float)),
    "heuristic_TEST_spearman_dependency":spearman(S_ecs[test], np.array([y_dep[i] for i in test],float))}

# ---------- Major 12: low-mutation conditional incremental value ----------
# bin genes by mutation frequency tier
mut_tier=np.array([mut_n(g) for g in U],float)
order=np.argsort(mut_tier); n=len(order); hi=order[int(0.75*n):]; lo=order[:int(0.75*n)]
def del_auprc(sub):
    e=np.array([S_ecs[i] for i in sub],float); m=np.array([S_mut[i] for i in sub],float)
    yb=np.array([y_labels[i] for i in sub],int)
    return auprc(yb,e)-auprc(yb,m), auprc(yb,e), auprc(yb,m)
lo_d,lo_e,lo_m=del_auprc(lo); hi_d,hi_e,hi_m=del_auprc(hi)
lowmut={"n_genes":len(lo),"AUPRC_ECS":lo_e,"AUPRC_MutFreq":lo_m,"Delta_AUPRC":lo_d}
highmut={"n_genes":len(hi),"AUPRC_ECS":hi_e,"AUPRC_MutFreq":hi_m,"Delta_AUPRC":hi_d}
print(f"\nMajor12 low-mutation bin (n={len(lo)}): DeltaAUPRC(ECS-MutFreq)={lo_d:+.3f}", flush=True)
print(f"Major12 high-mutation bin (n={len(hi)}): DeltaAUPRC(ECS-MutFreq)={hi_d:+.3f}", flush=True)

# ---------- Major 20/21: novel discovery endpoint ----------
cancer_set=set(KNOWN)|set(g for g,v in mut.items() if isinstance(v,dict))
ranked_ecs=sorted(range(len(U)),key=lambda i:-S_ecs[i])
novel_top=[U[i] for i in ranked_ecs if U[i] not in cancer_set][:20]
novel_endpoint=[{"gene":g,"ECS":round(float(S_ecs[U.index(g)]),4),
                 "DepMap_essential":gt[g]["essential"],
                 "DepMap_dependency":round(gt[g]["dependency_score"],4)} for g in novel_top]
print(f"\nMajor20/21 novel discovery endpoint (top20 outside driver/cancer set):", novel_top[:10],"...", flush=True)

# ---------- Major 10: colorectal zero-shot transfer ----------
crc_section={}
try:
    crc_ecs=load("crc_ecs_v9.json")
    try:
        crc_scores=load("crc_ecs_scores.json")
    except FileNotFoundError:
        crc_scores={r["gene"]:r["ECS"] for r in crc_ecs["ranking"]}
    crc_mut=load("crc_mutfreq_genomewide.json")
    cmvals=[v["freq_pct"] for g,v in crc_mut.items() if isinstance(v,dict) and v.get("freq_pct") is not None]
    cmmax=max(cmvals) if cmvals else 1.0
    def crc_mut_n(g):
        v=crc_mut.get(g)
        return (v["freq_pct"]/cmmax) if isinstance(v,dict) and v.get("freq_pct") is not None else 0.0
    gt_crc=load("depmap_crc_dependency.json")
    Uc=[g for g in gt_crc if g in string_syms]
    yc=np.array([1 if gt_crc[g]["essential"] else 0 for g in Uc],float)
    sc_ecs=np.array([crc_scores.get(g,0.0) for g in Uc],float)
    sc_mut=np.array([crc_mut_n(g) for g in Uc],float)
    sc_net=np.array([net_n(g) for g in Uc],float)
    sc_imp=np.array([imp_n(g) for g in Uc],float)
    sc_mean=np.array([np.mean([x for x in (crc_mut_n(g),net_n(g),imp_n(g)) if x>0]) if any(x>0 for x in (crc_mut_n(g),net_n(g),imp_n(g))) else 0.0 for g in Uc],float)
    crc_scorers={"ECS_proposed":sc_ecs,"MutationFreq":sc_mut,"STRING_centrality":sc_net,
                 "IMPC":sc_imp,"SimpleMean":sc_mean}
    crc_bench={n:full_metrics(yc,s) for n,s in crc_scorers.items()}
    crc_section={"n_genes":len(Uc),"frozen_params_from":"paad_ecs_v9.json",
                 "note":"D1 swapped to TCGA-COAD+READ; D2/D3 reused; NO re-tuning (true zero-shot)",
                 "crc_benchmark":crc_bench,
                 "failure_analysis":{"ECS_vs_MutFreq_AUROC_gap":round(float(crc_bench['ECS_proposed']['auroc']-crc_bench['MutationFreq']['auroc']),4)}}
    print(f"\nMajor10 CRC zero-shot transfer: ECS AUROC={crc_bench['ECS_proposed']['auroc']:.3f} "
          f"AUPRC={crc_bench['ECS_proposed']['auprc']:.3f}  (MutFreq AUROC={crc_bench['MutationFreq']['auroc']:.3f})", flush=True)
except FileNotFoundError as e:
    crc_section={"status":"skipped","reason":"crc_ecs_v9.json / crc_mutfreq not ready yet"}

out={
  "meta":{"ground_truth":"DepMap/CRISPR gene effect (PDAC cell lines)",
          "essential_rule":"top quartile -mean(gene_effect)",
          "benchmark_universe_size":len(U),
          "scorers":list(scorers.keys()),
          "known_positives":sorted(KNOWN),
          "limitations":"OpenTargets & legacy v6 baselines not re-run (require Open Targets platform / prior-version scores)."},
  "depmap_pdac_benchmark":bench,
  "low_mutation_conditional": {"low_mut_bin":lowmut,"high_mut_bin":highmut,
       "interpretation":"If Delta_AUPRC is larger in low-mutation bin, ECS adds value where mutation frequency fails."},
  "novel_discovery_endpoint":novel_endpoint,
  "colorectal_transfer":crc_section,
}
with open(os.path.join(OUT,"benchmark_v9.json"),"w") as f:
    json.dump(out,f,indent=2,ensure_ascii=False)
print("\nWROTE benchmark_v9.json", flush=True)
