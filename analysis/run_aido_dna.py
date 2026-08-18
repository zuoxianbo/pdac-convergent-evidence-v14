#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AIDO.DNA-300M zero-shot promoter-representation analysis for the PDAC convergent-evidence framework.

QUESTION (independent of ECS construction)
  Do the promoter sequence representations of a genome foundation model place ECS-specific
  candidates closer to the *independent* clinical-target reference (E6) than matched controls?

DESIGN (all gene groups derive from real project data; none is hand-picked to favour a result)
  reference   E6  : 35 independent PDAC clinical targets (ClinicalTrials.gov-derived, never used
                    in ECS construction)                       -> defines the reference centroid
  test        ECS : ECS-specific re-ranked candidates (benchmark_v14.json ecs_specific.conv_table),
                    with any E6 member removed                 -> disjoint from the reference
  control 1   HK  : canonical housekeeping genes               -> expected to be unrelated
  control 2   RND : random genes drawn from the scored universe (seed 42), excluding all above
                    -> empirical null

MODEL
  AIDO.DNA-300M (GenBio AI; RNABert architecture, 24 layers, hidden 1024, RoPE, 4,000-token context),
  loaded from local weights via ModelGenerator. Used strictly ZERO-SHOT: no PDAC fine-tuning,
  no label information, no gradient updates. Embedding = mean of final-layer token states over
  real (non-special, non-pad) nucleotides, L2-normalised.

READOUT
  cos(gene promoter embedding, E6 centroid). Group differences by two-sided Mann-Whitney U.
  GC content is reported and tested as a confounder (Spearman), because GC is the dominant
  compositional signal in promoters.

HONESTY
  Promoters are real Ensembl GRCh38 sequences. If the runtime cannot be provisioned the script
  still saves the real sequences and writes status='pending' with the exact reason.
  No embedding, similarity or P-value is ever synthesised.

Outputs: results/aido_dna_pdac.json , data/aido_pdac_promoters.json
"""
import json, os, random, sys, time, types, urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SKILL_SCRIPTS = "/Users/zuoxianbo/.workbuddy/skills/zuoxb-virtual-cell-platform/scripts"
ENSEMBL = "https://rest.ensembl.org"

PROMOTER_UP = 2000          # bp upstream of TSS (standard proximal promoter)
N_RANDOM = 30
SEED = 42
HOUSEKEEPING = ["GAPDH", "ACTB", "RPL13A", "TBP", "PPIA", "B2M"]

# stubbing these would corrupt the numerics -> never allowed
DENY_STUB = {"torch", "numpy", "scipy", "transformers", "lightning", "lightning_utilities",
             "safetensors", "huggingface_hub", "lmdb", "filelock"}

PROM_CACHE = os.path.join(ROOT, "data/aido_pdac_promoters.json")
OUT = os.path.join(ROOT, "results/aido_dna_pdac.json")


# ---------------------------------------------------------------- Ensembl
def _rest(path, tries=4):
    url = f"{ENSEMBL}{path}?content-type=application/json"
    req = urllib.request.Request(url, headers={"User-Agent": "aido-pdac"})
    last = None
    for _ in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return json.loads(r.read().decode("utf-8", "replace"))
        except Exception as e:
            last = e
            time.sleep(1.5)
    raise RuntimeError(f"{path}: {last}")


def fetch_promoter(symbol):
    loc = _rest(f"/lookup/symbol/homo_sapiens/{symbol}")
    chrom, strand = loc["seq_region_name"], loc["strand"]
    if strand == 1:
        p1, p2 = max(1, loc["start"] - PROMOTER_UP), loc["start"]
    else:
        p1, p2 = loc["end"], loc["end"] + PROMOTER_UP
    seq = _rest(f"/sequence/region/homo_sapiens/{chrom}:{p1}..{p2}:{strand}").get("seq", "").upper()
    if len(seq) < PROMOTER_UP * 0.9:
        raise RuntimeError(f"short sequence ({len(seq)}bp)")
    return {"symbol": symbol, "chr": str(chrom), "strand": strand,
            "region": f"{chrom}:{p1}..{p2}:{strand}", "seq": seq, "len": len(seq),
            "gc": round((seq.count("G") + seq.count("C")) / len(seq), 4)}


# ---------------------------------------------------------------- gene groups
def build_groups():
    B = json.load(open(os.path.join(ROOT, "results/benchmark_v14.json")))
    e6 = json.load(open(os.path.join(ROOT, "data/e6_clinical_validation.json")))
    ev = json.load(open(os.path.join(ROOT, "data/evidence_layers_v11.json")))

    e6_genes = list(dict.fromkeys(e6["genes"]))
    ecs = [r["gene"] for r in B["ecs_specific"]["conv_table"]]
    ecs = [g for g in dict.fromkeys(ecs) if g not in set(e6_genes)]

    universe = [g for g in ev["genes"] if isinstance(g, str)]
    used = set(e6_genes) | set(ecs) | set(HOUSEKEEPING)
    pool = sorted(g for g in universe if g not in used and g.isalnum())
    rnd = random.Random(SEED).sample(pool, min(N_RANDOM, len(pool)))

    return {"E6_clinical_reference": e6_genes, "ECS_specific": ecs,
            "housekeeping_control": list(HOUSEKEEPING), "random_background": rnd}, len(universe)


# ---------------------------------------------------------------- runtime
def import_backbones(budget=20):
    """Import modelgenerator.backbones, stubbing only missing *data-loading* side deps.

    ModelGenerator's package __init__ chain pulls cell/spatial data utilities
    (bionty, pyfaidx, tiledbsoma, wandb, ...) that are irrelevant to a DNA backbone
    forward pass. Numerical libraries are never stubbed (DENY_STUB).
    """
    stubbed = []

    class _Dummy:
        def __call__(self, *a, **k): return _Dummy()
        def __getattr__(self, n): return _Dummy()

    for _ in range(budget):
        try:
            import modelgenerator.backbones as bb
            return bb, stubbed
        except ModuleNotFoundError as e:
            name = getattr(e, "name", None)
            if not name or name.split(".")[0] in DENY_STUB or name.startswith("modelgenerator"):
                raise
            m = types.ModuleType(name)
            m.__path__ = []
            m.__getattr__ = lambda n, _d=_Dummy: _d()
            sys.modules[name] = m
            stubbed.append(name)
    raise RuntimeError(f"stub budget exhausted; stubbed={stubbed}")


def load_model():
    sys.path.insert(0, SKILL_SCRIPTS)
    import torch
    torch.set_num_threads(os.cpu_count() or 8)
    bb, stubbed = import_backbones()
    from modelgenerator.backbones.base import LegacyAdapterType
    from model_aido import variant_dir
    d = variant_dir("AIDO.DNA-300M")
    m = bb.gb_dna_300m(LegacyAdapterType.MASKED_LM, None)
    m.model_path = str(d)
    m.setup()
    m.eval()
    dev = "cpu"
    if torch.backends.mps.is_available():
        try:
            m.to("mps"); dev = "mps"
        except Exception:
            m.to("cpu"); dev = "cpu"
    return m, dev, stubbed, str(d)


def embed(m, dev, seq):
    import torch
    tok = m.tokenize([seq])
    ids = torch.tensor(tok["input_ids"])
    am = torch.tensor(tok["attention_mask"])
    stm = torch.tensor(tok["special_tokens_mask"])
    keep = am[0].bool() & (~stm[0].bool())          # real nucleotides only
    with torch.no_grad():
        out = m.forward(ids.to(dev), am.to(dev))
    h = out.last_hidden_state[0].float().cpu()
    v = h[keep].mean(0)
    return (v / v.norm()).numpy(), int(keep.sum())


# ---------------------------------------------------------------- main
def main():
    groups, n_universe = build_groups()
    wanted = [g for gl in groups.values() for g in gl]
    print(f"groups: " + ", ".join(f"{k}={len(v)}" for k, v in groups.items()), flush=True)

    proms = {}
    if os.path.exists(PROM_CACHE):
        try:
            cached = json.load(open(PROM_CACHE))
            proms = {g: p for g, p in cached.items() if isinstance(p, dict) and p.get("seq")}
            print(f"cache: reused {len(proms)} promoters", flush=True)
        except Exception:
            proms = {}

    failed = {}
    for g in wanted:
        if g in proms:
            continue
        try:
            proms[g] = fetch_promoter(g)
        except Exception as e:
            failed[g] = str(e)[:120]
    json.dump(proms, open(PROM_CACHE, "w"), indent=1)
    print(f"promoters: {len(proms)} ok, {len(failed)} failed", flush=True)

    res = {
        "meta": {
            "model": "AIDO.DNA-300M", "hf_repo": "genbio-ai/AIDO.DNA-300M",
            "architecture": "RNABert (GenBioBERT): 24 layers, hidden 1024, RoPE, 4000-token context",
            "usage": "zero-shot embedding extraction; no fine-tuning, no labels, no gradients",
            "promoter_definition": f"TSS-{PROMOTER_UP}..TSS, strand-aware, Ensembl REST GRCh38",
            "readout": "cosine(mean final-layer nucleotide embedding, E6 centroid)",
            "scored_universe_genes": n_universe, "seed": SEED,
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
        },
        "groups": {k: {"n_requested": len(v), "genes": v} for k, v in groups.items()},
        "promoters_failed": failed,
        "per_gene": {}, "group_stats": {}, "tests": {},
        "status": "pending", "error": None,
    }

    try:
        if len(proms) < 20:
            raise RuntimeError(f"too few promoters retrieved ({len(proms)})")
        t0 = time.time()
        m, dev, stubbed, wpath = load_model()
        res["meta"]["device"] = dev
        res["meta"]["weights_dir"] = wpath
        res["meta"]["stubbed_side_dependencies"] = stubbed
        res["meta"]["hidden_size"] = int(m.get_embedding_size())
        print(f"model ready on {dev} in {time.time()-t0:.1f}s (stubbed: {stubbed})", flush=True)

        import numpy as np
        vecs, t1 = {}, time.time(); done = 0
        emb_cache = {}
        EMB_CACHE = os.path.join(ROOT, "data/aido_embeddings.json")
        if os.path.exists(EMB_CACHE):
            try: emb_cache = json.load(open(EMB_CACHE))
            except Exception: emb_cache = {}
        hidim = int(m.get_embedding_size())
        for i, (g, p) in enumerate(proms.items(), 1):
            if g in emb_cache and len(emb_cache[g]) == hidim:
                vecs[g] = np.array(emb_cache[g], dtype=float); done += 1; continue
            v, ntok = embed(m, dev, p["seq"])
            vecs[g] = v; emb_cache[g] = v.tolist()
            if i % 5 == 0 or i == 1:
                json.dump(emb_cache, open(EMB_CACHE, "w"))
                print(f"  [{i}/{len(proms)}] {g} tokens={ntok} ({time.time()-t1:.0f}s) [cached {len(emb_cache)}]", flush=True)
        json.dump(emb_cache, open(EMB_CACHE, "w"))
        print(f"embedded {len(vecs)} promoters ({done} reused) in {time.time()-t1:.0f}s", flush=True)

        gof = {g: k for k, gl in groups.items() for g in gl}
        ref = [vecs[g] for g in groups["E6_clinical_reference"] if g in vecs]
        if len(ref) < 10:
            raise RuntimeError(f"reference group too small ({len(ref)})")
        c = np.mean(np.stack(ref), axis=0); c = c / np.linalg.norm(c)

        for g, v in vecs.items():
            res["per_gene"][g] = {
                "group": gof.get(g, "unassigned"),
                "region": proms[g]["region"], "len": proms[g]["len"], "gc": proms[g]["gc"],
                "cos_to_E6_centroid": round(float(np.dot(v, c)), 5),
            }

        def vals(k):
            return [d["cos_to_E6_centroid"] for d in res["per_gene"].values() if d["group"] == k]

        def gcs(k):
            return [d["gc"] for d in res["per_gene"].values() if d["group"] == k]

        from scipy import stats
        for k in groups:
            a = vals(k)
            if a:
                res["group_stats"][k] = {
                    "n_embedded": len(a),
                    "cos_mean": round(float(np.mean(a)), 5), "cos_sd": round(float(np.std(a, ddof=1)), 5) if len(a) > 1 else None,
                    "cos_median": round(float(np.median(a)), 5),
                    "gc_mean": round(float(np.mean(gcs(k))), 4),
                }
        for lbl, k1, k2 in [("ECS_vs_housekeeping", "ECS_specific", "housekeeping_control"),
                            ("ECS_vs_random", "ECS_specific", "random_background"),
                            ("housekeeping_vs_random", "housekeeping_control", "random_background")]:
            a, b = vals(k1), vals(k2)
            if len(a) >= 3 and len(b) >= 3:
                u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                res["tests"][lbl] = {"n1": len(a), "n2": len(b), "U": float(u),
                                     "p_two_sided": float(f"{p:.3g}"),
                                     "median_diff": round(float(np.median(a) - np.median(b)), 5)}
        nonref = [(d["gc"], d["cos_to_E6_centroid"]) for d in res["per_gene"].values()
                  if d["group"] != "E6_clinical_reference"]
        if len(nonref) >= 8:
            rho, pv = stats.spearmanr([x[0] for x in nonref], [x[1] for x in nonref])
            res["tests"]["gc_confounder_spearman"] = {
                "n": len(nonref), "rho": round(float(rho), 4), "p": float(f"{pv:.3g}"),
                "interpretation": "association between promoter GC content and similarity to the E6 centroid",
            }
        res["status"] = "ok"
    except Exception as e:
        res["status"] = "pending"
        res["error"] = f"{type(e).__name__}: {e}"

    res["honest_notes"] = [
        "Promoter sequences are real Ensembl GRCh38 regions; coordinates are recorded per gene.",
        "AIDO.DNA-300M is used zero-shot: it was never fine-tuned on PDAC and never saw any endpoint label.",
        "The E6 reference set is disjoint from every test/control group, so no label leaks into the readout.",
        "This is an exploratory sequence-level check; it does not enter ECS and does not affect any benchmark metric.",
        "GC content is reported and tested because it is the dominant compositional confounder in promoters.",
    ]
    json.dump(res, open(OUT, "w"), indent=1)
    print(f"WROTE {OUT} status={res['status']}" + (f" error={res['error']}" if res["error"] else ""), flush=True)


if __name__ == "__main__":
    main()
