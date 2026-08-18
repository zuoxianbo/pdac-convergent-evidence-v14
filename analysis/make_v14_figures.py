# -*- coding: utf-8 -*-
"""V14 publication figures (6 panels). Reads benchmark_v14.json."""
import json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Circle, FancyArrowPatch, FancyBboxPatch
from matplotlib_venn import venn3

import os
ROOT = os.environ.get("V14_ROOT", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data"))
RES = os.environ.get("V14_RES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
FIG = os.path.join(RES, "figures")
import os; os.makedirs(FIG, exist_ok=True)
B = json.load(open(f"{RES}/benchmark_v14.json"))

plt.rcParams.update({"font.size": 9.5, "axes.titlesize": 10.5, "axes.titleweight": "bold",
                     "figure.dpi": 150, "axes.grid": False, "font.family": "DejaVu Sans"})

C_ECS = "#1f6f54"; C_STR = "#2c5f8a"; C_ACC = "#c0392b"; C_GREY = "#7f8c8d"

def f(x, n=3):
    try: return f"{x:.{n}f}"
    except Exception: return str(x)

# ----------------------------------------------------------------- FIG 1: Evidence Convergence Framework (conceptual)
fig = plt.figure(figsize=(13.2, 5.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, 13.2); ax.set_ylim(0, 5.0)
ax.text(6.6, 4.78, "Evidence Convergence Framework: when does multimodal integration add information?",
        ha="center", fontsize=13, weight="bold", color="#222")
ax.text(6.6, 4.42, "Evidence is not information; orthogonal evidence is.",
        ha="center", fontsize=10.5, style="italic", color=C_ECS)
# left branch: topology
def box(x, y, w, h, fc, ec, title, body, fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.06",
                                fc=fc, ec=ec, lw=1.5))
    ax.text(x + w/2, y + h - 0.34, title, ha="center", va="center", fontsize=9.6, weight="bold", color=ec)
    ax.text(x + w/2, y + h/2 - 0.12, body, ha="center", va="center", fontsize=fs, color="#222", wrap=True)
box(0.4, 1.0, 3.9, 1.5, "#eef3f8", C_STR, "Single strong predictor",
    "STRING network topology\nencodes what the network already\n'knows' about a gene", 8.4)
ax.annotate("", xy=(2.35, 0.92), xytext=(2.35, 2.5), arrowprops=dict(arrowstyle="-|>", color=C_STR, lw=2))
ax.text(0.3, 0.62, "Topology-encoded tasks", fontsize=9.6, weight="bold", color=C_STR)
for i, (lab, v) in enumerate([("E1 pan-dependency", "fusion adds little"),
                               ("E4 genetic disease", "(OT genetics overlaps)"),
                               ("E6 clinical targets", "(already network-central)")]):
    yy = 0.30 - i*0.0
ax.text(0.3, 0.34, "E1 pan-dependency   E4 genetic disease   E6 clinical targets\n— fusion adds little beyond topology —",
        fontsize=8.6, color=C_STR)
# center operator
ax.add_patch(Circle((6.6, 2.5), 0.62, fc="#fff", ec=C_ECS, lw=2))
ax.text(6.6, 2.5, "ECS =\nD×(1+α·PHI)", ha="center", va="center", fontsize=9.2, weight="bold", color=C_ECS)
ax.annotate("", xy=(4.35, 2.5), xytext=(2.3, 2.5), arrowprops=dict(arrowstyle="-|>", color="#555", lw=1.8, ls="--"))
ax.text(3.3, 2.72, "STRING (D)", fontsize=8.2, color="#555")
# right branch: orthogonal evidence
box(8.9, 1.0, 3.9, 1.5, "#eaf4ef", C_ECS, "Orthogonal evidence integrated",
    "genetics · druggability · tissue\nspecificity · prognosis · driver\n(convergent support, PHI)", 8.4)
ax.annotate("", xy=(6.6, 1.88), xytext=(8.9, 2.1), arrowprops=dict(arrowstyle="-|>", color=C_ECS, lw=2))
ax.text(8.6, 2.66, "Conjunctive task", fontsize=9.6, weight="bold", color=C_ECS)
ax.text(8.2, 1.78, "E3 actionable target\n— fusion ADDS information\n(ΔAUROC = 0.090, DeLong p<0.001)",
        fontsize=8.8, color=C_ECS, weight="bold")
# bottom principle
ax.add_patch(FancyBboxPatch((0.4, 0.0), 12.4, 0.22, boxstyle="round,pad=0.01", fc="#f4f1e8", ec="#b7a96a", lw=1))
ax.text(6.6, 0.11, "Central thesis: multimodal evidence is worth integrating ONLY when the target property is itself conjunctive and not already encoded by a single strong predictor axis.",
        ha="center", va="center", fontsize=8.4, color="#6b5d1f")
plt.savefig(f"{FIG}/fig1_framework.png", bbox_inches="tight"); plt.close()
print("fig1 done")

# ----------------------------------------------------------------- FIG 2: Benchmark landscape (AUROC + AUPRC)
scorers = ["STRING_centrality", "ECS_proposed", "SimpleMean", "MutationFreq", "IMPC", "GeneticConstraint", "Druggability", "DepMap_oracle"]
ep_order = ["E1_pan_dependency", "E2_selective_dependency", "E3_actionable_target", "E4_genetic_disease", "E5_crc_transfer", "E6_clinical_validation"]
ep_lbl = ["E1\npan-dep", "E2\nselective", "E3\nactionable", "E4\ngenetic", "E5\nCRC", "E6\nclinical"]
A = np.array([[B["robustness_auroc"][s][e][0] if B["robustness_auroc"][s][e] else np.nan for e in ep_order] for s in scorers])
P = np.array([[B["robustness_auprc"][s][e][0] if B["robustness_auprc"][s][e] else np.nan for e in ep_order] for s in scorers])
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.4))
for axd, M, title in [(axes[0], A, "AUROC"), (axes[1], P, "AUPRC")]:
    im = axd.imshow(M, cmap="viridis", aspect="auto", vmin=0.45, vmax=1.0)
    axd.set_xticks(range(len(ep_order))); axd.set_xticklabels(ep_lbl, fontsize=8.5)
    axd.set_yticks(range(len(scorers))); axd.set_yticklabels([s.replace("_", " ") for s in scorers], fontsize=8.5)
    for i in range(len(scorers)):
        for j in range(len(ep_order)):
            v = M[i, j]
            if not math.isnan(v):
                axd.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=7.6,
                         color="white" if v > 0.72 else "black")
    axd.set_title(title, fontsize=11)
    # mark best per endpoint
    for j in range(len(ep_order)):
        col = M[:, j]; col = col[~np.isnan(col)]
        if len(col):
            best_i = int(np.nanargmax(M[:, j]))
            axd.add_patch(plt.Rectangle((j-0.5, best_i-0.5), 1, 1, fill=False, ec="red", lw=2.2))
    fig.colorbar(im, ax=axd, fraction=0.046, pad=0.02)
fig.suptitle("Fig 2. Benchmark landscape: 8 scorers x 6 endpoints", fontsize=12, weight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.96]); plt.savefig(f"{FIG}/fig2_landscape.png", bbox_inches="tight"); plt.close()
print("fig2 done")

# ----------------------------------------------------------------- FIG 3: Component attribution & incremental gain (E3)
ca = B["component_attribution"]
models = ["M1_STRING", "M2_STRING_druggability", "M3_STRING_genetics", "M4_STRING_drug_genetics",
          "M5_STRING_drug_genetics_tissue", "M6_Full_ECS"]
mlbl = ["STRING", "+druggability", "+genetics", "+drug+genetics", "+tissue", "Full ECS"]
vals = [ca[m]["auroc"] for m in models]
deltas = [ca[m]["delta_vs_STRING"] for m in models]
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.0))
ax = axes[0]
bars = ax.bar(range(len(models)), vals, color=[C_STR] + [C_ECS]*(len(models)-1), edgecolor="white")
ax.axhline(B["string_base_auroc"]["E3_actionable_target"], ls="--", color=C_STR, lw=1.5, label=f"STRING base = {f(B['string_base_auroc']['E3_actionable_target'])}")
for i, v in enumerate(vals):
    ax.text(i, v+0.006, f"{v:.3f}", ha="center", fontsize=8.4, weight="bold")
ax.set_xticks(range(len(models))); ax.set_xticklabels(mlbl, rotation=20, ha="right", fontsize=8.2)
ax.set_ylim(0.6, 0.86); ax.set_ylabel("AUROC (E3 actionable target)")
ax.set_title("Nested component attribution: the E3 gain is entirely druggability-driven (\u0394=0.087)", fontsize=10.5)
ax.legend(fontsize=8.5, loc="lower right")
ax2 = axes[1]
inc = [deltas[0]] + [ca[models[i]]["delta_vs_previous"] for i in range(1, len(models))]
colors = [C_STR] + [C_ECS if x > 0.01 else C_GREY for x in inc[1:]]
ax2.bar(range(len(models)), inc, color=colors, edgecolor="white")
for i, x in enumerate(inc):
    ax2.text(i, x + (0.002 if x >= 0 else -0.004), f"{x:+.3f}", ha="center", fontsize=8.0)
ax2.axhline(0, color="k", lw=1)
ax2.set_xticks(range(len(models))); ax2.set_xticklabels(mlbl, rotation=20, ha="right", fontsize=8.2)
ax2.set_ylabel("Incremental ΔAUROC vs previous step")
ax2.set_title("Incremental gain per added component", fontsize=10.5)
plt.tight_layout(); plt.savefig(f"{FIG}/fig3_component_attribution.png", bbox_inches="tight"); plt.close()
print("fig3 done")

# ----------------------------------------------------------------- FIG 4: ECS-specific target discovery (4 panels)
es = B["ecs_specific"]
ALL = set(json.load(open(f"{ROOT}/evidence_layers_v11.json"))["genes"])
ev = json.load(open(f"{ROOT}/evidence_layers_v11.json")); LAY = ev["layers"]
def gv(g, L):
    r = LAY.get(L, {}).get(g); return r["norm"] if (r and r.get("present")) else None
STRING = {g: gv(g, "string_centrality") for g in ALL}
# E3 actionable target (essential intersect druggable)
dpd = json.load(open(f"{ROOT}/depmap_pdac_dependency.json"))
E1 = {g: (bool(d.get("essential")) if isinstance(d, dict) else False) for g, d in dpd.items()}
DRUG = {g: gv(g, "druggability") for g in ALL}
E3 = {g: (E1.get(g, False) and DRUG.get(g) is not None) for g in ALL}
ALPHA=0.60; DEPW={"string_centrality":0.80,"mutation_freq":0.10,"impc_animal_ko":0.10}
SUP=["cancer_driver","ot_genetics_pdac","druggability","hpa_pdac_prognostic","hpa_rna_tissue_spec"]
def ecs(g):
    num=den=0.0
    for L,w in DEPW.items():
        v=gv(g,L)
        if v is not None: num+=w*v; den+=w
    if den==0: return None
    D=num/den
    vs=[gv(g,L) for L in SUP if gv(g,L) is not None]
    return D*(1+ALPHA*float(np.mean(vs))) if vs else D
ECS={g:ecs(g) for g in ALL}
gs=[g for g in ALL if ECS[g] is not None and STRING[g] is not None]
esc_rank={g:i for i,g in enumerate(sorted(gs,key=lambda x:-ECS[x]))}
str_rank={g:i for i,g in enumerate(sorted(gs,key=lambda x:-STRING[x]))}
N=len(gs)
xs=[str_rank[g] for g in gs]; ys=[esc_rank[g] for g in gs]
# ECS-specific actionable (top100 ECS outside top500 STRING)
top100_ecs=set(sorted(gs,key=lambda x:-ECS[x])[:100]); top500_str=set(sorted(gs,key=lambda x:-STRING[x])[:500])
ecs_spec_act=[g for g in (top100_ecs-top500_str) if E3[g]]
fig, axes = plt.subplots(2, 2, figsize=(13.6, 9.2))
# 4A rank scatter
axA=axes[0,0]
axA.scatter(xs, ys, s=3, c="#bcc6cf", alpha=0.35, edgecolors="none")
axA.scatter([str_rank[g] for g in ecs_spec_act],[esc_rank[g] for g in ecs_spec_act], s=22, c=C_ACC, edgecolors="none", label=f"ECS-specific actionable (n={len(ecs_spec_act)})")
axA.plot([0,N],[0,N], ls="--", color="#888", lw=1)
axA.set_xlabel("STRING rank (1 = best)"); axA.set_ylabel("ECS rank (1 = best)")
axA.set_title("4A. ECS vs STRING re-ranking", fontsize=10.5)
axA.invert_yaxis(); axA.legend(fontsize=8, loc="lower right")
# 4B top ECS-specific genes by delta
deltas_sorted=sorted(ecs_spec_act, key=lambda x:-(ECS[x]-STRING[x]))[:22]
gn=[g for g in deltas_sorted]; dv=[ECS[g]-STRING[g] for g in deltas_sorted]
axB=axes[0,1]
axB.barh(range(len(gn)), dv, color=[C_ACC if E3[g] else C_GREY for g in gn])
axB.set_yticks(range(len(gn))); axB.set_yticklabels(gn, fontsize=8.5)
axB.invert_yaxis(); axB.set_xlabel("ECS - STRING score")
axB.set_title("4B. Top ECS-specific genes (STRING low, ECS high)", fontsize=10.5)
# 4C evidence heatmap for top ECS-specific genes
axC=axes[1,0]
cols=["string_centrality","ot_genetics_pdac","druggability","hpa_rna_tissue_spec","hpa_pdac_prognostic","cancer_driver"]
colnames=["STRING","genetics","drug","tissue","progn","driver"]
Hg=ecs_spec_act[:18]
H=np.array([[gv(g,c) if gv(g,c) is not None else np.nan for c in cols] for g in Hg])
im=axC.imshow(H, cmap="viridis", aspect="auto", vmin=0, vmax=1)
axC.set_xticks(range(len(cols))); axC.set_xticklabels(colnames, fontsize=8.5)
axC.set_yticks(range(len(Hg))); axC.set_yticklabels(Hg, fontsize=8.5)
axC.set_title("4C. Evidence profile of ECS-specific targets", fontsize=10.5)
fig.colorbar(im, ax=axC, fraction=0.046, pad=0.02)
# 4D enrichment curve
axD=axes[1,1]
ks=[10,20,50,100,200]
for name,c in [("ECS",C_ECS),("STRING",C_STR),("SimpleMean",C_GREY)]:
    ec=es["enrichment_curve"][name]
    axD.plot(ks,[ec[str(k)]["fold"] for k in ks], "o-", color=c, label=name)
axD.axhline(1.0, ls="--", color="#999", lw=1)
axD.set_xlabel("Top-K genes"); axD.set_ylabel("Actionable-fold enrichment")
axD.set_title("4D. Actionable enrichment curve", fontsize=10.5)
axD.legend(fontsize=8.5)
fig.suptitle("Fig 4. ECS-specific targets: convergence re-ranks STRING-buried actables", fontsize=12, weight="bold")
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig(f"{FIG}/fig4_ecs_specific.png", bbox_inches="tight"); plt.close()
print("fig4 done")

# ----------------------------------------------------------------- FIG 5: Robustness controls
fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.4))
# 5A alpha sensitivity (E3 & E1)
asens=B["alpha_sensitivity"]
alphas=sorted([float(k[1:]) for k in asens.keys()])
e3v=[asens[f"a{a}"]["E3_actionable_target"] for a in alphas]
e1v=[asens[f"a{a}"]["E1_pan_dependency"] for a in alphas]
axA=axes[0]
axA.plot(alphas, e3v, "o-", color=C_ECS, label="E3 actionable")
axA.plot(alphas, e1v, "s-", color=C_STR, label="E1 pan-dependency")
axA.axvline(0.6, ls="--", color="#999", lw=1)
axA.set_xlabel("α (convergence weight)"); axA.set_ylabel("AUROC")
axA.set_title("5A. Hyper-parameter sensitivity", fontsize=10.5); axA.legend(fontsize=8.5)
# 5B permutation
perm=B["permutation"]
axB=axes[1]
# null band (summary only) + observed; draw a light null-density proxy from CI
nm=perm["null_mean"]; nci=perm["null_ci"]; obs=perm["observed"]
xs_=np.linspace(nci[0]-0.02, nci[1]+0.02, 200)
from scipy.stats import norm
null_density=norm.pdf(xs_, nm, (nci[1]-nci[0])/3.92)
axB.fill_between(xs_, null_density, color="#dfe6ee", label="permutation null (95% band)")
axB.axvline(nm, color=C_GREY, lw=1.5, label=f"null mean {f(nm)}")
axB.axvspan(nci[0], nci[1], color="#dfe6ee", label="permutation null 95%")
axB.axvline(nm, color=C_GREY, lw=1.5, label=f"null mean {f(nm)}")
axB.axvline(obs, color=C_ECS, lw=2.2, label=f"observed {f(obs)}")
axB.set_xlabel("E3 AUROC under layer permutation"); axB.set_title(f"5B. Permutation control (p={perm['p_ge_observed']:.3f})", fontsize=10.5)
axB.legend(fontsize=8)
# 5C annotation bias
ab=B["annotation_bias"]
axC=axes[2]
cats=["ECS","STRING"]; raw=[ab["ECS_proposed"]["raw_auroc"], ab["STRING_centrality"]["raw_auroc"]]
res=[ab["ECS_proposed"]["residualized_auroc"], ab["STRING_centrality"]["residualized_auroc"]]
x=np.arange(2); w=0.35
axC.bar(x-w/2, raw, w, color=C_ECS, label="raw")
axC.bar(x+w/2, res, w, color=C_STR, label="residualised on annotation")
for i,v in enumerate(raw): axC.text(i-w/2, v+0.005, f"{v:.3f}", ha="center", fontsize=8)
for i,v in enumerate(res): axC.text(i+w/2, v+0.005, f"{v:.3f}", ha="center", fontsize=8)
axC.set_xticks(x); axC.set_xticklabels(cats, fontsize=9)
axC.set_ylim(0.6, 0.86); axC.set_ylabel("AUROC (E3)")
axC.set_title("5C. Annotation-bias control", fontsize=10.5); axC.legend(fontsize=8.5)
plt.tight_layout(); plt.savefig(f"{FIG}/fig5_robustness.png", bbox_inches="tight"); plt.close()
print("fig5 done")

# ----------------------------------------------------------------- FIG 6: Task-dependence & evidence complementarity
fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2))
# 6A per-endpoint delta, colored by role
cs=B["crc_shift"]
eps=["E1_pan_dependency","E2_selective_dependency","E3_actionable_target","E4_genetic_disease","E5_crc_transfer","E6_clinical_validation"]
elbl=["E1\npan-dep","E2\nselective","E3\nactionable","E4\ngenetic","E5\nCRC","E6\nclinical"]
deltas=[cs[e]["delta"] for e in eps]
rolecol={"negative_control":C_GREY,"positive":C_ECS,"weak":C_STR}
cols=[rolecol[cs[e]["role"]] for e in eps]
axA=axes[0]
axA.bar(range(len(eps)), deltas, color=cols, edgecolor="white")
for i,v in enumerate(deltas): axA.text(i, v+(0.004 if v>=0 else -0.008), f"{v:+.3f}", ha="center", fontsize=8.2)
axA.axhline(0, color="k", lw=1)
axA.set_xticks(range(len(eps))); axA.set_xticklabels(elbl, fontsize=8.4)
axA.set_ylabel("ECS - STRING ΔAUROC")
axA.set_title("6A. Integration helps only on the conjunctive task (E3)", fontsize=10.5)
from matplotlib.patches import Patch as _P
axA.legend(handles=[_P(color=C_ECS,label="E3: convergence helps"),_P(color=C_GREY,label="E1/E4/E6: negative controls"),_P(color=C_STR,label="E2: weak")], fontsize=8, loc="upper right")
# 6B complementarity scatter: |corr| vs incremental gain, per support layer
comp=B["evidence_complementarity"]
sup=SUP
xc=[abs(comp["support_corr_with_string"][L]) for L in sup]
yg=[comp["support_incremental_gain"][L] for L in sup]
axB=axes[1]
axB.scatter(xc, yg, s=70, c=C_ECS, edgecolors="k", zorder=3)
for i,L in enumerate(sup):
    axB.annotate(L.replace("_","\n"), (xc[i],yg[i]), textcoords="offset points", xytext=(6,4), fontsize=7.6)
axB.axvline(0.5, ls="--", color="#bbb", lw=1); axB.axhline(0, ls="--", color="#bbb", lw=1)
axB.set_xlabel("|correlation| with STRING (redundancy)")
axB.set_ylabel("Incremental E3 gain (STRING + layer)")
axB.set_title("6B. Evidence complementarity is concentrated in one layer", fontsize=10.5)
axB.text(0.02, 0.95, "druggability:\nthe only layer\nthat adds signal", transform=axB.transAxes, fontsize=7.8, color=C_ECS, va="top")
axB.text(0.98, 0.05, "other layers:\n~0 incremental gain", transform=axB.transAxes, fontsize=7.8, color=C_GREY, ha="right", va="bottom")
plt.tight_layout(); plt.savefig(f"{FIG}/fig6_task_dependence.png", bbox_inches="tight"); plt.close()
print("fig6 done")
print("ALL FIGURES DONE")
