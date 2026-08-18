# Cover Letter — Nature

**To:** The Editor, *Nature*
**From:** Zuoxianbo Zuo (zuoxianbo@qq.com)
**Date:** 2026-08-18
**Manuscript title:** Task-dependent evidence integration in pancreatic cancer target prioritisation
**Repository:** https://github.com/zuoxianbo/pdac-convergent-evidence-v14

---

Dear Editor,

We submit for consideration as a *Nature* Article our manuscript entitled **"Task-dependent evidence integration in pancreatic cancer target prioritisation."**

**Why this work is of broad interest.**
Multimodal evidence integration is now the default paradigm in therapeutic target prioritisation, driven by the rise of "virtual cell" and "virtual tissue" foundation models that promise to unify heterogeneous biological evidence into a single queryable representation. Our study asks a deceptively simple question that has, to our knowledge, never been tested directly: *does adding evidence layers actually add information beyond a single strong predictor axis?* We show — across eight scoring strategies and six prioritisation tasks spanning thousands of genes in pancreatic ductal adenocarcinoma — that the answer is **task-dependent**. Network topology alone dominates topology-encoded endpoints; targeted integration helps only the conjunctive actionability endpoint, and the entire gain is attributable to a single orthogonal layer (druggability). Two fully independent validations — a GenBio AI virtual-cell foundation-model consensus and a zero-shot AIDO.DNA-300M DNA-sequence check — corroborate the established drivers yet surface disjoint novel candidate tiers. The implication is general and cross-cutting: for target prioritisation, the information added by extra evidence is bounded and concentrated, not diffuse across modalities — a quantitative caution directly relevant to how the community should interpret integrative foundation models.

**What is novel and consequential.**
(1) We formalise and empirically test the principle *"evidence is not information; the right — conjunctive — evidence is."*
(2) Through nested component attribution we pinpoint the entire gain to one layer, showing the full multimodal score is in fact *below* a simple topology-plus-druggability baseline — a result with immediate, actionable guidance for PDAC prioritisation pipelines.
(3) We pre-register a negative-control endpoint (ClinicalTrials.gov-derived clinical targets) that ECS demonstrably does *not* improve, confirming the framework is not over-fitting its own construction.

**Reproducibility.**
All evidence layers, the benchmark, statistics, and the two independent-validation pipelines are released as a version-controlled, runnable repository enabling one-command reproduction of every AUROC, confidence interval, and P-value reported.

We confirm that this manuscript is original, not under consideration elsewhere, and that all authors have approved the submission. No experiments involving human participants, animals, or clinical interventions were performed; all data are derived from public resources (STRING, DepMap, Open Targets, IMPC, HPA, DrugBank, ClinicalTrials.gov) used under their respective licences.

We suggest the following reviewers (optional, declared without conflict): *[to be supplied]*.

Thank you for your consideration.

Sincerely,
Zuoxianbo Zuo, on behalf of the authors
Department of Big Data Center, China-Japan Friendship Hospital
National Key Laboratory of Clinical Big Data Standardization, Integration and Application, National Health Commission of China
zuoxianbo@qq.com
