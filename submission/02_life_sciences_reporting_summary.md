# Life Sciences Reporting Summary — Completed for *Nature*

> 说明：本稿为**计算/网络医学（computational & network-medicine）研究**，无湿实验、无细胞系/动物/植物/临床受试者操作。
> 凡不适用项均标注 **N/A（计算研究）**，并给出对应计算等价说明，以满足 Nature 报告规范。

---

## Experimental design

**1. Sample size**
- 主分析涵盖 **{n_genes} 个基因**（PDAC 评分宇宙，详见 `benchmark_v14.json` 的 `meta.n_genes`）。
- 六个端点阳性集规模：E1/E2/E3/E4/E5 由公开组学阈值定义；E6 临床验证集 n = 35（ClinicalTrials.gov 映射，详见 Methods）。
- *计算等价说明*：样本量由证据层与端点定义确定，非统计功效估算；所有 AUROC/AUPRC 均附 1,000–2,000 次 bootstrap 置信区间。

**2. Data exclusions**
- **无事后排除**。缺失证据值按该层最小值插补（Methods 已述）。
- E4（Open Targets PDAC top-500）被显式标记为**非独立端点**（因其属于 ECS 构成层），用作构造性冗余检验，而非外部验证。E6 从未进入 ECS 构建、权重、超参或候选筛选。

**3. Replication**
- 所有关键结论经多重独立重采样重复：DeLong 配对检验、1,000 次 bootstrap、1,000 次支持层置换（permutation）、超参 α∈{0,0.2,0.4,0.6,0.8,1.0} 与 D 权重扫描。
- 两路**完全独立**验证（GenBio 虚拟细胞共识 + AIDO.DNA-300M 零样本序列检查）对 ECS 自体构建无泄漏。

**4. Randomization**
- **N/A（计算研究）**。无随机分组；置换检验以固定随机种子（scripts 内记录）对支持层做 1,000 次随机打乱作为零分布。

**5. Blinding**
- **N/A（计算研究）**。端点标签在 ECS 构建前即固定，E6 临床端点独立于 ECS（互不泄漏）。

---

## Materials & reagents

| 项目 | 状态 | 说明 |
|---|---|---|
| Antibodies | N/A | 无免疫实验 |
| Eukaryotic cell lines | N/A | 无细胞培养 |
| Palaeontology and archaeology | N/A | — |
| Animals and dual use research of concern | N/A | 无动物实验 |
| Plants | N/A | 无植物材料 |
| Clinical data (human research participants) | N/A（聚合公开数据） | 仅使用 ClinicalTrials.gov 公开试验登记元数据（干预→药物→靶点映射），不含个体受试者数据，无需伦理审批/知情同意；符合公开数据使用条款 |

---

## Dual use research of concern (DURC)
- **否**。本稿为靶点优先级排序的方法学/基准研究，不涉及可造成广泛危害的双用途实验。

---

## Reporting for specific studies
- **Bioinformatics / computational modelling**：已遵循 — 数据来源与版本（STRING v11、DepMap 26Q1、Open Targets、IMPC、HPA、DrugBank、ClinicalTrials.gov）、代码开源（MIT，GitHub 仓库）、随机种子与全部过程参数在 scripts 中记录、全部数值指标可一键复现。

---

## 提交配套文件
- 完整 Reporting Summary 字段已落入本文件；投稿系统内对应表单可逐项勾选/粘贴。
- 关联材料：`PDAC_convergent_evidence_Nature.docx`、`01_cover_letter.md`、`03_submission_checklist.md`、GitHub 仓库 `zuoxianbo/pdac-convergent-evidence-v14`。
