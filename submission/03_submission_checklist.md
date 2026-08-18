# Nature 投稿材料清单与自检表（Submission Checklist）

> 稿件：《Task-dependent evidence integration in pancreatic cancer target prioritisation》
> 仓库：https://github.com/zuoxianbo/pdac-convergent-evidence-v14

## ✅ 已就绪（本包内含）

- [x] **Nature 格式稿件** `PDAC_convergent_evidence_Nature.docx`
  - 标题 + 作者单位 + 摘要（128 词 ≤200） + 正文（~2,319 词，6–8 页区间）
  - Introduction / Results / Discussion / Methods / Data & Code availability
  - References（26 条 ≤50） + Figure legends（6 条，<300 词/条） + Additional Information
- [x] **6 幅图** `results/figures/fig1_framework.png` … `fig6_task_dependence.png`
- [x] **投稿信** `submission/01_cover_letter.md`
- [x] **Life Sciences Reporting Summary** `submission/02_life_sciences_reporting_summary.md`
- [x] **中文投稿指南与流程** `submission/00_Nature投稿指南与流程_中文.md`
- [x] **代码与数据** GitHub 仓库（含 AIDO.DNA-300M 零样本 + GenBio 虚拟细胞交叉验证脚本及输出）
- [x] **End-matter 声明**（作者贡献 / 竞争利益 / 资助 / 致谢 / 通讯）已写入稿件末尾

## ⏳ 提交前需您补充

- [ ] **Funding 资助号**：填入稿件 Additional Information 与 Cover letter
- [ ] **全体作者名单 + 各自贡献**：当前仅列通讯作者
- [ ] **ORCID**：每位作者在 ScholarOne 录入
- [ ] **矢量图**：终稿将 6 幅 .png 替换为 .pdf/.eps
- [ ] **Supplementary Information 打包**：E6 药物–靶点全表、证据层版本、超参敏感性全表（可放仓库 `data/`、`results/` 或单独 .zip）
- [ ] **推荐/回避审稿人**（Cover letter 末段，可选）
- [ ] **Zenodo DOI**：接收后归档（不影响初投）
- [ ] **系统声明勾选**：无抄袭、无重复发表、无数据篡改

## 🔁 可选增强（不阻塞投稿）

- [ ] **AIDO 真实推理回填**：当前 AIDO.DNA-300M 权重因分片下载零洞损坏，§2.7 以"框架已发布、嵌入待跑"诚实陈述。权重修复后重跑 `analysis/run_aido_dna.py` → 真实余弦相似度 + Mann–Whitney P 值回填 §2.7，再推送更新。
- [ ] **湿实验/功能验证章节**：若后续补 CRISPR  viability / 药物敏感性数据，可补强 Discussion 的 limitation。

## 提交流程提醒

1. （建议）先做 **Presubmission enquiry**（附 fully referenced summary + 参考列表 + cover letter 草稿）。
2. 若编辑有意 → **Initial submission** 走 ScholarOne：上传稿件、图、SI、Reporting Summary、Cover letter，录入 ORCID。
3. 技术审查 → 编辑评估 → 送审（2–3 人，4–8 周）→ 修改（point-by-point response）→ 接收 → 校样 → 在线发表。
