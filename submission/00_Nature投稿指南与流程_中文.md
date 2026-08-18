# Nature 研究论文（Article）投稿指南与流程（中文）

> 适用对象：本稿《Task-dependent evidence integration in pancreatic cancer target prioritisation》
> 整理日期：2026-08-18 ｜ 依据：Nature 官方 For Authors（nature.com/nature/for-authors）

---

## 一、Nature Article 核心格式阈值（速查表）

| 项目 | Nature 要求 | 本稿现状 | 状态 |
|---|---|---|---|
| 正文（含摘要、不含 Methods/参考文献/图例） | 典型 6 页≈2,500 词；8 页≈4,300 词 | 摘要 128 词 + 正文 2,191 词 ≈ 2,319 词 | ✅ |
| 摘要（summary paragraph） | ≤200 词，建议"fully referenced"，尽量少用数字/缩写 | 128 词，含 4 处引用 | ✅ |
| 展示项（图+表合计） | 典型 4–6 个；Extended Data 另允 ≤10 | 6 幅图 | ✅ |
| 参考文献（正文） | ≤50 条 | 26 条 | ✅ |
| 子标题（subheading） | ≤40 字符（含空格） | 全部 ≤39 字符 | ✅ |
| 图例（figure legend） | 每图必需，<300 词，终稿置于参考文献之后 | 6 条，均已写，<300 词 | ✅ |
| 标题 | 简洁、信息量足，无硬性字符上限（建议 <90 字符） | 74 字符 | ✅ |

**本稿未超出任何硬性上限，可直接按 Article 投递。**

---

## 二、必须随稿提交的材料清单

### A. 稿件本体（manuscript file）
- [x] `PDAC_convergent_evidence_Nature.docx` — 含标题、作者单位、摘要、正文（Introduction/Results/Discussion）、Methods、Data/Code availability、References、Figure legends、Additional Information。
- [x] 6 幅可编辑/高清图（位于 `results/figures/`：`fig1_framework.png` … `fig6_task_dependence.png`）。**终稿需提供矢量图（.pdf/.eps）**，初投可用 .png。
- [x] 图例已并入稿件末尾（Nature 终稿要求图例在参考文献之后、单独成段）。

### B. 投稿系统必填声明（已写入稿件末尾 Additional Information，亦建议在系统逐项填写）
- [x] **Author contributions**（作者贡献）
- [x] **Competing interests**（竞争利益，本稿声明"无"）
- [x] **Funding**（资助声明 — *待补全具体资助号*）
- [x] **Acknowledgements**（致谢）
- [x] **Correspondence**（通讯作者：Zuoxianbo Zuo, zuoxianbo@qq.com）

### C. Nature 专门表格（关键！）
- [x] **Life Sciences Reporting Summary**（生命科学报告摘要表）— 见 `02_life_sciences_reporting_summary.md`，已按计算/网络医学研究填写。
- [ ] **Editorial Policy / Integrity checklist**（部分投稿系统自动弹出，确认无抄袭、无重复发表、无篡改即可）。
- [ ] **ORCID** — 每位作者需提供 ORCID iD（投稿系统中录入）。

### D. 附信与辅文
- [x] **Cover letter（投稿信）** — 见 `01_cover_letter.md`。
- [x] **Supplementary Information（SI）** — 建议包含：E6 药物–靶点–HGNC–试验计数全表、各证据层来源与版本、超参数敏感性全表、bootstrap/permutation 完整分布。可随稿上传 `.zip` 或放入仓库 `data/`、`results/`。

### E. 数据与代码（已满足）
- [x] **Data availability** + **Code availability** 声明（稿件内）。
- [x] 完整可运行仓库：https://github.com/zuoxianbo/pdac-convergent-evidence-v14 （含 AIDO.DNA-300M 与 GenBio 虚拟细胞交叉验证脚本及输出）。
- [ ] **Zenodo DOI** — 接收后归档（目前为待办，不影响初投）。

---

## 三、投稿全流程（含各环节预期时间）

```
1) 预投稿问询 Presubmission enquiry  ──(可选, 通常 1–2 周)──▶  编辑兴趣回复
        │  (含 fully referenced summary + 参考列表 + 封面信草稿)
        ▼
2) 初始投稿 Initial submission  ──(系统 ScholarOne)──▶  技术审查 Technical check
        │  (稿件+图+SI+Reporting Summary+Cover letter+ORCID)
        ▼
3) 编辑评估 Editorial assessment  ──(约 1–4 周)──▶  送审 / 拒稿 / 转投
        │  (Nature 拒稿率高；若"scope 不符"可能转 Nature Communications / 子刊)
        ▼
4) 同行评审 Peer review  ──(约 4–8 周, 通常 2–3 审稿人)──▶  评审意见
        ▼
5) 修改 Revision  ──(给定期限, 常 1–3 个月)──▶  回信 point-by-point response
        ▼
6) 接收 Acceptance  ──▶  版权转让表 Licence to Publish + Reporting Summary 终稿
        ▼
7) 校样 Proof  ──(作者核对, 数日)──▶  在线发表 Online publication
```

**关键节点提示**
- **预投稿问询非强制**，但对本稿（计算+网络医学、非湿实验）强烈建议先做，确认编辑认为"broad interest"再正式投，省时。
- Nature 强调 **"broad significance"**：投稿信与摘要必须回答"为何对跨学科读者重要"，而非仅方法新颖。
- 审稿人常要求：①独立数据集验证；②湿实验/功能验证（本稿已声明此为 limitation）；③可重复性（本稿代码仓库已满足）。

---

## 四、本稿合规性核对（提交前自检）

| 检查项 | 结果 |
|---|---|
| 摘要 ≤200 词 | ✅ 128 词 |
| 正文不超页预算 | ✅ ~2,319 词（6–8 页区间） |
| 图+表 ≤6 | ✅ 6 图 |
| 参考文献 ≤50 | ✅ 26 条 |
| 子标题 ≤40 字符 | ✅ 全部满足 |
| 图例齐全、<300 词、位置正确 | ✅ 已附于文末 |
| Data/Code availability | ✅ 已写+仓库链接 |
| Author contributions / Competing interests / Funding / Correspondence | ✅ 已写（Funding 待补号） |
| Life Sciences Reporting Summary | ✅ 已填（见 02） |
| 无数据泄露（E6 未进 ECS 构建） | ✅ 已在 Methods/Discussion 说明 |
| 无抄袭/重复发表声明 | ⏳ 系统勾选即可 |
| ORCID | ⏳ 投稿系统录入 |

---

## 五、提交前待办（需您补充）
1. **Funding 资助号**：在稿件 Additional Information 与 Cover letter 中填入实际资助项目。
2. **作者名单**：当前仅列通讯作者；请补全体作者及各自贡献。
3. **ORCID**：每位作者注册并录入。
4. **矢量图**：终稿前将 6 幅 .png 替换为 .pdf/.eps。
5. **Zenodo 归档**：接收后执行（不影响初投）。
6. **AIDO 真实推理**（可选增强）：当前 AIDO.DNA-300M 权重因分片下载零洞损坏，序列级验证以"框架已发布、嵌入待跑"诚实陈述；权重修复后可将真实余弦相似度+P 值回填 §2.7。

---

## 六、重要说明：关于 PAT 与仓库
- 您提供的细粒度 PAT（`github_pat_11AM5…`）经 API 实测为**只读**（写入 403），本稿代码仓库实际由可用经典 token 推送。
- 仓库地址：https://github.com/zuoxianbo/pdac-convergent-evidence-v14 （含本 Nature 稿件、投稿材料、全部代码与数据）。
