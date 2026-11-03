# 论文初稿说明

- `manuscript.md`：英文初稿（JER 格式：Highlights、摘要、正文、表、图注、参考文献）。
- `figures/`：图 1–4，由 `make_figures.py` 生成；图 5–7 直接用 `analysis/results/figures/fig_b1..b3`。
  重新出图：`MPLBACKEND=Agg python docs/manuscript/make_figures.py`（先运行 `validation/build_inputs.py` 和 `validation/run_validation.py`）。
- 所有数字都来自仓库里已提交的结果（`validation/README.md` 第 4 节、`analysis/README.md`），没有新跑分析，模型仍是 `4e09ac5`。

## 投稿前必须核实（正文中标为 [VERIFY]）

1. **剂量的摄入量**：Fan 定的是用国际标准（IAEA）。我没能在线读到 IAEA SRS-19 的表 XIV / XVIII（成人默认摄入量），
   正文暂时沿用分析里的 100 kg/年作比较用。需要拿到 SRS-19 原文后替换，并重算表 5 的年剂量（只是乘一个系数）。
2. **三篇需要读全文的前序工作**：Keum et al. 2005（Health Phys.）、Melintescu 2023（FST）、Galeriu et al. 1995（Fusion Technol.）。
   确认它们具体做了什么，引言和讨论里的对比措辞要据此调整。
3. **“A versatile model…”(2005, Radioprotection) 的作者**：不同来源写法不同（Barescut 等 / Melintescu & Galeriu），以期刊页面为准。
4. **参考文献的卷、页、DOI**：大部分来自检索结果的摘要页，没有逐条核对；Raskob 1993（UFOTRI 描述）、de Wit 2019、SALib 要补全。
5. **TECDOC-1991 公式编号**（OBT 生成用的“式 43”、叶–气交换公式）和 EMRAS TC14 假想情景的出处章节。
6. **Highlights 字数**：Elsevier 要求每条不超过 85 个字符。
7. **AI 使用声明**：按 Elsevier 现行政策核对措辞；作者、CRediT、基金、利益冲突由 Fan 和导师补。
8. **知网补查**：国内聚变堆氚事故食物链剂量方面的工作（本次检索覆盖不到）。

## 还缺的东西

- 表 1（基准汇总）、表 3（逐条结果）：可以从 `validation/benchmarks/` 和 `validation/results/` 直接生成，下一步做。
- 补充材料 S1（求解器对比）、S2（65 条基准全表、情景卡）。
- Zenodo 存档和 DOI，仓库加 LICENSE。
- 篇幅：正文约 6000 词，JER 一般不严格限字，投稿前再看指南。
