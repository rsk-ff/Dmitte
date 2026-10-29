# 新颖性调研：Dmitte 论文与已有工作的关系

日期：2026-10-29。目的：在起草正文之前，确认本文与已有工作是否重复，以及能主张哪些新意。
方法：网络检索（英文公开来源、出版社摘要页、IAEA 报告、项目报告），检索词涵盖作物动态氚/OBT 模型、生育期依赖、释放时间、
敏感性分析，以及 UFOTRI、SOLVEG、CERES、CROPTRIT、EMRAS/MODARIA、Galeriu、Melintescu、Strack、Diabaté、KAERI 等。
**局限**：我无法访问 Scopus、Web of Science 和 CNKI，也打不开部分付费全文；第 5 节列出需要 Fan 通过学校图书馆补查的内容。

---

## 1 结论

- **“用作物生长模型驱动氚迁移 / OBT 生成”本身不是新的。** 至少有四条已有工作这样做：
  Barescut、Melintescu、Galeriu 等 2005 年用 WOFOST 的光合作用驱动 OBT 生成；罗马尼亚 IFIN-HH 的 CROPTRIT（冬小麦，生长模块改自 WOFOST/ORYZA）；
  日本 JAEA 的 SOLVEG（光合驱动的碳水化合物子模型）；韩国 KAERI 的动态库室模型（含随时间变化的植物生长方程，按不同生育期暴露做过验证）。
  规划里原来的题目和论点 1 需要改，不能把“耦合 WOFOST”当卖点。
- **没有找到与本文核心结果重复的工作。** 下面三点在检索到的文献中都没有出现，可以作为主要新意：
  1. 在**同一模型结构**内，把时间变化逐项冻结（只保留日变化 / 整季常数 / 冻结在释放时刻），定量分解“生长 vs 日变化”对收获期剂量的贡献；
  2. 用全局敏感性分析（Morris + Sobol）说明参数不确定性造成的范围（1.2–2 倍）远小于释放日期造成的范围（约 100 倍）；
  3. 用公开文献整理的 65 条基准、按实验活动划分并冻结测试集的验证流程，代码和基准全部公开。
- **撞车风险：中等，集中在 CROPTRIT。** 它与本文最接近（作物生长驱动、冬小麦、事故释放剂量、报告中提到季节和昼夜释放的差异），
  目前只查到项目报告和会议摘要，没有查到期刊论文。如果他们先发表释放时间分析，本文第 1 点的新意会被削弱，但第 2、3 点不受影响。
  另外，这个圈子很小，审稿人很可能就是 Galeriu、Melintescu、Raskob、Strack 这些人，必须充分引用他们的工作。

---

## 2 最接近的已有工作

| 工作 | 做了什么 | 与本文的关系 |
|---|---|---|
| Barescut, Gariel, Péres, Melintescu, Galeriu (2005), *Radioprotection* 40(S1): S437–S442，“A versatile model for tritium transfer from atmosphere to plant and soil” | 过程模型；冠层阻力用光合驱动的 Jarvis-Calvet 方法；**WOFOST 提供光合速率，用于冠层阻力和 OBT 生成**；小麦预测在 2 倍以内，另测水稻、大豆 | 直接说明“WOFOST 驱动 OBT”已有先例。本文差别：开源、逐小时矩阵指数求解、植物水平衡按构造闭合、马铃薯；系统的基准验证；释放日期与基线分解 |
| CROPTRIT（IFIN-HH，Melintescu 等；OPTRAS 项目报告 2017–2019；Tritium 2022 会议摘要） | 运行型动态模型，Cernavoda 电站；冬小麦为主，另有牧草、玉米、葡萄；生长模块改自 WOFOST/ORYZA，干物质和 OBT 逐日计算、按器官分配；加入夜间 OBT；报告称做过不确定性和敏感性分析（未见定量结果）；白天晴空释放剂量秋季最高，夜间或降雨释放高约一个数量级 | **最接近，撞车风险主要在这里。** 尚未查到期刊论文。本文差别：公开代码与基准、冻结测试集、两种作物、在同一结构内分解动态 vs 常数、Sobol 方差分解 |
| SOLVEG（Ota & Nagai 2011, *J. Environ. Radioact.* 102: 813–823，及 JAEA 系列） | 研究型、复杂的土壤–植被–大气模型；叶片中光合驱动的碳水化合物子模型生成 OBT，生成速率随植物类型和生育期变化 | 机制更细，但计算量大、参数多，不适合大量情景和全局敏感性分析。本文定位为“简单到能做全局 SA、又保留生长驱动”的模型 |
| KAERI 动态库室模型（Keum, Lee, Kang, Jun, Choi, Lee；韩国核学会会议论文；相关期刊文章见 PubMed 16340607，约 2006 年，需核实刊名） | 大气–土壤–植物库室，叶菜、根菜、谷物/块茎分 HTO 和 OBT；**含随时间变化的植物生长方程，结果依赖暴露时的生育期**；用水稻、大豆、白菜、萝卜在不同生育期暴露的实验验证 | “预测结果依赖生育期”已有先例，要引用。本文差别：用通用作物模型而不是经验生长方程；马铃薯/小麦；定量到剂量并与常数基线比较；全局 SA |
| Choi 等 (2005), *J. Environ. Radioact.*（S0265931X05001396），白菜和萝卜在不同生育期急性暴露 HTO | 实验：不同生育期暴露的氚水平 | 实验证据，可作为“生育期依赖”的支持文献 |
| Diabaté & Strack (1997), *J. Environ. Radioact.* 36: 157–175；Strack, Diabaté, Raskob (2005), *Fusion Sci. Technol.* 48: 767–770；Strack 等 (2011) | 小麦舱室/田间箱实验，籽粒 TLI 随灌浆阶段变化；Strack 2011 给出由叶、穗 HTO 积分推算籽粒 OBT 的经验关系（夜间权重 0.2） | 本文的测试集数据（B53–B55）就来自这里。实验层面早已发现灌浆期最敏感，本文是在模型中复现并推广到整季和剂量 |
| UFOTRI / PLANT-OBT（Raskob，KIT） | 运行型模型；OBT 依赖叶片 HTO、植物类型和发育期，夜间生成单独标定；PLANT-OBT 需要大量按作物标定的参数 | 发育期通过参数体现，不由作物模型驱动 |
| Melintescu & Galeriu，“Uncertainty of current understanding regarding OBT formation in plants”（arXiv:1609.05052）；Galeriu 等 (2013) *J. Environ. Radioact.* 118: 40–56（短时暴露 OBT 实验综述）；Melintescu 等 (2015) *Fusion Sci. Technol.* 67: 479–482；Galeriu & Melintescu (2017) *Fusion Sci. Technol.* 71: 595–599（夜间 OBT） | 综述 EMRAS II 小麦情景中 SOLVEG、CERES、FDMH+、UFOTRI 的表现：收获时叶片 OBT 没有模型满足统计标准；“现有动态模型的不确定性不允许定量评估”；呼吁更多跨作物、跨环境的检验；没有对动态模型做正式敏感性分析 | **这正是本文的切入点**：系统的多数据集验证 + 动态模型的全局敏感性分析。引言可以直接引用他们提出的缺口 |
| IAEA-TECDOC-1738 (2014, EMRAS II)、IAEA-TECDOC-1991 (2022, MODARIA)、EMRAS TC14 | 模型比对与情景；TECDOC-1991 式 43 是本文 OBT 生成的基础 | 本文的基准、判定标准（2 倍）、假想情景都来自这里 |

---

## 3 本文可以主张的新意（按强弱排列）

1. **释放时间依赖性的定量分解（主卖点）。** 已有工作只是定性说明“灌浆期最敏感”“夜间释放更高”，或者用一个模型报告了剂量。
   本文在同一模型结构内构造三种非动态基线，结果是：同样的 1 h 释放，收获期剂量随释放日期相差约 100 倍（5–95 %），
   常数迁移率只有约 2 倍；最坏情况常数模型低估 5 倍（小麦）到 40 倍（马铃薯），籽粒/块茎形成前高估约一个数量级；
   差别几乎全部来自作物生长，日变化只贡献 1.0–1.3 倍。这为“运行型模型能否用常数迁移率”给出了直接的量化答案。
2. **参数不确定性与释放时间的对比（Sobol）。** 参数不确定性只造成 1.2–2 倍的范围；在校准集上标定的参数对剂量几乎没有影响（ST ≤ 0.05）；
   小麦由同位素歧视因子 ID_OBT 主导，马铃薯由最小气孔阻力主导。Melintescu & Galeriu 明确指出动态模型缺少这类定量分析。
3. **验证流程。** 65 条公开文献基准、实验元数据、情景卡、按实验活动划分并冻结的测试集（只运行一次：13 项中 11 项在 2 倍以内），
   加上留一交叉验证和测试后补充分析。这在氚模型文献中很少见，而且全部公开可复现。
4. **开源、可移植的实现。** PCSE/WOFOST 有现成的作物参数库，换作物或地点只换参数文件；逐小时矩阵指数精确求解。单独作为新意偏弱，作为第 1–3 点的支撑。

不能主张的：耦合作物模型本身；“首次”考虑生育期依赖；夜间 OBT 生成（已有专门工作，而且本文这部分恰是已知局限）。

---

## 4 对规划的修改建议

- **题目**：去掉以“coupling a crop growth model”为中心的说法，改为突出释放时间。例如：
  *Release timing dominates the ingestion dose from acute HTO releases to potato and wheat: a growth-driven compartment model with benchmark validation and global sensitivity analysis*
- **引言**：按“已有动态模型（Barescut 2005、CROPTRIT、SOLVEG、KAERI、UFOTRI）→ EMRAS II 指出的缺口（验证不足、不确定性无法定量）→ 本文三个问题”的顺序写，
  并在第一段就承认作物驱动的思路已有先例。
- **讨论**：单独一段与 CROPTRIT 和 KAERI 模型比较；与 Strack 2011 的经验关系比较（本文能否复现其夜间权重 0.2）。
- **时间**：CROPTRIT 有可能先发表类似的释放时间分析，宜尽快投稿。

---

## 5 需要补查的内容（Fan 通过学校图书馆）

我检索不到的部分，建议投稿前补一次，重点看 2015 年以后：

1. **Scopus / Web of Science**，检索式示例：
   `(tritium OR HTO OR "organically bound tritium") AND (crop OR wheat OR potato OR rice) AND (model*) AND ("growth stage" OR "release time" OR "crop growth" OR WOFOST OR "sensitivity analysis")`
   重点期刊：*J. Environ. Radioact.*、*Fusion Sci. Technol.*、*Fusion Eng. Des.*、*Health Phys.*、*Radioprotection*、*J. Radioanal. Nucl. Chem.*、*Sci. Total Environ.*
2. **CNKI**（中文）：“氚 作物 OBT 模型”“有机结合氚 小麦 动态”“聚变堆 氚 事故释放 农作物 剂量”。国内聚变方向（如合肥的聚变堆安全团队）有事故氚释放的动态评价工作，
   检索到的英文条目“Dynamic evaluation of environmental impact due to tritium accidental release from the fusion reactor”没能打开全文，需要确认它是否涉及作物生长驱动。
3. **CROPTRIT 是否已有期刊论文**：查 Melintescu A. 2019 年以后的发表记录。
4. **KAERI 模型的期刊文章**：核实 PubMed 16340607 的刊名、年份和具体做法。
5. **IRSN 的 TOCATTA / VATO** 动态模型：检索到的是草地生态系统，似乎不涉及大田作物，建议确认。

---

## 来源

- Barescut et al. 2005（摘要页）：https://www.cambridge.org/core/journals/radioprotection/article/versatile-model-for-tritium-transfer-from-atmosphere-to-plant-and-soil/18013F5628F81E3EDFBFC2658F9DCE2C
- CROPTRIT（OPTRAS 项目报告）：https://www.nipne.ro/proiecte/pn3/optras/results1.pdf ，https://www.nipne.ro/proiecte/pn3/optras/results2019.pdf
- Tritium 2022 会议摘要集：https://www.tritium2022.ro/wp-content/uploads/2022/10/Book-Of-Abstracts.pdf
- Melintescu & Galeriu，OBT 形成不确定性综述：https://arxiv.org/pdf/1609.05052
- Melintescu et al. 2015，*Fusion Sci. Technol.* 67：https://wx1.ans.org/pubs/journals/fst/a_36954
- Galeriu & Melintescu 2017，*Fusion Sci. Technol.* 71：https://wx1.ans.org/pubs/journals/fst/a_40470
- KAERI 动态模型（会议论文）：https://www.kns.org/files/pre_paper/18/8금동권.pdf ；期刊条目：https://pubmed.ncbi.nlm.nih.gov/16340607
- Choi et al. 2005，*J. Environ. Radioact.*：https://www.sciencedirect.com/science/article/abs/pii/S0265931X05001396
- IAEA TECDOC-1991：https://iaea.org/publications/15015/harmonization-and-intercomparison-of-models-for-accidental-tritium-releases-to-the-atmosphere
- IAEA TECDOC-1738：https://iaea.org/publications/10617/transfer-of-tritium-in-the-environment-after-accidental-releases-from-nuclear-facilities
