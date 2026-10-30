# 氚（HTO/OBT）作物迁移模型文献筛选清单

检索日期：2026-10-29/30｜数据：OpenAlex、Crossref（去重后 1640 条；原始 `results.csv`、`search_log.json` 在 Fan 电脑
`D:\Work\Workspace\【期刊投稿】Dmitte论文\literature_search\`）。汇总结论见 `novelty_review.md` 第 6 节。

## 检索说明与局限
- OpenAlex：23 个指定检索式 + 19 个定向检索式（TOCATTA、Sobol/Morris、MODARIA、CFETR/ITER 等），每式 ≤50 条；Crossref：23 式，每式 ≤50 条。
  Semantic Scholar 没有 API key，几乎全部 HTTP 429，只成功 1 个检索式。知网连不上。Scopus、WoS 需要订阅，未使用。
- 只看标题和摘要。不少 JER、FST 早期论文没有摘要，标“?”的必须看全文核实。
- 标注：✔ 是；◐ 部分/间接；✘ 摘要中未见；? 无法判断。
  (a) 耦合作物生长模型；(b) 量化释放日期/生育期对剂量或收获期 OBT 的影响；(c) 比较动态与恒定迁移率；(d) 全局敏感性分析（Sobol/Morris）；(e) 冻结测试集或系统的多数据集基准验证。

## 一、2015 年及以后  [a/b/c/d/e]
1. Melintescu, Galeriu, Diabaté, Strack (2015) Preparatory steps for a robust dynamic model for OBT dynamics in agricultural crops. *Fusion Sci. Technol.* doi:10.13182/fst14-t59。以干物质增长和光合驱动 OBT 的模型框架（CROPTRIT 前身）。[◐/?/✘/✘/✘]
2. Galeriu, Melintescu (2015) Progresses in tritium accident modelling in the frame of IAEA EMRAS II. *FST* doi:10.13182/fst14-t26。EMRAS II 进展与模型比对。[?/?/?/✘/◐]
3. Melintescu, Galeriu (2017) Uncertainty of current understanding regarding OBT formation in plants. *JER* doi:10.1016/j.jenvrad.2016.11.026。机理不确定性综述（定性）。[✘/?/?/✘/✘]
4. Galeriu, Melintescu (2017) Relevance of night production of OBT in crops. *FST* doi:10.1080/15361055.2016.1273705。含夜间过程的 OBT 概念模型。[◐/◐/✘/✘/✘]
5. Le Dizès, Aulagnier, Maro 等 (2017) The VATO project: dynamic transfer model of tritium in grassland ecosystem. *JER* doi:10.1016/j.jenvrad.2016.11.028。TOCATTA-χ 草地小时步长模型，La Hague 单站点验证。[◐/?/?/✘/✘]
6. Le Dizès, Maro, Rozet, Hébert (2015) Modeling and validating tritium transfer in a grassland ecosystem. *FST* doi:10.13182/fst14-t51。同上早期版本。[◐/?/?/✘/✘]
7. Renard H. (2017) Hourly modelisation of tritium transfers in a grassland ecosystem. 博士论文。[◐/?/?/?/✘]
8. Ota, Kwamena, Mihok, Korolevych (2017) Role of soil-to-leaf tritium transfer in controlling leaf tritium dynamics. *JER* doi:10.1016/j.jenvrad.2017.09.002。CTEM-CLASS-TT 与 SOLVEG-II 对比。[◐/✘/✘/✘/◐]
9. Yokoyama 等 (2019) Development of field estimation technique and improvement of environmental tritium behavior model. *Plasma Fusion Res.* doi:10.1585/pfr.14.3405099。[?/?/?/✘/✘]
10. Melintescu (2023) Radiological impact assessment of acute tritium releases — a soil dynamic model. *FST* doi:10.1080/15361055.2023.2223723。CROPTRIT 升级（土壤子模型）。[✔?/?/?/✘/✘]
11. Melintescu, Patryl (2025) Interception and uptake by plant leaves of tritium from precipitation. *JER* doi:10.1016/j.jenvrad.2025.107677。[?/?/?/?/?]
12. Satoh, Tani, Nagai 等 (2025) Necessity of considering non-photosynthetic metabolism to reproduce the behavior of OBT in a leafy vegetable. *J. Hazard. Mater.* doi:10.1016/j.jhazmat.2025.139518。[?/◐/?/?/?]
13. Priebe, Rowan, Saurette, Melintescu 等 (2026) The enhanced Craig-Gordon deuterium model as a tool for estimating free and organically bound tritium in plants. *JER*（预印本 doi:10.2139/ssrn.6377409）。机理模型与恒定转移因子比较。[✘/✘/◐/✘/✘]
14. Murphy J. (2024) Evaluating postulated tritium releases on a large DOE site using UFOTRI. OSTI。气象抽样的概率剂量，单因素敏感性。[◐/◐/✘/✘/✘]
15. Clavier 等 (2023) Comparison of tritium dose calculations from MACCS, UFOTRI, and ETMOD. doi:10.2172/2293570；Garcia (2023) doi:10.2172/2463045。[✘/?/◐/✘/◐]
16. Park, Kim (2024) Evaluation of tritium ingestion dose in the Republic of Korea using methodologies from three organizations. *Nucl. Eng. Technol.* doi:10.1016/j.net.2024.11.025。稳态方法比较。[✘/✘/✘/✘/✘]
17. Raskob, Hasemann, Schichtel, Ottenburger (2026) Dose to the population from unit releases of tritium and activation products in EU DEMO. KITopen doi:10.5445/ir/1000195173。[?/?/?/?/✘]
18. Cui, Li, Cao, Cao (2025) Dose assessment of gaseous tritium release from D-T fusion devices based on steady-state plume model. *Fusion Eng. Des.* doi:10.1016/j.fusengdes.2025.114965。[✘/✘/✘/?/✘]
19. 福田ほのか 等 (2025) UFOTRI 用于六所村 JA-DEMO 候选厂址。IRDB。[◐/?/✘/✘/✘]
20. Wu, Liu, Kearfott, Sun (2020) Evaluation of public dose from FHR tritium release with consideration of meteorological uncertainties. *Sci. Total Environ.* doi:10.1016/j.scitotenv.2019.136085。[✘/◐/✘/?/✘]
21. Lebel, Brown, Vucicevic, Korolevych (2022) Emergency response implications of tritium in CANDU severe accidents. *Ann. Nucl. Energy*。[✘/✘/✘/✘/✘]
22. Nie, Chen, Wang 等 (2023) 多形态氚环境行为与辐射评价体系研究进展. 中国科学: 技术科学 doi:10.1360/sst-2022-0356。中文综述。[✘/✘/✘/✘/✘]
23. Blanchard A. (2021) A posteriori verification and validation of a tritium dispersion and consequence model (UFOTRI). OSTI。[◐/?/✘/✘/◐]
24. Baburajan 等 (2025) TFWT and OBT in terrestrial plants at Tarapur … *J. Radioanal. Nucl. Chem.*。实测与稳态模型比较。[✘/✘/✘/✘/✘]
25. Polivkina 等 (2023–2024) 辣椒、茄子、生菜、番茄、豆类在舱室和田间不同生育期短时 HTO 暴露（NNC RK Bull.、Radiats. Gig.、PLoS ONE）。可作验证数据。[✘/◐/✘/✘/—]
26. Shen, Liu (2016) Soybean HTO day/night. *Nucl. Sci. Tech.*；Shen, Guo (2019) Maize HTO day/night. *Appl. Radiat. Isot.*。[✘/◐/✘/✘/—]

## 二、2015 年以前的关键前序工作
27. ★ Keum, Lee, Kang, Jun, Choi, Lee (2005) Prediction of tritium level in agricultural plants after short term exposure to HTO vapor and its comparison with experimental results. *Health Phys.* doi:10.1097/01.hp.0000174811.24386.cc。随时间变化的生长方程代替恒定产量，不同生育期暴露的预测与验证。[✔ 简化生长方程/✔/◐ 与恒定产量模型对比/✘/◐]
28. Choi 等 (2005) Chinese cabbage and radish acutely exposed at different growth stages. *JER* doi:10.1016/j.jenvrad.2005.04.004；Choi 等 (2002) rice. *JER*。实验。[✘/✔/✘/✘/—]
29. ★ Melintescu, Galeriu (2005) A versatile model for tritium transfer from atmosphere to plant and soil. *Radioprotection* doi:10.1051/radiopro:2005s1-064。WOFOST 光合驱动冠层阻力和 OBT。[✔/◐/?/✘/✘]
30. Melintescu, Galeriu, Marica (2002) Using WOFOST crop model for data base derivation of tritium and terrestrial food chain modules in RODOS. *Radioprotection* doi:10.1051/radiopro/2002154。[✔/◐/✘/✘/✘]
31. ★ Galeriu, Davis, Chouhan, Raskob (1995) Uncertainty and sensitivity analysis for the environmental tritium code UFOTRI. *Fusion Technol.* doi:10.13182/fst95-a30511。不确定性 < 1 个数量级，主导参数 LAI、最小气孔阻力、产量、初始土壤水分。[◐/✘/✘/◐ 抽样法/✘]
32. Raskob, Barry (1997) Importance and variability in processes relevant to environmental tritium ingestion dose models. *JER* doi:10.1016/s0265-931x(96)00090-2。[?/?/?/?/✘]
33. Barry 等 (1999) Intercomparison of model predictions of tritium concentrations in soil and foods following acute airborne HTO exposure. *JER* doi:10.1016/s0265-931x(98)00054-x。[?/?/?/✘/◐]
34. Galeriu, Melintescu, Strack, Atarashi-Andoh, Kim (2013) An overview of OBT experiments in plants following a short atmospheric HTO exposure. *JER* doi:10.1016/j.jenvrad.2012.11.005。实验汇编。[—/◐/—/—/◐]
35. Strack, Diabaté, Müller, Raskob (1995) OBT formation and translocation in crop plants. *Fusion Technol.*；Diabaté, Strack (1997) *JER*。Plant-OBT 模型。[◐/◐/✘/✘/✘]
36. Balashov, Golubev 等 (2011) Modeling the build-up of OBT in crops after acute tritium exposure in greenhouse's environment. *FST*。[◐/◐/✘/✘/✘]
37. Ota, Nagai (2011) SOLVEG-II. *JER*；Ota, Nagai (2011) nighttime wet deposition. *Radioprotection*。[◐/◐/✘/✘/✘]
38. Korolevych, Kim (2011) OBT in tomato and potato (CTEM+CLASS). *FST*；Korolevych, Kim, Davis (2014) OBT/HTO ratio in produce. *JER*。[✔ 物候/✘/✘/✘/◐]
39. Raskob (2002) UFOTRI rice. *FST*；Raskob, Velarde, Perlado (2005) ITER Vandellós. *FST*。全年气象概率剂量。[◐/◐/✘/✘/✘]
40. Galeriu, Melintescu 等 (2009) Accidental release of tritiated water — toward a better radiological assessment. *Radioprotection*。[◐/?/✘/✘/✘]

未检出、需人工补查：IAEA TECDOC-1738、MODARIA 报告（已在本文基准中使用）；TOCATTA-χ 模型描述论文；知网中 CFETR/聚变堆氚事故食物链剂量文献。

## 三、与本文三项主张的重叠
- **(1) 释放日期 × 动态/常数分解**：没有找到在同一模型内冻结时间变化、并定量给出约 100 倍和最坏 5–40 倍低估的论文。
  最接近：Keum 2005（#27，生长方程 + 生育期 + 与恒定产量对比）、Melintescu/Galeriu 的 WOFOST/CROPTRIT 系列（#1、#10、#29、#30）。须正面引用并说明差异。
- **(2) Sobol：参数不确定性 ≪ 释放时间**：没有找到对氚作物模型做 Sobol/Morris 的论文。必须引用 Galeriu 等 1995（#31）的 UFOTRI 不确定性分析。
- **(3) 65 条基准 + 冻结测试集 + 开源**：没有找到。已有验证多为单站点、单组实验或多模型情景比对；Galeriu 2013（#34）的汇编可能与本文基准来源重叠。
