# 模型验证：公开文献数据 × Dmitte HTO 迁移模型

本目录用公开文献中的 HTO 暴露实验结果和公开气象/土壤/作物数据，对 `dmitte` 的空气–土壤–作物
氚迁移模型（`dmitte/transfer_equotions.py` + `dmitte/calc_para.py` + WOFOST）做外部验证。

```
validation/
├── benchmarks/literature_benchmarks.csv   # 整理的文献基准（21 条，含出处链接）
├── build_inputs.py                         # 从公开数据构建 ./data（气象、土壤、作物参数、农事）
├── run_validation.py                       # 虚拟实验：1 h 急性暴露（昼/夜）+ 慢性暴露
├── analyze.py                              # 与文献基准对比，输出表格和图
└── results/                                # 运行结果（CSV、图）
```

## 1 结论速览

| 结论 | 依据 |
|---|---|
| 数值求解本身可靠 | 本验证用的逐小时矩阵指数传播与仓库 `solve_con()` 最大相对差 2×10⁻⁴（马铃薯 default）、5×10⁻⁵（baomi）、8×10⁻⁷（UFOTRI） |
| 暴露后叶片 TFWT 的衰减幅度（baomi、UFOTRI 参数）落在文献范围 | 马铃薯：暴露结束到收获下降 660 倍（baomi）/ 5 400 倍（UFOTRI），文献（水稻）600–95 000 倍 |
| 慢性释放下块茎 OBT/空气 HTO：仅 UFOTRI 常数版接近文献 | 1.15（UFOTRI） vs 0.93±0.21；baomi 5.3，default 225 |
| **昼夜差异完全没有表现出来** | 夜间/白天叶片吸收比 1.0–1.7（文献 ≈0.25–0.4）；OBT 生成夜/昼比 1–4（文献 0.1–0.33） |
| **叶–气交换速率偏离文献几个数量级** | 白天 default 0.012 h⁻¹、baomi 1.2 h⁻¹（玉米 0.10–0.21）；夜间 10⁻⁸–10⁻⁶ h⁻¹（文献只比白天低 2–10 倍） |
| **块茎 TLI 偏高** | 0.2–0.3 %（文献）vs 1.2 %（baomi）、3.7 %（UFOTRI）、93 %（default） |
| **谷物（CEREAL）分支在 default/UFOTRI 下数值发散** | `kbh_fo` 恒为负、`kfh_bh` 达 6×10⁴ h⁻¹，库存发散到 10²⁸⁹，`solve_con()` 直接报错 |
| **OBT 库存会变负** | `TROBT ∝ (PGASS − PMRES)` 在衰老期为负；凡是用 `TROBT` 的分支（马铃薯 default/baomi、小麦全部三种）叶片 OBT 都会出现负值 |

总的判断：模型框架和求解器没问题。但植物交换和 OBT 生成这两块过程的参数化，和公开实验相比在机理上有系统偏差
（第 5 节）。剂量计算目前依赖的作物 OBT（`dose_adult_potato` 里的 `As[harvest_time, -1]`）在 default
参数下可能高估 1–2 个数量级，baomi 下约高估 4–5 倍，用于论文结论前建议先修正第 5 节的问题。

## 2 数据

### 2.1 文献基准（`benchmarks/literature_benchmarks.csv`）

本环境无法访问 IAEA、ScienceDirect、OSTI 等站点下载全文，因此**基准值取自论文摘要、会议摘要和二次引用**。
每条都在 `evidence` 列里注明了来源类型，`doi_or_url` 列给出链接，正式使用前建议对照全文核实。

| ID | 量 | 作物 / 条件 | 文献值 | 出处 |
|---|---|---|---|---|
| B1 | 暴露结束时叶片 TFWT / 空气水汽 HTO | 白菜、萝卜，1 h 室外箱暴露 | 10–50 % | Choi et al. 2005, JER 84:79 ([PubMed](https://pubmed.ncbi.nlm.nih.gov/15936121/)) |
| B2 / B3 | 同上，白天 / 夜间 | 水稻 | ≈100 % / 30–40 % | Choi et al.（[ETDEWEB 21051446](https://www.osti.gov/etdeweb/biblio/21051446)） |
| B4 | 夜间 / 白天叶片 HTO 吸收 | 小麦，实验室 | ≈1/4 | Diabaté & Strack 1997（[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0265931X97849855)） |
| B5 / B6 | 叶–气交换速率；昼/夜比 | 综述 | < 1 h⁻¹；2–10 | Galeriu et al. 2013, JER 118:40（[doi](https://doi.org/10.1016/j.jenvrad.2012.11.005)） |
| B7 / B8 | 叶片 HTO 吸收速率常数 白天 / 夜间 | 玉米，1 h 田间暴露 | 0.10–0.21 / 0.035–0.13 h⁻¹ | [PubMed 31581057](https://pubmed.ncbi.nlm.nih.gov/31581057/) |
| B9–B11 | TFWT 从暴露结束到收获的下降倍数 | 水稻 / 萝卜根 / 白菜叶 | 600–95 000 / ≤1.3×10⁴ / ≤1.1×10⁶ | Choi et al. 2002 JER 58:67；Choi et al. 2005 |
| B12 | TLI = 收获时可食部 OBT / 暴露结束时叶片 TFWT | 马铃薯块茎 | 0.2–0.3 % | EMRAS II WG7 报告（IAEA，[pdf](https://www-ns.iaea.org/downloads/rw/projects/emras/emras-two/third-technical-meeting/wgroup-seven/presentation-5th-wg7-obt-in-night-time.pdf)） |
| B13 / B14 | TLI | 萝卜 / 白菜 | 0.1–0.3 % / 0.1–0.5 % | Choi et al.（[ETDEWEB 20305115](https://www.osti.gov/etdeweb/biblio/20305115)） |
| B15 | 夜间 / 白天 OBT 生成（同等叶片 HTO） | 综述 | 0.1–0.33 | Galeriu et al. 2013 |
| B16 / B17 | 暴露结束 OBT / 空气 HTO；24 h 后 OBT 占植物总氚 | 菜豆 | 0.2 %；2–4 % | [ETDEWEB 410495](https://www.osti.gov/etdeweb/biblio/410495) |
| B18 / B19 | 慢性释放：果实/块茎 OBT / 空气 HTO；TFWT / 空气 HTO | 加拿大 CRL 菜园 2008–2011 | 0.93±0.21；1.20±0.63 | Korolevych & Kim 2013, JER 118（[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0265931X12002937)） |
| B20 | 叶片 OBT/HTO（慢性） | EMRAS 情景 | ≈0.7 | IAEA EMRAS 氚与 C-14 工作组报告（TECDOC-1678） |
| B21 | 土壤 HTO 再释放 | 田间 D₂O 示踪 | 持续数周（定性） | JER 2004（[ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0265931X03001693)） |

其他相关但未取到数值的资料：IAEA-TECDOC-1738（EMRAS II 氚事故释放）、IAEA-TECDOC-1991（MODARIA 氚模型比对，2022）、
EMRAS「Potato Scenario」最终报告（2008）。网络允许时，这些报告里的逐时序列数据是下一步最值得补充的。

### 2.2 模型输入（`build_inputs.py`，全部公开、固定到 commit）

| 输入 | 来源 |
|---|---|
| 气象 | 荷兰 Wageningen Haarweg 站 2004–2008 逐日数据（Wageningen 大学气象组），取自 [ajwdewit/pcse_notebooks](https://github.com/ajwdewit/pcse_notebooks) `data/meteo/nl1.xlsx` |
| 土壤 | `ec2.soil`（EC2-medium，与原研究同名同参数：SMW 0.099、SMFCF 0.272、SM0 0.39），同上仓库 |
| 作物 | WOFOST 7.2 作物参数 [ajwdewit/WOFOST_crop_parameters@wofost72](https://github.com/ajwdewit/WOFOST_crop_parameters/tree/wofost72)：`Potato_701`、`Winter_wheat_102` |
| 农事 | 生成：马铃薯 5 月 1 日出苗，冬小麦前一年 10 月 20 日播种 |
| `meteo_usefor_cttm.xlsx` 中的派生量 | RH = 水汽压 / 日均温饱和水汽压；气压取标准大气（站点海拔 7 m）；5/10/20 cm 土壤含水量 = 同年马铃薯 WOFOST 根区含水量 SM（无作物期取田间持水量 0.272） |

## 3 方法

1. **虚拟实验**，按文献实验流程设计：
   - *急性*：植物在空气水汽 HTO 浓度恒定的箱中暴露 1 h，然后回到洁净空气直到收获。整个生长季每周做一次，
     每次分 11:00（白天）和 23:00（夜间）两组。马铃薯 2004–2008 共 5 季、小麦 2005–2008 共 4 季。
   - *慢性*：整个生长季空气 HTO 恒定（对应 Korolevych & Kim 的连续释放菜园）。
2. **模型原样使用**：迁移率来自 `transfer_rates`（default，`iaea-case1-HTO.py` 在用）、`transfer_rates_baomi`
   （`potato_HTO.py`、`cereals_HTO.py` 在用）、`transfer_rates_UFORTI`（UFOTRI 常数）；ODE 右端就是
   `transfer_equotions()` 本身。系统对 y 线性、迁移率逐小时分段常数，所以每小时的传播子 = 生成矩阵的矩阵指数，
   生成矩阵从 `transfer_equotions()` 逐列取出。唯一的改动是实验施加的边界条件：空气库室在暴露小时内固定为暴露浓度，
   其余时间为 0（洁净空气）。
3. **昼夜强迫**：仓库的 `meteodata()` 把日数据线性插值到小时，没有昼夜变化。为检验昼夜基准，`run_validation.py`
   用日数据生成逐时强迫（温度在 TMIN/TMAX 间余弦变化、辐射按太阳高度角分配、RH 由逐时温度计算），通过替换
   `calc_para.meteodata` 注入。同时也用仓库原始的逐日插值强迫跑了一遍（结果见 CSV 中 `forcing=daily`），两者除昼夜指标外结论一致。
4. **浓度换算**：叶片（地上部）TFWT = `A_bh / plant_w`；块茎/籽粒 TFWT = `A_fh / friut_w`；OBT 以燃烧水计，
   = OBT 库存 / (干物质 × 0.556 L/kg)（淀粉/纤维素含氢 6.2 %）；空气水汽量 = `ML × a_h`。
5. **判定**：取模型中位数（全部年份 × 暴露日期，逐时强迫），落在文献范围内为 ✅ 符合；与范围相差 3 倍以内为 ⚠️ 接近；
   超过 3 倍为 ❌ 偏离；≤0 或发散（|log₁₀| > 30）为 ⛔ 非物理。只给中心值的基准按 ×/÷1.5 作为容差。

## 4 结果

完整表见 `results/benchmark_comparison.csv`（含 5–95 % 分位）。中位数与判定：

| 基准 | 指标 | 作物 | UFOTRI | baomi | default | 文献 |
|---|---|---|---|---|---|---|
| B1 | 暴露结束地上部 TFWT/空气，% | 马铃薯 | 1.3 ❌ | 4.8 ⚠️ | 6.8 ⚠️ | 10–50 |
|  |  | 小麦 | 25 ✅ | 52 ⚠️ | 117 ⚠️ | 10–50 |
| B2 | 同上（白天），% | 马铃薯 | 1.3 ❌ | 3.5 ❌ | 6.9 ❌ | ≈100 |
| B3 | 同上（夜间），% | 马铃薯 | 1.3 ❌ | 6.7 ❌ | 6.8 ❌ | 30–40 |
| B4 | 夜/昼吸收比 | 马铃薯 | 1.0 ⚠️ | 1.7 ❌ | 1.0 ⚠️ | ≈0.25 |
|  |  | 小麦 | 1.0 ⚠️ | 23 ❌ | 1.1 ⚠️ | ≈0.25 |
| B7 | 白天叶→气交换 kbh_a2，h⁻¹ | 马铃薯 | 0.35 ⚠️ | 1.2 ❌ | 0.012 ❌ | 0.10–0.21 |
|  |  | 小麦 | 0.35 ⚠️ | 24 ❌ | 0.24 ⚠️ | 0.10–0.21 |
| B8 | 夜间叶→气交换，h⁻¹ | 马铃薯 | 0.35 ⚠️ | 5×10⁻⁶ ❌ | 5×10⁻⁸ ❌ | 0.035–0.13 |
| B6 | 昼/夜交换比 | 马铃薯 | 1 ⚠️ | 3×10⁵ ❌ | 3×10⁵ ❌ | 2–10 |
| B9 | 叶片 TFWT 暴露→收获下降倍数 | 马铃薯 | 5 400 ✅ | 660 ✅ | 5.5 ❌ | 600–95 000 |
|  |  | 小麦 | 发散 ⛔ | 1 200 ✅ | 发散 ⛔ | 600–95 000 |
| B10 | 块茎 TFWT 峰值→收获下降倍数 | 马铃薯 | 83 ✅ | 12 ✅ | 1.4 ✅ | ≤1.3×10⁴ |
| B12 | TLI 块茎，% | 马铃薯 | 3.7 ❌ | 1.2 ❌ | 93 ❌ | 0.2–0.3 |
| B15 | 夜/昼 OBT 生成 | 马铃薯 | 1.0 ❌ | 4.0 ❌ | 1.0 ❌ | 0.1–0.33 |
| B16 | 暴露结束 OBT/空气，% | 马铃薯 | 0.006 ❌ | 0.17 ✅ | 0.32 ⚠️ | ≈0.2 |
| B17 | 24 h 后 OBT 占植物总氚，% | 马铃薯 | 7.2 ⚠️ | 46 ❌ | 27 ❌ | 2–4 |
| B18 | 慢性：块茎/籽粒 OBT / 空气 | 马铃薯 | 1.15 ⚠️ | 5.3 ❌ | 225 ❌ | 0.93±0.21 |
|  |  | 小麦 | 发散 ⛔ | 1.4 ⚠️ | 发散 ⛔ | 0.93±0.21 |
| B19 | 慢性：块茎/籽粒 TFWT / 空气 | 马铃薯 | 0.18 ❌ | 8.5 ❌ | 13 ❌ | 1.20±0.63 |
| B20 | 慢性：叶片 OBT/TFWT | 马铃薯 | 4.1 ❌ | 负值 ⛔ | 负值 ⛔ | ≈0.7 |

![叶片吸收](results/figures/fig1_leaf_uptake.png)
![交换速率](results/figures/fig2_exchange_rate.png)
![TLI](results/figures/fig3_TLI_potato.png)
![慢性](results/figures/fig4_chronic.png)

## 5 发现的问题（按影响排序）

1. **叶→气交换用的是蒸腾量，不是交换通量**（`transfer_equotions.py` 中 `kbh_a2 = plantfx.ETRM / plantfx.plant_w`）。
   吸收端是沉积速度 `VDPF = 1/(RAM+RB+RC1)` 乘空气水汽量，相当于 g·ρ_v·C_air；释放端却是蒸腾 E = g·(ρ_sat − ρ_v)。
   稳态时 C_leaf/C_air ≈ ρ_v/(ρ_sat − ρ_v) = RH/(1 − RH)，RH = 0.75 时为 3，而物理上应 ≤ RH/α（≈0.7）。
   慢性情景中地上部 TFWT/空气达到 5–65（default/baomi）就是这个原因。夜间 ETRM→0，交换几乎停止
   （10⁻⁸–10⁻⁶ h⁻¹），这和夜间只低 2–10 倍的观测不符，也导致 default 下叶片 HTO 长期滞留、TLI 被高估约 300 倍。
   建议把释放项改为 g·ρ_sat(T_leaf)/W_leaf（与 UFOTRI/CTEM 的做法一致），蒸腾只作为根系吸水通道。
2. **单位不一致**：`calc_ETRM` 把 `RC1 × 100` 当成 s/m 用（注释按 s/cm 理解），`calc_VDPF` 直接把 `RC1` 当 s/m 用。
   baomi 版本里的 `× 100` 看起来是在补偿这一点，但补偿后白天交换速率又高出 6 倍（马铃薯）到 100 倍（小麦）。
3. **气孔阻力没有昼夜变化**：`calc_st_res` 的条件 `... | (ISTRE == 1)` 在默认 `ISTRE=1` 时恒为真，
   RC1 恒等于 RC/LAI，所以吸收（`ka2_bh`）白天夜里一样，夜/昼吸收比 ≈1（文献 ≈0.25–0.4）。
4. **OBT 生成**：`calc_transOBT` 用逐日 `PGASS − PMRES` 插值到小时，没有昼夜变化（夜/昼比 ≈1，文献 0.1–0.33）。
   衰老期呼吸大于同化时为负，`kbh_bo < 0`，OBT 库存变负（凡用 `TROBT` 的分支都会出现，只有 UFOTRI 的 ROOT_VEG 常数分支不受影响）。建议至少 `np.maximum(…, 0)`，
   并用逐时辐射调制。
5. **CEREAL 分支**（`transfer_rates`、`transfer_rates_UFORTI`、`transfer_rates_Korea`）：
   - `kbh_fo = np.log(2/(HWZ/2)) * …`，log 的自变量 < 1，恒为负（ROOT_VEG 分支写的是 `np.log(2)/(HWZ/2)/24`，疑为笔误）；
   - `kbh_fh = np.log2(10)/2 = 1.66 h⁻¹`（疑为 `np.log(2)/…`），`kfh_bh = kbh_fh·plant_w/friut_w` 在灌浆初期达 6×10⁴ h⁻¹，
     籽粒为 0 时为 inf；
   - 结果：库存发散到 10²⁸⁹，`solve_con()` 中 RK45 失败（`AttributeError: 'list' object has no attribute 'squeeze'`）。
     `iaea-case1-HTO.py` 里的 `np.abs(k_array)` 掩盖了负号问题，但没有解决发散。
6. **ROOT_VEG default 分支**：`kfh_bh = kbh_fh·plant_w/friut_w` 在块茎形成前为 inf（每季 265 个非有限值）；
   `ks2_s3 = Va_b2/150` 在向上水流时为负，应改为反向迁移项。
7. **干物质比例**：`cereals_HTO.py` 用籽粒干物质比例 0.86 计算整株含水量（`plant_w = tagp·(1−dm)/dm`），
   整株含水被低估约 10 倍，地上部 TFWT 超过空气浓度（>100 %，非物理）。马铃薯的 `plant_w` 用 TAGP（含块茎）计算，
   把“叶片”HTO 稀释到了块茎水里，这是马铃薯 B1–B3 偏低的一部分原因。
8. **已修复的可复现性问题**（本次提交）：
   - PCSE 5.5 默认从 `WOFOST_crop_parameters` 的 `master` 分支下载作物参数，而该分支现已清空，新克隆无法运行。
     现改为优先读 `./data/crop/wofost72`，否则用 `wofost72` 分支。
   - 原代码依赖一份手工改过的 PCSE 配置（输出 `PGASS`、`PMRES`），标准 PCSE 不输出这两个变量，`Calc_plant` 会 KeyError。
     现在随仓库提供 `dmitte/conf/Wofost72_WLP_FD_dmitte.conf`，并通过 `Engine(config=...)` 使用。

## 6 局限

- 文献值来自摘要/二次引用，没有逐时序列；不同作物（白菜、萝卜、水稻、玉米、菜豆）的基准被用来约束马铃薯和小麦，
  只适合判断量级和定性行为（昼夜比、衰减幅度），不适合精细标定。
- 气象是荷兰 Wageningen，不是实验所在地（韩国大田、德国卡尔斯鲁厄、加拿大 Chalk River）；逐时强迫是由日数据合成的。
- 土壤三层含水量都用 WOFOST 根区 SM，没有分层观测。
- 大气扩散（高斯烟团）部分没有验证：可用的公开示踪实验（如 Prairie Grass）在本环境无法下载。

## 7 复现

```bash
pip install -r requirements-dev.txt          # Python 3.12
python validation/build_inputs.py           # 下载公开数据，生成 ./data（约 1 MB）
python validation/run_validation.py         # 约 1 分钟，写 validation/results/*.csv
python validation/analyze.py                # 对比表 + 图
```
