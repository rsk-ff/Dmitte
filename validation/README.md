# 模型验证：公开文献数据 × Dmitte HTO 迁移模型

本目录用公开文献中的 HTO 暴露实验结果和公开气象/土壤/作物数据，对 `dmitte` 的空气–土壤–作物
氚迁移模型（`dmitte/transfer_equotions.py` + `dmitte/calc_para.py` + WOFOST）做外部验证。

```
validation/
├── benchmarks/                             # 文献基准 v2：46 条标量基准 + 17 个实验数据集（见 benchmarks/README.md）
├── build_inputs.py                         # 从公开数据构建 ./data（气象、土壤、作物参数、农事）
├── run_validation.py                       # 虚拟实验：1 h 急性暴露（昼/夜）+ 慢性暴露
├── analyze.py                              # 与文献基准对比，输出表格和图
├── exp_attribution.py                      # 归因实验：马铃薯慢性指标偏高的原因（第 5.3 节）
└── results/                                # 运行结果（CSV、图）
```

## 1 结论速览

本节和第 4 节是**基准 v3**（第 2.1 节）在当前模型上的结果，**只含校准集**。测试集（15 条）已冻结，
在模型修改完成前不运行、不查看（第 3 节第 6 条）。v3 改为“本作物数据优先”，所以和之前版本的判定不能直接比较；
之前的结果见本文件的 git 历史。

| 结论 | 依据（校准集） |
|---|---|
| 迁移率全部有限且非负，没有非物理结果 | 三个变体、两种作物的负值和非有限值都为 0；判定中没有 ⛔ |
| 小麦白天吸收和白天交换速率在文献范围内 | 1 h 暴露结束叶片 TFWT/空气 58 %（田间箱 54–92 %，B49）；白天交换 1.3 h⁻¹（0.85–1.66，B5） |
| **小麦夜间吸收和夜间交换都偏高** | 夜间叶片 TFWT/空气 67 %（田间箱 17–19 %，B50）；夜间交换 0.91 h⁻¹（0.24–0.38，B29）；夜/昼吸收比 1.2（0.25，B4）。气孔阻力没有昼夜变化（review #15） |
| **马铃薯叶–气交换偏慢约 8 倍** | 白天 0.064 h⁻¹，马铃薯温室拟合值 0.5–1.1 h⁻¹（B56）；温室和舱室强制通风，田间速率应偏低，所以实际差距小于 8 倍 |
| **叶片 OBT 生成偏低，籽粒 TLI 偏低** | 小麦白天暴露结束叶片 OBT/TFWT 0.28 %（0.91–1.59 %，B22）；籽粒 TLI 0.065 %（0.18–0.32 %，B26） |
| **块茎 TLI 和慢性块茎 OBT 偏高** | 根茎代理（萝卜）TLI 6.6 %（0.11–0.33 %，B13）；慢性块茎 OBT/空气 16（0.72–1.14，B18），原因见第 5.3 节 |
| 均匀比活度检验不通过 | 植物水净流入/流出 +0.5～+0.9，baomi 果实水 +21～+165（第 4.1 节） |

总的判断：叶–气交换的白天行为对小麦已经合理，昼夜差、OBT 生成和向可食部的转运仍有系统偏差。
需要修的过程依次是：植物水平衡（根系吸水 = 蒸腾，review #10）、气孔昼夜变化（#15）、OBT 按干物质增长生成（#17）。

## 2 数据

### 2.1 文献基准（`benchmarks/`，v3）

详见 [benchmarks/README.md](benchmarks/README.md)。v3 在 v2 的基础上补读了 5 篇原文
（Diabaté & Strack 1997、Strack et al. 2005、Korolevych & Kim 2013、Shen & Guo 2019、Balashov et al. 2011）：

- `literature_benchmarks.csv`：65 条标量基准，每条带实验类型（舱室/田间箱/温室/田间）、示踪剂、暴露时长、C3/C4 和数据性质（实测/拟合/汇编/参数）；
- `datasets/`：19 个逐次实验数据集，新增小麦舱室逐时序列（空气 HTO、温度、湿度、光照和各器官 TFWT/OBT）和玉米 1 h 暴露；
- `scenario_cards.csv`：15 张情景卡，列出每组实验的暴露条件和缺失项，供用实测条件驱动模型；
- `split.csv`：校准集 50 条、测试集 15 条（小麦舱室、马铃薯舱室、Pickering 慢性）。

马铃薯的直接氚数据只有 FZK 舱室实验（TECDOC-1738 表 24，在测试集）和 Balashov 温室拟合速率常数（在校准集）。

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
   - *急性*：植物在空气水汽 HTO 浓度恒定的箱中暴露 1 h 或 2 h（田间箱实验为 1 h，舱室实验为 2 h），然后回到洁净空气直到收获。
     整个生长季每周做一次，每次分 11:00（白天）和 23:00（夜间）两组。马铃薯 2004–2008 共 5 季、小麦 2005–2008 共 4 季。
     每次暴露记录开花后天数和所处的灌浆阶段（WOFOST 贮藏器官日增量不低于季节最大值一半的时段为线性灌浆期），
     用于对应 Diabaté & Strack (1997) 按灌浆阶段给出的 TLI。
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
5. **判定**：取模型中位数（全部年份 × 暴露日期，逐时强迫），落在文献范围内为 ✅ 符合；与范围相差 2 倍以内为 ⚠️ 接近
   （TECDOC-1738 §7.9：模型与实验之间的差异至少为 2 倍，实验误差至少 50 %）；超过 2 倍为 ❌ 偏离；
   ≤0 或发散（|log₁₀| > 30）为 ⛔ 非物理。只给中心值的基准按 ×/÷1.5 作为容差。
   每项比较按同作物、同暴露时长对应；没有本作物数据时用代理作物，并在结果中标明。
6. **校准集与测试集**：`benchmarks/split.csv` 把 15 条基准和 4 个数据集划为测试集，按实验活动整体划分：
   FZK 小麦舱室实验（Diabaté & Strack 1997）、FZK 马铃薯舱室实验（TECDOC-1738 表 24）、Pickering 慢性释放。
   除 B12 外它们此前都没有参与过模型对比（B12 在 v1–v2 报告中出现过，属于已知泄漏）。
   `analyze.py` 默认只比较和绘制校准集；测试集的指纹（sha256）记录在 `analyze.py` 中，`tests/test_smoke.py` 检查它没有被改动。
   模型修改全部完成后，用 `python validation/analyze.py --include-test` 运行一次测试集。

## 4 结果

完整表见 `results/benchmark_comparison.csv`（含 5–95 % 分位）。下表为校准集各项的模型中位数与判定。
“匹配”列：本作物 = 同一作物的实验；代理 = 其他作物（马铃薯用萝卜、小麦用水稻）；通用 = 多作物汇总值。

| 基准 | 作物 | 文献指标 | 实验作物与类型 | 匹配 | UFOTRI | baomi | default | 文献 |
|---|---|---|---|---|---|---|---|---|
| B1 | 马铃薯 | leaf TFWT at end of 1 h exposure / mean air-moisture HTO | radish；field box | 代理 | 1.3 ❌ | 6.5 ❌ | 6.5 ❌ | 22–46 % |
| B10 | 马铃薯 | TFWT decrease factor from end of exposure to harvest | radish root；field box | 代理 | 83 ✅ | 6.8 ✅ | 8.7 ✅ | 0–1.3×10⁴ - |
| B13 | 马铃薯 | TLI = OBT in root at harvest / leaf TFWT at end of exposure | radish；field box | 代理 | 3.7 ❌ | 7.1 ❌ | 6.6 ❌ | 0.11–0.33 % |
| B56 | 马铃薯 | leaf TFWT fast loss rate constant k1 after daytime exposure | potato；greenhouse | 本作物 | 0.35 ⚠️ | 0.064 ❌ | 0.064 ❌ | 0.5–1.1 h-1 |
| B57 | 马铃薯 | leaf TFWT fast loss rate constant k1 after night exposure | potato；greenhouse | 本作物 | 0.35 ✅ | 0.041 ❌ | 0.041 ❌ | 0.19–0.42 h-1 |
| B58 | 马铃薯 | leaf HTO uptake rate constant (day) | potato；greenhouse | 本作物 | 0.35 ❌ | 0.064 ❌ | 0.064 ❌ | 0.93–2.1 h-1 |
| B6 | 马铃薯 | day/night ratio of leaf HTO uptake rate | most crops (not young lettuce)；review | 通用 | 1 ❌ | 1.5 ❌ | 1.5 ❌ | 3–5 - |
| B15 | 马铃薯 | night/day ratio of OBT formation at equal leaf HTO | wheat, rice, potato, bean, Chlorella；review | 通用 | 1 ⚠️ | 1.1 ❌ | 1.1 ❌ | 0.1–0.5 - |
| B18 | 马铃薯 | fruit/tuber OBT (Bq/L comb. water) / air HTO (Bq/L) averaged at pod filling | tomato, pepper, potato (pooled)；field | 本作物 | 1.1 ✅ | 16 ❌ | 16 ❌ | 0.72–1.1 - |
| B19 | 马铃薯 | fruit/tuber HTO / air HTO averaged 3 d before harvest | tomato, pepper, potato (pooled)；field | 本作物 | 0.18 ❌ | 12 ❌ | 6.8 ❌ | 0.57–1.8 - |
| B20 | 马铃薯 | plant OBT / plant TFWT at equilibrium (isotopic discrimination factor) | various；review | 通用 | 4.1 ❌ | 12 ❌ | 11 ❌ | 0.64–1.3 - |
| B4 | 小麦 | night/day ratio of leaf HTO uptake | wheat；field box | 本作物 | 0.99 ❌ | 1.2 ❌ | 1.2 ❌ | 0.17–0.38 - |
| B49 | 小麦 | leaf TFWT / air-moisture HTO at end of 1 h exposure (day) | wheat；field box | 本作物 | 25 ❌ | 59 ✅ | 58 ✅ | 54–92 % |
| B50 | 小麦 | leaf TFWT / air-moisture HTO at end of 1 h exposure (night) | wheat；field box | 本作物 | 25 ⚠️ | 70 ❌ | 67 ❌ | 17–19 % |
| B5 | 小麦 | leaf TFWT loss rate constant within 1 h after daytime exposure | wheat；field box | 本作物 | 0.35 ❌ | 1.3 ✅ | 1.3 ✅ | 0.85–1.7 h-1 |
| B29 | 小麦 | leaf TFWT loss rate constant within 1 h after night exposure | wheat；field box | 本作物 | 0.35 ✅ | 0.91 ❌ | 0.91 ❌ | 0.24–0.38 h-1 |
| B59 | 小麦 | leaf TFWT fast loss rate constant k1 after daytime exposure | wheat；greenhouse | 本作物 | 0.35 ⚠️ | 1.3 ⚠️ | 1.3 ⚠️ | 0.55–1.2 h-1 |
| B60 | 小麦 | leaf TFWT fast loss rate constant k1 after night exposure | wheat；greenhouse | 本作物 | 0.35 ⚠️ | 0.91 ❌ | 0.91 ❌ | 0.1–0.22 h-1 |
| B6 | 小麦 | day/night ratio of leaf HTO uptake rate | most crops (not young lettuce)；review | 通用 | 1 ❌ | 1.5 ❌ | 1.5 ❌ | 3–5 - |
| B22 | 小麦 | leaf OBT (comb. water) / leaf TFWT at end of 1 h exposure (day) | wheat；field box | 本作物 | 0.23 ❌ | 0.28 ❌ | 0.28 ❌ | 0.91–1.6 % |
| B23 | 小麦 | leaf OBT (comb. water) / leaf TFWT at end of 1 h exposure (night) | wheat；field box | 本作物 | 0.24 ⚠️ | 0.26 ⚠️ | 0.27 ⚠️ | 0.33–0.43 % |
| B15 | 小麦 | night/day ratio of OBT formation at equal leaf HTO | wheat, rice, potato, bean, Chlorella；review | 通用 | 1 ❌ | 1.3 ❌ | 1.3 ❌ | 0.1–0.5 - |
| B26 | 小麦 | TLI grain (OBT at harvest / leaf TFWT at end of 1 h exposure) | wheat；field box | 本作物 | 0.12 ⚠️ | 0.025 ❌ | 0.065 ❌ | 0.18–0.32 % |
| B9 | 小麦 | TFWT decrease factor from end of exposure to harvest | rice；field box | 代理 | 3 387 ✅ | 8 493 ✅ | 2 960 ✅ | 600–9.5×10⁴ - |
| B20 | 小麦 | plant OBT / plant TFWT at equilibrium (isotopic discrimination factor) | various；review | 通用 | 0.49 ⚠️ | 0.04 ❌ | 0.031 ❌ | 0.64–1.3 - |

判定数（校准集）：

| | ✅ 符合 | ⚠️ 接近 | ❌ 偏离 |
|---|---|---|---|
| 马铃薯 UFOTRI / baomi / default | 3 / 1 / 1 | 2 / 0 / 0 | 6 / 10 / 10 |
| 小麦 UFOTRI / baomi / default | 2 / 3 / 3 | 6 / 2 / 2 | 6 / 9 / 9 |

读表时注意：
- 温室（B56–B60）和舱室实验强制通风，叶–气速率比田间偏高；田间箱实验箱内温度在暴露期间升高约 10 °C（TECDOC-1738 §7.9）。
- UFOTRI 的叶–气交换是常数 0.35 h⁻¹，昼夜相同，所以它在夜间指标上“符合”只是巧合。
- 玉米（B7、B8、B64）是 C4 作物，不再用于判定马铃薯和小麦，只作 C4 参照保留在基准表中。

### 4.1 均匀比活度检验

`water_balance_residual()`：所有水库室取同一比活度、不计衰变和 OBT、去掉叶片释放中的分馏，计算 dA/dt ÷ A，
再除以该库室的流出速率 k_out。水量闭合时应 ≈ 0（只剩水量的相对变化 dlnW/dt，相对 k_out 很小）；
正值表示迁移率隐含的水流入多于流出，会造成虚假的比活度富集，负值则是虚假亏损。下表为逐时中位数（`diagnostics.csv` 中的 `wb_*`）：

| 作物 | 变体 | 土壤 1 | 土壤 2 | 土壤 3 | 植物水 | 果实水 |
|---|---|---|---|---|---|---|
| 马铃薯 | UFOTRI | −0.07 | −0.24 | −0.35 | −0.86 | +0.15 |
|  | baomi | −0.35 | −0.20 | −0.65 | +0.47 | +21 |
|  | default | −0.78 | +0.32 | +1.8 | +0.78 | +15 |
| 小麦 | UFOTRI | −0.23 | +0.08 | +0.15 | +0.38 | 0 |
|  | baomi | −0.36 | −0.11 | −0.65 | +0.90 | +165 |
|  | default | −0.97 | +0.41 | +23 | +0.66 | 0 |

- 植物水 +0.5～+0.9：根系吸水按“大气向植物的沉积通量”计算（`ka2_bh·atm_h·9`），约为蒸腾的 5 倍；对称交换的总释放已包含蒸腾，多出的水没有出口。
- 果实水：baomi 的 `kfh_bh = kbh_fh·plant_w/np.maximum(friut_w, 20)` 用 20 kg/m² 作下限，实际块茎/籽粒水只有 1–4 kg/m²；default 马铃薯是根系向块茎的吸水没有对应流出。
- 土壤：第 1 层蒸发、根吸水和下渗之和大于沉积和上行通量；default 第 3 层的大正值来自达西式（review #7）。

![叶片吸收](results/figures/fig1_leaf_uptake.png)
![交换速率](results/figures/fig2_exchange_rate.png)
![TLI](results/figures/fig3_TLI_potato.png)
![慢性](results/figures/fig4_chronic.png)

## 5 发现的问题

### 5.1 已修正

**第一步**（单位与索引，上一次提交；回归测试 `tests/test_rate_fixes.py`）：

| 问题 | 修正 | 影响的函数 |
|---|---|---|
| OBT 食入剂量系数写成 4.2e-7 mSv/Bq | ICRP-72 成人值 4.2e-8（HTO 1.8e-8 不变），集中到 `para_constant.DOSE_COEF_OBT/HTO` | 所有剂量脚本，OBT 剂量降为原来的 1/10 |
| Penman–Monteith：`CP = 1.013`、`GAM = 0.67` | `CP = 1013` J/(kg·K)、`GAM = 0.067` kPa/K | `calc_ESOIL`、`calc_ETRM` |
| `ks1_a2 = ESOIL / BODW1`（分母是体积含水率） | 改为第 1 层水量 `soil1w` | default、Korea、baomi、morris |
| `ks3_bh`、`ks3_fh` 分母用了 `soil1w` | 改为 `soil3w` | CEREAL/ROOT_VEG 的全部变体 |
| `kbh_fo = np.log(2/(HWZ/2))`（恒为负、缺 `/24`） | `np.log(2)/(HWZ/2)/24` | default、UFOTRI、Korea 的 CEREAL；TUBER_VEG |
| `HWZ[airfx.sta-1]`（用稳定度索引作物参数，`sta=6` 时 IndexError） | `HWZ[airfx.ICOMP-1]` | 同上 |
| `np.log2(10)` | `np.log(2)` | default、UFOTRI、Korea 的 CEREAL `kbh_fh`；Korea `kbo_bh` |

**第二步**（守恒与非负，本次提交；测试 `tests/test_conservation.py`）：

| 问题 | 修正 |
|---|---|
| `kbh_so` 从植物库室流出后没有去向（质量泄漏） | 流入土壤第 1 层（`dA_s1_dt` 加 `+ kbh_so·A_bh`）。现在 ODE 只通过三个外部汇损失活度：`ka1_a1`、`ka2_a2`（扩散出混合层）和 `ks3_s3`（深层渗漏） |
| `kbh_bo = TROBT` 在衰老期为负 | `calc_transOBT` 中 `TROBT = max(…, 0)`；OBT 分解由 `kbo_bh` 描述 |
| 达西式 `ks1_s2 = Va_b1/100`、`ks2_s3 = Va_b2/150` 在向上水流时为负 | `_split_darcy()` 拆成两个非负速率：负值部分按相同水通量改为下层 → 上层（加到 `ks2_s1`、`ks3_s2`） |
| 果实/块茎形成前 `friut_w = 0`，`kfh_bh = kbh_fh·plant_w/friut_w` 为 inf | `_assemble()`：果实水为 0 的时刻不向果实迁移，根系吸水全部进入植物水；植物水为 0（出苗前）的时刻进出植物的迁移率都为 0 |
| `iaea-case1-HTO.py` 用 `np.abs(k_array)` 掩盖负值 | 去掉 `abs()`；`solve_con()` 调用 `check_rates()`，遇到负值或非有限值直接报错并列出迁移率名 |
| 没有水量闭合检验 | `water_balance_residual()`：均匀比活度检验，结果写入 `diagnostics.csv`（`wb_*`、`wbabs_*`），见第 4.1 节 |

**review #8：叶–气 HTO 对称交换**（本次提交）：

- `leaf_exchange_rate()`：`kbh_a2 = VDPF · ρ_sat / (ALF_HTO · plant_w)`，与吸收项 `ka2_bh = VDPF / ML` 共用总导度
  `VDPF = 1/(RAM+RB+RC1)`。吸收通量 g·ρ_v·C_air，释放通量 g·ρ_sat·C_leaf/α，稳态 C_leaf = α·RH·C_air，夜间照常交换。
  ρ_sat 用 `calc_row()` 的饱和水汽密度（叶温取气温），`ALF_HTO = 1.1`。
- 用于 default、Korea、baomi、morris；UFOTRI 是文献常数，不变。
- baomi/morris 的 `kbh_a2 = ETRM/plant_w × 100` 一并去掉。
- 蒸腾 `k_tr = ETRM/plant_w` 只保留在液流路径上（default、Korea ROOT_VEG 的 `kbh_fh`，以及 TUBER_VEG 分支的根系吸水；
  TUBER_VEG 不在 `Plant_dict` 中，`Calc_air` 无法构建，目前是死代码）。
- `dmitte_morris.py` 扰动的 `ETRM` 不再进入 baomi ROOT_VEG 的任何迁移率，该参数的 Morris 指标会变为 0；
  `VDPF`、`plant_w` 的扰动现在同时作用于吸收和释放。

### 5.2 仍未解决（按影响排序）

1. **植物水不守恒**（review #10、#11）：对称交换的总释放 g·ρ_sat 已包含蒸腾，进入植物水的其他水流就必须等于蒸腾。
   但根系吸水按“大气向植物的沉积通量”计算（`ka2_bh·atm_h·9`），马铃薯季节平均 0.22 kg m⁻² h⁻¹，是蒸腾 ETRM（0.089）的 2.5 倍；
   `ka2_bh` 里的降雨项（平均 0.09 kg m⁻² h⁻¹）把雨水整体放进植物水，同一份雨水在 `ka2_s1` 中又进入土壤一次。
   均匀比活度检验中植物水 +0.5～+0.9 就是这两项。建议根系吸水改为 `TRA · f_root,i / soil_iw`（取 WOFOST 的 TRA），
   降雨截留单独成项并有对应流出。
2. **OBT 生成过快，且与干物质增长无关**（review #17）：24 h 后 OBT 占植物总氚 36–90 %（文献 2–4 %）；马铃薯慢性叶片 OBT/TFWT 11–12
   （文献 0.64–1.3），块茎 OBT/空气 16。`TROBT` 和 `kbh_fo` 都是一阶速率，没有按干物质增长换算。
3. **果实水下限**：baomi/morris 的 `kfh_bh = kbh_fh·plant_w/np.maximum(friut_w, 20)` 中 20 kg/m² 远大于实际块茎/籽粒水（1–4 kg/m²），
   果实水流出偏小约 20 倍，慢性块茎 TFWT/空气 12（文献 0.57–1.83）。有了 `_assemble()` 的处理，这个下限可以去掉。
4. **气孔阻力单位与昼夜变化**（review #15）：`calc_st_res` 的条件 `... | (ISTRE == 1)` 恒为真，RC1 昼夜不变，夜/昼吸收比 ≈1
   （文献 0.23–0.25）；昼/夜交换比 1.5（文献 3–5）只来自温度。`calc_ETRM` 把 `RC1 × 100` 当 s/m，`calc_VDPF` 直接把 `RC1` 当 s/m。
   按 UFOTRI 的写法 RCMIN 单位是 s/cm，白天应为 `RCMIN·RPHO`（`RPHO` 已算出但没有用），夜间为 `15·RCMIN`。
5. **干物质比例**（review #9）：`cereals_HTO.py` 用籽粒干物质比例 0.86 计算整株含水量，整株含水被低估约 10 倍；
   小麦叶–气交换因此偏快，夜间尤其明显（夜间 0.91 h⁻¹，田间箱 0.24–0.38、温室 0.15 h⁻¹）。马铃薯 `plant_w` 用 TAGP（含块茎），
   把叶片 HTO 稀释到块茎水里，是马铃薯暴露结束叶片 TFWT（6.5 %，文献 22–46 %）偏低的一部分原因。
6. **土壤–大气交换也是不对称的**：沉积用 `VDSO`、再释放用 `ESOIL`（净蒸发），和第 8 项改之前的叶片一样；
   土壤层间达西式（review #7）量纲不自洽，default 第 3 层水量闭合误差 +1.8～+23。
7. **`solve_con` 用 RK45**（review #19）：小麦 default/UFOTRI 籽粒形成初期 `kfh_bh` 达 3 900 h⁻¹，籽粒 HTO 与矩阵指数解相差 57–90 %；
   本验证用逐小时矩阵指数，不受影响。
8. **已修复的可复现性问题**：PCSE 作物参数改从 `wofost72` 分支或 `./data/crop/wofost72` 读取；随仓库提供输出 `PGASS`、`PMRES` 的
   PCSE 配置（`dmitte/conf/Wofost72_WLP_FD_dmitte.conf`）。

### 5.3 马铃薯慢性指标为什么变差

`exp_attribution.py`（结果 `results/attribution_stdout.txt`）在 baomi 上分别去掉果实水下限、把 OBT 生成放慢 10 倍、把根系吸水缩到等于蒸腾：

| 变体 | 块茎 OBT/空气 | 块茎 TFWT/空气 | 叶片 TFWT/空气 | TLI，% | 24 h OBT 占比，% |
|---|---|---|---|---|---|
| baomi（本轮） | 15.5 | 12.3 | 4.6 | 6.9 | 34 |
| 去掉 `friut_w` 下限 | 15.8 | 3.6 | 4.6 | 6.6 | 35 |
| `kbh_bo / 10` | 15.2 | 11.7 | 4.2 | 6.5 | 5.4 |
| 根系吸水 = 蒸腾 | 13.3 | 10.2 | 3.6 | 7.0 | 37 |
| 三项同时 | 13.0 | 2.2 | 3.0 | 6.3 | 5.8 |
| 文献 | 0.72–1.14 | 0.57–1.83 | α·RH ≈ 0.85（理论值） | 0.2–0.3 | 2–4 |

- 去掉果实水下限解决块茎 TFWT；放慢 OBT 生成解决 24 h OBT 占比；
- 但叶片 TFWT 仍是空气的 3 倍（对称交换单独作用时应为 α·RH ≈ 0.85）：慢性情景中叶片 HTO 只有 43–47 % 来自水汽交换，
  其余来自降雨截留（14–17 %）、根系（14–17 %）和块茎水回流（13–20 %），这些水都进入植物水而没有对应的流出；
- 块茎 OBT/空气在三项同时修改后仍为 13，说明 `kbh_fo`（块茎 OBT 形成）也需要按块茎干物质增长重写（第 5.2 节 2）。

所以下一步应把植物水的进出放到同一个水量平衡里（根系吸水 = 蒸腾，降雨截留有出口），再按干物质增长改写 OBT 生成；
之后再看是否还需要任何经验系数。

## 6 局限

- 部分基准仍来自其他作物（水稻、大豆、萝卜等），只适合判断量级和定性行为（昼夜比、衰减幅度），不适合精细标定。
- 气象是荷兰 Wageningen，不是实验所在地（韩国大田、德国卡尔斯鲁厄、加拿大 Chalk River）；逐时强迫是由日数据合成的，夜间 RH 常被截到 100 %。
- 土壤三层含水量都用 WOFOST 根区 SM，没有分层观测。
- 大气扩散（高斯烟团）部分没有验证：可用的公开示踪实验（如 Prairie Grass）在本环境无法下载。

## 7 复现

```bash
pip install -r requirements-dev.txt          # Python 3.12
python validation/build_inputs.py           # 下载公开数据，生成 ./data（约 1 MB）
python validation/run_validation.py         # 约 1 分钟，写 validation/results/*.csv
python validation/analyze.py                # 对比表 + 图
python validation/exp_attribution.py        # 归因实验（第 5.3 节）
python validation/analyze.py --include-test # 只在模型修改全部完成后运行一次：评估冻结的测试集
python -m pytest tests                      # 含第一、二步修正和对称交换的回归测试
```
