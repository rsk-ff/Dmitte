# Dmitte 交接说明（2026-10-29）

给接手的会话：先读本文件，再读 `validation/README.md`（完整的验证记录）和 `validation/benchmarks/README.md`（基准说明）。
用户用中文交流，偏好简洁、结论先行；论文方向是“基于作物生长模型的动态氚迁移模型 + 敏感性分析（导师要求）”。

## 1 仓库与环境

- 仓库 `rsk-ff/Dmitte`，工作分支 `claude/adoring-turing-wlqwpu`（全部改动都在这里，未开 PR）。最新提交见第 6 节。
- Python 3.12（pcse 5.5.6 在 3.13 下装不上）：`pip install -r requirements-dev.txt`。
- `./data` 不入库：先运行 `python validation/build_inputs.py`（下载 Wageningen 气象、土壤、WOFOST 7.2 作物参数，约 1 MB）。
- 画图脚本需要 `MPLBACKEND=Agg`。
- **换行符**：`dmitte/calc_para.py`、`dmitte/transfer_equotions.py`、`dmitte/para_constant.py`、根目录 `*_HTO.py` 是 CRLF，
  前两个还混有少量 LF 行。编辑时保持原有换行符（不要整文件转换），否则 diff 会变成全文件。
- 提交信息不写模型名；按会话要求在末尾加 Co-Authored-By / Claude-Session 两行。

## 2 模型现状（`dmitte/`）

- 11 个库室，顺序见 `transfer_equotions.A_NAMES`：a1 s0 a2 s1 s2 s3 bh bo fh fo ba（`ba` 为叶片新鲜同化物 OBT，放在最后，前 10 个编号不变）。
  迁移率列顺序见 `K_NAMES`（33 个）。按名称取：`A_INDEX['fo']`、`K_INDEX['kbh_a2']`。
- 只有一个迁移率函数 `transfer_rates(char, plant_type, airfx, soilfx, plantfx, percolation='constant', kbh_fh=ln2/80)`；
  default/baomi/Korea/morris/UFOTRI 变体已合并或删除。`solve_con()` 用逐小时矩阵指数（精确，旧 RK45 在刚性时刻误差可达 100 %）。
- 关键过程：
  - 叶–气 HTO 对称交换 `kbh_a2 = VDPF·ρsat/(α·W_leaf)`，`ALF_HTO = 1.1`；
  - 植物水 = 活叶片水（TWLV × LAI/LAI_max × 0.8/0.2）；根系吸水 = 净蒸腾 + kbh_so·W + d(W_bh+W_fh)/dt（`plant_water_routes`），植物水和果实水按构造闭合；降雨只进土壤；
  - 风速廓线按 10 m 测量高度外推（原来当成 100 m，u* 低 5 倍）；
  - 气孔：白天 `RCMIN·RPHO`，夜间及 RH > 98 % 为 `RC_NIGHT·RCMIN`，`RGES = RAM + RB + 100·RC1`；
  - OBT（TECDOC-1991 式 43 的增长形式）：叶片水 → `ba` 按总光合（WOFOST PGASS）；`ba` 以 `TAU_ASSIM` 流出，按当天干物质分配到 bo / fo / 回叶片水；
    结构 OBT 不周转（`KBO_BH = 0`）；日内按辐射分配，夜间为白天均值的 `OBT_DARK` 倍。
- 参数（`para_constant.py`）及来源：`DMF_LEAF 0.2`、`ID_OBT 0.8`、`COMB_WATER 0.556`、`COMB_CH2O 0.6`（文献/化学计量）；
  `RC_NIGHT 4`（校准集标定，原 15）、`TAU_ASSIM 6 h`（校准集标定）、`OBT_DARK 0.3`（取基准 B15 中值，B15 因此不独立）、
  `PERCOLATION 'constant'`、`KBH_FH ln2/80`（基准无法区分，留给敏感性分析）。
- 剂量系数已改为 ICRP-72（`DOSE_COEF_HTO 1.8e-8`、`DOSE_COEF_OBT 4.2e-8` mSv/Bq）。

## 3 验证体系（`validation/`）

运行顺序（约 3 分钟）：

```bash
python validation/run_validation.py        # 通用虚拟实验（Wageningen 天气，马铃薯 2004–08、小麦 2005–08）
python validation/scenario_field_box.py    # 小麦田间箱 14 次暴露按实测条件逐次驱动
python validation/analyze.py               # 只评校准集；结果 results/benchmark_comparison.csv
python -m pytest tests                     # 62 个测试
# 标定脚本：calibrate_night_resistance.py（约 15 分钟）、calibrate_assimilate_residence.py（约 2 分钟）
```

- 基准 v3：`benchmarks/literature_benchmarks.csv` 65 条、`datasets/` 20 个、`scenario_cards.csv` 15 张、`split.csv`（校准 50 / 测试 15）。
  测试集指纹在 `analyze.py`（`FROZEN_TEST_FINGERPRINT`），`tests/test_smoke.py` 检查。
- 判定：落在文献范围内 ✅，2 倍以内 ⚠️（TECDOC-1738 §7.9），更远 ❌。
- 田间箱 8 条基准（B4 B5 B22 B23 B26 B29 B49 B50）用实测条件驱动的结果评分，箱内空气动力学阻力取 15 s/m（8–30 s/m 作范围，原文未报告）。

## 4 结果

- 校准集（当前模型）：马铃薯 ✅7 ⚠️4 ❌0，小麦 ✅6 ⚠️7 ❌1；25 项中 24 项在 2 倍以内。修正前为马铃薯 1/0/10、小麦 3/2/9。
- **测试集已于 2026-10-29 用 commit `4e09ac5` 运行一次**（`validation/README.md` 第 4.6 节）：13 项有模型对应，✅3 ⚠️8 ❌2，11/13 在 2 倍以内。
  ❌ 为舱室强光小麦叶片吸收 B46（低 2.1 倍，通用实验未驱动舱室条件）和 Pickering 块根慢性 OBT/空气 B34（高 2.3 倍，与 Chalk River B18 矛盾）。
  B36、B37 冻结时就没有模型检查，未补加。籽粒 TLI 随灌浆阶段的变化（B53–B55）在测试集上成立，是论文“动态迁移率”的主要证据。
- **测试集已经用掉**：之后改模型，测试结果不再是独立检验；如需再看测试数据，必须作为明确标注的补充分析。

## 5 已知问题（按影响）

1. 叶片 OBT 短期过程：模型在暴露后 1–3 h 继续生成、4 h 偏高 1.6 倍，实测在暴露结束时已形成；夜间暴露结束时 OBT 低约 2.6 倍（B23 ❌）。需要比总光合更快的标记机制。
2. 气孔在黎明/黄昏的响应：模型只随瞬时光照，黎明开得太小、黄昏开得太大（田间箱暴露后损失速率黎明 0.32 vs 0.69–1.04 h⁻¹，黄昏 0.37 vs 0.06–0.18）。
3. 土壤水不闭合（达西式量纲不自洽；土壤沉积/再释放不对称）。
4. B15 指标是 24 h 积分，与基准定义（相同叶片 HTO 下的生成比）不同。
5. 近似：活叶按 LAI/LAI_max 折算；籽粒含水全程按 0.86 干物质比例；`figures/`、`dose_HTO*` 等是修正前的旧结果（原始站点数据不在仓库）。
6. `pcse_morris.py` 已改为调用 `dmitte.calc_para`，以前的 Morris 结果来自修正前的模型，需要重做。

## 6 提交记录（本分支，最新在上）

`fbd5091` 运行测试集 · `4e09ac5` 叶片同化物库、KBO_BH → 0 · `76b47eb` 田间箱用驱动结果评分、重新标定 RC_NIGHT ·
`d4c06a4` 田间箱按实测条件驱动 · `bab578f` 标定 RC_NIGHT · `aeda958` 合并迁移率变体、删除 UFOTRI ·
`5ffc134` 植物水平衡、气孔昼夜、OBT 按生长 · `bad8dc1` 基准 v3 与冻结测试集 · `56de19c` 质量守恒与对称交换 · `0faccbf` 单位与索引错误

## 7 下一步（用户尚未选定）

- A. 补充分析：用 `WHEAT_CH_HOURLY_DAY/NIGHT` 情景卡（Diabaté & Strack 1997 逐时数据，`datasets/wheat_DE_chamber_timeseries_1993.csv`，测试数据）
  驱动舱室实验，明确标注为测试后的补充分析，说明 B46 的偏差有多少来自实验条件。我的建议是先做这一项。
- B. 论文分析：动态 vs 常数迁移率（用本模型冻结时间变化的版本作基线，不用 UFOTRI）、释放日期扫描、EMRAS TC14 第 9 章假想情景
  （仓库 Case 1/2/3 对应稳定度 1/4/6）、Morris → Sobol 敏感性分析（候选参数：RC_NIGHT 3–6、TAU_ASSIM、KBH_FH、PERCOLATION、箱内阻力等）。
