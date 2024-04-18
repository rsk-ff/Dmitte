import numpy as np

# ------------------- 库室模型用 -----------------------
RHO = 1.2923       # 空气密度，单位 kg/m^3
DISSOZ = 1.0       # HTO分布耗散系数
CO2DRY = 1.5       # co2到干物质的转化因子
Fsoil = 3.4E4      # 土层之间的质量流量
mixLayerH = [1600, 1200, 800, 560, 320, 200]                  # 混合层高度
LAMBDA_T_H = np.log(2) / (12.43 * 365 * 24)   # 氚衰变常数, 单位 h-1
LAMBDA_T_D = np.log(2) / (12.43 * 365)        # 氚衰变常数, 单位 d-1

# 大气参数[对应A-F大气稳定度]
ISK = [0.07, 0.13, 0.21, 0.34, 0.44, 0.44]                    # 风廓线幂指数
A0 = [-0.1135, -0.0385, -0.0081, 0, 0.0081, 0.0385]
B0 = [-0.1025, -0.1710, -0.3045, -0.5030, -0.3045, -0.1710]

# --------------------- 土壤植物模型--------------------
# SMW: soil moisture content at wilting point [cm3/cm3]
# SMFCF: soil moisture content at field capacity [cm3/cm3]
# SM0: soil moisture content at saturation [cm3/cm3]
# K0: hydraulic conductivity of saturated soil [cm day-1]
Soil_dict = {'SOIL1': {'SMW': 0.04, 'SWFCF': 0.11, 'SM0': 0.39, 'K0': 99.7},
             'SOIL2': {'SMW': 0.09, 'SWFCF': 0.27, 'SM0': 0.39, 'K0': 23.9},
             'SOIL3': {'SMW': 0.10, 'SWFCF': 0.30, 'SM0': 0.41, 'K0': 25.6}}

# 作物参数[对应区域1、区域2、土豆、小麦]
HWZ = [60, 60, 90, 60, 50]             # 土豆的生长期为90天,小麦生长期为60天；
LT = [90, 303, 455, 670]               # [第一年植物生长期的第一天，最后一天，第二年植物生长期的第一天，最后一天]
RCMIN = [2.0, 2.0, 1.5, 2.0]           # 按叶面指数归一化的最小气孔扩散阻滞因子
Plant_dict = {'LEAFY_PLANT': {'ICOMP': 1, 'TMAXL': 45.0, 'TMINL': 0, 'TOPTL': 25.0, 'VPD': 0.2, 'RCPAR': 20.0},
              'GRASS':     {'ICOMP': 2, 'TMAXL': 45.0, 'TMINL': 0, 'TOPTL': 25.0, 'VPD': 0.2, 'RCPAR': 20.0},
              'ROOT_VEG':  {'ICOMP': 3, 'TMAXL': 45.0, 'TMINL': 0, 'TOPTL': 25.0, 'VPD': 0.2, 'RCPAR': 50.0},
              'CEREAL':    {'ICOMP': 4, 'TMAXL': 45.0, 'TMINL': 0, 'TOPTL': 25.0, 'VPD': 0.2, 'RCPAR': 30.0}}

# 动物参数
breathRate = 130.  # 奶牛的呼吸速率
NCOW = 250.        # 在每平方公里放牧牛的数量
WWCOW = 350.       # 牛的无机部分(水)的质量, KG/牛
WOCOW = 150.       # 牛的有机部分的质量, KG/牛
M_cHTO = 38.9      # 一头奶牛的库室中的HTO中氢的质量，kg/ km2;
M_cOBT = 11.5      # 一头奶牛的库室中的OBT中氢的质量，kg/ km2;
MgOBT = 11800      # 牧草有机库室中氢的质量，kg/ km2;
M_mHTO = 0.097     # 牛奶HTO中氢的质量，kg/ km2;



