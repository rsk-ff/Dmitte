# %%
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze

from plot_config import *


def model_function(params):
    """
    RELFEU (float): 相对湿度 (0-1)
    TEMP (float): 温度 (摄氏度)
    PAR (float): 光合有效辐射 (W/m^2)
    LEAFA (float): 叶面积指数 (m^2/m^2)
    """
    RATMOS, RSTOMA, TEMP, r_h, ea, PAR, LAI = params
    # RATMOS = RAM + RB             # 大气阻滞因子之和,单位 m/s
    # RSTOMA = RC1 * 100.                # 植物冠层气孔阻力因子,单位 m/s
    es = 0.6108 * np.exp(17.27 * TEMP / (TEMP + 237.3))  # 饱和蒸汽压,kpa,Tetens公式
    EFEUDE = (es - ea * r_h)         # 空气实际与饱和的蒸汽压之差,N m-2
    DELE = 4098 * es / (TEMP + 237.3)**2    # 饱和水汽压的温度依赖,单位 Pa
    
    GAM = 0.67                            # 大气干湿度压力系数
    RHO = 1.2923                          # 空气密度,单位 kg/m^3
    CP = 1.013                            # 比热容,单位 J/(kg K)
    LE = 2.45E+6                          # 水的蒸发潜热,单位 J/kg
    
    QSTR = PAR * (1- np.exp(-0.398 * LAI))
    ETRM = ((DELE * QSTR) + (RHO * CP * EFEUDE / RATMOS)) / (DELE + GAM * (1. + RSTOMA / RATMOS)) / LE

    return ETRM

problem = {
    'num_vars': 7, # Number of parameters
    'names': ['RA', 'RST', 'TEMP', 'RH', 'EA', 'PAR', 'LAI'],
    'bounds': [[40, 250], [5, 160], [-5, 30], [0.5, 0.9], [0.1, 3], [5, 400], [0.1, 4]]
}
# Generate samples
param_values = sample(problem, N=1000, num_levels=4, optimal_trajectories=None)

Y = np.array([model_function(params) for params in param_values])

Si = analyze(problem, param_values, Y, print_to_console=True, num_levels=4, num_resamples=100)
Si.plot()

# %%
