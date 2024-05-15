#%%
import numpy as np
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze

from dmitte.guassian import gaussian_puff_model
from dmitte import run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con, transfer_rates_UFORTI, transfer_rates_baomi
from dmitte.para_constant import LAMBDA_T_H

from plot_config import *
import tqdm

fig_path = os.path.join(FIG_PATH, 'compartments')

# %% ---------------------------01 计算大气库室初始浓度----------------------------


def gaussian_puff_model(x, y, z, t, Q_total, UU, HEG, sigmax, sigmay, sigmaz, t_release, puff_num: int = 600, humidity=50, rain=0, temperature=25):
    '''
    高斯烟团模型计算瞬时浓度

    Parameters:
    ---
    x, y, z: 目标位置
    t: 释放后时间点, s
    Q_total: 泄漏源瞬时泄漏的泄漏量, Bq
    t_release: 释放时间, s
    UU: 环境风速, m/s
    HEG: 烟团的有效高度, m
    stability: 大气稳定度, 1 ~ 6 对应 A ~ F
    puff_num: 烟团离散数量
    humidity: 相对湿度, %
    rain: 降雨量, mm
    temperature: 温度, ℃

    Return:
    ---
    con : t时刻空间(x,y,z)上的浓度, Bq/m3
    '''
    
    def single_gaussian(Q_single, t_single, t_rain):
        '''
        计算单个烟团的行为

        Parameters:
        ---
        Q_single: 单个烟团的释放量
        t_single: 该烟团释放后的时间
        t_rain: 烟团被雨冲刷的时间

        Return:
        ---
        result: 单个烟团导致的浓度
        '''
        np.seterr(divide='ignore', invalid='ignore')  # 消除被除数为0的警告
        
        result = Q_single / ((2 * np.pi) ** 1.5 * sigmax * sigmay * sigmaz) * \
                 np.exp(-0.5 * ((x - UU * t_single) ** 2) / (sigmax ** 2) - 0.5 * y ** 2 / sigmay ** 2) * \
                 (np.exp(-0.5 * ((z - HEG) ** 2) / sigmaz ** 2) + np.exp(-0.5 * ((z + HEG) ** 2) / sigmaz ** 2))
        
        result[t_single < 0] = 0

        # 湿沉积
        saturation = 10 ** (8.07131 - 1730.63 / (233.426 + temperature)) * 133.322           # 饱和蒸汽压, Pa
        absolute_humidity = humidity / 100 * saturation / (temperature + 273) / 461.5 * 1000 # 绝对湿度, g/m3
        k_wet = rain * 1e3 / (3600 * absolute_humidity * 1000)                               # 迁移率, s-1
        
        t_single[t_single > t_rain] = t_rain
        f_wet = np.exp(-k_wet * t_single)

        return result * f_wet


    con = np.zeros_like(x)
    for i in range(puff_num):
        con += single_gaussian(Q_total / puff_num, t - t_release / puff_num * i, t_release - i * (t_release / puff_num)) 
    
    return con
#%%
def model_function(params):
    UU, sigmay, sigmaz, RH,TEMP, Rain = params
    downwind = 1e3
    hight = np.linspace(0, 1000, 101)
    t_serial = np.linspace(0, 7200, 121)

    dd, hh, tt = np.meshgrid(downwind, hight, t_serial)
    Q_total=3.7e15
    sigmax = sigmay

    con_1 = gaussian_puff_model(dd,0,hh,tt,Q_total,UU,20,sigmax,sigmay,sigmaz,3600,600,RH,Rain,TEMP)
    init_con_1 = np.mean(con_1[:, :, :], axis=(0, 2)) * 1000 * 120 / 60
    return init_con_1


problem = {
    'num_vars': 6 ,
    'names' : ['UU','sigmay','sigmaz','RH','TEMP','Rain'],
    'bounds' :[[2,5],[35,250],[10,780],[50,90],[-5,30],[0.001,20]]
}

# %% time
# Generate samples
param_values = sample(problem, N=1000, num_levels=4, optimal_trajectories=None)

# Run model (evaluate the output of the model at each set of parameters)
Y = np.array([model_function(param) for param in param_values])

# Perform Morris analysis
Si = analyze(problem, param_values, Y, print_to_console=False, num_levels=4, num_resamples=100)

# Print results
print("Elementary Effects (mu):", Si['mu'])
print("Standard deviations (sigma):", Si['sigma'])
print("Means of absolute deviations (mu_star):", Si['mu_star'])

# %%
