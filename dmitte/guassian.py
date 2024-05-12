import numpy as np


def calc_sigmaxyz(stability, x):
    if stability == 1:
        sigmax = sigmay = 0.22 * x * (1 + 0.0001 * x) ** (-0.5)
        sigmaz = 0.2 * x
    elif stability == 2:
        sigmax = sigmay = 0.16 * x * (1 + 0.0001 * x) ** (-0.5)
        sigmaz = 0.12 * x
    elif stability == 3:
        sigmax = sigmay = 0.11 * x * (1 + 0.0001 * x) ** (-1 / 2)
        sigmaz = 0.08 * x * (1 + 0.0002 * x) ** (-1 / 2)
    elif stability == 4:
        sigmax = sigmay = 0.08 * x * (1 + 0.0001 * x) ** (-1 / 2)
        sigmaz = 0.06 * x * (1 + 0.0015 * x) ** (-1 / 2)
    elif stability == 5:
        sigmax = sigmay = 0.06 * x * (1 + 0.0001 * x) ** (-1 / 2)
        sigmaz = 0.03 * x * (1 + 0.0003 * x) ** (-1)
    elif stability == 6:
        sigmax = sigmay = 0.04 * x * (1 + 0.0001 * x) ** (-1 / 2)
        sigmaz = 0.016 * x * (1 + 0.0003 * x) ** (-1)
    else:
        raise ValueError('Invalid stability class')
    return sigmax, sigmay, sigmaz


def gaussian_puff_model(x, y, z, t, Q_total, UU, HEG, stability, t_release, puff_num: int = 600, humidity=50, rain=0, temperature=25):
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

    # 使用Briggs经验公式 计算大气稳定度：
    sigmax, sigmay, sigmaz = calc_sigmaxyz(stability, x)
   
    con = np.zeros_like(x)
    for i in range(puff_num):
        con += single_gaussian(Q_total / puff_num, t - t_release / puff_num * i, t_release - i * (t_release / puff_num)) 
    
    return con
