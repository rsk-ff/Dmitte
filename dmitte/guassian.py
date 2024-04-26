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
    sigma = [sigmax, sigmay, sigmaz]
    return sigma

def gaussian_puff_model(x, y, z, t, Q_total, UU, HEG, stability, t_release, puff_num = 600):
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

    Return:
    ---
    con : t时刻空间(x,y,z)上的浓度, Bq/m3
    '''
    # 使用Briggs经验公式 计算大气稳定度：
    sigma = calc_sigmaxyz(stability, x)
    sigmax = sigma[0]
    sigmay = sigma[1]
    sigmaz = sigma[2]

    def single_gaussian(Q_single, t_single):
        '''
        计算单个烟团的行为

        Parameters:
        ---
        Q_single: 单个烟团的释放量
        t_single: 该烟团释放后的时间

        Return:
        ---
        result: 单个烟团导致的浓度
        '''
        np.seterr(divide='ignore', invalid='ignore')  # 消除被除数为0的警告
        exp0 = Q_single / ((2 * np.pi) ** 1.5 * sigmax * sigmay * sigmaz)
        exp1 = np.exp(-0.5 * ((x - UU * t_single) ** 2) / (sigmax ** 2))
        exp2 = np.exp(- y ** 2 / sigmay ** 2)
        exp3 = np.exp(-0.5 * ((z - HEG) ** 2) / sigmaz ** 2)
        exp4 = np.exp(-0.5 * ((z + HEG) ** 2) / sigmaz ** 2)
        result = exp0 * exp1 * exp2 * (exp3 + exp4)
        
        return result

      # 离散烟团数量
    release_puff_num = min(puff_num, int(t / t_release * puff_num) + 1) # t 时刻已经释放的烟团数
    
    con = np.zeros_like(x)
    for i in range(int(release_puff_num)):
        con += single_gaussian(Q_total / puff_num, t - t_release / release_puff_num * i) 
    
    return con
