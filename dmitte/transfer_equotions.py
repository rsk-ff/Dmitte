# %%
import numpy as np
from functools import partial
from scipy.integrate import solve_ivp

from dmitte import para_constant
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant


def transfer_rates(char: str, 
                   plant_type: str, 
                   airfx: Calc_air, 
                   soilfx: Calc_soil, 
                   plantfx: Calc_plant) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG' or 'TUBER_VEG'

    Return:
    ---
    k_array: 迁移率矩阵
    '''
    ka2_s1 = (soilfx.VDSO / airfx.ML * 3600) + (airfx.rainfall / 9) / airfx.atm_h
    ks3_s2 = ks3_s3 = 3.42E-4
    ks2_s1 = ks3_s3 * soilfx.soil3w / soilfx.soil2w
    ks1_s2 = soilfx.Va_b1 / 100
    ks2_s3 = soilfx.Va_b2 / 150
    kbh_so = 3E-4
    ka2_bh = (plantfx.VDPF / airfx.ML * 3600) + (airfx.rainfall /9) / airfx.atm_h
    kbh_a2 = plantfx.ETRM / (plantfx.plant_w)
    kbo_bh = (np.log(2) / 10) /24
    
    if char.upper() == 'HT':
        ka1_a1 = 0.693 
        ka2_a2 = 0 
        ka1_a2 = 4E-3
        ka1_s0 = 4.5e-4 / airfx.ML * 3600 
        ks0_s1 = 0.2
        ks0_s2 = 0.02
        ks0_s3 = 0.02
        ks0_a1 = ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (0.5+0.1) 
        
    elif char.upper() == 'HTO':
        ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
        ka2_a2 = 0.693
        ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (1+0.1) 

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))


    if plant_type.upper() == 'LEAFY_PLANT':  
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil3w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = 0
        kfh_bh = 0
        kbh_fo = 0
        kbh_bo = 0

    elif plant_type.upper() == 'CEREAL':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = (np.log2(10) / 2) 
        kfh_bh = kbh_fh * plantfx.plant_w / plantfx.friut_w        
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'ROOT_VEG':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (5/6)
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (5/6)
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (5/6)
        ks1_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (1/6)
        ks2_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (1/6)
        ks3_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (1/6)
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = kbh_fh * (plantfx.plant_w / plantfx.friut_w ) 
        kbh_fo = np.log(2)/(para_constant.HWZ[airfx.sta-1]/2) * (plantfx.friut_wh / plantfx.plant_wh) / 24
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'TUBER_VEG':
        ks1_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (5/6)
        ks2_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (5/6)
        ks3_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (5/6)
        ks1_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (1/6)
        ks2_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (1/6)
        ks3_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (1/6)
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = (kbh_a2 /2) * (plantfx.plant_w / plantfx.friut_w )  * (2/6)
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    else:
        raise ValueError("Unsupported plant type: '{}'.".format(plant_type))


    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, 
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, 
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, 
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh]

    n_t_step = len(airfx.U10)

    k_array = np.array([np.full(n_t_step, k, dtype=np.float64) if np.isscalar(k) else k for k in k_list])

    return k_array.T


def transfer_equotions(t, y, transfer_rates, LAMBDA_T: float):
    '''
    氚转移微分方程组
    
    Parameters:
    ---
    y, t: 形式需要, y 为待求解的变量向量, t 表示时间
    transfer_rates: d-1 or h-1, 迁移率列表或向量
    LAMBDA_T: d-1 or h-1, 氚的衰变常数, 单位与迁移率保持一致
    '''
    A_a1, A_s0, A_a2, A_s1, A_s2, A_s3, A_bh, A_bo, A_fh, A_fo = y
    ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, ka2_s1, \
        ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, ks3_bh, ka2_a2, \
        kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, kbh_fh, kbh_fo, ks1_fh, \
        ks2_fh, ks3_fh = transfer_rates

    
    dA_a1_dt = ks0_a1 * A_s0 - (ka1_a1 + ka1_a2 + ka1_s0 + LAMBDA_T) * A_a1
    dA_s0_dt = ka1_s0 * A_a1 - (ks0_a1 + ks0_s1 + ks0_s2 + ks0_s3 + LAMBDA_T) * A_s0
    dA_a2_dt = ka1_a2 * A_a1 + ks1_a2 * A_s1 + kbh_a2 * A_bh - (ka2_a2 + ka2_bh + ka2_s1 + LAMBDA_T) * A_a2
    dA_s1_dt = ks0_s1 * A_s0 + ka2_s1 * A_a2 + ks2_s1 * A_s2 - (ks1_a2 + ks1_bh + ks1_fh + ks1_s2 + LAMBDA_T) * A_s1
    dA_s2_dt = ks0_s2 * A_s0 + ks1_s2 * A_s1 + ks3_s2 * A_s3 - (ks2_s1 + ks2_bh + ks2_fh + ks2_s3 + LAMBDA_T) * A_s2
    dA_s3_dt = ks0_s3 * A_s0 + ks2_s3 * A_s2 - (ks3_s3 + ks3_bh + ks3_fh + ks3_s2 + LAMBDA_T) * A_s3
    dA_bh_dt = ka2_bh * A_a2 + ks1_bh * A_s1 + ks2_bh * A_s2 + ks3_bh * A_s3 + kbo_bh * A_bo + kfh_bh * A_fh - (kbh_a2 + kbh_so + kbh_bo + kbh_fh + kbh_fo + LAMBDA_T) * A_bh
    dA_bo_dt = kbh_bo * A_bh - (kbo_bh + LAMBDA_T) * A_bo
    dA_fh_dt = ks1_fh * A_s1 + ks2_fh * A_s2 + ks3_fh * A_s3 + kbh_fh * A_bh - (kfh_bh + LAMBDA_T) * A_fh
    dA_fo_dt = kbh_fo * A_bh - LAMBDA_T * A_fo

    return dA_a1_dt, dA_s0_dt, dA_a2_dt, dA_s1_dt, dA_s2_dt, dA_s3_dt, dA_bh_dt, dA_bo_dt, dA_fh_dt, dA_fo_dt


def solve_con(char, init_con, k_array, LAMBDA_T):
    '''
    计算各库室浓度

    Parameters:
    ---
    char: 'HT' or 'HTO'
    init_con: 初始浓度, Bq/m2
    k_array: 迁移率矩阵
    LAMBDA_T: 氚的衰变常数

    Return:
    ---
    As: 库室浓度矩阵
    '''
    def next_con(Ai, ki):
        '''
        定义单步求解过程, 请勿外部调用
        '''
        equotions_with_params = partial(transfer_equotions, transfer_rates=ki, LAMBDA_T=LAMBDA_T)
        sol = solve_ivp(equotions_with_params, t_span=(0, 1), y0=Ai, t_eval=[1])
        return sol.y.squeeze()

    k_array[0, 0] = 0
    k_array[0, 17] = 0

    n_t_step = len(k_array)
    As = np.zeros((n_t_step + 1, 10))

    if char == 'HT':
        As[0] = [init_con, 0, 0, 0, 0, 0, 0, 0, 0, 0]

        As[1] = next_con(As[0], k_array[0])

        As[1, 0] = 0
        As[1, 2] = 0
        
        for i in range(1, n_t_step):
            As[i + 1] = next_con(As[i], k_array[i])

    elif char == 'HTO':
        As[0] = [0, 0, init_con, 0, 0, 0, 0, 0, 0, 0]

        As[1] = next_con(As[0], k_array[0])

        As[1, 0] = 0
        As[1, 2] = 0
        
        for i in range(1, n_t_step):
            As[i + 1] = next_con(As[i], k_array[i])

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))

    return As


def transfer_rates_UFORTI(char: str, 
                   plant_type: str, 
                   airfx: Calc_air, 
                   soilfx: Calc_soil, 
                   plantfx: Calc_plant) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG' or 'TUBER_VEG'

    Return:
    ---
    k_array: 迁移率矩阵
    '''
    ka2_s1 = 0.68 / 24
    ks3_s2 = ks3_s3 = 8.2e-3 / 24
    ks2_s1 = 1.2e-2 / 24
    ks1_s2 = 0.15 / 24
    ks2_s3 = 5.0e-2 / 24
    kbh_so = 0
    ka2_bh = 0.205 / 24
    kbh_a2 = 8.3 / 24
    kbo_bh = 6.9e-2 / 24
    
    if char.upper() == 'HT':
        ka1_a1 = 0.693 
        ka2_a2 = 0 
        ka1_a2 = 4E-3
        ka1_s0 = 4.5e-4 / airfx.ML * 3600 
        ks0_s1 = 0.2
        ks0_s2 = 0.02
        ks0_s3 = 0.02
        ks0_a1 = ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (0.5+0.1) 
        
    elif char.upper() == 'HTO':
        ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
        ka2_a2 = 0.693
        ks1_a2 = 0.27 / 24

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))


    if plant_type.upper() == 'LEAFY_PLANT':  
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil3w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = 0
        kfh_bh = 0
        kbh_fo = 0
        kbh_bo = 0

    elif plant_type.upper() == 'CEREAL':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = (np.log2(10) / 2) 
        kfh_bh = kbh_fh * plantfx.plant_w / plantfx.friut_w        
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'ROOT_VEG':
        ks1_bh = 4.4e-2 / 24
        ks2_bh = 4.4e-2 / 24
        ks3_bh = 2.9e-2 / 24
        ks1_fh = 8.8e-3 / 24
        ks2_fh = 8.8e-3 / 24
        ks3_fh = 5.8e-3 / 24
        kbh_fh = 0.2 / 24
        kfh_bh = 0.47 / 24
        kbh_fo = 2.0e-2 / 24
        kbh_bo = 1.2e-2 / 24

    elif plant_type.upper() == 'TUBER_VEG':
        ks1_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (5/6)
        ks2_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (5/6)
        ks3_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (5/6)
        ks1_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (1/6)
        ks2_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (1/6)
        ks3_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (1/6)
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = (kbh_a2 /2) * (plantfx.plant_w / plantfx.friut_w )  * (2/6)
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    else:
        raise ValueError("Unsupported plant type: '{}'.".format(plant_type))


    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, 
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, 
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, 
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh]

    n_t_step = len(airfx.U10)

    k_array = np.array([np.full(n_t_step, k, dtype=np.float64) if np.isscalar(k) else k for k in k_list])

    return k_array.T


def transfer_rates_Korea(char: str, 
                   plant_type: str, 
                   airfx: Calc_air, 
                   soilfx: Calc_soil, 
                   plantfx: Calc_plant) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG' or 'TUBER_VEG'

    Return:
    ---
    k_array: 迁移率矩阵
    '''
    ka2_s1 = (soilfx.VDSO / airfx.ML * 3600) + (airfx.rainfall / 9) / airfx.atm_h
    ks3_s2 = ks3_s3 = 3.42E-4
    ks2_s1 = ks3_s3 * soilfx.soil3w / soilfx.soil2w
    ks1_s2 = soilfx.Va_b1 / 100
    ks2_s3 = soilfx.Va_b2 / 150
    kbh_so = 3E-4
    ka2_bh = (plantfx.VDPF / airfx.ML * 3600) + (airfx.rainfall /9) / airfx.atm_h
    kbh_a2 = plantfx.ETRM / (plantfx.plant_w)
    kbo_bh = (np.log2(10) / 10) /24
    
    if char.upper() == 'HT':
        ka1_a1 = 0.693 
        ka2_a2 = 0 
        ka1_a2 = 4E-3
        ka1_s0 = 4.5e-4 / airfx.ML * 3600 
        ks0_s1 = 0.2
        ks0_s2 = 0.02
        ks0_s3 = 0.02
        ks0_a1 = ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (0.5+0.1) 
        
    elif char.upper() == 'HTO':
        ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
        ka2_a2 = 0.693
        ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (1+0.1) 

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))


    if plant_type.upper() == 'LEAFY_PLANT':  
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil3w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = 0
        kfh_bh = 0
        kbh_fo = 0
        kbh_bo = 0

    elif plant_type.upper() == 'CEREAL':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = (np.log2(10) / 2) 
        kfh_bh = kbh_fh * plantfx.plant_w / plantfx.friut_w        
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'ROOT_VEG':
        ks1_bh = 0
        ks2_bh = 0
        ks3_bh = 0
        ks1_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = kbh_fh * (plantfx.plant_w / plantfx.friut_w ) 
        kbh_fo = np.log(2)/(para_constant.HWZ[airfx.sta-1]/2) * (plantfx.friut_wh / plantfx.plant_wh) / 24
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'TUBER_VEG':
        ks1_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (5/6)
        ks2_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (5/6)
        ks3_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (5/6)
        ks1_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (1/6)
        ks2_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (1/6)
        ks3_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (1/6)
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = (kbh_a2 /2) * (plantfx.plant_w / plantfx.friut_w )  * (2/6)
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    else:
        raise ValueError("Unsupported plant type: '{}'.".format(plant_type))


    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, 
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, 
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, 
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh]

    n_t_step = len(airfx.U10)

    k_array = np.array([np.full(n_t_step, k, dtype=np.float64) if np.isscalar(k) else k for k in k_list])

    return k_array.T


def transfer_rates_baomi(char: str, 
                   plant_type: str, 
                   airfx: Calc_air, 
                   soilfx: Calc_soil, 
                   plantfx: Calc_plant) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG' or 'TUBER_VEG'

    Return:
    ---
    k_array: 迁移率矩阵
    '''
    ka2_s1 = (soilfx.VDSO / airfx.ML * 3600) + (airfx.rainfall / 9) / airfx.atm_h
    ks3_s2 = ks3_s3 = 3.42E-4
    ks2_s1 = ks3_s3 * soilfx.soil3w / soilfx.soil2w
    kbh_so = 3E-4
    ka2_bh = (plantfx.VDPF / airfx.ML * 3600) + (airfx.rainfall /9) / airfx.atm_h
    kbh_a2 = plantfx.ETRM / (plantfx.plant_w) * 100
    kbo_bh = (np.log(2) / 10) /24
    
    if char.upper() == 'HT':
        ka1_a1 = 0.693 
        ka2_a2 = 0 
        ka1_a2 = 4E-3
        ka1_s0 = 4.5e-4 / airfx.ML * 3600 
        ks0_s1 = 0.2
        ks0_s2 = 0.02
        ks0_s3 = 0.02
        ks0_a1 = ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (0.5+0.1) 
        
    elif char.upper() == 'HTO':
        ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
        ka2_a2 = 0.693
        ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (1+0.1) 

    else:
        raise ValueError("Unsupported char: '{}'.".format(char))


    if plant_type.upper() == 'LEAFY_PLANT':  
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil3w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = 0
        kfh_bh = 0
        kbh_fo = 0
        kbh_bo = 0

    elif plant_type.upper() == 'CEREAL':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4
        ks1_fh = 0
        ks2_fh = 0
        ks3_fh = 0
        kbh_fh = (np.log(2) / 2) / 40
        kfh_bh = kbh_fh * plantfx.plant_w / np.maximum(plantfx.friut_w, 20)        
        kbh_fo = np.log(2) / 60 * (plantfx.friut_wh / plantfx.plant_wh) / 24
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'ROOT_VEG':
        ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (5/6)
        ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (5/6)
        ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (5/6)
        ks1_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (1/6)
        ks2_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (1/6)
        ks3_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (1/6)
        kbh_fh = np.log(2) / 2 / 40
        kfh_bh = kbh_fh * (plantfx.plant_w / np.maximum(plantfx.friut_w, 20) ) 
        kbh_fo = np.log(2)/ 45 * (plantfx.friut_wh / plantfx.plant_wh) / 24
        kbh_bo = plantfx.TROBT

    elif plant_type.upper() == 'TUBER_VEG':
        ks1_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (5/6)
        ks2_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (5/6)
        ks3_bh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (5/6)
        ks1_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil1w) * 0.2 * (1/6)
        ks2_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil2w) * 0.4 * (1/6)
        ks3_fh = ((kbh_a2/2) * plantfx.plant_w  / soilfx.soil3w) * 0.4 * (1/6)
        kbh_fh = (kbh_a2 /2) * (1/6) 
        kfh_bh = (kbh_a2 /2) * (plantfx.plant_w / plantfx.friut_w )  * (2/6)
        kbh_fo = np.log(2/(para_constant.HWZ[airfx.sta-1]/2)) * (plantfx.friut_wh / plantfx.plant_wh)
        kbh_bo = plantfx.TROBT

    else:
        raise ValueError("Unsupported plant type: '{}'.".format(plant_type))

    
    # ks1_s2 = (ka2_s1 * airfx.atm_w + ks2_s1 * soilfx.soil2w - ks1_bh * soilfx.soil1w - ks1_a2 * soilfx.soil1w) / soilfx.soil1w
    ks1_s2 = 0.181 * 9 / soilfx.soil1w / 24
    ks2_s3 = 0.072 * 9 / soilfx.soil2w / 24

    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, 
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, 
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, 
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh]

    n_t_step = len(airfx.U10)

    k_array = np.array([np.full(n_t_step, k, dtype=np.float64) if np.isscalar(k) else k for k in k_list])

    return k_array.T


def transfer_rates_morris(char: str, 
                   plant_type: str, 
                   airfx: Calc_air, 
                   soilfx: Calc_soil, 
                   plantfx: Calc_plant) -> np.ndarray:
    '''
    计算迁移率矩阵

    Parameters:
    ---
    char: 'HT' or 'HTO', HT 或者 HTO 释放
    airfx: 大气相关参数类
    soilfx: 土壤参数计算类
    plantfx: 植物相关参数类
    plant_type: 植物类型, 'LEAFY_PLANT' or 'CEREAL' or 'ROOT_VEG' or 'TUBER_VEG'

    Return:
    ---
    k_array: 迁移率矩阵
    '''
    ka2_s1 = (soilfx.VDSO / airfx.ML * 3600) + (airfx.rainfall / 9) / airfx.atm_h
    ks3_s2 = ks3_s3 = 3.42E-4
    ks2_s1 = ks3_s3 * soilfx.soil3w / soilfx.soil2w
    kbh_so = 3E-4
    ka2_bh = (plantfx.VDPF / airfx.ML * 3600) + (airfx.rainfall /9) / airfx.atm_h
    kbh_a2 = plantfx.ETRM / (plantfx.plant_w) * 100
    kbo_bh = (np.log(2) / 10) /24
            
    ka1_a1 = ka1_s0 =ka1_a2 =ks0_a1= ks0_s1 =ks0_s2 = ks0_s3=0
    ka2_a2 = 0.693
    ks1_a2 = soilfx.ESOIL / soilfx.BODW1 * (1+0.1) 

    ks1_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (5/6)
    ks2_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (5/6)
    ks3_bh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (5/6)
    ks1_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.2* (1/6)
    ks2_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil2w) * 0.4* (1/6)
    ks3_fh = (ka2_bh * airfx.atm_h * 9 / soilfx.soil1w) * 0.4* (1/6)
    kbh_fh = np.log(2) / 2 / 40
    kfh_bh = kbh_fh * (plantfx.plant_w / np.maximum(plantfx.friut_w, 20) ) 
    kbh_fo = np.log(2)/ 45 * (plantfx.friut_wh / plantfx.plant_wh) / 24
    kbh_bo = plantfx.TROBT
   
    ks1_s2 = 0.181 * 9 / soilfx.soil1w / 24
    ks2_s3 = 0.072 * 9 / soilfx.soil2w / 24

    k_list = [ka1_a1, ks0_a1, ka1_s0, ka1_a2, ks0_s1, ks0_s2, ks0_s3, ks1_a2, 
              ka2_s1, ks2_s1, ks1_s2, ks3_s2, ks2_s3, ks3_s3, ks1_bh, ks2_bh, 
              ks3_bh, ka2_a2, kbh_a2, ka2_bh, kbh_so, kbh_bo, kbo_bh, kfh_bh, 
              kbh_fh, kbh_fo, ks1_fh, ks2_fh, ks3_fh]

    n_t_step = len(airfx.U10)

    k_array = np.array([np.full(n_t_step, k, dtype=np.float64) if np.isscalar(k) else k for k in k_list])

    return k_array.T
