# %%
from dmitte.guassian import gaussian_puff_model
from dmitte import run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con, transfer_rates_UFORTI
from dmitte.para_constant import LAMBDA_T_H

from plot_config import *

# %% ------------------------02 计算迁移率---------------------------
plantmodel = ['sweetpotato', 'Sweetpotato_VanHeemst_1988', 'ec2.soil', 'tomato.agro']
plant_type = 'CEREAL'

start_date = '2021-06-20'
end_date = '2021-09-20'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 1
HEG = 20
x = 1000

char='HTO'
drymatter = 0.06
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates(char, plant_type, airfx, soilfx, plantfx)
k_array_ufotri = transfer_rates_UFORTI(char, plant_type, airfx, soilfx, plantfx)

columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
k_df = pd.DataFrame(k_array, columns=columns_k)
k_ufotri_df = pd.DataFrame(k_array_ufotri, columns=columns_k)

k_df.to_csv('data/output/transfer_rates.csv', index=False)

# %%
As = solve_con('HTO', 9.78728318e+08, k_array, LAMBDA_T_H)

# 产量 & 剂量
harvest_time = 28 * 24
yield_potato = pcse_results['TWSO'][harvest_time] / 10000 / drymatter
print(f'The yeild of potato is {yield_potato} kg/m2.')

dose_adult_potato = (As[harvest_time, -1] / 3) * (100 / 1000 * 2 * 30) * 4.2e-7 + \
                    (As[harvest_time, -2] / 3) * (100 / 1000 * 2 * 30) * 1.8e-8
print(f'The dose for a adult from potato consumption is {dose_adult_potato} mSv')
# %%
