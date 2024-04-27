# %%
from dmitte import run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con
from dmitte.para_constant import LAMBDA_T_H

import numpy as np
import pandas as pd

# %% 01 计算迁移率
plantmodel = ['sugarbeet', 'Sugarbeet_601', 'ec2.soil', 'sugarbeet_calendar.agro']
plant_type = 'LEAFY_PLANT'

start_date = '2021-04-20'
end_date = '2021-09-15'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 3
HEG = 30
x =500

char='HTO'
drymatter = 0.34
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates(char, plant_type, airfx, soilfx, plantfx)

# %% 02 求解库室浓度
inits = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
As = solve_con(inits, np.abs(k_array), LAMBDA_T_H)

# %% 如果需要可以把迁移率和浓度结果输出为 CSV 表格查看
# 转换为 Dataframe
columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
columns_A = ['A_a1', 'A_s0', 'A_a2', 'A_s1', 'A_s2', 'A_s3', 'A_bh', 'A_bo', 'A_fh', 'A_fo']
k_df = pd.DataFrame(k_array, columns=columns_k)
A_df = pd.DataFrame(As, columns=columns_A)

# 输出到 CSV 文件
k_df.to_csv('data/output/transfer_rates.csv', index=False)
A_df.to_csv('data/output/T_cons.csv', index=False)
# %%
