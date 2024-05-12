# %%
from dmitte.guassian import gaussian_puff_model
from dmitte import run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con, transfer_rates_UFORTI, transfer_rates_baomi
from dmitte.para_constant import LAMBDA_T_H

from plot_config import *

# %% ---------------------------01 计算大气库室初始浓度----------------------------
downwind = 1e3
hight = np.linspace(0, 1000, 101)
t_serial = np.linspace(0, 7200, 121)

dd, hh, tt = np.meshgrid(downwind, hight, t_serial)

con_1 = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=2, stability=1, HEG=20, t_release=3600)
con_2 = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=5, stability=4, HEG=20, t_release=3600)
con_3 = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=2, stability=6, HEG=20, t_release=3600)

# 大气库室初始氚浓度 Bq/m2
init_con_1 = np.mean(con_1[:, :, :], axis=(0, 2)) * 1000 * 120 / 60
init_con_2 = np.mean(con_2[:, :, :], axis=(0, 2)) * 1000 * 120 / 60
init_con_3 = np.mean(con_3[:, :, :], axis=(0, 2)) * 1000 * 120 / 60
print('大气库室初始氚浓度 Bq·s·m-3:\n', init_con_1,'\n', init_con_2, '\n', init_con_3, '\n')

# %% ------------------------02 计算迁移率---------------------------
plantmodel = ['wheat', 'Winter_wheat_107', 'ec2.soil', 'wheat.agro']
plant_type = 'CEREAL'

start_date = '2021-06-20'
end_date = '2021-08-15'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 1
HEG = 20
x = 1000

char='HTO'
drymatter = 0.86
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates_baomi(char, plant_type, airfx, soilfx, plantfx)
k_array_ufotri = transfer_rates_UFORTI(char, plant_type, airfx, soilfx, plantfx)

columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
k_df = pd.DataFrame(k_array, columns=columns_k)
k_ufotri_df = pd.DataFrame(k_array_ufotri, columns=columns_k)

k_df.to_csv('data/output/transfer_rates.csv', index=False)

# %%
for k_name in columns_k:
    fig, ax = plt.subplots()
    ax.plot(k_df[k_name] / k_ufotri_df[k_name])
    ax.set_title(k_name)
    fig.savefig('figures/k_compare/'+k_name+'.png', dpi=600, format='png')

# %% --------------------------Case 1-------------------------------------
As = solve_con('HTO', 9.78728318e+08, k_array, LAMBDA_T_H)

# 产量 & 剂量
harvest_time = 49 * 24
yield_potato = pcse_results['TWSO'][harvest_time] / 10000 / drymatter
print(f'The yeild of potato is {yield_potato} kg/m2.')

# dose_adult_potato = (As[harvest_time, -2] / yield_potato) * (430 / 1000 * 12 * 30) * 1.8e-8

dose_adult_potato = (As[harvest_time, -1] / yield_potato) * (430 / 1000 * 12 * 30) * 4.2e-7 + \
                    (As[harvest_time, -2] / yield_potato) * (430 / 1000 * 12 * 30) * 1.8e-8
print(f'The dose for a adult from potato consumption is {dose_adult_potato} mSv')

# %%
figname = '土壤库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.set_yscale('log')
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 3:6], label=['Soil 1', 'Soil 2', 'Soil 3'])
ax.legend(fontsize="small")
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '土豆库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, -4:], label=['Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '所有库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 2:], label=['Air', 'Soil 1', 'Soil 2', 'Soil 3', 'Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(ncol=2, fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %% 如果需要可以把迁移率和浓度结果输出为 CSV 表格查看
# 转换为 Dataframe
columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
columns_A = ['A_a1', 'A_s0', 'A_a2', 'A_s1', 'A_s2', 'A_s3', 'A_bh', 'A_bo', 'A_fh', 'A_fo']
k_df = pd.DataFrame(k_array, columns=columns_k)
# A_df = pd.DataFrame(As, columns=columns_A)

# 输出到 CSV 文件
k_df.to_csv('data/output/transfer_rates.csv', index=False)
# A_df.to_csv('data/output/T_cons.csv', index=False)
# %%



# %% ------------------------Case 2---------------------------
plantmodel = ['wheat', 'Winter_wheat_107', 'ec2.soil', 'wheat.agro']
plant_type = 'CEREAL'

start_date = '2021-06-20'
end_date = '2021-08-15'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 4
HEG = 20
x = 1000

char='HTO'
drymatter = 0.86
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates_baomi(char, plant_type, airfx, soilfx, plantfx)
k_array_ufotri = transfer_rates_UFORTI(char, plant_type, airfx, soilfx, plantfx)

columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
k_df = pd.DataFrame(k_array, columns=columns_k)
k_ufotri_df = pd.DataFrame(k_array_ufotri, columns=columns_k)

k_df.to_csv('data/output/transfer_rates.csv', index=False)

# %%
for k_name in columns_k:
    fig, ax = plt.subplots()
    ax.plot(k_df[k_name] / k_ufotri_df[k_name])
    ax.set_title(k_name)
    fig.savefig('figures/k_compare/'+k_name+'.png', dpi=600, format='png')

# %% --------------------------Case 1-------------------------------------
As = solve_con('HTO', 1.15224004e+09, k_array, LAMBDA_T_H)

# 产量 & 剂量
harvest_time = 49 * 24
yield_potato = pcse_results['TWSO'][harvest_time] / 10000 / drymatter
print(f'The yeild of potato is {yield_potato} kg/m2.')

dose_adult_potato = (As[harvest_time, -1] / yield_potato) * (430 / 1000 * 12 * 30) * 4.2e-7 + \
                    (As[harvest_time, -2] / yield_potato) * (430 / 1000 * 12 * 30) * 1.8e-8
print(f'The dose for a adult from potato consumption is {dose_adult_potato} mSv')

# %%
figname = '土壤库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.set_yscale('log')
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 3:6], label=['Soil 1', 'Soil 2', 'Soil 3'])
ax.legend(fontsize="small")
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '土豆库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, -4:], label=['Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '所有库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 2:], label=['Air', 'Soil 1', 'Soil 2', 'Soil 3', 'Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(ncol=2, fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))




# %% ------------------------Case 3---------------------------
plantmodel = ['wheat', 'Winter_wheat_107', 'ec2.soil', 'wheat.agro']
plant_type = 'CEREAL'

start_date = '2021-06-20'
end_date = '2021-08-15'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 6
HEG = 20
x = 1000

char='HTO'
drymatter = 0.86
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates_baomi(char, plant_type, airfx, soilfx, plantfx)
# k_array_ufotri = transfer_rates_UFORTI(char, plant_type, airfx, soilfx, plantfx)

columns_k = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2',
              'ka2_s1', 'ks2_s1', 'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 
              'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh', 'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 
              'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
k_df = pd.DataFrame(k_array, columns=columns_k)
# k_ufotri_df = pd.DataFrame(k_array_ufotri, columns=columns_k)

k_df.to_csv('data/output/transfer_rates.csv', index=False)

# %%
for k_name in columns_k:
    fig, ax = plt.subplots()
    ax.plot(k_df[k_name] / k_ufotri_df[k_name])
    ax.set_title(k_name)
    fig.savefig('figures/k_compare/'+k_name+'.png', dpi=600, format='png')

# %% ---------------------------------------------------------------
As = solve_con('HTO', 5.7351502e+09, k_array, LAMBDA_T_H)

# 产量 & 剂量
harvest_time = 49 * 24
yield_potato = pcse_results['TWSO'][harvest_time] / 10000 / drymatter
print(f'The yeild of potato is {yield_potato} kg/m2.')

dose_adult_potato = (As[harvest_time, -1] / yield_potato) * (430 / 1000 * 12 * 30) * 4.2e-7 + \
                    (As[harvest_time, -2] / yield_potato) * (430 / 1000 * 12 * 30) * 1.8e-8
print(f'The dose for a adult from potato consumption is {dose_adult_potato} mSv')

# %%
figname = '土壤库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.set_yscale('log')
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 3:6], label=['Soil 1', 'Soil 2', 'Soil 3'])
ax.legend(fontsize="small")
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '土豆库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, -4:], label=['Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))

# %%
figname = '所有库室.svg'
fig, ax = plt.subplots(figsize=(4.5, 3))
ax.plot(np.arange(len(As[1:])) / 24, As[1:, 2:], label=['Air', 'Soil 1', 'Soil 2', 'Soil 3', 'Body HTO', 'Body OBT', 'Root HTO', 'Root OBT'])
ax.legend(ncol=2, fontsize="small")
ax.set_yscale('log')
# ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=100))
ax.set_xlabel('时间 / d')
ax.set_ylabel(r'氚浓度$\;/\;\mathrm{(Bq/m^{2})}$', fontdict={'family': 'SimSun'})
# fig.savefig(os.path.join(FIG_PATH, figname))