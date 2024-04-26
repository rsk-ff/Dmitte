# %%
from dmitte.guassian import gaussian_puff_model
from dmitte import run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con
from dmitte.para_constant import LAMBDA_T_H

from plot_config import *


# %% -----------------------00 大气扩散----------------------
downwind = np.array([1e3, 3e3, 10e3, 30e3])
crosswind = 0
hight = np.linspace(0, 1000, 101)

dd, hh = np.meshgrid(downwind, hight)

t_serial = np.linspace(0, 25200, 421)

con_T = np.empty((len(t_serial), len(hight), len(downwind)))
for i, t in enumerate(t_serial):
    con_T[i] = gaussian_puff_model(dd, 0, hh, t, Q_total=3.7e15, UU=2, stability=1, HEG=20, t_release=3600)

# 不考虑库室时的空气积分浓度 Bq·s·m-3
int_con_T = np.mean(con_T[:, :21, :], axis=(0, 1)) * 25200
# 下风向各距离处的大气库室氚浓度（随时间变化）
con_T_1 = np.mean(con_T[:, :, 0], axis=1) * 1000
con_T_3 = np.mean(con_T[:, :, 1], axis=1) * 1000
con_T_10 = np.mean(con_T[:, :, 2], axis=1) * 1000
con_T_30 = np.mean(con_T[:, :, 3], axis=1) * 1000

# %%
figname = '高斯.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
titles = ['(a) 1 km', '(b) 3 km', '(c) 10 km', '(d) 30 km']

for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_T[:, :, i], axis=1) * 1000)
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])

# fig.savefig(os.path.join(FIG_PATH, figname))
# 释放期间大气库室氚浓度 Bq/m2
init_con = np.mean(con_T[:, :, :], axis=(0, 1)) * 1000 * 420 / 60

# %% ------------------------01 土豆---------------------------
plantmodel = ['potato', 'Potato_702', 'ec2.soil', 'potato.agro']
plant_type = 'ROOT_VEG'

start_date = '2021-06-20'
end_date = '2021-09-20'
date = [start_date, end_date]    # 【开始日期，结束日期】
stability = 1
HEG = 20
x = 1000

char='HTO'
drymatter = 0.21
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)

k_array = transfer_rates(char, plant_type, airfx, soilfx, plantfx)

# %%
As = solve_con('HTO', 6.15066907e+08, np.abs(k_array), LAMBDA_T_H)

# %%
harvest_time = 2160
yield_potato = pcse_results['TWSO'][harvest_time] / 10000 / drymatter
print(f'The yeild of potato is {yield_potato} kg/m2.')

dose_adult_potato = (As[harvest_time, -1] / yield_potato) * (200 / 1000 * 8 * 30) * 4.2e-7 + \
                    (As[harvest_time, -2] / yield_potato) * (200 / 1000 * 8 * 30) * 1.8e-8
print(f'The dose for a adult from potato consumption is {dose_adult_potato} mSv')
# %%
dose_adult_air = int_con_T[0]
print(f'The dose for a adult from potato consumption is {dose_adult_air} mSv')
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

# %%
