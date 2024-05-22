# %%
import numpy as np
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze

from dmitte import para_constant, run_wofost
from dmitte.calc_para import Calc_air,Calc_soil,Calc_plant
from dmitte.transfer_equotions import transfer_rates, solve_con, transfer_rates_UFORTI, transfer_rates_baomi
from dmitte.para_constant import LAMBDA_T_H
from openpyxl import load_workbook
from joblib import Parallel, delayed
from tqdm import tqdm
from adjustText import adjust_text

import uuid

from plot_config import *

def model_function(params):
    airfx.ML *= params[0]
    airfx.atm_h *= params[1]
    airfx.rainfall *= params[2]

    soilfx.VDSO *= params[3]
    soilfx.soil1w *= params[4]
    soilfx.soil2w *= params[5]
    soilfx.soil3w *= params[6]
    soilfx.ESOIL *= params[7]

    plantfx.VDPF *= params[8]
    plantfx.ETRM *= params[9]
    plantfx.plant_w *= params[10]
    plantfx.friut_w *= params[11]

    k_array = transfer_rates_baomi(char, plant_type, airfx, soilfx, plantfx)

    As = solve_con('HTO', 9.78728318e+08, k_array, LAMBDA_T_H)

    airfx.ML /= params[0]
    airfx.atm_h /= params[1]
    airfx.rainfall /= params[2]

    soilfx.VDSO /= params[3]
    soilfx.soil1w /= params[4]
    soilfx.soil2w /= params[5]
    soilfx.soil3w /= params[6]
    soilfx.ESOIL /= params[7]

    plantfx.VDPF /= params[8]
    plantfx.ETRM /= params[9]
    plantfx.plant_w /= params[10]
    plantfx.friut_w /= params[11]

    return As[2160]

def evaluate_model(params):
    return model_function(params)
#%%
problem = {
    'num_vars': 12, # Number of parameters
    'names': ['ML','A_h','Rain','VDSO','S1_w','S2_w','S3_w','ESOIL','VDPF','ETRM','P_w','F_w'],
    'bounds': [[0.9, 1.1]] * 12
}

param_values = sample(problem, N=1000, num_levels=4, optimal_trajectories=None)
print(param_values.shape)

plantmodel = ['potato', 'Potato_702', 'ec2.soil', 'potato.agro']
plant_type = 'ROOT_VEG'

start_date = '2021-06-20'
end_date = '2021-09-20'
date = [start_date, end_date]
stability = 1
HEG = 20
x = 1000

char='HTO'
drymatter = 0.21
pcse_results = run_wofost.plantModel(plant_type, plantmodel, date)

airfx =  Calc_air(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
soilfx = Calc_soil(pcse_results, start_date, end_date, plant_type, stability, HEG, x)
plantfx = Calc_plant(pcse_results, start_date, end_date, plant_type, stability, HEG, x,drymatter)
# %%
num_cores = -1
Y = Parallel(n_jobs=num_cores)(delayed(evaluate_model)(params) for params in tqdm(param_values))


#%%
figname = 'figures/morris/2.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 2], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
ax.set_xlim(-2, 50)
ax.set_ylim(-0.2, 4.2)
fig.savefig(figname, format='svg')
# %%
figname = 'figures/morris/3.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 3], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
ax.set_xlim(-0.4, 7.2)
ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')

# %%
figname = 'figures/morris/4.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 4], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')

# %%
figname = 'figures/morris/5.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 5], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')
# %%
# %%
figname = 'figures/morris/6.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 6], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')
# %%
# %%
figname = 'figures/morris/7.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 7], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')
# %%
# %%
figname = 'figures/morris/8.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 8], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')
# %%
# %%
figname = 'figures/morris/9.svg'
Si = analyze(problem, param_values, np.array(Y)[:, 9], print_to_console=True, num_levels=4, num_resamples=100)
df = Si.to_df()
df: pd.DataFrame
fig, ax = plt.subplots(figsize=(4, 3))
ax.scatter(df['mu_star'][1:-1], df['sigma'][1:-1], s=10)
labels = ['Ah','Rain','VDSO','S1w','S2w','S3w','ESOIL','VDPF','ETRM','Pw']
texts = []
for j, label in enumerate(labels[:]):
    texts.append(plt.text(df['mu_star'][j+1], df['sigma'][j+1], label, fontsize=6))
adjust_text(texts, arrowprops=dict(arrowstyle='->', color='black', shrinkA=0, shrinkB=0, lw=0.5))
ax.set_xlabel(r'$\mu*\;/\;\mathrm{(Bq/m^2)}$')
ax.set_ylabel(r'$\sigma\;/\;\mathrm{(Bq/m^2)}$')
# ax.set_xlim(-0.4, 7.2)
# ax.set_ylim(-0.05, 1.4)
fig.savefig(figname, format='svg')
# %%
