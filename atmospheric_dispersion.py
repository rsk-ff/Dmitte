# %%
import pandas as pd

from dmitte.guassian import gaussian_puff_model, calc_sigmaxyz

from plot_config import *

fig_path = os.path.join(FIG_PATH, 'air_dis')
# %% -----------------------Gaussian----------------------
downwind = np.array([1e3, 3e3, 10e3, 30e3])
hight = np.linspace(0, 1000, 101)
t_serial = np.linspace(0, 25200, 421)

dd, hh, tt = np.meshgrid(downwind, hight, t_serial)

con_1 = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=2, stability=1, HEG=20, t_release=3600)
con_2_w = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=5, stability=4, HEG=20, t_release=3600, humidity=90, rain=15, temperature=20)
con_2_o = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=5, stability=4, HEG=20, t_release=3600, humidity=90, rain=0, temperature=20)
con_3 = gaussian_puff_model(dd, 0, hh, tt, Q_total=3.7e15, UU=2, stability=6, HEG=20, t_release=3600)

# 空气积分浓度 Bq·s·m-3
int_con_1 = np.mean(con_1[:2, :, :], axis=(0, 2)) * 25200
int_con_2_o = np.mean(con_2_o[:2, :, :], axis=(0, 2)) * 25200
int_con_2_w = np.mean(con_2_w[:2, :, :], axis=(0, 2)) * 25200
int_con_3 = np.mean(con_3[:2, :, :], axis=(0, 2)) * 25200
print('空气积分浓度 Bq·s·m-3:\n', int_con_1,'\n', int_con_2_w, '\n', int_con_3, '\n')

# %% ----------------------Plot_Gauss-------------------------------
titles = ['(a) 1 km', '(b) 3 km', '(c) 10 km', '(d) 30 km']

figname = 'Guass1.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_1[:, i, :], axis=0) * 1000)
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])
fig.savefig(os.path.join(fig_path, figname))

figname = 'Guass2.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_2_o[:, i, :], axis=0) * 1000, label='无修正')
    ax.plot(t_serial / 3600, np.mean(con_2_w[:, i, :], axis=0) * 1000, label='湿沉积')
    ax.legend(fontsize='small')
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])
fig.savefig(os.path.join(fig_path, figname))

figname = 'Guass3.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_3[:, i, :], axis=0) * 1000)
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])
fig.savefig(os.path.join(fig_path, figname))

# %%
figname = 'Guass123.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_1[:, i, :], axis=0) * 1000, label='Case 1')
    ax.plot(t_serial / 3600, np.mean(con_2_o[:, i, :], axis=0) * 1000, label='Case 2')
    ax.plot(t_serial / 3600, np.mean(con_3[:, i, :], axis=0) * 1000, label='Case 3')
    ax.legend()
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])
fig.savefig(os.path.join(fig_path, figname))

figname = 'Guass123w.svg'
fig, axs = plt.subplots(2, 2, figsize=(6, 4), sharex=True, constrained_layout=True)
for i, ax in enumerate(axs.flat):
    ax: axes.Axes
    ax.plot(t_serial / 3600, np.mean(con_1[:, i, :], axis=0) * 1000, label='Case 1')
    ax.plot(t_serial / 3600, np.mean(con_2_w[:, i, :], axis=0) * 1000, label='Case 2')
    ax.plot(t_serial / 3600, np.mean(con_3[:, i, :], axis=0) * 1000, label='Case 3')
    ax.legend()
    if i == 2 or i == 3:
        ax.set_xlabel('Time after release / h')
    if i == 0 or i == 2:
        ax.set_ylabel(r'Tritium / ($\rm{Bq/m^{2}}$)')
    ax.set_title(titles[i])
fig.savefig(os.path.join(fig_path, figname))

# %% --------------------Plot_Gauss_compare-------------------------------
figname = 'Guass_compare1.svg'
data_to_compare = pd.read_csv(os.path.join(fig_path, 'con1.csv'), header=0, index_col=0)
marker_list = ['D', 'o', '^', 'v', 's', 'p', '*', '+', 'x']
labels = data_to_compare.index
fig, ax = plt.subplots(figsize=(4, 2.5), constrained_layout=True)
ax.set_yscale('log')
ax.set_ylim(1e5, 1e12)
ax.set_xlim(-2.5, 35)
for i, row in enumerate(np.array(data_to_compare)):
    ax.plot([1, 3, 10, 30], row, c='gray', alpha=0.5, marker=marker_list[i], linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label=labels[i])
ax.plot([1, 3, 10, 30], int_con_1, marker='x', linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label='Dmitte')
ax.set_xlabel('下风向距离 / km')
ax.set_ylabel(r'空气积分浓度 / $\mathrm{Bq\cdot s\cdot m^{-3}}$', fontdict={'family': 'SimSun'})
ax.legend(fontsize='x-small', ncol=3, loc='lower left')
ax.get_yaxis().set_major_locator(LogLocator(base=10.0, numticks=10))
ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
fig.savefig(os.path.join(fig_path, figname))
print((int_con_1-np.array(data_to_compare)[1]) / np.array(data_to_compare)[1] * 100)

figname = 'Guass_compare2.svg'
data_to_compare = pd.read_csv(os.path.join(fig_path, 'con2.csv'), header=0, index_col=0)
marker_list = ['D', 'o', '^', 'v', 's', 'p', '*', '+', 'x']
labels = data_to_compare.index
fig, ax = plt.subplots(figsize=(4, 2.5), constrained_layout=True)
ax.set_yscale('log')
ax.set_ylim(1e5, 1e12)
ax.set_xlim(-2.5, 35)
for i, row in enumerate(np.array(data_to_compare)):
    ax.plot([1, 3, 10, 30], row, c='gray', alpha=0.5, marker=marker_list[i], linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label=labels[i])
ax.plot([1, 3, 10, 30], int_con_2_w, marker='x', linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label='Dmitte')
ax.set_xlabel('下风向距离 / km')
ax.set_ylabel(r'空气积分浓度 / $\mathrm{Bq\cdot s\cdot m^{-3}}$', fontdict={'family': 'SimSun'})
ax.legend(fontsize='x-small', ncol=3, loc='lower left')
ax.get_yaxis().set_major_locator(LogLocator(base=10.0, numticks=10))
ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
fig.savefig(os.path.join(fig_path, figname))
print((int_con_2_o-np.array(data_to_compare)[1]) / np.array(data_to_compare)[1] * 100)
print((int_con_2_w-np.array(data_to_compare)[1]) / np.array(data_to_compare)[1] * 100)

figname = 'Guass_compare3.svg'
data_to_compare = pd.read_csv(os.path.join(fig_path, 'con3.csv'), header=0, index_col=0)
marker_list = ['D', 'o', '^', 'v', 's', 'p', '*', '+', 'x']
labels = data_to_compare.index
fig, ax = plt.subplots(figsize=(4, 2.5), constrained_layout=True)
ax.set_yscale('log')
ax.set_ylim(1e5, 1e12)
ax.set_xlim(-2.5, 35)
for i, row in enumerate(np.array(data_to_compare)):
    ax.plot([1, 3, 10, 30], row, c='gray', alpha=0.5, marker=marker_list[i], linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label=labels[i])
ax.plot([1, 3, 10, 30], int_con_3, marker='x', linewidth = 0.5, markersize=4, fillstyle='none', markeredgewidth=0.8, label='Dmitte')
ax.set_xlabel('下风向距离 / km')
ax.set_ylabel(r'空气积分浓度 / $\mathrm{Bq\cdot s\cdot m^{-3}}$', fontdict={'family': 'SimSun'})
ax.legend(fontsize='x-small', ncol=3, loc='lower left')
ax.get_yaxis().set_major_locator(LogLocator(base=10.0, numticks=10))
ax.get_yaxis().set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1, numticks=10))
fig.savefig(os.path.join(fig_path, figname))
print((int_con_3-np.array(data_to_compare)[1]) / np.array(data_to_compare)[1] * 100)

# %% --------------------Plot_sigma_compare-------------------------------
print('sigma 值(m):\n', 
      calc_sigmaxyz(1, downwind)[1:], '\n',
      calc_sigmaxyz(4, downwind)[1:], '\n',
      calc_sigmaxyz(6, downwind)[1:], '\n',
)
# %%
