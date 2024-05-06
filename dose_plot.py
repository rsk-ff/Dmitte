# %%
from plot_config import *

fig_path = os.path.join(FIG_PATH, 'dose_HTO')

# %% HTO
fig_name = os.path.join(fig_path, 'Case1.svg')
normalize = np.array([16604148410, 7231123455, 13601099170, 69367370914, 3365703174, 13156294122, 2113049872, 4391362390, 1.39249055e+10]) / 6e9
dose1 = pd.read_csv(os.path.join(fig_path, 'dose1_HTO.csv'), index_col=0, header=0).div(normalize,  axis=0)
colors = [
    '#fb6a4a',  # 中等饱和度红色
    '#74c476',  # 中等饱和度绿色
    '#fd8d3c',  # 中等饱和度橙色
    '#6baed6',  # 中等饱和度蓝色
    '#9e9ac8',  # 中等饱和度紫色
    '#fddc6c'   # 中等饱和度黄色
]
headers = ['Inhalation', 'Potatoes', 'Cereals']
fig, ax = plt.subplots(figsize=(5, 3.5))
bottom = 0
for i, header in enumerate(headers):
    ax.bar(dose1.index, dose1[header], bottom=bottom, label=header, color=colors[i], edgecolor='black', linewidth=0.3)
    bottom += dose1[header]
ax.tick_params(axis='x', which='both', length=0, labelsize='small', labelrotation=0)
ax.tick_params(axis='y', labelsize='small')
ax.set_ylabel('剂量 / mSv')
ax.legend(fontsize='small')
fig.savefig(fig_name, format='svg')

# %%
fig_name = os.path.join(fig_path, 'Case2.svg')
normalize = np.array([57718744060, 54014491019, 70426584293, 1.19726e+11, 48899048840, 19320049307, 70426584293, 88830001774, 6.62318372e+10]) / 3e10
dose1 = pd.read_csv(os.path.join(fig_path, 'dose2_HTO.csv'), index_col=0, header=0).div(normalize,  axis=0)
colors = [
    '#fb6a4a',  # 中等饱和度红色
    '#74c476',  # 中等饱和度绿色
    '#fd8d3c',  # 中等饱和度橙色
    '#6baed6',  # 中等饱和度蓝色
    '#9e9ac8',  # 中等饱和度紫色
    '#fddc6c'   # 中等饱和度黄色
]
headers = ['Inhalation', 'Potatoes', 'Cereals']
fig, ax = plt.subplots(figsize=(5, 3.5))
bottom = 0
for i, header in enumerate(headers):
    ax.bar(dose1.index, dose1[header], bottom=bottom, label=header, color=colors[i], edgecolor='black', linewidth=0.3)
    bottom += dose1[header]
ax.tick_params(axis='x', which='both', length=0, labelsize='small', labelrotation=0)
ax.tick_params(axis='y', labelsize='small')
ax.set_ylabel('剂量 / mSv')
ax.legend(fontsize='small')
fig.savefig(fig_name, format='svg')

# %%
# %%
fig_name = os.path.join(fig_path, 'Case3.svg')
normalize = np.array([1.89193E+11, 4.08695E+11, 3.37113E+11, 69959754185, 3.71182E+11, 1.0964E+11, 3.59459E+11, 4.22023E+11, 4.08080453e+11]) / 3e11
dose1 = pd.read_csv(os.path.join(fig_path, 'dose3_HTO.csv'), index_col=0, header=0).div(normalize,  axis=0)
colors = [
    '#fb6a4a',  # 中等饱和度红色
    '#74c476',  # 中等饱和度绿色
    '#fd8d3c',  # 中等饱和度橙色
    '#6baed6',  # 中等饱和度蓝色
    '#9e9ac8',  # 中等饱和度紫色
    '#fddc6c'   # 中等饱和度黄色
]
headers = ['Inhalation', 'Potatoes', 'Cereals']
fig, ax = plt.subplots(figsize=(5, 3.5))
bottom = 0
for i, header in enumerate(headers):
    ax.bar(dose1.index, dose1[header], bottom=bottom, label=header, color=colors[i], edgecolor='black', linewidth=0.3)
    bottom += dose1[header]
ax.tick_params(axis='x', which='both', length=0, labelsize='small', labelrotation=0)
ax.tick_params(axis='y', labelsize='small')
ax.set_ylabel('剂量 / mSv')
ax.legend(fontsize='small')
fig.savefig(fig_name, format='svg')

# %%
