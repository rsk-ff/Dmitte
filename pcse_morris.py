# %%
import numpy as np
from SALib.sample.morris import sample
from SALib.analyze.morris import analyze

from dmitte import calc_para, para_constant, run_wofost
from dmitte.transfer_equotions import transfer_rates, solve_con
from dmitte.para_constant import LAMBDA_T_H
from openpyxl import load_workbook
from joblib import Parallel, delayed
from tqdm import tqdm
import uuid

from plot_config import *

def meteodata(start_date,end_date, params): 
    df = pd.read_excel('./data/meteo/meteo_usefor_cttm.xlsx')
    df['TIMESTAMP'] = pd.to_datetime(df['TIMESTAMP'])
    df.set_index('TIMESTAMP', inplace=True)

    full_range = pd.date_range(start=start_date, end=end_date, freq='D')
    df = df.reindex(full_range)

    df.fillna(method='bfill', inplace=True)  # 先使用后一天的数据填充
    df.fillna(method='ffill', inplace=True)  # 然后使用前一天的数据填充

    # 得到指定日期范围的数据
    filtered_df = df.loc[start_date:end_date]   
    # 使用线性插值填充缺失值
    hourly_df = filtered_df.resample('H').interpolate(method='linear')

    names = ['WS10m_avg', 'Ta_Avg', 'DR_Avg', 'Rain_Tot', 'RH_Avg', 'Pvapor_Avg', 'VWC_5cm_Avg', 'VWC_10cm_Avg', 'VWC_20cm_Avg']
    for i, name in enumerate(names):
        hourly_df[name] = float(params[i])

    return hourly_df

# 气孔、土壤、植物过程直接用 dmitte.calc_para, 只把气象换成 Morris 样本值 (见 model_function)


from pcse.fileinput import CABOFileReader, YAMLCropDataProvider,YAMLAgroManagementReader, ExcelWeatherDataProvider
from pcse.util import WOFOST72SiteDataProvider
from pcse.base import ParameterProvider
from pcse.models import Wofost72_WLP_FD
from datetime import datetime, timedelta

def plantModel(plant_type, plantmodel, date, params):

    file_name = f"result_{uuid.uuid4()}.xlsx"
    for row in range(13, 1195):
        ws[f'B{row}'].value = float(params[2]) * 86.4
        ws[f'C{row}'].value = float(params[1])
        ws[f'D{row}'].value = float(params[1])
        ws[f'E{row}'].value = float(params[5])
        ws[f'F{row}'].value = float(params[0])
        ws[f'G{row}'].value = float(params[3])
    
    wb.save(os.path.join(data_dir, 'meteo', file_name))

    # 设置气象参数
    weatherfile = os.path.join(data_dir, 'meteo', file_name)
    wdp = ExcelWeatherDataProvider(weatherfile)
    
    if plant_type.upper() == 'GRASS':
        # 设置植物生理参数
        cropfile = os.path.join(data_dir, 'crop', 'Grass.crop')
        cropd = CABOFileReader(cropfile)
        # 设置作物生长收获相关参数
        agromanagement_file = os.path.join (data_dir, 'agro', "Grass.agro")
        agromanagement = YAMLAgroManagementReader(agromanagement_file)
    else:
        cropd = YAMLCropDataProvider()
        # plantmodel_EXP = [crop, available_varieties, soilfile_name, agrofile_name]
        cropd.set_active_crop(plantmodel[0], plantmodel[1])
        agromanagement_file = os.path.join(data_dir, 'agro', plantmodel[3])
        agromanagement = YAMLAgroManagementReader(agromanagement_file)

# 初始化参数，运行wofost
    parameters = ParameterProvider(cropdata=cropd, soildata=soild, sitedata=sited)
    wofsim = Wofost72_WLP_FD(parameters, wdp, agromanagement)

    start_date = datetime.strptime(date[0], '%Y-%m-%d').date()
    end_date = datetime.strptime(date[1], '%Y-%m-%d').date()
    wofsim.run_till(end_date)
    # wofsim.run_till_terminate()
    # C:\Users\auror\miniconda3\envs\cttm\Lib\site-packages\pcse\conf\Wofost72_WLP_FD.conf
    df_results = pd.DataFrame(wofsim.get_output())
    
    # 输出模拟日期内所需的参数值
    first_date = df_results.iloc[0]['day']
    delta = first_date - start_date
    day_delta = delta.days
    df_results = df_results.set_index("day")
    df_results = df_results.iloc[abs(day_delta):]
    
    # 使用线性插值将每日的数据转变为每小时数据
    df_results.index = pd.to_datetime(df_results.index)
    df_resampled = df_results.resample('H').interpolate(method='linear')

    # df_results.to_csv(f'./data/pcse_output/{plant_type}_pcse_results_daily.cvs')
    # df_resampled.to_csv(f'./data/pcse_output/{plant_type}_pcse_results_hourly.cvs')
    os.remove(weatherfile)
    return df_resampled

def model_function(params):
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
    pcse_results = plantModel(plant_type, plantmodel, date, params)

    calc_para.meteodata = lambda s, e: meteodata(s, e, params)    # 气象参数取 Morris 样本值
    args = (pcse_results, start_date, end_date, plant_type, stability, HEG, x)
    airfx = calc_para.Calc_air(*args)
    soilfx = calc_para.Calc_soil(*args)
    plantfx = calc_para.Calc_plant(*args, drymatter)

    k_array = transfer_rates(char, plant_type, airfx, soilfx, plantfx)
    As = solve_con('HTO', 9.78728318e+08, k_array, LAMBDA_T_H)

    return As[2160]

#%%
# Problem definition
problem = {
    'num_vars': 9, # Number of parameters
    'names': ['WS10m_avg', 'Ta_Avg', 'DR_Avg', 'Rain_Tot', 'RH_Avg', 'Pvapor_Avg', 'VWC_5cm_Avg', 'VWC_10cm_Avg', 'VWC_20cm_Avg'],
    'bounds': [[0.1, 10], [-5, 38], [20, 300], [0, 5], [50, 95], [0.1, 3], [0.1, 0.4], [0.1, 0.4], [0.1, 0.4]]
}
# Generate samples
param_values = sample(problem, N=1000, num_levels=4, optimal_trajectories=None)
print(param_values.shape)

data_dir = os.path.join(os.getcwd(), "./data")
# 设置土壤相关参数
soilfile = os.path.join(data_dir, 'soil', 'ec2.soil')

soild = CABOFileReader(soilfile)

# 设置站点参数
sited = WOFOST72SiteDataProvider(WAV=10)

pcsefile = os.path.join(data_dir, 'meteo', 'meteo_usefor_pcse.xlsx')
wb = load_workbook(pcsefile)

ws = wb.active
def evaluate_model(params):
    return model_function(params)

# %%
num_cores = -1
Y = Parallel(n_jobs=num_cores)(delayed(evaluate_model)(params) for params in tqdm(param_values))


# %%
