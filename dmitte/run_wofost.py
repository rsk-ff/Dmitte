import os
import pandas as pd

from pcse.fileinput import CABOFileReader, YAMLCropDataProvider,YAMLAgroManagementReader, ExcelWeatherDataProvider
from pcse.util import WOFOST72SiteDataProvider
from pcse.base import ParameterProvider
from pcse.models import Wofost72_WLP_FD
from datetime import datetime, timedelta


#%%
# def plantModel(plant_type, plantmodel, date, meteo_file):
def plantModel(plant_type, plantmodel, date):
    data_dir = os.path.join(os.getcwd(), "./data")
    # 设置土壤相关参数
    soilfile = os.path.join(data_dir, 'soil', plantmodel[2])
    soild = CABOFileReader(soilfile)

    # 设置站点参数
    sited = WOFOST72SiteDataProvider(WAV=10)

    # 设置气象参数
    # weatherfile = os.path.join(data_dir, 'meteo', meteo_file)
    weatherfile = os.path.join(data_dir, 'meteo', 'meteo_usefor_pcse.xlsx')
    wdp = ExcelWeatherDataProvider(weatherfile)
    df = pd.DataFrame (wdp.export ())

    if plant_type.upper() == 'GRASS':
        # 设置植物生理参数
        cropfile = os.path.join(data_dir, 'crop', 'Grass.crop')
        cropd = CABOFileReader(cropfile)
        # 设置作物生长收获相关参数
        agromanagement_file = os.path.join (data_dir, 'agro', "Grass.agro")
        agromanagement = YAMLAgroManagementReader (agromanagement_file)
        # https://github.com/ajwdewit/WOFOST_crop_parameters
        print(agromanagement)
    else:
        cropd = YAMLCropDataProvider()
        # plantmodel_EXP = [crop, available_varieties, soilfile_name, agrofile_name]
        cropd.set_active_crop(plantmodel[0], plantmodel[1])
        agromanagement_file = os.path.join(data_dir, 'agro', plantmodel[3])
        agromanagement = YAMLAgroManagementReader(agromanagement_file)

# 初始化参数，运行wofost
    parameters = ParameterProvider (cropdata=cropd, soildata=soild, sitedata=sited)
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

    return df_resampled

#%%


if __name__ == '__main__':
    plantmodel = ['sugarbeet', 'Sugarbeet_601', 'ec2.soil', 'sugarbeet_calendar.agro']
    plant_type = 'GREEN_VEG'
    date = ['2021-04-20', '2021-08-15']    # 【开始日期，结束日期】
    # meteo_file = 'meteo_usefor_pcse.xlsx'
    pcse_results = plantModel(plant_type, plantmodel, date)
