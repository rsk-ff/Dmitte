"""Assemble a complete ./data input set for the model from public data.

Everything comes from public GitHub repositories, pinned to a commit so the
inputs are reproducible:

* Weather: Wageningen (Haarweg) station 2004-2008, daily, from the PCSE example
  notebooks (Meteorology and Air Quality Group, Wageningen University).
* Soil:    EC2-medium CABO soil file (ec2.soil), same as the original runs.
* Crops:   WOFOST 7.2 crop parameters (potato, wheat, ...).

The tritium model additionally needs relative humidity, air pressure and soil
water content at 5/10/20 cm (meteo_usefor_cttm.xlsx). Those are derived:
RH from vapour pressure and mean temperature, pressure = standard atmosphere
(the station is 7 m a.s.l.), soil water content = WOFOST root-zone soil
moisture SM from the potato run of the same year (field capacity 0.272 when no
crop is in the field). See validation/README.md.

Run from the repository root:  python validation/build_inputs.py
"""
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')

PCSE_NOTEBOOKS = ('https://raw.githubusercontent.com/ajwdewit/pcse_notebooks/'
                  '4719f3b931137135b79993ecbdf46c29cd3f53eb/data/')
WOFOST72_PARAMS = ('https://raw.githubusercontent.com/ajwdewit/WOFOST_crop_parameters/'
                   'f0a6491f23685998fa2172b397ff959a3b5ea738/')
CROP_FILES = ['crops.yaml', 'barley.yaml', 'cassava.yaml', 'chickpea.yaml', 'cotton.yaml', 'cowpea.yaml',
              'fababean.yaml', 'groundnut.yaml', 'maize.yaml', 'millet.yaml', 'mungbean.yaml',
              'pigeonpea.yaml', 'potato.yaml', 'rapeseed.yaml', 'rice.yaml', 'seed_onion.yaml',
              'sorghum.yaml', 'soybean.yaml', 'sugarbeet.yaml', 'sugarcane.yaml', 'sunflower.yaml',
              'sweetpotato.yaml', 'tobacco.yaml', 'wheat.yaml']

YEARS = range(2004, 2009)
FIELD_CAPACITY = 0.272          # SMFCF of ec2.soil
STATION_PRESSURE_HPA = 1013.25

POTATO = ('potato', 'Potato_701')
WHEAT = ('wheat', 'Winter_wheat_102')

POTATO_AGRO = """Version: 1.0
AgroManagement:
- {y}-04-01:
    CropCalendar:
        crop_name: potato
        variety_name: {variety}
        crop_start_date: {y}-05-01
        crop_start_type: emergence
        crop_end_date: {y}-09-30
        crop_end_type: harvest
        max_duration: 300
    TimedEvents: null
    StateEvents: null
"""

WHEAT_AGRO = """Version: 1.0
AgroManagement:
- {y0}-10-01:
    CropCalendar:
        crop_name: wheat
        variety_name: {variety}
        crop_start_date: {y0}-10-20
        crop_start_type: sowing
        crop_end_date: {y}-08-15
        crop_end_type: harvest
        max_duration: 400
    TimedEvents: null
    StateEvents: null
"""


def download(url, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.exists(dest):
        return
    print('downloading', url)
    with urllib.request.urlopen(url, timeout=60) as r, open(dest, 'wb') as f:
        f.write(r.read())


def saturation_vapour_pressure_kpa(t_c):
    return 0.6108 * np.exp(17.27 * t_c / (t_c + 237.3))


def read_station_daily(path):
    """Daily Wageningen data from a PCSE ExcelWeatherDataProvider sheet."""
    raw = pd.read_excel(path, header=None, skiprows=12)
    raw = raw.iloc[:, :7]
    raw.columns = ['DAY', 'IRRAD', 'TMIN', 'TMAX', 'VAP', 'WIND', 'RAIN']
    raw['DAY'] = pd.to_datetime(raw['DAY'])
    df = raw.set_index('DAY').astype(float).replace(-999.0, np.nan)
    return df.interpolate(limit_direction='both')


def potato_soil_moisture(year):
    """Daily WOFOST root-zone soil moisture (SM) for the potato season of `year`."""
    sys.path.insert(0, ROOT)
    from dmitte import run_wofost
    res = run_wofost.plantModel('ROOT_VEG', [*POTATO, 'ec2.soil', f'potato_{year}.agro'],
                                [f'{year}-04-01', f'{year}-09-30'])
    return res['SM'].resample('D').mean()


def main():
    os.chdir(ROOT)
    download(PCSE_NOTEBOOKS + 'soil/ec2.soil', os.path.join(DATA, 'soil', 'ec2.soil'))
    download(PCSE_NOTEBOOKS + 'meteo/nl1.xlsx', os.path.join(DATA, 'meteo', 'meteo_usefor_pcse.xlsx'))
    for name in CROP_FILES:
        download(WOFOST72_PARAMS + name, os.path.join(DATA, 'crop', 'wofost72', name))

    agro_dir = os.path.join(DATA, 'agro')
    os.makedirs(agro_dir, exist_ok=True)
    for y in YEARS:
        with open(os.path.join(agro_dir, f'potato_{y}.agro'), 'w') as f:
            f.write(POTATO_AGRO.format(y=y, variety=POTATO[1]))
        if y > YEARS[0]:
            with open(os.path.join(agro_dir, f'wheat_{y}.agro'), 'w') as f:
                f.write(WHEAT_AGRO.format(y0=y - 1, y=y, variety=WHEAT[1]))

    daily = read_station_daily(os.path.join(DATA, 'meteo', 'meteo_usefor_pcse.xlsx'))
    sm = pd.concat([potato_soil_moisture(y) for y in YEARS])
    sm = sm.reindex(daily.index).fillna(FIELD_CAPACITY)

    t_mean = (daily['TMIN'] + daily['TMAX']) / 2
    rh = np.clip(daily['VAP'] / saturation_vapour_pressure_kpa(t_mean), 0.0, 1.0) * 100
    cttm = pd.DataFrame({
        'TIMESTAMP': daily.index,
        'WS10m_avg': daily['WIND'].values,
        'Ta_Avg': t_mean.values,
        'DR_Avg': (daily['IRRAD'] / 86.4).values,      # kJ m-2 d-1 -> W m-2 (daily mean)
        'Rain_Tot': daily['RAIN'].values,               # mm d-1; calc_para divides by 24
        'RH_Avg': rh.values,
        'Pvapor_Avg': daily['VAP'].values,              # kPa
        'P_Avg': STATION_PRESSURE_HPA,
        'VWC_5cm_Avg': sm.values,
        'VWC_10cm_Avg': sm.values,
        'VWC_20cm_Avg': sm.values,
        'TMIN': daily['TMIN'].values,                   # extra columns, used for the diurnal forcing
        'TMAX': daily['TMAX'].values,
        'IRRAD': daily['IRRAD'].values,
    })
    out = os.path.join(DATA, 'meteo', 'meteo_usefor_cttm.xlsx')
    cttm.to_excel(out, index=False)
    os.makedirs(os.path.join(DATA, 'output'), exist_ok=True)
    print('wrote', out, cttm.shape)


if __name__ == '__main__':
    main()
