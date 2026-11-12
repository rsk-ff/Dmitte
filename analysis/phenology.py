"""WOFOST runs with perturbed phenology for the sensitivity analysis (analysis/sensitivity_phenology.py).

run(crop, year, tsum1_scale, tsum2_scale, start_shift_days) repeats dmitte.run_wofost.plantModel (same files,
same PCSE configuration, same daily-to-hourly interpolation) with the temperature sums of the variety scaled
and the crop start date (sowing for wheat, emergence for potato) shifted. With (1, 1, 0) it returns exactly
the plantModel output (tests/test_analysis.py). The model package dmitte/ is not modified.
"""
import datetime as dt
import os

import pandas as pd
import yaml
from pcse.base import ParameterProvider
from pcse.engine import Engine
from pcse.fileinput import CABOFileReader, ExcelWeatherDataProvider, YAMLCropDataProvider
from pcse.util import WOFOST72SiteDataProvider

import common as c
from dmitte.run_wofost import WOFOST_CONF

DATA = os.path.join(c.ROOT, 'data')
_WDP = {}


def run(crop, year, tsum1_scale=1.0, tsum2_scale=1.0, start_shift_days=0):
    cfg = c.CROPS[crop]
    start, end = (f'{year}-{md}' for md in cfg['window'])
    cropd = YAMLCropDataProvider(fpath=os.path.join(DATA, 'crop', 'wofost72'))
    cropd.set_active_crop(*cfg['crop'])
    with open(os.path.join(DATA, 'agro', cfg['agro'].format(y=year))) as f:
        agro = yaml.safe_load(f)['AgroManagement']
    cal = list(agro[0].values())[0]['CropCalendar']
    cal['crop_start_date'] = cal['crop_start_date'] + dt.timedelta(days=int(round(start_shift_days)))
    if 'wdp' not in _WDP:
        _WDP['wdp'] = ExcelWeatherDataProvider(os.path.join(DATA, 'meteo', 'meteo_usefor_pcse.xlsx'))
    params = ParameterProvider(cropdata=cropd, soildata=CABOFileReader(os.path.join(DATA, 'soil', 'ec2.soil')),
                               sitedata=WOFOST72SiteDataProvider(WAV=10))
    sim = Engine(params, _WDP['wdp'], agro, config=WOFOST_CONF)
    # the crop parameters are loaded when the crop starts; override after reading the nominal values
    nominal = {n: cropd[n] for n in ('TSUM1', 'TSUM2')}
    params.set_override('TSUM1', nominal['TSUM1'] * tsum1_scale, check=False)
    params.set_override('TSUM2', nominal['TSUM2'] * tsum2_scale, check=False)
    sim.run_till(dt.datetime.strptime(end, '%Y-%m-%d').date())
    df = pd.DataFrame(sim.get_output()).set_index('day')
    df.index = pd.to_datetime(df.index)
    df = df.loc[pd.Timestamp(start):pd.Timestamp(end)]
    return df.resample('h').interpolate(method='linear')
