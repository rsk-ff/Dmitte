import sys, os.path

import yaml
import numpy as np
import pandas as pd
from SALib.sample import saltelli
from SALib.analyze import sobol
import matplotlib.pyplot as plt
from IPython.display import display
pd.set_option("display.max_rows", None)
pd.set_option("display.max_colwidth", 250)

import pcse
from pcse.models import Wofost72_PP
from pcse.base import ParameterProvider
from pcse.db import NASAPowerWeatherDataProvider
from pcse.fileinput import YAMLCropDataProvider
from pcse.util import WOFOST72SiteDataProvider, DummySoilDataProvider
from progressbar import printProgressBar

print("This notebook was built with:")
print("python version: %s " % sys.version)
print("PCSE version: %s" %  pcse.__version__)


df = pd.read_excel("ScalarParametersOfWofost-Potential.xlsx")
# Define location, crop and season
latitude = 52.0
longitude = 5.0
crop_name = 'sugarbeet'
variety_name = 'Sugarbeet_601'
campaign_start_date = '2021-01-01'
emergence_date = "2021-03-31"
harvest_date = "2021-10-20"
max_duration = 300

# Here we define the agromanagement for sugar beet
agro_yaml = """
- {start}:
    CropCalendar:
        crop_name: {cname}
        variety_name: {vname}
        crop_start_date: {startdate}
        crop_start_type: emergence
        crop_end_date: {enddate}
        crop_end_type: harvest
        max_duration: {maxdur}
    TimedEvents: null
    StateEvents: null
""".format(cname=crop_name, vname=variety_name, 
           start=campaign_start_date, startdate=emergence_date, 
           enddate=harvest_date, maxdur=max_duration)
agro = yaml.safe_load(agro_yaml)
print(agro_yaml)

# Weather data for Netherlands
wdp = NASAPowerWeatherDataProvider(latitude=latitude, longitude=longitude)

# Parameter sets for crop, soil and site
# Standard crop parameter library
cropd = YAMLCropDataProvider()
# We don't need soil for potential production, so we use dummy values
soild = DummySoilDataProvider()
# Some site parameters
sited = WOFOST72SiteDataProvider(WAV=50, CO2=360.)

# Retrieve all parameters in the form of a single object. 
# In order to see all parameters for the selected crop already, we
# synchronise data provider cropd with the crop/variety: 
firstkey = list(agro[0])[0]
cropcalendar = agro[0][firstkey]['CropCalendar'] 
cropd.set_active_crop(cropcalendar['crop_name'], cropcalendar['variety_name'])
params = ParameterProvider(cropdata=cropd, sitedata=sited, soildata=soild)
# For each scalar parameter, determine a sensible interval 
problem_yaml = """
    num_vars: 5
    names: 
    - TSUM1
    - TSUM2
    - SPAN
    - Q10
    - TDWI
    bounds:
    - [500, 800]
    - [1200, 1600]
    - [28, 37]
    - [1.8, 2.2]
    - [0.4, 0.6]
"""
problem = yaml.safe_load(problem_yaml)

calc_second_order = True
nsamples = 50
paramsets = saltelli.sample(problem, nsamples, calc_second_order=calc_second_order)
print("We are going to do %s simulations" % len(paramsets))