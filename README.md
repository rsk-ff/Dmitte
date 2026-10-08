# Dmitte: Tritium Transfer

Compartment model of tritium (HT/HTO) transfer through air, soil and crops, coupled
with the WOFOST crop model (via [PCSE](https://pcse.readthedocs.io)) and a Gaussian
puff dispersion model, plus Morris/Sobol sensitivity analyses (SALib).

## Setup

Python 3.12 is what CI uses.

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt    # or requirements-dev.txt to also get pytest
```

Two pins matter:

- `pcse<6`: `dmitte/run_wofost.py` imports `WOFOST72SiteDataProvider` from `pcse.util`, which PCSE 6 moved.
- `pandas<3`: `dmitte/calc_para.py` uses `fillna(method=...)`, which pandas 3 removed.

## Input data (`./data`, not in git)

`data/` is gitignored, so a fresh clone needs it copied in. All paths are relative to
the **current working directory**, so run scripts from the repo root.

```
data/
├── meteo/
│   ├── meteo_usefor_cttm.xlsx   # hourly/daily station data for the tritium model
│   └── meteo_usefor_pcse.xlsx   # weather in PCSE ExcelWeatherDataProvider format
├── soil/
│   └── ec2.soil                 # CABO soil file (PCSE CABOFileReader)
├── crop/
│   └── Grass.crop               # CABO crop file, only for plant_type 'GRASS'
├── agro/                        # PCSE YAML agromanagement files
│   ├── potato.agro
│   ├── wheat.agro
│   ├── tomato.agro
│   ├── sugarbeet_calendar.agro
│   └── Grass.agro
└── output/                      # created by you; *_HTO.py scripts write CSVs here
```

`meteo_usefor_cttm.xlsx` is read by `dmitte.calc_para.meteodata()` and needs these columns:

| Column | Meaning |
|---|---|
| `TIMESTAMP` | date/time (parsed with `pd.to_datetime`) |
| `WS10m_avg` | 10 m wind speed, m/s |
| `Ta_Avg` | air temperature, °C |
| `DR_Avg` | solar radiation, W/m² |
| `Rain_Tot` | rainfall, mm |
| `RH_Avg` | relative humidity, % |
| `Pvapor_Avg` | vapour pressure, kPa |
| `P_Avg` | air pressure |
| `VWC_5cm_Avg`, `VWC_10cm_Avg`, `VWC_20cm_Avg` | soil volumetric water content at 5/10/20 cm |

Which soil and agro file a run uses comes from the `plantmodel` list in each script, e.g.
`['potato', 'Potato_702', 'ec2.soil', 'potato.agro']` = crop, variety, soil file, agro file.

Crop parameters for non-grass crops come from `YAMLCropDataProvider()`, which downloads
the [WOFOST crop parameters](https://github.com/ajwdewit/WOFOST_crop_parameters) on first
use, so the first run needs internet access.

`wofost_morris.py` reads `ScalarParametersOfWofost-Potential.xlsx` from the repo root
(committed). `atmospheric_dispersion.py` and `dose_plot.py` read reference CSVs from `figures/`.

## Running

```bash
python potato_HTO.py        # also cereals_HTO.py, tomatoes_HTO.py, iaea-case1-HTO.py
python dmitte_morris.py     # Morris sensitivity analysis
```

## Tests

```bash
pip install -r requirements-dev.txt
python -m pytest tests
```

The smoke tests in `tests/` run without `./data`; the one test that needs it is skipped
when the folder is missing. CI (`.github/workflows/tests.yml`) runs them on every push to
`dev` and on pull requests.
