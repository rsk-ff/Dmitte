"""Validate the Dmitte HTO transfer model against published experiments.

Virtual experiments that mirror the published protocols (see
validation/benchmarks/literature_benchmarks.csv):

* acute: crops exposed for 1 h to air of constant HTO concentration (outdoor
  exposure box), then kept in clean air until harvest. Exposures are repeated
  weekly through the season, once at 11:00 (day) and once at 23:00 (night).
* chronic: constant air HTO from the start of the season to harvest
  (garden plot near a continuous source).

The model is used as is: transfer rates come from dmitte.transfer_equotions.transfer_rates
(with its default options) and the ODE
right-hand side is dmitte.transfer_equotions.transfer_equotions itself. Since
the system is linear with rates constant over each hour, the hourly step is the
matrix exponential of the generator, which is read off transfer_equotions
column by column. check_against_solve_con() confirms this reproduces
solve_con().

The only change in the experiments is the boundary condition the experiments
impose on the air compartment: held at the exposure concentration during the
exposure hour and at zero (clean air) otherwise.

Run from the repository root after validation/build_inputs.py:
    python validation/run_validation.py
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.linalg import expm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
warnings.filterwarnings('ignore')

from dmitte import calc_para, run_wofost  # noqa: E402
from dmitte.calc_para import Calc_air, Calc_plant, Calc_soil  # noqa: E402
from dmitte.para_constant import ALF_HTO, COMB_WATER, LAMBDA_T_H  # noqa: E402
from dmitte.transfer_equotions import (A_INDEX, K_INDEX, K_NAMES, solve_con, transfer_equotions, transfer_rates,  # noqa: E402
                                       water_balance_residual)

OUT = os.path.join(ROOT, 'validation', 'results')
LAT = 51.97                   # Wageningen
STABILITY, HEG, X = 4, 20, 1000   # irrelevant for clamped-air experiments, required by Calc_air
EXPOSURE_HOURS = (1, 2)           # field-box experiments are 1 h, chamber experiments (potato, wheat) 2 h

# compartment order of transfer_equotions (ba: leaf assimilate OBT)
A_A2, A_S1, A_BH, A_BO, A_FH, A_FO, A_BA = (A_INDEX[c] for c in ('a2', 's1', 'bh', 'bo', 'fh', 'fo', 'ba'))
N_COMP = len(A_INDEX)

VARIANTS = {'model': transfer_rates}    # name -> rate function; add entries to compare model options

CROPS = {
    'potato': dict(plant_type='ROOT_VEG', crop=('potato', 'Potato_701'), agro='potato_{y}.agro',
                   years=range(2004, 2009), window=('05-10', '09-05'),
                   exposures=('05-17', '08-23'), drymatter=0.21),   # drymatter as in potato_HTO.py
    'wheat': dict(plant_type='CEREAL', crop=('wheat', 'Winter_wheat_102'), agro='wheat_{y}.agro',
                  years=range(2005, 2009), window=('04-10', '08-05'),
                  exposures=('04-17', '07-17'), drymatter=0.86),    # drymatter as in cereals_HTO.py
}


# --------------------------------------------------------------------------- forcing
def diurnal_meteo(start, end, sm_hourly):
    """Hourly forcing with a diurnal cycle built from the daily station data.

    Temperature: cosine between TMIN and TMAX (max at 15:00). Radiation: daily
    total distributed with the sine of solar elevation. Vapour pressure and wind:
    daily value. RH: vapour pressure / saturation at the hourly temperature.
    Soil water: WOFOST SM of the same run.
    """
    daily = pd.read_excel('./data/meteo/meteo_usefor_cttm.xlsx')
    daily['TIMESTAMP'] = pd.to_datetime(daily['TIMESTAMP'])
    daily = daily.set_index('TIMESTAMP')
    idx = pd.date_range(start, end, freq='h')
    d = daily.reindex(idx.normalize()).set_index(idx)
    hour = idx.hour.values + 0.5
    doy = idx.dayofyear.values

    t_mean = (d['TMIN'] + d['TMAX']).values / 2
    amp = (d['TMAX'] - d['TMIN']).values / 2
    temp = t_mean + amp * np.cos(2 * np.pi * (hour - 15) / 24)

    decl = np.radians(23.45) * np.sin(2 * np.pi * (284 + doy) / 365)
    lat = np.radians(LAT)
    sin_elev = np.sin(lat) * np.sin(decl) + np.cos(lat) * np.cos(decl) * np.cos(np.radians(15 * (hour - 12)))
    sin_elev = np.clip(sin_elev, 0, None)
    day_sum = pd.Series(sin_elev, index=idx).groupby(idx.normalize()).transform('sum').values
    rad = d['IRRAD'].values * sin_elev / np.where(day_sum > 0, day_sum, 1) / 3.6     # kJ/m2/h -> W/m2

    vap = d['Pvapor_Avg'].values
    es = 0.6108 * np.exp(17.27 * temp / (temp + 237.3))
    rh = np.clip(vap / es, 0, 1) * 100
    sm = sm_hourly.reindex(idx).interpolate(limit_direction='both').values
    return pd.DataFrame({'WS10m_avg': d['WS10m_avg'].values, 'Ta_Avg': temp, 'DR_Avg': rad,
                         'Rain_Tot': d['Rain_Tot'].values, 'RH_Avg': rh, 'Pvapor_Avg': vap,
                         'P_Avg': d['P_Avg'].values, 'VWC_5cm_Avg': sm, 'VWC_10cm_Avg': sm,
                         'VWC_20cm_Avg': sm}, index=idx)


def model_state(crop, year, forcing):
    """Run WOFOST + calc_para for one crop season; forcing is 'diurnal' or 'daily' (stock)."""
    cfg = CROPS[crop]
    start, end = (f'{year}-{md}' for md in cfg['window'])
    plantmodel = [*cfg['crop'], 'ec2.soil', cfg['agro'].format(y=year)]
    pr = run_wofost.plantModel(cfg['plant_type'], plantmodel, [start, end])
    stock = calc_para.meteodata
    if forcing == 'diurnal':
        calc_para.meteodata = lambda s, e: diurnal_meteo(s, e, pr['SM'])
    try:
        args = (pr, start, end, cfg['plant_type'], STABILITY, HEG, X)
        air, soil, plant = Calc_air(*args), Calc_soil(*args), Calc_plant(*args, cfg['drymatter'])
    finally:
        calc_para.meteodata = stock
    return pr, air, soil, plant


# --------------------------------------------------------------------------- propagators
def generators(k_array):
    """Hourly 10x10 generators of transfer_equotions (linear in y)."""
    n = len(k_array)
    M = np.zeros((n, N_COMP, N_COMP))
    eye = np.eye(N_COMP)
    for i in range(n):
        for j in range(N_COMP):
            M[i, :, j] = transfer_equotions(0, eye[j], k_array[i], LAMBDA_T_H)
    return M


def clean_rates(k_array):
    """Model rates with non-finite entries set to 0 (they are counted and reported)."""
    k = np.array(k_array, dtype=float)
    bad = ~np.isfinite(k)
    k[bad] = 0.0
    return k, int(bad.sum()), int((k < 0).sum())


def check_against_solve_con(k_array):
    """Max relative difference between expm propagation and solve_con (model's own protocol)."""
    k = np.array(k_array, dtype=float)
    k0 = k.copy()
    k0[0, 0] = k0[0, 17] = 0
    As_ref = solve_con('HTO', 1.0, k, LAMBDA_T_H)
    y = np.zeros(N_COMP)
    y[A_A2] = 1.0
    out = [y]
    Ms = generators(k0)
    for i, M in enumerate(Ms):
        y = expm(M) @ y
        if i == 0:
            y[0] = y[A_A2] = 0
        out.append(y)
    As = np.array(out)
    scale = np.abs(As_ref).max(axis=0) + 1e-30
    return float(np.max(np.abs(As - As_ref) / scale))


# --------------------------------------------------------------------------- experiments
def season_arrays(air, soil, plant):
    """Water / dry-matter pools (kg m-2) used to turn inventories (Bq m-2) into concentrations."""
    w_air = air.ML * air.a_h / 1000
    dm_body = np.maximum(plant.tagp - plant.twso, 0) / 1e4
    dm_edible = plant.twso / 1e4
    return dict(w_air=w_air, w_plant=plant.plant_w, w_fruit=plant.friut_w,
                comb_body=dm_body * COMB_WATER, comb_edible=dm_edible * COMB_WATER,
                comb_leaf=plant.leaf_dm * COMB_WATER)


def leaf_obt_conc(y, pools, i):
    """Leaf OBT per combustion water: assimilate OBT in living leaves plus structural OBT of the shoot."""
    assim = y[A_BA] / pools['comb_leaf'][i] if pools['comb_leaf'][i] > 0 else 0.0
    return assim + safe_div(y[A_BO], pools['comb_body'][i])


def run_clamped(P, w_air, i0, i1):
    """Propagate with air held at concentration 1 (Bq per kg air moisture) for hours [i0, i1)."""
    n = len(P)
    ys = np.zeros((n + 1, N_COMP))
    y = np.zeros(N_COMP)
    for i in range(i0, n):
        y[A_A2] = w_air[i] if i < i1 else 0.0
        y = P[i] @ y
        y[A_A2] = 0.0
        ys[i + 1] = y
    return ys


def water_balance_summary(k, air, soil, plant, variant):
    """Uniform-specific-activity test per water compartment (median over hours, relative to outflow rate).

    wb_<c>:     median of residual / k_out (signed; > 0 spurious enrichment, < 0 spurious depletion)
    wbabs_<c>:  median of |residual - dlnW/dt| / k_out (part not explained by the storage change)
    """
    res = water_balance_residual(k, air, soil, plant, alpha=ALF_HTO)
    out = {}
    for c, r in res.items():
        with np.errstate(divide='ignore', invalid='ignore'):
            ok = np.isfinite(r['residual']) & (r['k_out'] > 0)
            rel = r['residual'][ok] / r['k_out'][ok]
            unexplained = np.abs(r['residual'][ok] - r['dlnW_dt'][ok]) / r['k_out'][ok]
        out[f'wb_{c}'] = float(np.median(rel)) if ok.any() else np.nan
        out[f'wbabs_{c}'] = float(np.median(unexplained)) if ok.any() else np.nan
    return out


def safe_div(a, b):
    return a / b if b > 0 else np.nan


def grain_phases(pr):
    """Anthesis time and the linear phase of storage-organ growth (WOFOST TWSO) of one season.

    The linear phase is where the daily TWSO increment is at least half its seasonal maximum
    (the "linear phase of grain growth" of Diabate & Strack 1997, Table 5).
    """
    dvs, twso = pr['DVS'], pr['TWSO'].fillna(0)
    anthesis = dvs.index[np.argmax(dvs.values >= 1)] if (dvs.values >= 1).any() else None
    daily = twso.resample('D').last().diff()
    linear = daily.index[daily.values >= 0.5 * np.nanmax(daily.values)] if daily.notna().any() else []
    start = linear[0] if len(linear) else None
    end = linear[-1] + pd.Timedelta(days=1) if len(linear) else None
    return anthesis, start, end


def phase_of(t, start, end):
    if start is None:
        return ''
    return 'before' if t < start else ('during' if t < end else 'after')


def acute_metrics(ys, i0, i1, pools, k):
    """Metrics of an exposure over hours [i0, i1); h0 is the end of exposure."""
    h0 = i1
    hv = len(ys) - 1
    c_leaf = lambda i: safe_div(ys[i, A_BH], pools['w_plant'][min(i, hv - 1)])
    c_fruit_tfwt = lambda i: safe_div(ys[i, A_FH], pools['w_fruit'][min(i, hv - 1)])
    obt_total_conc = safe_div(ys[h0, A_BO] + ys[h0, A_FO] + ys[h0, A_BA],
                              pools['comb_body'][i0] + pools['comb_edible'][i0])
    obt_leaf_conc = leaf_obt_conc(ys[h0], pools, min(h0, hv - 1))
    plant_t = ys[:, [A_BH, A_BO, A_FH, A_FO, A_BA]]
    i24 = min(h0 + 24, hv)
    return dict(
        R0_leaf=100 * c_leaf(h0),
        k_exchange=k[i0, K_INDEX['kbh_a2']],                    # leaf -> air
        k_out_leaf=k[i0, [K_INDEX[n] for n in ('kbh_a2', 'kbh_so', 'kbh_bo', 'kbh_fh', 'kbh_fo', 'kbh_ba')]].sum(),
        tfwt_drop_leaf=safe_div(c_leaf(h0), c_leaf(hv)),
        # edible-part TFWT keeps rising after h0 while it is fed from the shoot, so use its peak
        tfwt_drop_edible=safe_div(np.nanmax([c_fruit_tfwt(i) for i in range(h0, hv + 1)]), c_fruit_tfwt(hv)),
        TLI_edible=100 * safe_div(safe_div(ys[hv, A_FO], pools['comb_edible'][hv - 1]), c_leaf(h0)),
        obt_end_over_air=100 * obt_total_conc,
        obt_end_over_leaf=100 * safe_div(obt_leaf_conc, c_leaf(h0)),
        obt_share_24h=100 * safe_div(plant_t[i24, [1, 3, 4]].sum(), plant_t[i24].sum()),
        obt_formed_24h_per_leaf_h0=safe_div(ys[i24, A_BO] + ys[i24, A_FO] + ys[i24, A_BA], ys[h0, A_BH]),
        min_inventory=float(ys.min()),
    )


def chronic_metrics(ys, pools):
    hv = len(ys) - 1
    return dict(
        chronic_OBT_over_air=safe_div(ys[hv, A_FO], pools['comb_edible'][hv - 1]),
        chronic_TFWT_over_air=safe_div(ys[hv, A_FH], pools['w_fruit'][hv - 1]),
        chronic_leaf_TFWT_over_air=safe_div(ys[hv, A_BH], pools['w_plant'][hv - 1]),
        chronic_leaf_OBT_over_TFWT=safe_div(leaf_obt_conc(ys[hv], pools, hv - 1),
                                            safe_div(ys[hv, A_BH], pools['w_plant'][hv - 1])),
    )


def main():
    os.makedirs(OUT, exist_ok=True)
    acute_rows, chronic_rows, diag_rows = [], [], []
    for crop, cfg in CROPS.items():
        for year in cfg['years']:
            for forcing in ('diurnal', 'daily'):
                pr, air, soil, plant = model_state(crop, year, forcing)
                pools = season_arrays(air, soil, plant)
                anthesis, lin_start, lin_end = grain_phases(pr)
                idx = pd.date_range(f'{year}-{cfg["window"][0]}', f'{year}-{cfg["window"][1]}', freq='h')
                for vname, vfun in VARIANTS.items():
                    k_raw = vfun('HTO', cfg['plant_type'], air, soil, plant)
                    k, n_bad, n_neg = clean_rates(k_raw)
                    M = generators(k)
                    Mc = M.copy()
                    Mc[:, A_A2, :] = 0          # air is a boundary condition in the experiments
                    P = expm(Mc)
                    diag = dict(crop=crop, year=year, forcing=forcing, variant=vname,
                                nonfinite_rates=n_bad, negative_rates=n_neg,
                                negative_rate_names=','.join(sorted({c for c, m in zip(K_NAMES, (k < 0).any(axis=0)) if m})),
                                **water_balance_summary(k, air, soil, plant, vname))
                    if forcing == 'diurnal' and year == cfg['years'][0]:
                        try:
                            diag['expm_vs_solve_con_maxrel'] = check_against_solve_con(k)
                        except Exception as exc:     # solve_con itself fails on some rate sets
                            diag['expm_vs_solve_con_maxrel'] = np.nan
                            diag['solve_con_error'] = f'{type(exc).__name__}: {exc}'
                    diag_rows.append(diag)

                    days = pd.date_range(f'{year}-{cfg["exposures"][0]}', f'{year}-{cfg["exposures"][1]}', freq='7D')
                    for day in days:
                        for label, hour in (('day', 11), ('night', 23)):
                            for duration in EXPOSURE_HOURS:
                                if forcing == 'daily' and (label == 'night' or duration != 1):
                                    continue
                                t0 = day + pd.Timedelta(hours=hour)
                                i0 = idx.get_loc(t0)
                                ys = run_clamped(P, pools['w_air'], i0, i0 + duration)
                                m = acute_metrics(ys, i0, i0 + duration, pools, k)
                                daa = (t0 - anthesis) / pd.Timedelta(days=1) if anthesis is not None else np.nan
                                acute_rows.append(dict(crop=crop, year=year, forcing=forcing, variant=vname,
                                                       exposure=day.strftime('%m-%d'), time=label,
                                                       duration_h=duration, DVS=float(pr['DVS'].iloc[i0]),
                                                       days_after_anthesis=daa,
                                                       grain_phase=phase_of(t0, lin_start, lin_end), **m))
                    ys = run_clamped(P, pools['w_air'], 0, len(P))
                    chronic_rows.append(dict(crop=crop, year=year, forcing=forcing, variant=vname,
                                             **chronic_metrics(ys, pools)))
                print(crop, year, forcing, 'done', flush=True)

    acute = pd.DataFrame(acute_rows)
    chronic = pd.DataFrame(chronic_rows)
    diag = pd.DataFrame(diag_rows)
    acute.to_csv(os.path.join(OUT, 'acute_exposures.csv'), index=False)
    chronic.to_csv(os.path.join(OUT, 'chronic_exposure.csv'), index=False)
    diag.to_csv(os.path.join(OUT, 'diagnostics.csv'), index=False)
    print(diag.groupby(['crop', 'variant', 'forcing'])[['nonfinite_rates', 'negative_rates']].sum())


if __name__ == '__main__':
    main()
