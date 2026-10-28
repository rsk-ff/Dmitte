"""Why potato chronic OBT/TFWT is high in baomi after the symmetric leaf exchange.

Two checks on transfer_rates_baomi, potato 2004-2008, diurnal forcing:
1. variants: drop the friut_w floor in kfh_bh, slow OBT formation (kbh_bo / 10), scale root uptake
   (ks*_bh + ks*_fh) so that it equals transpiration ETRM, or all three; chronic and acute
   (07-05 11:00) metrics as in run_validation.py
2. HTO budget of the plant-water compartment in the 2004 chronic run: season-mean share of each inflow

Run from the repository root after validation/build_inputs.py:
    python validation/exp_attribution.py
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

import run_validation as rv  # noqa: E402
from dmitte.transfer_equotions import K_INDEX as K, transfer_rates_baomi  # noqa: E402


ROOT_RATES = [f'ks{n}_{c}' for n in '123' for c in ('bh', 'fh')]


def variants(soil, plant, base):
    water = {'1': soil.soil1w, '2': soil.soil2w, '3': soil.soil3w}
    uptake = sum(base[:, K[name]] * water[name[2]] for name in ROOT_RATES)
    root_eq_etrm = base.copy()
    with np.errstate(divide='ignore', invalid='ignore'):
        scale = np.where(uptake > 0, plant.ETRM / uptake, 0)
    for name in ROOT_RATES:
        root_eq_etrm[:, K[name]] *= scale
    no_floor = base.copy()
    with np.errstate(divide='ignore', invalid='ignore'):
        no_floor[:, K['kfh_bh']] = np.where(plant.friut_w > 0, base[:, K['kbh_fh']] * plant.plant_w / plant.friut_w, 0)
    slow_obt = base.copy()
    slow_obt[:, K['kbh_bo']] /= 10
    all_three = no_floor.copy()
    all_three[:, K['kbh_bo']] /= 10
    for name in ROOT_RATES:
        all_three[:, K[name]] *= scale
    return {'baomi': base, 'no friut_w floor': no_floor, 'kbh_bo / 10': slow_obt,
            'root uptake = ETRM': root_eq_etrm, 'all three': all_three}


def leaf_inflow_shares(k, air, plant, pools):
    """Season-mean share of each HTO inflow into plant water under chronic exposure (air held at 1)."""
    ys = rv.run_clamped(propagators(k), pools['w_air'], 0, len(k))[:-1].copy()
    ys[:, rv.A_A2] = pools['w_air'][:len(k)]
    inflow = {'vapour (VDPF/ML)': plant.VDPF / air.ML * 3600 * ys[:, rv.A_A2],
              'rain interception': air.rainfall / 9 / air.atm_h * ys[:, rv.A_A2],
              'roots (ks*_bh)': sum(k[:, K[f'ks{n}_bh']] * ys[:, j] for n, j in zip('123', (3, 4, 5))),
              'tuber water (kfh_bh)': k[:, K['kfh_bh']] * ys[:, rv.A_FH],
              'OBT turnover (kbo_bh)': k[:, K['kbo_bh']] * ys[:, rv.A_BO]}
    total = sum(inflow.values())
    with np.errstate(divide='ignore', invalid='ignore'):
        return {n: float(np.nanmean(v / total)) for n, v in inflow.items()}


def propagators(k):
    M = rv.generators(k)
    M[:, rv.A_A2, :] = 0
    return expm(M)


def main():
    rows = []
    for year in rv.CROPS['potato']['years']:
        pr, air, soil, plant = rv.model_state('potato', year, 'diurnal')
        pools = rv.season_arrays(air, soil, plant)
        idx = pd.date_range(f'{year}-05-10', f'{year}-09-05', freq='h')
        i0 = idx.get_loc(pd.Timestamp(f'{year}-07-05 11:00'))
        base = transfer_rates_baomi('HTO', 'ROOT_VEG', air, soil, plant)
        for name, k in variants(soil, plant, base).items():
            P = propagators(k)
            chronic = rv.chronic_metrics(rv.run_clamped(P, pools['w_air'], 0, len(P)), pools)
            acute = rv.acute_metrics(rv.run_clamped(P, pools['w_air'], i0, i0 + 1), i0, i0 + 1, pools, k)
            rows.append(dict(year=year, variant=name, chronic_OBT_over_air=chronic['chronic_OBT_over_air'],
                             chronic_TFWT_over_air=chronic['chronic_TFWT_over_air'],
                             chronic_leaf_TFWT_over_air=chronic['chronic_leaf_TFWT_over_air'],
                             TLI_edible=acute['TLI_edible'], obt_share_24h=acute['obt_share_24h']))
        if year == 2004:
            budget = {name: leaf_inflow_shares(k, air, plant, pools)
                      for name, k in variants(soil, plant, base).items() if name in ('baomi', 'all three')}
            roots = sum(base[:, K[name]] * {'1': soil.soil1w, '2': soil.soil2w, '3': soil.soil3w}[name[2]]
                        for name in ROOT_RATES)
            water = dict(root_uptake=float(np.mean(roots)), ETRM=float(np.mean(plant.ETRM)),
                         rain=float(np.mean(air.rainfall)),
                         gross_vapour_release=float(np.mean(plant.VDPF * plant.ROW / 1000 * 3600)))
    df = pd.DataFrame(rows)
    print(df.groupby('variant', sort=False).median(numeric_only=True).drop(columns='year').round(3).to_string())
    print('literature: OBT/air 0.72-1.14, TFWT/air 0.57-1.83, TLI 0.2-0.3 %, OBT share at 24 h 2-4 %')
    print('\nplant-water HTO inflow shares, 2004 chronic run, season mean:')
    print(pd.DataFrame(budget).T.round(2).to_string())
    print('mean water fluxes into / out of plant water, kg m-2 h-1:', {n: round(v, 3) for n, v in water.items()})


if __name__ == '__main__':
    main()
