"""Calibrate the residence time TAU_ASSIM of fresh assimilate OBT in leaves.

Leaf water HTO enters the assimilate pool (ba) in proportion to gross assimilation; the pool empties
with residence time TAU_ASSIM into shoot OBT, fruit OBT and (respiration) back to leaf water
(dmitte.calc_para.calc_transOBT). TAU_ASSIM only shapes how leaf OBT evolves after an exposure, so it
is fitted to the shape of the measured time course in the wheat field box (TECDOC-1738 Fig. 27,
benchmarks/datasets/wheat_DE_field_box_1h.csv, calibration data): leaf OBT 2 h, 4 h and 1 d after the
start of exposure relative to its value at the end of exposure (1 h). The level at the end of
exposure (benchmarks B22, B23) is not fitted and stays a check of the formation rate.

For each value on a grid the 14 exposures are rerun with their measured conditions
(scenario_field_box.py, box RAM + RB of analyze.BOX_CASE); per exposure the model median over the
WOFOST seasons is compared. Objective: sum over exposures and times of
(log10 model ratio - log10 observed ratio)^2.

Run from the repository root after validation/build_inputs.py:
    python validation/calibrate_assimilate_residence.py
writes validation/results/calibration_tau_assim.csv.
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

import analyze as an  # noqa: E402
import scenario_field_box as fb  # noqa: E402
from dmitte import para_constant  # noqa: E402

GRID = [2, 4, 6, 9, 12, 18, 24, 36, 48]
TIMES = ['obt_2h', 'obt_4h', 'obt_1d']


def main():
    exps = pd.read_csv(os.path.join(fb.DATA, 'wheat_DE_field_box_1h.csv'))
    rows = []
    default = para_constant.TAU_ASSIM
    try:
        for tau in GRID:
            para_constant.TAU_ASSIM = tau
            res = fb.run({an.BOX_CASE: fb.RA_CASES[an.BOX_CASE]})
            m = res.groupby('experiment').median(numeric_only=True).reset_index().merge(exps, on='experiment')
            for _, r in m.iterrows():
                for t in ['obt_1h'] + TIMES:
                    rows.append(dict(tau=tau, experiment=r.experiment, start_hour=r.start_hour, time=t,
                                     model=r[t], observed=r[fb.OBS[t]], model_1h=r.obt_1h, observed_1h=r[fb.OBS['obt_1h']]))
            print(f'TAU_ASSIM = {tau} h done', flush=True)
    finally:
        para_constant.TAU_ASSIM = default
    d = pd.DataFrame(rows)
    d.to_csv(os.path.join(fb.OUT, 'calibration_tau_assim.csv'), index=False)
    report(d)


def report(d):
    shape = d[d.time.isin(TIMES)].dropna(subset=['observed']).copy()
    shape['miss'] = np.log10(shape.model / shape.model_1h) - np.log10(shape.observed / shape.observed_1h)
    obj = shape.groupby('tau').miss.apply(lambda x: float((x ** 2).sum()))
    print('\nobjective (shape of the leaf OBT time course):')
    print(obj.round(3).to_string())
    print(f'best TAU_ASSIM = {obj.idxmin()} h')
    ratio = d.assign(r=d.model / d.observed).groupby(['tau', 'time']).r.apply(lambda x: float(np.exp(np.log(x).mean())))
    print('\nmodel / observed leaf OBT, geometric mean over exposures:')
    print(ratio.unstack()[['obt_1h'] + TIMES].round(2).to_string())
    obs = d[d.tau == d.tau.iloc[0]].groupby('time')['observed'].median()
    mod = d.groupby(['tau', 'time'])['model'].median().unstack()
    print('\nmedian leaf OBT / leaf TFWT at end of exposure, %:')
    print(pd.concat([mod, obs.to_frame('observed').T])[['obt_1h'] + TIMES].round(3).to_string())


if __name__ == '__main__':
    main()
