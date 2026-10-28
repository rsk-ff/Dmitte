"""Calibrate the night stomatal resistance factor RC_NIGHT on the calibration benchmarks.

Night (and RH > 98 %) stomatal resistance is RC_NIGHT * RCMIN (dmitte.calc_para.calc_st_res).
For each value on a grid the virtual experiments of run_validation.py are rerun and compared
with the calibration benchmarks (analyze.py). The objective is the sum of squared log10 factors
by which the model median misses the night targets (0 when inside the published range):

    wheat:  B29, B60 (night exchange rate), B50 (night 1 h uptake), B4 (night/day uptake), B6 (day/night exchange)
    potato: B57 (night exchange rate), B6

The value is chosen on both crops together; each crop is also calibrated alone and checked on
the other (leave-one-crop-out). All other calibration benchmarks are reported for every value,
so that changes outside the targets are visible. The frozen test benchmarks are not used.

Run from the repository root after validation/build_inputs.py:
    python validation/calibrate_night_resistance.py
writes validation/results/calibration_rc_night.csv.
"""
import functools
import os
import sys
import tempfile
import warnings

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

import analyze as an  # noqa: E402
import run_validation as rv  # noqa: E402
from dmitte import para_constant  # noqa: E402

GRID = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 15]
TARGETS = {'wheat': ['B29', 'B60', 'B50', 'B4', 'B6'], 'potato': ['B57', 'B6']}


def objective(comp, crops):
    t = comp[comp.apply(lambda r: r.crop in crops and r.benchmark in TARGETS[r.crop], axis=1)]
    return float((np.log10(t.factor_off.astype(float)) ** 2).sum())


_plant_model = rv.run_wofost.plantModel


@functools.lru_cache(maxsize=None)
def _wofost(plant_type, plantmodel, date):
    return _plant_model(plant_type, list(plantmodel), list(date))


def cached_plant_model(plant_type, plantmodel, date):
    '''WOFOST does not depend on RC_NIGHT: run each season once.'''
    return _wofost(plant_type, tuple(plantmodel), tuple(date)).copy()


def main():
    rv.run_wofost.plantModel = cached_plant_model
    bench = pd.read_csv(os.path.join(an.BENCH_DIR, 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index
    test_ids, _ = an.load_split()
    calib = [c for c in an.CHECKS if c['benchmark'] not in test_ids]
    assert not {b for v in TARGETS.values() for b in v} & test_ids

    rows = []
    default = para_constant.RC_NIGHT
    try:
        for f in GRID:
            para_constant.RC_NIGHT = f
            with tempfile.TemporaryDirectory() as out:
                rv.OUT = out
                rv.main()
                acute = pd.read_csv(os.path.join(out, 'acute_exposures.csv'), dtype={'exposure': str})
                chronic = pd.read_csv(os.path.join(out, 'chronic_exposure.csv'))
            comp = an.comparison(acute, chronic, bench, calib)
            comp['rc_night'] = f
            rows.append(comp)
            print(f'RC_NIGHT = {f}: objective both {objective(comp, ("wheat", "potato")):.3f}', flush=True)
    finally:
        para_constant.RC_NIGHT = default
    res = pd.concat(rows, ignore_index=True)
    res['target'] = res.apply(lambda r: r.benchmark in TARGETS[r.crop], axis=1)
    res.to_csv(os.path.join(rv.ROOT, 'validation', 'results', 'calibration_rc_night.csv'), index=False)
    report(res)


def report(res):
    by_f = {f: g for f, g in res.groupby('rc_night')}
    scores = pd.DataFrame({crops: {f: objective(g, crops) for f, g in by_f.items()}
                           for crops in (('wheat', 'potato'), ('wheat',), ('potato',))})
    scores.columns = ['both', 'wheat', 'potato']
    print('\nobjective (sum of squared log10 misses on the night targets):')
    print(scores.round(3).to_string())
    for crop, other in (('wheat', 'potato'), ('potato', 'wheat')):
        best = scores[crop].idxmin()
        print(f'calibrated on {crop} alone: RC_NIGHT = {best}; {other} objective there {scores.loc[best, other]:.3f} '
              f'(at its own best {scores[other].min():.3f})')
    print(f'calibrated on both: RC_NIGHT = {scores["both"].idxmin()}')
    t = res.pivot_table(index=['crop', 'benchmark', 'metric', 'target'], columns='rc_night', values='model_median',
                        sort=False)
    print('\nmodel medians by RC_NIGHT (all calibration benchmarks):')
    print(t.round(3).to_string())
    v = res.groupby(['rc_night', 'crop']).verdict.value_counts().unstack(fill_value=0)
    print('\nverdicts by RC_NIGHT:')
    print(v.to_string())


if __name__ == '__main__':
    main()
