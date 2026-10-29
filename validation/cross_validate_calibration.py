"""Leave-one-group-out cross-validation of the two fitted parameters, RC_NIGHT and TAU_ASSIM.

Supplementary analysis for the paper: how stable the calibration is, and how well the calibration
targets are predicted when their own data are left out of the fit. Only calibration data are used;
the frozen test set is not touched. The scheme below was fixed before the run.

RC_NIGHT (validation/README section 4.2). The night targets come from three experiment groups; each
is left out in turn, RC_NIGHT is refitted on the others with the same grid and objective as
calibrate_night_resistance.py, and the left-out targets are scored at the refitted value:

    field box   wheat B4, B29, B50 (driven field-box exposures, TECDOC-1738 section 7)
    greenhouse  potato B57, wheat B60 (Balashov et al. 2011)
    review      B6 day/night exchange ratio, both crops
    and, as before, each crop's targets as a group (leave-one-crop-out).

TAU_ASSIM (section 4.5). It is fitted to the shape of the leaf OBT time course in the 14 field-box
exposures, the only time course in the calibration set, so a whole-campaign fold cannot refit it.
Instead each exposure is left out in turn (leave-one-exposure-out): TAU_ASSIM is refitted on the
other 13 with the grid and objective of calibrate_assimilate_residence.py, and the left-out
exposure's shape (2 h, 4 h, 1 d relative to 1 h) is predicted. A second split leaves out each time-of-day
class (dawn, day, dusk, night; scenario_field_box.time_class). The TAU_ASSIM grid is run at
RC_NIGHT = 4 and at every other RC_NIGHT value an RC_NIGHT fold selects.

OBT_DARK is not fitted but taken as the B15 central value, so B15 has no out-of-sample estimate and is
not part of this analysis; it is a sensitivity-analysis parameter.

Metrics: for RC_NIGHT, the log10 factor by which each left-out target misses its published range
(0 inside) at the fold value and at the full-data value, and the verdict counts over all calibration
benchmarks at the fold value. For TAU_ASSIM, the RMS of log10(model/observed shape ratio) of the
left-out exposures, against the in-sample RMS of the full-data fit.

Run from the repository root after validation/build_inputs.py (about 10 minutes):
    python validation/cross_validate_calibration.py
writes validation/results/cv_rc_night_grid.csv, cv_tau_assim_grid.csv, cv_rc_night.csv, cv_tau_assim.csv
and cv_stdout.txt.
"""
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
import calibrate_assimilate_residence as car  # noqa: E402
import calibrate_night_resistance as cnr  # noqa: E402
import run_validation as rv  # noqa: E402
import scenario_field_box as fb  # noqa: E402
from dmitte import para_constant  # noqa: E402

OUT = os.path.join(HERE, 'results')
RC_FOLDS = {
    'field box': {'wheat': ['B4', 'B29', 'B50']},
    'greenhouse': {'potato': ['B57'], 'wheat': ['B60']},
    'review': {'wheat': ['B6'], 'potato': ['B6']},
    'wheat targets': {'wheat': cnr.TARGETS['wheat']},
    'potato targets': {'potato': cnr.TARGETS['potato']},
}


def held_out(r, fold):
    return r.benchmark in fold.get(r.crop, [])


def miss_sq(comp, keep):
    t = comp[comp.apply(lambda r: r.benchmark in cnr.TARGETS[r.crop] and keep(r), axis=1)]
    return float((np.log10(t.factor_off.astype(float)) ** 2).sum())


def rc_night_grid():
    rv.run_wofost.plantModel = cnr.cached_plant_model
    bench = pd.read_csv(os.path.join(an.BENCH_DIR, 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index
    test_ids, _ = an.load_split()
    calib = [c for c in an.CHECKS if c['benchmark'] not in test_ids]
    rows = []
    default = para_constant.RC_NIGHT
    try:
        for f in cnr.GRID:
            para_constant.RC_NIGHT = f
            with tempfile.TemporaryDirectory() as out:
                rv.OUT = out
                rv.main()
                acute = pd.read_csv(os.path.join(out, 'acute_exposures.csv'), dtype={'exposure': str})
                chronic = pd.read_csv(os.path.join(out, 'chronic_exposure.csv'))
            scenario = fb.run({case: fb.RA_CASES[case] for case in (an.BOX_CASE, 'RAM+RB 8 s/m', 'RAM+RB 30 s/m')})
            comp = an.comparison(acute, chronic, bench, calib, scenario)
            comp['rc_night'] = f
            rows.append(comp)
            print(f'RC_NIGHT = {f} done', flush=True)
    finally:
        para_constant.RC_NIGHT = default
        rv.OUT = os.path.join(rv.ROOT, 'validation', 'results')
    res = pd.concat(rows, ignore_index=True)
    res['target'] = res.apply(lambda r: r.benchmark in cnr.TARGETS[r.crop], axis=1)
    return res


def rc_night_cv(grid):
    by_f = dict(tuple(grid.groupby('rc_night')))
    full = min(by_f, key=lambda f: miss_sq(by_f[f], lambda r: True))
    rows = []
    for name, fold in [('none (full data)', {})] + list(RC_FOLDS.items()):
        f = min(by_f, key=lambda v: miss_sq(by_f[v], lambda r: not held_out(r, fold)))
        g, g_full = by_f[f], by_f[full]
        out = g[g.apply(lambda r: held_out(r, fold), axis=1)]
        verdicts = g.verdict.value_counts()
        for _, r in (out.iterrows() if len(out) else [(None, None)]):
            row = dict(fold=name, rc_night=f, rc_night_full=full,
                       consistent=int(verdicts.get('consistent', 0)), close=int(verdicts.get('close', 0)),
                       off=int(verdicts.get('off', 0)))
            if r is not None:
                rf = g_full[(g_full.benchmark == r.benchmark) & (g_full.crop == r.crop)].iloc[0]
                row.update(benchmark=r.benchmark, crop=r.crop, bench_low=r.bench_low, bench_high=r.bench_high,
                           model_fold=r.model_median, log10_miss_fold=float(np.log10(r.factor_off)),
                           model_full=rf.model_median, log10_miss_full=float(np.log10(rf.factor_off)))
            rows.append(row)
    return pd.DataFrame(rows)


def tau_grid(rc_values):
    exps = pd.read_csv(os.path.join(fb.DATA, 'wheat_DE_field_box_1h.csv'))
    rows = []
    default_rc, default_tau = para_constant.RC_NIGHT, para_constant.TAU_ASSIM
    try:
        for rc in rc_values:
            para_constant.RC_NIGHT = rc
            for tau in car.GRID:
                para_constant.TAU_ASSIM = tau
                res = fb.run({an.BOX_CASE: fb.RA_CASES[an.BOX_CASE]})
                m = res.groupby('experiment').median(numeric_only=True).reset_index().merge(exps, on='experiment')
                for _, r in m.iterrows():
                    for t in car.TIMES:
                        rows.append(dict(rc_night=rc, tau=tau, experiment=r.experiment, start_hour=r.start_hour,
                                         time_class=fb.time_class(r.start_hour), time=t,
                                         miss=np.log10(r[t] / r.obt_1h) - np.log10(r[fb.OBS[t]] / r[fb.OBS['obt_1h']])))
                print(f'RC_NIGHT = {rc}, TAU_ASSIM = {tau} h done', flush=True)
    finally:
        para_constant.RC_NIGHT, para_constant.TAU_ASSIM = default_rc, default_tau
    return pd.DataFrame(rows).dropna(subset=['miss'])


def tau_cv(grid):
    rows = []
    for rc, g in grid.groupby('rc_night'):
        obj = g.groupby('tau').miss.apply(lambda x: float((x ** 2).sum()))
        full = obj.idxmin()
        for split, key in (('exposure', 'experiment'), ('time of day', 'time_class')):
            for left, out in g.groupby(key):
                fit = g[g[key] != left].groupby('tau').miss.apply(lambda x: float((x ** 2).sum()))
                tau = fit.idxmin()
                rows.append(dict(rc_night=rc, split=split, left_out=left, n_points=len(out) // len(car.GRID),
                                 tau_full=full, tau_fold=tau,
                                 rms_out_of_sample=float(np.sqrt((out[out.tau == tau].miss ** 2).mean())),
                                 rms_in_sample=float(np.sqrt((out[out.tau == full].miss ** 2).mean()))))
    return pd.DataFrame(rows)


def report(rc, tau):
    print('\n=== RC_NIGHT: leave-one-group-out ===')
    print(rc.round(3).to_string(index=False))
    o = rc.dropna(subset=['benchmark'])
    print('\nleft-out targets: RMS log10 miss at fold value '
          f'{np.sqrt((o.log10_miss_fold ** 2).mean()):.3f}, at full-data value {np.sqrt((o.log10_miss_full ** 2).mean()):.3f}; '
          f'within a factor 2: {(o.log10_miss_fold < np.log10(2)).sum()} / {len(o)} (full {(o.log10_miss_full < np.log10(2)).sum()})')
    print('\n=== TAU_ASSIM: leave-one-out ===')
    print(tau.round(3).to_string(index=False))
    for (r, split), g in tau.groupby(['rc_night', 'split']):
        w = g.n_points
        print(f'RC_NIGHT {r}, by {split}: TAU_ASSIM folds {sorted(g.tau_fold.unique().tolist())}; pooled RMS '
              f'out-of-sample {np.sqrt((g.rms_out_of_sample ** 2 * w).sum() / w.sum()):.3f}, '
              f'in-sample {np.sqrt((g.rms_in_sample ** 2 * w).sum() / w.sum()):.3f}')


def main():
    rc_grid = rc_night_grid()
    rc_grid.to_csv(os.path.join(OUT, 'cv_rc_night_grid.csv'), index=False)
    rc = rc_night_cv(rc_grid)
    rc.to_csv(os.path.join(OUT, 'cv_rc_night.csv'), index=False)
    rc_values = sorted({para_constant.RC_NIGHT, *rc.rc_night.unique().tolist()})
    tg = tau_grid(rc_values)
    tg.to_csv(os.path.join(OUT, 'cv_tau_assim_grid.csv'), index=False)
    tau = tau_cv(tg)
    tau.to_csv(os.path.join(OUT, 'cv_tau_assim.csv'), index=False)
    report(rc, tau)


if __name__ == '__main__':
    main()
