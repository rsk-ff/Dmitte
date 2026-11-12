"""Sobol indices of the harvest OBT dose for releases on fixed calendar dates, with crop phenology uncertain.

analysis/sensitivity.py varies the tritium-transfer parameters with the crop calendar fixed, and looks at
releases placed relative to the simulated crop stage. In an assessment the release has a calendar date and the
crop stage on that date is itself uncertain. Here the eight tritium parameters kept by the Morris screening of
sensitivity.py are varied together with three phenology parameters of WOFOST:

    TSUM1_SCALE  0.9-1.1   [1]  factor on the temperature sum emergence -> anthesis / tuber initiation
    TSUM2_SCALE  0.9-1.1   [1]  factor on the temperature sum anthesis -> maturity
    START_SHIFT  -10-10 d  [0]  shift of the crop start date (sowing of wheat, emergence of potato)

(uniform; plausibility ranges that move anthesis of wheat by about +-5 d, see the printed table).

Outputs (per crop; log10 of the geometric mean over seasons of the OBT part of the dose per kg eaten, per unit
time-integrated air-moisture HTO; 12:00 releases):

    obt_mean    mean over releases at every hour while the crop stands
    obt_onset   release on the calendar date where the linear phase of storage-organ growth starts in the nominal run
    obt_mid     release on the nominal mid date of the linear phase
    obt_end     release on the nominal end date of the linear phase
    dose_mid    total dose (OBT + HTO in organ water) for the release of obt_mid

The same sample with the phenology parameters held at their nominal values gives the output spread due to the
tritium parameters alone. The model package dmitte/ is not modified (phenology.py reruns WOFOST with overrides).

Run from the repository root after validation/build_inputs.py (about 1 h on 4 cores):
    OMP_NUM_THREADS=1 python analysis/sensitivity_phenology.py [--n 256] [--jobs 4]
writes analysis/results/sobol_phenology.csv, phenology_shift.csv and figures/fig_b4_sensitivity_phenology.png.
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from SALib.analyze import sobol as sobol_an  # noqa: E402
from SALib.sample import sobol as sobol_sa  # noqa: E402

import common as c  # noqa: E402
import phenology  # noqa: E402
import sensitivity as sa  # noqa: E402
from dmitte import calc_para  # noqa: E402
from dmitte.calc_para import Calc_air, Calc_plant, Calc_soil  # noqa: E402

TRITIUM = [n for n in sa.PARAMS if n not in ('ALF_HTO', 'KBO_BH')]      # kept by the Morris screening
PHENO = {'TSUM1_SCALE': (0.9, 1.1, 1.0), 'TSUM2_SCALE': (0.9, 1.1, 1.0), 'START_SHIFT': (-10.0, 10.0, 0.0)}
NAMES = TRITIUM + list(PHENO)
BOUNDS = [sa.PARAMS[n][:2] for n in TRITIUM] + [v[:2] for v in PHENO.values()]
NOMINAL = [sa.PARAMS[n][2] for n in TRITIUM] + [v[2] for v in PHENO.values()]
OUTPUTS = ('obt_mean', 'obt_onset', 'obt_mid', 'obt_end', 'dose_mid')
FLOOR = 1e-18                                                           # mSv per kg; outputs below it are counted
_DATES = {}


def nominal_dates(crop):
    """Per season: calendar dates of the start, middle and end of the nominal linear phase (12:00)."""
    if crop not in _DATES:
        out = {}
        for year in c.CROPS[crop]['years']:
            _, start, end = c.rv.grain_phases(phenology.run(crop, year))
            noon = pd.Timedelta(hours=12)
            out[year] = dict(onset=start.normalize() + noon, mid=(start + (end - start) / 2).normalize() + noon,
                             end=(end - pd.Timedelta(days=1)).normalize() + noon)
        _DATES[crop] = out
    return _DATES[crop]


def model(crop, values):
    p = dict(zip(NAMES, values))
    full = {n: v[2] for n, v in sa.PARAMS.items()}
    full.update({n: p[n] for n in TRITIUM})
    rates = sa.apply([full[n] for n in sa.PARAMS])
    cfg = c.CROPS[crop]
    res = []
    try:
        for year, dates in nominal_dates(crop).items():
            start, end = (f'{year}-{md}' for md in cfg['window'])
            pr = phenology.run(crop, year, p['TSUM1_SCALE'], p['TSUM2_SCALE'], p['START_SHIFT'])
            forcing = c.rv.diurnal_meteo(start, end, pr['SM'])
            stock = calc_para.meteodata
            calc_para.meteodata = lambda s, e: forcing
            try:
                args = (pr, start, end, cfg['plant_type'], c.rv.STABILITY, c.rv.HEG, c.rv.X)
                state = (pr, Calc_air(*args), Calc_soil(*args), Calc_plant(*args, cfg['drymatter']))
                s = c.Season(crop, year, rates=rates, state=state)
            finally:
                calc_para.meteodata = stock
            names, R = s.release_response('dynamic')
            R = pd.DataFrame(R, columns=names, index=s.index[:s.n])
            dvs = pr['DVS'].reindex(R.index).values
            standing = R[np.isfinite(dvs) & (dvs >= 0) & (dvs < 2)]
            anth = s.anthesis
            res.append(dict(obt_mean=standing.dose_obt_per_kg.mean(),
                            obt_onset=R.dose_obt_per_kg.get(dates['onset'], np.nan),
                            obt_mid=R.dose_obt_per_kg.get(dates['mid'], np.nan),
                            obt_end=R.dose_obt_per_kg.get(dates['end'], np.nan),
                            dose_mid=R.dose_per_kg.get(dates['mid'], np.nan),
                            anthesis_doy=anth.dayofyear if anth is not None else np.nan,
                            lin_start_doy=s.lin_start.dayofyear, lin_end_doy=s.lin_end.dayofyear))
    finally:
        sa.apply([sa.PARAMS[n][2] for n in sa.PARAMS])
    r = pd.DataFrame(res)
    y = [float(np.log10(np.exp(np.log(r[o].clip(lower=FLOOR)).mean()))) for o in OUTPUTS]
    return y + [float(r[k].mean()) for k in ('anthesis_doy', 'lin_start_doy', 'lin_end_doy')]


def evaluate(crop, X, jobs):
    return np.array(Parallel(n_jobs=jobs, batch_size=8)(delayed(model)(crop, x) for x in X))


def spread(y):
    return float(10 ** (np.quantile(y, .95) - np.quantile(y, .05)))


def run(n, jobs, seed=3):
    pb = dict(num_vars=len(NAMES), names=NAMES, bounds=BOUNDS)
    X = sobol_sa.sample(pb, N=n, calc_second_order=False, seed=seed)
    A = X[::len(NAMES) + 2].copy()                    # the base sample: a plain random sample of all parameters
    A_fixed = A.copy()
    A_fixed[:, len(TRITIUM):] = [v[2] for v in PHENO.values()]
    rows, shifts = [], []
    for crop in c.CROPS:
        nom = model(crop, NOMINAL)
        Y = evaluate(crop, X, jobs)
        Yf = evaluate(crop, A_fixed, jobs)
        for k, o in enumerate(OUTPUTS):
            a = sobol_an.analyze(pb, Y[:, k], calc_second_order=False, seed=seed)
            for i, nm in enumerate(NAMES):
                rows.append(dict(crop=crop, output=o, parameter=nm, S1=a['S1'][i], S1_conf=a['S1_conf'][i],
                                 ST=a['ST'][i], ST_conf=a['ST_conf'][i]))
            ya = Y[::len(NAMES) + 2, k]
            rows.append(dict(crop=crop, output=o, parameter='(spread, all parameters)', ST=spread(ya)))
            rows.append(dict(crop=crop, output=o, parameter='(spread, tritium parameters only)', ST=spread(Yf[:, k])))
            rows.append(dict(crop=crop, output=o, parameter='(share at floor)', ST=float((Y[:, k] <= np.log10(FLOOR) + 1e-9).mean())))
            rows.append(dict(crop=crop, output=o, parameter='(nominal, log10)', ST=nom[k]))
        ya = Y[::len(NAMES) + 2]
        shifts.append(dict(crop=crop, **{f'{nm}_{q}': float(np.quantile(ya[:, len(OUTPUTS) + j], q / 100) - nom[len(OUTPUTS) + j])
                                         for j, nm in enumerate(('anthesis', 'lin_start', 'lin_end')) for q in (5, 95)}))
        print('sobol', crop, len(X) + len(A_fixed), 'runs done', flush=True)
    return pd.DataFrame(rows), pd.DataFrame(shifts)


def figure(s):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    colors = dict(zip(OUTPUTS, ['#52514e', '#eb6834', '#2a78d6', '#1baf7a', '#8a4fd3']))
    label = dict(obt_mean='OBT dose, season mean', obt_onset='OBT dose, release at nominal onset of filling',
                 obt_mid='OBT dose, release in mid filling', obt_end='OBT dose, release at end of filling',
                 dose_mid='total dose, release in mid filling')
    w = 0.16
    for ax, crop in zip(axes, c.CROPS):
        g = s[(s.crop == crop) & ~s.parameter.str.startswith('(')]
        for k, o in enumerate(OUTPUTS):
            h = g[g.output == o].set_index('parameter').reindex(NAMES)
            x = np.arange(len(NAMES)) + (k - 2) * w
            ax.bar(x, h.ST, w, color=colors[o], label=label[o])
            ax.errorbar(x, h.ST, yerr=h.ST_conf, fmt='none', ecolor='#0b0b0b', elinewidth=0.5)
        ax.axvline(len(TRITIUM) - 0.5, color='#b0aea8', linewidth=0.8, linestyle='--')
        ax.set_xticks(np.arange(len(NAMES)))
        ax.set_xticklabels(NAMES, rotation=35, ha='right', fontsize=8)
        ax.set_title(f'{crop}: Sobol total indices, 12:00 releases on fixed calendar dates', fontsize=10)
        for sp in ('top', 'right'):
            ax.spines[sp].set_visible(False)
    axes[0].set_ylabel('Sobol total index ST')
    axes[1].legend(fontsize=7.5, frameon=False, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(c.FIG, 'fig_b4_sensitivity_phenology.png'), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=256)
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    os.makedirs(c.FIG, exist_ok=True)
    pd.set_option('display.width', 220)
    s, shifts = run(a.n, a.jobs)
    s.to_csv(os.path.join(c.OUT, 'sobol_phenology.csv'), index=False)
    shifts.to_csv(os.path.join(c.OUT, 'phenology_shift.csv'), index=False)
    print('\nshift of the simulated dates over the sample (days, 5 % and 95 % relative to the nominal run):')
    print(shifts.round(1).to_string(index=False))
    print('\nSobol total index ST (log10 output); rows in brackets: 95 %/5 % ratio of the output over the sample,')
    print('share of runs at the numerical floor, and the nominal output:')
    print(s.pivot_table(index='parameter', columns=['crop', 'output'], values='ST', sort=False).round(3).to_string())
    figure(s)


if __name__ == '__main__':
    main()
