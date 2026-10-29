"""Morris screening followed by Sobol indices of the harvest response to a one-hour release (thesis analysis B3).

Model outputs (per crop, geometric mean over the WOFOST seasons of run_validation.py; log10 of the value):

    dose_mean        dose per kg eaten, mean over 1 h releases at every hour between emergence and maturity
    dose_fill_noon   dose per kg eaten, 1 h release at 12:00 in mid grain / tuber filling
    dose_fill_night  the same release at 00:00
    obt_mean         edible OBT per combustion water, mean over all release hours

(all per unit time-integrated air-moisture HTO, see common.Season.release_response). Mid filling is the middle
of the linear phase of storage-organ growth (validation/run_validation.grain_phases).

Parameters and ranges (uniform; nominal value in brackets):

    RC_NIGHT     3-6            [4]      night stomatal resistance factor (calibrated; handoff range)
    TAU_ASSIM    2-12 h         [6]      residence time of fresh leaf assimilate (calibrated)
    KBH_FH       half-life 20-200 h [80] plant water -> fruit water exchange (not identifiable from benchmarks)
    PERCOLATION  constant / darcy [constant]  soil percolation scheme (categorical: lower / upper half of [0, 1])
    OBT_DARK     0.1-0.5        [0.3]    night OBT formation relative to the daytime mean (TECDOC-1738 Sec. 7.1-7.2)
    ID_OBT       0.6-1.0        [0.8]    isotope discrimination of OBT formation (CSA N288.1 0.64-1.3; B20)
    ALF_HTO      1.08-1.12      [1.1]    liquid-vapour HTO fractionation
    DMF_LEAF     0.15-0.25      [0.2]    leaf dry-matter fraction
    RCMIN_SCALE  0.5-1.5        [1]      factor on the crop's minimum stomatal resistance RCMIN
    KBO_BH       0-ln2/10 d     [0]      turnover of structural OBT (set to 0 in validation section 4.5)

The chamber/box aerodynamic resistance of the validation is not a parameter of the field scenario and is not
included. Sobol indices are computed for the parameters whose Morris mu* reaches 10 % of the largest mu* of
any output of either crop; the others are held at their nominal values.

Run from the repository root after validation/build_inputs.py (about 1 h on 4 cores):
    python analysis/sensitivity.py [--morris-only] [--r 20] [--n 256] [--jobs 4]
writes analysis/results/morris.csv, sobol.csv and figures/fig_b3_sensitivity.png.
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from SALib.analyze import morris as morris_an, sobol as sobol_an  # noqa: E402
from SALib.sample import morris as morris_sa, sobol as sobol_sa  # noqa: E402

import common as c  # noqa: E402
from dmitte import calc_para, para_constant  # noqa: E402
from dmitte.calc_para import Calc_air, Calc_plant, Calc_soil  # noqa: E402

LN2 = np.log(2)
PARAMS = {
    'RC_NIGHT': (3.0, 6.0, 4.0),
    'TAU_ASSIM': (2.0, 12.0, 6.0),
    'KBH_FH_HALF_LIFE_H': (20.0, 200.0, 80.0),
    'PERCOLATION': (0.0, 1.0, 0.25),
    'OBT_DARK': (0.1, 0.5, 0.3),
    'ID_OBT': (0.6, 1.0, 0.8),
    'ALF_HTO': (1.08, 1.12, 1.1),
    'DMF_LEAF': (0.15, 0.25, 0.2),
    'RCMIN_SCALE': (0.5, 1.5, 1.0),
    'KBO_BH': (0.0, LN2 / 240, 0.0),
}
OUTPUTS = ('dose_mean', 'dose_fill_noon', 'dose_fill_night', 'obt_mean')
_CACHE = {}


def _seasons_inputs(crop):
    """WOFOST runs and site forcing of a crop, cached per worker process."""
    if crop not in _CACHE:
        cfg = c.CROPS[crop]
        out = []
        for year in cfg['years']:
            start, end = (f'{year}-{md}' for md in cfg['window'])
            pr = c.rv.run_wofost.plantModel(cfg['plant_type'], [*cfg['crop'], 'ec2.soil', cfg['agro'].format(y=year)],
                                            [start, end])
            out.append((year, start, end, pr, c.rv.diurnal_meteo(start, end, pr['SM'])))
        _CACHE[crop] = out
    return _CACHE[crop]


def apply(values):
    """Set the model parameters; returns the transfer_rates keyword arguments."""
    p = dict(zip(PARAMS, values))
    for name in ('RC_NIGHT', 'TAU_ASSIM', 'OBT_DARK', 'ID_OBT', 'ALF_HTO', 'DMF_LEAF', 'KBO_BH'):
        setattr(para_constant, name, float(p[name]))
    para_constant.RCMIN = [r * p['RCMIN_SCALE'] for r in NOMINAL_RCMIN]
    return dict(kbh_fh=LN2 / p['KBH_FH_HALF_LIFE_H'], percolation='darcy' if p['PERCOLATION'] >= 0.5 else 'constant')


NOMINAL_RCMIN = list(para_constant.RCMIN)
NOMINAL_CONST = {n: getattr(para_constant, n) for n in ('RC_NIGHT', 'TAU_ASSIM', 'OBT_DARK', 'ID_OBT', 'ALF_HTO',
                                                        'DMF_LEAF', 'KBO_BH')}


def model(crop, values):
    rates = apply(values)
    cfg = c.CROPS[crop]
    res = []
    try:
        for year, start, end, pr, forcing in _seasons_inputs(crop):
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
            mid = s.lin_start + (s.lin_end - s.lin_start) / 2
            day = mid.normalize()
            res.append(dict(dose_mean=standing.dose_per_kg.mean(),
                            dose_fill_noon=R.dose_per_kg.get(day + pd.Timedelta(hours=12), np.nan),
                            dose_fill_night=R.dose_per_kg.get(day + pd.Timedelta(hours=24), np.nan),
                            obt_mean=standing.edible_OBT.mean()))
    finally:
        apply([PARAMS[n][2] for n in PARAMS])
    r = pd.DataFrame(res)
    return [float(np.log10(np.exp(np.log(r[o]).mean()))) for o in OUTPUTS]


def evaluate(crop, X, jobs):
    Y = Parallel(n_jobs=jobs, batch_size=8)(delayed(model)(crop, x) for x in X)
    return np.array(Y)


def problem(names):
    return dict(num_vars=len(names), names=list(names), bounds=[PARAMS[n][:2] for n in names])


def run_morris(r, jobs, seed=1):
    pb = problem(PARAMS)
    X = morris_sa.sample(pb, N=r, num_levels=4, seed=seed)
    rows = []
    for crop in c.CROPS:
        Y = evaluate(crop, X, jobs)
        for k, o in enumerate(OUTPUTS):
            a = morris_an.analyze(pb, X, Y[:, k], num_levels=4, seed=seed)
            for i, n in enumerate(pb['names']):
                rows.append(dict(crop=crop, output=o, parameter=n, mu=a['mu'][i], mu_star=a['mu_star'][i],
                                 sigma=a['sigma'][i], mu_star_conf=a['mu_star_conf'][i]))
        print('morris', crop, len(X), 'runs done', flush=True)
    return pd.DataFrame(rows)


def screened(m, share=0.1):
    top = m.groupby(['crop', 'output']).mu_star.transform('max')
    keep = m[m.mu_star >= share * top].parameter.unique()
    return [n for n in PARAMS if n in keep]


def run_sobol(names, n, jobs, seed=2):
    pb = problem(names)
    Xs = sobol_sa.sample(pb, N=n, calc_second_order=False, seed=seed)
    fixed = {k: v[2] for k, v in PARAMS.items()}
    X = np.array([[dict(fixed, **dict(zip(names, x)))[k] for k in PARAMS] for x in Xs])
    rows = []
    for crop in c.CROPS:
        Y = evaluate(crop, X, jobs)
        for k, o in enumerate(OUTPUTS):
            a = sobol_an.analyze(pb, Y[:, k], calc_second_order=False, seed=seed)
            for i, nm in enumerate(names):
                rows.append(dict(crop=crop, output=o, parameter=nm, S1=a['S1'][i], S1_conf=a['S1_conf'][i],
                                 ST=a['ST'][i], ST_conf=a['ST_conf'][i]))
            rows.append(dict(crop=crop, output=o, parameter='(output spread)', S1=np.nan, S1_conf=np.nan,
                             ST=float(10 ** (np.quantile(Y[:, k], .95) - np.quantile(Y[:, k], .05))), ST_conf=np.nan))
        print('sobol', crop, len(X), 'runs done', flush=True)
    return pd.DataFrame(rows)


def figure(m, s):
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    colors = dict(zip(OUTPUTS, ['#2a78d6', '#eb6834', '#1baf7a', '#52514e']))
    for col, crop in enumerate(c.CROPS):
        ax = axes[0, col]
        g = m[m.crop == crop]
        for o in OUTPUTS:
            h = g[g.output == o]
            ax.scatter(h.mu_star, h.sigma, s=24, color=colors[o], label=o, alpha=0.85, linewidths=0)
            for _, r in h[h.mu_star >= 0.1 * g.mu_star.max()].iterrows():
                ax.annotate(r.parameter, (r.mu_star, r.sigma), fontsize=7, color='#52514e',
                            xytext=(3, 2), textcoords='offset points')
        ax.set_xlabel('Morris mu* (log10 output)')
        ax.set_ylabel('sigma')
        ax.set_title(f'{crop}: Morris screening', fontsize=10)
        ax = axes[1, col]
        if s is not None:
            g = s[(s.crop == crop) & (s.parameter != '(output spread)')]
            names = list(dict.fromkeys(g.parameter))
            w = 0.2
            for k, o in enumerate(OUTPUTS):
                h = g[g.output == o].set_index('parameter').reindex(names)
                ax.bar(np.arange(len(names)) + (k - 1.5) * w, h.ST, w, color=colors[o], alpha=0.85, label=o)
                ax.errorbar(np.arange(len(names)) + (k - 1.5) * w, h.ST, yerr=h.ST_conf, fmt='none',
                            ecolor='#0b0b0b', elinewidth=0.6)
            ax.set_xticks(np.arange(len(names)))
            ax.set_xticklabels(names, rotation=30, ha='right', fontsize=8)
            ax.set_ylabel('Sobol total index ST')
            ax.set_title(f'{crop}: Sobol total indices (screened parameters)', fontsize=10)
        for a in (axes[0, col], ax):
            for sp in ('top', 'right'):
                a.spines[sp].set_visible(False)
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(c.FIG, 'fig_b3_sensitivity.png'), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--morris-only', action='store_true')
    ap.add_argument('--r', type=int, default=20)
    ap.add_argument('--n', type=int, default=256)
    ap.add_argument('--jobs', type=int, default=4)
    a = ap.parse_args()
    os.makedirs(c.FIG, exist_ok=True)
    pd.set_option('display.width', 200)
    nominal = {crop: model(crop, [PARAMS[n][2] for n in PARAMS]) for crop in c.CROPS}
    print('nominal outputs (log10):', {k: np.round(v, 3).tolist() for k, v in nominal.items()})
    m = run_morris(a.r, a.jobs)
    m.to_csv(os.path.join(c.OUT, 'morris.csv'), index=False)
    print('\nMorris mu* (log10 output; 1 = a factor 10 per elementary effect):')
    print(m.pivot_table(index='parameter', columns=['crop', 'output'], values='mu_star').round(3).to_string())
    s = None
    if not a.morris_only:
        names = screened(m)
        print('\nscreened for Sobol:', names)
        s = run_sobol(names, a.n, a.jobs)
        s.to_csv(os.path.join(c.OUT, 'sobol.csv'), index=False)
        print('\nSobol first-order S1 / total ST (log10 output); "(output spread)" = 95 %/5 % ratio of the output:')
        print(s.pivot_table(index='parameter', columns=['crop', 'output'], values=['S1', 'ST']).round(3).to_string())
    figure(m, s)


if __name__ == '__main__':
    main()
