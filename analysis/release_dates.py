"""Dynamic versus constant transfer rates, and the release-date scan (thesis analysis B1).

For every hour of every WOFOST season (potato 2004-2008, wheat 2005-2008) a one-hour exposure to air-moisture
HTO of 1 Bq/kg is applied and the harvest endpoints of the edible part are computed (common.Season.
release_response): OBT per combustion water, TFWT, total tritium per fresh weight and adult ingestion dose
per kg eaten. Three rate variants of the same model are compared (common.VARIANTS):

* dynamic: the hourly rates of the crop-model-driven model;
* diurnal_only: crop growth frozen, i.e. for every hour of day the season-mean rates of that hour;
* constant: season-mean rates (no diurnal or seasonal variation), the analogue of a static model with fixed
  transfer parameters, but with this model's own parameter values, so that only the time variation differs;
* constant_norm: the constant response scaled, per season and endpoint, to the season mean of the dynamic
  model. Because the model is linear, the season mean of the acute responses is the harvest response to a
  constant air concentration, so this is a constant-rate model calibrated to reproduce the chronic case; it
  differs from the dynamic model only in how the response is distributed over release dates;
* release_frozen: the rates of the release hour held constant to harvest (a static model fed with the
  conditions at the time of release); computed for the 00:00 and 12:00 releases only.

Season statistics use the 00:00 and 12:00 releases of every day, so that all variants see the same releases.

Run from the repository root after validation/build_inputs.py:
    python analysis/release_dates.py
writes analysis/results/release_response.csv (daily 00:00 and 12:00 releases between emergence and maturity), release_summary.csv,
release_dates_stdout.txt (via tee) and analysis/results/figures/fig_b1_release_dates.png.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import common as c  # noqa: E402

ENDPOINT = 'dose_per_kg'
ENDPOINTS = ('dose_per_kg', 'dose_obt_per_kg', 'dose_hto_per_kg', 'edible_OBT')
COLORS = {'dynamic': '#2a78d6', 'diurnal_only': '#eb6834', 'constant': '#52514e', 'release_frozen': '#1baf7a',
          'constant_norm': '#8a4fd3'}
ORGAN_EVENT = {'potato': 'tuber initiation', 'wheat': 'anthesis'}      # what WOFOST DVS = 1 is for the crop


def run():
    rows = []
    for s in c.seasons():
        dvs = s.pr['DVS'].reindex(s.index[:s.n]).values
        for v in c.VARIANTS:
            names, r = s.release_response(v)
            for j in range(s.n):
                rows.append(dict(crop=s.crop, year=s.year, variant=v, time=s.index[j], hour=s.index[j].hour,
                                 DVS=dvs[j], days_after_anthesis=s.days_after_anthesis(j),
                                 **dict(zip(names, r[j]))))
        for j in range(s.n):
            if s.index[j].hour in (0, 12):
                names, r = s.release_frozen_response(j)
                rows.append(dict(crop=s.crop, year=s.year, variant='release_frozen', time=s.index[j],
                                 hour=s.index[j].hour, DVS=dvs[j], days_after_anthesis=s.days_after_anthesis(j),
                                 **dict(zip(names, r[0]))))
        print(s.crop, s.year, 'done', flush=True)
    d = pd.DataFrame(rows)
    d = d[d.DVS.notna() & (d.DVS >= 0) & (d.DVS < 2)]          # releases while the crop stands (emergence to maturity)
    return pd.concat([d, constant_norm(d)], ignore_index=True)


def constant_norm(d):
    """The constant response scaled per season and endpoint to the season mean (all hours) of the dynamic model."""
    names = [n for n in d.columns if n.startswith(('edible_', 'dose_'))]
    out = []
    for (crop, year), g in d.groupby(['crop', 'year']):
        dyn, con = g[g.variant == 'dynamic'], g[g.variant == 'constant'].copy()
        for n in names:
            con[n] *= dyn[n].mean() / con[n].mean()
        out.append(con.assign(variant='constant_norm'))
    return pd.concat(out)


def summary(d, endpoint):
    """Per crop and variant: season statistics of the response to releases at any hour while the crop stands."""
    out = []
    for (crop, v), g in d[d.hour.isin([0, 12])].groupby(['crop', 'variant'], sort=False):
        per_year = []
        for year, gy in g.groupby('year'):
            x = gy[endpoint]
            noon = gy[gy.hour == 12].set_index(gy[gy.hour == 12].time.dt.normalize())[endpoint]
            midnight = gy[gy.hour == 0].set_index(gy[gy.hour == 0].time.dt.normalize())[endpoint]
            per_year.append(dict(
                season_mean=x.mean(), season_max=x.max(),
                max_over_mean=x.max() / x.mean(),
                max_over_median=x.max() / x.median(),
                p95_over_p5=np.quantile(x, .95) / np.quantile(x, .05),
                noon_over_midnight=(noon / midnight).median(),
                peak_DAA=float(gy.loc[x.idxmax(), 'days_after_anthesis'])))
        out.append(dict(crop=crop, endpoint=endpoint, variant=v, **pd.DataFrame(per_year).median().to_dict()))
    s = pd.DataFrame(out)
    dyn = s[s.variant == 'dynamic'].set_index('crop')
    for col in ('season_mean', 'season_max'):
        s[f'{col}_vs_dynamic'] = s.apply(lambda r: r[col] / dyn.loc[r.crop, col], axis=1)
    return s


def ratio_by_stage(d, endpoint=ENDPOINT):
    """Dynamic / constant ratio of the response by days after anthesis (10 d bins), median over years and hours."""
    p = d.pivot_table(index=['crop', 'year', 'time', 'days_after_anthesis', 'hour'], columns='variant',
                      values=endpoint).reset_index()
    p['bin'] = (np.floor(p.days_after_anthesis / 10) * 10).astype(int)
    p['dyn_over_const'] = p.dynamic / p.constant
    p['dyn_over_norm'] = p.dynamic / p.constant_norm
    p['dyn_over_diurnal'] = p.dynamic / p.diurnal_only
    p['dyn_over_release'] = p.dynamic / p.release_frozen
    p['tod'] = np.where(p.hour.between(8, 17), 'day 08-17', 'night 18-07')
    return (p.groupby(['crop', 'bin', 'tod']).agg(dyn_over_const=('dyn_over_const', 'median'),
                                                  dyn_over_norm=('dyn_over_norm', 'median'),
                                                  dyn_over_diurnal=('dyn_over_diurnal', 'median'),
                                                  dyn_over_release=('dyn_over_release', 'median'))
            .unstack('tod'))


def figure(d):
    """Top: OBT part of the dose (the validated pathway). Bottom: total dose, including HTO in organ water."""
    rows = (('dose_obt_per_kg', 'OBT dose per kg eaten'), ('dose_per_kg', 'total dose (OBT + HTO) per kg eaten'))
    order = ('dynamic', 'constant_norm', 'constant', 'diurnal_only', 'release_frozen')
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6), sharex='col')
    for col, crop in enumerate(c.CROPS):
        g = d[d.crop == crop].copy()
        g['day'] = np.floor(g.days_after_anthesis)
        for row, (endpoint, ylabel) in enumerate(rows):
            ax = axes[row, col]
            for v in order:
                for hour, ls in ((12, '-'), (0, ':')):
                    if v != 'dynamic' and hour == 0:
                        continue
                    h = g[(g.variant == v) & (g.hour == hour)]
                    q = h.groupby('day')[endpoint].quantile([.1, .5, .9]).unstack()
                    q = q.where(q > 0)
                    lab = f'{c.VARIANT_LABEL[v]}' + (', 12:00 release' if v == 'dynamic' and hour == 12 else
                                                     ', 00:00 release' if v == 'dynamic' else '')
                    ax.plot(q.index, q[.5], ls, color=COLORS[v], linewidth=1.8 if v == 'dynamic' else 1.2, label=lab)
                    if v == 'dynamic':
                        ax.fill_between(q.index, q[.1], q[.9], color=COLORS[v], alpha=0.12, linewidth=0)
            ax.set_yscale('log')
            ax.set_ylim(1e-14, 3e-10)
            ax.set_ylabel(f'{ylabel},\nmSv per (Bq h kg$^{{-1}}$)')
            if row == 0:
                ax.set_title(f'{crop}: 1 h release, median over seasons (band: 10-90 %)', fontsize=10)
            else:
                ax.set_xlabel(f'release date, days after {ORGAN_EVENT[crop]} (WOFOST DVS = 1)')
            for s in ('top', 'right'):
                ax.spines[s].set_visible(False)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, fontsize=8, frameon=False)
    fig.tight_layout(rect=(0, 0.09, 1, 1))
    fig.savefig(os.path.join(c.FIG, 'fig_b1_release_dates.png'), dpi=150)


def main():
    os.makedirs(c.FIG, exist_ok=True)
    d = run()
    d[d.hour.isin([0, 12])].to_csv(os.path.join(c.OUT, 'release_response.csv'), index=False)
    s = pd.concat([summary(d, e) for e in ENDPOINTS])
    s.to_csv(os.path.join(c.OUT, 'release_summary.csv'), index=False)
    pd.set_option('display.width', 220)
    print('\nseason statistics of the harvest response to a 1 h release at any hour between emergence and maturity')
    print('(median over seasons; dose in mSv per kg eaten, OBT in Bq per kg combustion water, both per Bq h kg-1):')
    print(s.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    print('\nratio of the dynamic response to the frozen baselines by release date (days after anthesis, 10 d bins;')
    print(' dyn_over_release uses the 00:00 and 12:00 releases only):')
    print(ratio_by_stage(d).round(2).to_string())
    print('\nthe same for the OBT part of the dose:')
    print(ratio_by_stage(d, 'dose_obt_per_kg').round(2).to_string())
    print('\nshare of OBT in the dose of the dynamic model (median over seasons): season mean, and at the dose maximum')
    g = d[(d.variant == 'dynamic') & d.hour.isin([0, 12])]
    for crop, gc in g.groupby('crop'):
        per = [(gy.dose_obt_per_kg.mean() / gy.dose_per_kg.mean(),
                gy.loc[gy.dose_per_kg.idxmax(), 'dose_obt_per_kg'] / gy.dose_per_kg.max()) for _, gy in gc.groupby('year')]
        print(f'  {crop}: {np.median([a for a, _ in per]):.2f}, {np.median([b for _, b in per]):.2f}')
    figure(d)


if __name__ == '__main__':
    main()
