"""Dynamic versus constant transfer rates, and the release-date scan (thesis analysis B1).

For every hour of every WOFOST season (potato 2004-2008, wheat 2005-2008) a one-hour exposure to air-moisture
HTO of 1 Bq/kg is applied and the harvest endpoints of the edible part are computed (common.Season.
release_response): OBT per combustion water, TFWT, total tritium per fresh weight and adult ingestion dose
per kg eaten. Three rate variants of the same model are compared (common.VARIANTS):

* dynamic: the hourly rates of the crop-model-driven model;
* diurnal_only: crop growth frozen, i.e. for every hour of day the season-mean rates of that hour;
* constant: season-mean rates (no diurnal or seasonal variation), the analogue of a static model with fixed
  transfer parameters, but with this model's own parameter values, so that only the time variation differs;
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
COLORS = {'dynamic': '#2a78d6', 'diurnal_only': '#eb6834', 'constant': '#52514e', 'release_frozen': '#1baf7a'}


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
    return d[d.DVS.notna() & (d.DVS >= 0) & (d.DVS < 2)]       # releases while the crop stands (emergence to maturity)


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


def ratio_by_stage(d):
    """Dynamic / constant ratio of the response by days after anthesis (10 d bins), median over years and hours."""
    p = d.pivot_table(index=['crop', 'year', 'time', 'days_after_anthesis', 'hour'], columns='variant',
                      values=ENDPOINT).reset_index()
    p['bin'] = (np.floor(p.days_after_anthesis / 10) * 10).astype(int)
    p['dyn_over_const'] = p.dynamic / p.constant
    p['dyn_over_diurnal'] = p.dynamic / p.diurnal_only
    p['dyn_over_release'] = p.dynamic / p.release_frozen
    p['tod'] = np.where(p.hour.between(8, 17), 'day 08-17', 'night 18-07')
    return (p.groupby(['crop', 'bin', 'tod']).agg(dyn_over_const=('dyn_over_const', 'median'),
                                                  dyn_over_diurnal=('dyn_over_diurnal', 'median'),
                                                  dyn_over_release=('dyn_over_release', 'median'))
            .unstack('tod'))


def figure(d):
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    for ax, crop in zip(axes, c.CROPS):
        g = d[d.crop == crop].copy()
        g['day'] = np.floor(g.days_after_anthesis)
        for v in c.ALL_VARIANTS:
            for hour, ls in ((12, '-'), (0, ':')):
                if v != 'dynamic' and hour == 0:
                    continue
                h = g[(g.variant == v) & (g.hour == hour)]
                q = h.groupby('day')[ENDPOINT].quantile([.1, .5, .9]).unstack()
                lab = f'{c.VARIANT_LABEL[v]}' + (', 12:00 release' if v == 'dynamic' and hour == 12 else
                                                 ', 00:00 release' if v == 'dynamic' else '')
                ax.plot(q.index, q[.5], ls, color=COLORS[v], linewidth=1.4, label=lab)
                if v == 'dynamic':
                    ax.fill_between(q.index, q[.1], q[.9], color=COLORS[v], alpha=0.12, linewidth=0)
        ax.set_yscale('log')
        ax.set_xlabel('release date, days after anthesis (WOFOST DVS = 1)')
        ax.set_ylabel('dose per kg eaten, mSv per (Bq h kg-1)')
        ax.set_title(f'{crop}: 1 h release, median over seasons (band: 10-90 %)', fontsize=10)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(os.path.join(c.FIG, 'fig_b1_release_dates.png'), dpi=150)


def main():
    os.makedirs(c.FIG, exist_ok=True)
    d = run()
    d[d.hour.isin([0, 12])].to_csv(os.path.join(c.OUT, 'release_response.csv'), index=False)
    s = pd.concat([summary(d, e) for e in (ENDPOINT, 'edible_OBT')])
    s.to_csv(os.path.join(c.OUT, 'release_summary.csv'), index=False)
    pd.set_option('display.width', 220)
    print('\nseason statistics of the harvest response to a 1 h release at any hour between emergence and maturity')
    print('(median over seasons; dose in mSv per kg eaten, OBT in Bq per kg combustion water, both per Bq h kg-1):')
    print(s.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    print('\nratio of the dynamic response to the frozen baselines by release date (days after anthesis, 10 d bins;')
    print(' dyn_over_release uses the 00:00 and 12:00 releases only):')
    print(ratio_by_stage(d).round(2).to_string())
    figure(d)


if __name__ == '__main__':
    main()
