"""Compare the virtual-experiment results with the literature benchmarks.

Reads validation/results/*.csv (from run_validation.py) and
validation/benchmarks/literature_benchmarks.csv, writes
validation/results/benchmark_comparison.csv and the figures in
validation/results/figures/.

Verdict on the model median (diurnal forcing, all years and exposure dates):
  within the published range           -> consistent
  within a factor 3 of the range       -> close
  further away                         -> off (factor given)
  median <= 0, not finite or diverged  -> unphysical (|log10| > 30)
Benchmarks that only give a central value get a factor-1.5 tolerance.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')
FIG = os.path.join(RES, 'figures')

VARIANT_ORDER = ['UFOTRI', 'baomi', 'default']
COLORS = {'UFOTRI': '#2a78d6', 'baomi': '#eb6834', 'default': '#1baf7a'}
BAND = '#c9c8c2'
INK, INK2 = '#0b0b0b', '#52514e'

# (benchmark id, crop, model metric, subset of acute rows, how to aggregate)
CHECKS = [
    ('B1', 'potato', 'R0_leaf', 'all'), ('B1', 'wheat', 'R0_leaf', 'all'),
    ('B2', 'potato', 'R0_leaf', 'day'), ('B2', 'wheat', 'R0_leaf', 'day'),
    ('B3', 'potato', 'R0_leaf', 'night'), ('B3', 'wheat', 'R0_leaf', 'night'),
    ('B4', 'potato', 'R0_night_over_day', 'pair'), ('B4', 'wheat', 'R0_night_over_day', 'pair'),
    ('B5', 'potato', 'k_exchange', 'day'), ('B5', 'wheat', 'k_exchange', 'day'),
    ('B7', 'potato', 'k_exchange', 'day'), ('B7', 'wheat', 'k_exchange', 'day'),
    ('B8', 'potato', 'k_exchange', 'night'), ('B8', 'wheat', 'k_exchange', 'night'),
    ('B6', 'potato', 'k_day_over_night', 'pair'), ('B6', 'wheat', 'k_day_over_night', 'pair'),
    ('B9', 'potato', 'tfwt_drop_leaf', 'all'), ('B9', 'wheat', 'tfwt_drop_leaf', 'all'),
    ('B10', 'potato', 'tfwt_drop_edible', 'all'),
    ('B12', 'potato', 'TLI_edible', 'all'), ('B13', 'potato', 'TLI_edible', 'all'),
    ('B15', 'potato', 'obt_night_over_day', 'pair'), ('B15', 'wheat', 'obt_night_over_day', 'pair'),
    ('B16', 'potato', 'obt_end_over_air', 'day'), ('B16', 'wheat', 'obt_end_over_air', 'day'),
    ('B17', 'potato', 'obt_share_24h', 'all'), ('B17', 'wheat', 'obt_share_24h', 'all'),
    ('B18', 'potato', 'chronic_OBT_over_air', 'chronic'), ('B18', 'wheat', 'chronic_OBT_over_air', 'chronic'),
    ('B19', 'potato', 'chronic_TFWT_over_air', 'chronic'), ('B19', 'wheat', 'chronic_TFWT_over_air', 'chronic'),
    ('B20', 'potato', 'chronic_leaf_OBT_over_TFWT', 'chronic'),
    ('B20', 'wheat', 'chronic_leaf_OBT_over_TFWT', 'chronic'),
]
ONE_SIDED = {'B10': (0, None)}   # B10: "up to 1.3e4" (upper bound)


def bench_range(row):
    lo, mid, hi = row['low'], row['central'], row['high']
    if row['id'] in ONE_SIDED:
        return 0.0, (hi if pd.notna(hi) else mid)
    if pd.notna(lo) and pd.notna(hi):
        return lo, hi
    return mid / 1.5, mid * 1.5


def verdict(median, lo, hi):
    if not np.isfinite(median) or median <= 0 or not 1e-30 < median < 1e30:
        return 'unphysical', np.nan
    if lo <= median <= hi:
        return 'consistent', 1.0
    factor = median / hi if median > hi else lo / median
    return ('close' if factor <= 3 else 'off'), factor


def paired(acute, metric):
    """Night/day ratios for the same crop, year, variant and exposure date."""
    d = acute[acute.forcing == 'diurnal'].pivot_table(index=['crop', 'year', 'variant', 'exposure'],
                                                       columns='time', values=metric)
    return d['night'] / d['day']


def model_values(acute, chronic, crop, variant, metric, subset):
    a = acute[(acute.forcing == 'diurnal') & (acute.crop == crop) & (acute.variant == variant)]
    if subset == 'chronic':
        c = chronic[(chronic.forcing == 'diurnal') & (chronic.crop == crop) & (chronic.variant == variant)]
        return c[metric].values
    if subset == 'pair':
        sub = acute[(acute.crop == crop) & (acute.variant == variant)]
        if metric == 'R0_night_over_day':
            return paired(sub, 'R0_leaf').values
        if metric == 'k_day_over_night':
            return 1 / paired(sub, 'k_exchange').values
        if metric == 'obt_night_over_day':
            return paired(sub, 'obt_formed_24h_per_leaf_h0').values
    if subset in ('day', 'night'):
        a = a[a.time == subset]
    return a[metric].values


def fmt(x):
    return f'{x:.3g}' if np.isfinite(x) else 'nan'


def comparison(acute, chronic, bench):
    rows = []
    for bid, crop, metric, subset in CHECKS:
        b = bench.loc[bid]
        lo, hi = bench_range(b)
        for variant in VARIANT_ORDER:
            v = np.asarray(model_values(acute, chronic, crop, variant, metric, subset), dtype=float)
            v = v[np.isfinite(v)]
            med = float(np.median(v)) if len(v) else np.nan
            p5, p95 = (np.percentile(v, [5, 95]) if len(v) else (np.nan, np.nan))
            verd, factor = verdict(med, lo, hi)
            rows.append(dict(benchmark=bid, crop=crop, variant=variant, metric=metric,
                             benchmark_quantity=b['quantity'], benchmark_crop=b['crop'],
                             bench_low=lo, bench_high=hi, unit=b['unit'],
                             model_median=med, model_p5=p5, model_p95=p95,
                             n=len(v), verdict=verd, factor_off=factor))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- figures
def style(ax):
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.yaxis.grid(True, color='#e6e5e0', linewidth=0.6)
    ax.set_axisbelow(True)


def band(ax, lo, hi, label):
    ax.axhspan(lo, hi, color=BAND, alpha=0.55, linewidth=0)
    ax.text(0.99, np.sqrt(lo * hi) if lo > 0 else hi, label, transform=ax.get_yaxis_transform(),
            ha='right', va='center', fontsize=7, color=INK2)


def strip(ax, groups, values_by_group, ylog=True):
    for gi, (variant, vals) in enumerate(values_by_group):
        vals = np.asarray(vals, float)
        ok = np.isfinite(vals) & (vals > 0) & (vals < 1e6) if ylog else np.isfinite(vals)
        if ylog and (~ok).any():
            ax.annotate(f'{(~ok).sum()}/{len(vals)}\nn/a', (gi, 0.02), xycoords=('data', 'axes fraction'),
                        ha='center', fontsize=6, color=INK2)
        vals = vals[ok]
        jitter = (np.random.default_rng(gi).random(len(vals)) - 0.5) * 0.35
        ax.scatter(gi + jitter, vals, s=9, color=COLORS[variant], alpha=0.55, linewidths=0)
        if len(vals):
            ax.plot([gi - 0.25, gi + 0.25], [np.median(vals)] * 2, color=INK, linewidth=2)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels(groups, fontsize=8, color=INK)
    if ylog:
        ax.set_yscale('log')
    style(ax)


def fig_leaf_uptake(acute):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for ax, crop in zip(axes, ['potato', 'wheat']):
        a = acute[(acute.forcing == 'diurnal') & (acute.crop == crop)]
        groups, vals = [], []
        for variant in VARIANT_ORDER:
            for t in ('day', 'night'):
                groups.append(f'{variant}\n{t}')
                vals.append((variant, a[(a.variant == variant) & (a.time == t)]['R0_leaf']))
        strip(ax, groups, vals)
        band(ax, 10, 50, 'cabbage/radish 10–50 %')
        band(ax, 100 / 1.5, 150, 'rice, day ~100 %')
        band(ax, 30, 40, 'rice, night 30–40 %')
        ax.set_title(f'{crop}: shoot TFWT at end of 1 h exposure', fontsize=9, color=INK, loc='left')
        ax.set_ylabel('% of air-moisture HTO', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig1_leaf_uptake.png'), dpi=150)


def fig_exchange(acute):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for ax, crop in zip(axes, ['potato', 'wheat']):
        a = acute[(acute.forcing == 'diurnal') & (acute.crop == crop)]
        groups, vals = [], []
        for variant in VARIANT_ORDER:
            for t in ('day', 'night'):
                groups.append(f'{variant}\n{t}')
                vals.append((variant, a[(a.variant == variant) & (a.time == t)]['k_exchange']))
        strip(ax, groups, vals)
        band(ax, 0.10, 0.21, 'maize, day 0.10–0.21')
        band(ax, 0.035, 0.13, 'maize, night 0.035–0.13')
        ax.axhline(1.0, color=INK2, linewidth=0.8, linestyle='--')
        ax.set_title(f'{crop}: leaf→air HTO exchange rate kbh_a2', fontsize=9, color=INK, loc='left')
        ax.set_ylabel('h$^{-1}$', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig2_exchange_rate.png'), dpi=150)


def fig_tli(acute):
    fig, ax = plt.subplots(figsize=(6.5, 3.4), constrained_layout=True)
    a = acute[(acute.forcing == 'diurnal') & (acute.crop == 'potato') & (acute.time == 'day')]
    for variant in VARIANT_ORDER:
        g = a[a.variant == variant].groupby('exposure')['TLI_edible']
        med, lo, hi = g.median(), g.quantile(0.1), g.quantile(0.9)
        x = pd.to_datetime('2000-' + med.index)
        ax.fill_between(x, lo.clip(lower=1e-4), hi.clip(lower=1e-4), color=COLORS[variant], alpha=0.15, linewidth=0)
        ax.plot(x, med.clip(lower=1e-4), color=COLORS[variant], linewidth=2, label=variant)
        ax.text(x[-1], med.iloc[-1], f' {variant}', color=INK, fontsize=8, va='center')
    band(ax, 0.2, 0.3, 'potato tuber 0.2–0.3 %')
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b %d'))
    style(ax)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    ax.set_title('Potato: tuber OBT at harvest / shoot TFWT at end of exposure (TLI), daytime',
                 fontsize=9, color=INK, loc='left')
    ax.set_ylabel('TLI, %', fontsize=8, color=INK2)
    ax.set_xlabel('exposure date (median and 10–90 % over 2004–2008)', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig3_TLI_potato.png'), dpi=150)


def fig_chronic(chronic):
    metrics = [('chronic_OBT_over_air', 'edible OBT / air HTO', (0.72, 1.14)),
               ('chronic_TFWT_over_air', 'edible TFWT / air HTO', (0.57, 1.83)),
               ('chronic_leaf_OBT_over_TFWT', 'shoot OBT / shoot TFWT', (0.7 / 1.5, 0.7 * 1.5))]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), constrained_layout=True)
    c = chronic[chronic.forcing == 'diurnal']
    for ax, (m, title, (lo, hi)) in zip(axes, metrics):
        groups, vals = [], []
        for crop in ('potato', 'wheat'):
            for variant in VARIANT_ORDER:
                groups.append(f'{crop[:3]}\n{variant}')
                vals.append((variant, c[(c.crop == crop) & (c.variant == variant)][m]))
        strip(ax, groups, vals)
        band(ax, lo, hi, 'literature')
        ax.set_title(title, fontsize=9, color=INK, loc='left')
        ax.tick_params(axis='x', labelsize=7)
    axes[0].set_ylabel('ratio at harvest', fontsize=8, color=INK2)
    fig.supxlabel('n/a = negative or diverged (>1e6) runs, not plotted; bars = median of plotted runs', fontsize=7, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig4_chronic.png'), dpi=150)


def main():
    os.makedirs(FIG, exist_ok=True)
    acute = pd.read_csv(os.path.join(RES, 'acute_exposures.csv'), dtype={'exposure': str})
    chronic = pd.read_csv(os.path.join(RES, 'chronic_exposure.csv'))
    bench = pd.read_csv(os.path.join(HERE, 'benchmarks', 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index
    comp = comparison(acute, chronic, bench)
    comp.to_csv(os.path.join(RES, 'benchmark_comparison.csv'), index=False)
    show = comp.assign(model=comp.apply(lambda r: f"{fmt(r.model_median)} [{fmt(r.model_p5)}–{fmt(r.model_p95)}]", axis=1),
                       bench=comp.apply(lambda r: f"{fmt(r.bench_low)}–{fmt(r.bench_high)}", axis=1))
    print(show[['benchmark', 'crop', 'variant', 'metric', 'bench', 'model', 'verdict', 'factor_off']].to_string(index=False))
    print(comp.groupby(['crop', 'variant']).verdict.value_counts().unstack(fill_value=0))
    fig_leaf_uptake(acute)
    fig_exchange(acute)
    fig_tli(acute)
    fig_chronic(chronic)


if __name__ == '__main__':
    main()
