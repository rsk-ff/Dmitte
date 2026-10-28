"""Compare the virtual-experiment results with the literature benchmarks.

Reads validation/results/*.csv (from run_validation.py) and
validation/benchmarks/literature_benchmarks.csv, writes
validation/results/benchmark_comparison.csv and the figures in
validation/results/figures/.

The benchmarks are split into a calibration and a frozen test set
(validation/benchmarks/split.csv). By default only calibration benchmarks are
compared and drawn. The test set is evaluated once, after model development, with
    python validation/analyze.py --include-test
which also writes validation/results/benchmark_comparison_test.csv.

Each check compares a benchmark with the matching model metric for the same crop
where crop-specific data exist (match = crop), otherwise with a proxy crop
(match = proxy) or a crop-independent value (match = general), and uses the
exposure duration of the experiment (1 h field box, 2 h chamber).

Verdict on the model median (diurnal forcing, all years and exposure dates):
  within the published range           -> consistent
  within a factor 2 of the range       -> close  (model-data variability is at least a
                                          factor 2, IAEA-TECDOC-1738 Sec. 7.9)
  further away                         -> off (factor given)
  median <= 0, not finite or diverged  -> unphysical (|log10| > 30)
Benchmarks that only give a central value get a factor-1.5 tolerance.
"""
import argparse
import hashlib
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')
FIG = os.path.join(RES, 'figures')
BENCH_DIR = os.path.join(HERE, 'benchmarks')

VARIANT_ORDER = ['model']                 # variant names written by run_validation.py
COLORS = {'model': '#2a78d6'}
BAND = '#c9c8c2'
INK, INK2 = '#0b0b0b', '#52514e'

# sha256 of the test benchmarks and test datasets when the split was frozen (2026-10-28)
FROZEN_TEST_FINGERPRINT = '067765dc47f5b13370605917708c3e381386f9b6d1367f5c15839e20926b1dd8'


def check(bid, crop, metric, match, time='all', hours=1, phase=None, daa=None, subset=None):
    """One benchmark-model comparison.

    time: day / night / all (acute rows); hours: exposure duration; phase: grain-filling phase
    (before / during / after); daa: (min, max) days after anthesis; subset: 'pair' or 'chronic'.
    """
    return dict(benchmark=bid, crop=crop, metric=metric, match=match, time=time, hours=hours, phase=phase,
                daa=daa, subset=subset)


CHECKS = [
    # ---- potato
    check('B1', 'potato', 'R0_leaf', 'proxy'),                       # radish, 1 h field box
    check('B10', 'potato', 'tfwt_drop_edible', 'proxy'),             # radish root
    check('B13', 'potato', 'TLI_edible', 'proxy'),                   # radish root
    check('B56', 'potato', 'k_exchange', 'crop', time='day'),
    check('B57', 'potato', 'k_exchange', 'crop', time='night'),
    check('B58', 'potato', 'k_exchange', 'crop', time='day'),
    check('B6', 'potato', 'k_day_over_night', 'general', subset='pair'),
    check('B15', 'potato', 'obt_night_over_day', 'general', subset='pair'),
    check('B18', 'potato', 'chronic_OBT_over_air', 'crop', subset='chronic'),
    check('B19', 'potato', 'chronic_TFWT_over_air', 'crop', subset='chronic'),
    check('B20', 'potato', 'chronic_leaf_OBT_over_TFWT', 'general', subset='chronic'),
    check('B12', 'potato', 'TLI_edible', 'crop', hours=2, daa=(13, 45)),            # test
    check('B24', 'potato', 'R0_leaf', 'crop', time='day', hours=2, daa=(13, 45)),   # test
    check('B25', 'potato', 'R0_leaf', 'crop', time='night', hours=2, daa=(13, 45)),  # test
    check('B34', 'potato', 'chronic_OBT_over_air', 'crop', subset='chronic'),       # test (root crops)
    # ---- wheat
    check('B4', 'wheat', 'R0_night_over_day', 'crop', subset='pair'),
    check('B49', 'wheat', 'R0_leaf', 'crop', time='day'),
    check('B50', 'wheat', 'R0_leaf', 'crop', time='night'),
    check('B5', 'wheat', 'k_exchange', 'crop', time='day'),
    check('B29', 'wheat', 'k_exchange', 'crop', time='night'),
    check('B59', 'wheat', 'k_exchange', 'crop', time='day'),
    check('B60', 'wheat', 'k_exchange', 'crop', time='night'),
    check('B6', 'wheat', 'k_day_over_night', 'general', subset='pair'),
    check('B22', 'wheat', 'obt_end_over_leaf', 'crop', time='day'),
    check('B23', 'wheat', 'obt_end_over_leaf', 'crop', time='night'),
    check('B15', 'wheat', 'obt_night_over_day', 'general', subset='pair'),
    check('B26', 'wheat', 'TLI_edible', 'crop', daa=(12, 28)),
    check('B9', 'wheat', 'tfwt_drop_leaf', 'proxy'),                 # rice
    check('B20', 'wheat', 'chronic_leaf_OBT_over_TFWT', 'general', subset='chronic'),
    check('B46', 'wheat', 'R0_leaf', 'crop', time='day', hours=2),            # test
    check('B47', 'wheat', 'R0_leaf', 'crop', time='day', hours=2),            # test
    check('B48', 'wheat', 'R0_leaf', 'crop', time='night', hours=2),          # test
    check('B51', 'wheat', 'obt_end_over_leaf', 'crop', time='day', hours=2),  # test
    check('B52', 'wheat', 'obt_end_over_leaf', 'crop', time='night', hours=2),  # test
    check('B53', 'wheat', 'TLI_edible', 'crop', hours=2, phase='before'),     # test
    check('B54', 'wheat', 'TLI_edible', 'crop', hours=2, phase='during'),     # test
    check('B55', 'wheat', 'TLI_edible', 'crop', hours=2, phase='after'),      # test
    check('B35', 'wheat', 'chronic_OBT_over_air', 'crop', subset='chronic'),  # test (grain)
]
ONE_SIDED = {'B10': (0, None)}   # B10: "up to 1.3e4" (upper bound)


# --------------------------------------------------------------------------- split
def load_split():
    split = pd.read_csv(os.path.join(BENCH_DIR, 'split.csv'))
    test = split[split.split == 'test']
    ids = set(test.id[~test.id.str.startswith('datasets/')])
    files = sorted(test.id[test.id.str.startswith('datasets/')])
    return ids, files


def test_fingerprint(bench, test_ids, test_files):
    """sha256 over the test benchmark rows and the test dataset files."""
    h = hashlib.sha256()
    rows = bench.loc[sorted(test_ids)].drop(columns='id', errors='ignore').reset_index().astype(str).to_dict('records')
    h.update(json.dumps(rows, sort_keys=True).encode())
    for f in test_files:
        with open(os.path.join(BENCH_DIR, f), 'rb') as fh:
            h.update(fh.read())
    return h.hexdigest()


# --------------------------------------------------------------------------- comparison
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
    return ('close' if factor <= 2 else 'off'), factor


def paired(acute, metric):
    """Night/day ratios for the same crop, year, variant, exposure date and duration."""
    d = acute[acute.forcing == 'diurnal'].pivot_table(index=['crop', 'year', 'variant', 'exposure', 'duration_h'],
                                                       columns='time', values=metric)
    return d['night'] / d['day']


def select(acute, c, variant):
    a = acute[(acute.forcing == 'diurnal') & (acute.crop == c['crop']) & (acute.variant == variant)
              & (acute.duration_h == c['hours'])]
    if c['time'] in ('day', 'night'):
        a = a[a.time == c['time']]
    if c['phase']:
        a = a[a.grain_phase == c['phase']]
    if c['daa']:
        a = a[a.days_after_anthesis.between(*c['daa'])]
    return a


def model_values(acute, chronic, c, variant):
    if c['subset'] == 'chronic':
        ch = chronic[(chronic.forcing == 'diurnal') & (chronic.crop == c['crop']) & (chronic.variant == variant)]
        return ch[c['metric']].values
    if c['subset'] == 'pair':
        sub = select(acute, dict(c, time='all'), variant)
        if c['metric'] == 'R0_night_over_day':
            return paired(sub, 'R0_leaf').values
        if c['metric'] == 'k_day_over_night':
            return 1 / paired(sub, 'k_exchange').values
        if c['metric'] == 'obt_night_over_day':
            return paired(sub, 'obt_formed_24h_per_leaf_h0').values
    return select(acute, c, variant)[c['metric']].values


def fmt(x):
    return f'{x:.3g}' if np.isfinite(x) else 'nan'


def comparison(acute, chronic, bench, checks):
    rows = []
    for c in checks:
        b = bench.loc[c['benchmark']]
        lo, hi = bench_range(b)
        for variant in VARIANT_ORDER:
            v = np.asarray(model_values(acute, chronic, c, variant), dtype=float)
            v = v[np.isfinite(v)]
            med = float(np.median(v)) if len(v) else np.nan
            p5, p95 = (np.percentile(v, [5, 95]) if len(v) else (np.nan, np.nan))
            verd, factor = verdict(med, lo, hi)
            rows.append(dict(benchmark=c['benchmark'], crop=c['crop'], variant=variant, metric=c['metric'],
                             match=c['match'], exposure_h=c['hours'] if c['subset'] != 'chronic' else np.nan,
                             benchmark_quantity=b['quantity'], benchmark_crop=b['crop'], setting=b['setting'],
                             bench_low=lo, bench_high=hi, unit=b['unit'],
                             model_median=med, model_p5=p5, model_p95=p95,
                             n=len(v), verdict=verd, factor_off=factor))
    return pd.DataFrame(rows)


def report(comp, title):
    show = comp.assign(model=comp.apply(lambda r: f"{fmt(r.model_median)} [{fmt(r.model_p5)}–{fmt(r.model_p95)}]", axis=1),
                       bench=comp.apply(lambda r: f"{fmt(r.bench_low)}–{fmt(r.bench_high)}", axis=1))
    print(f'==== {title}')
    print(show[['benchmark', 'crop', 'match', 'variant', 'metric', 'bench', 'model', 'verdict', 'factor_off']]
          .to_string(index=False))
    print(comp.groupby(['crop', 'variant']).verdict.value_counts().unstack(fill_value=0))


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


def bands(ax, bench, allowed, items):
    """Draw literature bands for the listed benchmarks that are allowed (calibration unless --include-test)."""
    for bid, label in items:
        if bid in allowed:
            lo, hi = bench_range(bench.loc[bid])
            band(ax, lo, hi, f'{bid} {label}')


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


def fig_leaf_uptake(acute, bench, allowed):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    crop_bands = {'potato': [('B1', 'radish 1 h (proxy)'), ('B24', 'potato 2 h day'), ('B25', 'potato 2 h night')],
                  'wheat': [('B49', 'wheat 1 h day'), ('B50', 'wheat 1 h night')]}
    for ax, crop in zip(axes, ['potato', 'wheat']):
        a = acute[(acute.forcing == 'diurnal') & (acute.crop == crop) & (acute.duration_h == 1)]
        groups, vals = [], []
        for variant in VARIANT_ORDER:
            for t in ('day', 'night'):
                groups.append(f'{variant}\n{t}')
                vals.append((variant, a[(a.variant == variant) & (a.time == t)]['R0_leaf']))
        strip(ax, groups, vals)
        bands(ax, bench, allowed, crop_bands[crop])
        ax.set_title(f'{crop}: shoot TFWT at end of 1 h exposure', fontsize=9, color=INK, loc='left')
        ax.set_ylabel('% of air-moisture HTO', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig1_leaf_uptake.png'), dpi=150)


def fig_exchange(acute, bench, allowed):
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    crop_bands = {'potato': [('B56', 'potato day'), ('B57', 'potato night')],
                  'wheat': [('B5', 'wheat day'), ('B29', 'wheat night')]}
    for ax, crop in zip(axes, ['potato', 'wheat']):
        a = acute[(acute.forcing == 'diurnal') & (acute.crop == crop) & (acute.duration_h == 1)]
        groups, vals = [], []
        for variant in VARIANT_ORDER:
            for t in ('day', 'night'):
                groups.append(f'{variant}\n{t}')
                vals.append((variant, a[(a.variant == variant) & (a.time == t)]['k_exchange']))
        strip(ax, groups, vals)
        bands(ax, bench, allowed, crop_bands[crop])
        ax.axhline(1.0, color=INK2, linewidth=0.8, linestyle='--')
        ax.set_title(f'{crop}: leaf→air HTO exchange rate kbh_a2', fontsize=9, color=INK, loc='left')
        ax.set_ylabel('h$^{-1}$', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig2_exchange_rate.png'), dpi=150)


def fig_tli(acute, bench, allowed):
    fig, ax = plt.subplots(figsize=(6.5, 3.4), constrained_layout=True)
    a = acute[(acute.forcing == 'diurnal') & (acute.crop == 'potato') & (acute.time == 'day')
              & (acute.duration_h == 1)]
    for variant in VARIANT_ORDER:
        g = a[a.variant == variant].groupby('exposure')['TLI_edible']
        med, lo, hi = g.median(), g.quantile(0.1), g.quantile(0.9)
        x = pd.to_datetime('2000-' + med.index)
        ax.fill_between(x, lo.clip(lower=1e-4), hi.clip(lower=1e-4), color=COLORS[variant], alpha=0.15, linewidth=0)
        ax.plot(x, med.clip(lower=1e-4), color=COLORS[variant], linewidth=2, label=variant)
        ax.text(x[-1], med.iloc[-1], f' {variant}', color=INK, fontsize=8, va='center')
    bands(ax, bench, allowed, [('B13', 'radish root (proxy)'), ('B12', 'potato tuber, 2 h')])
    ax.set_yscale('log')
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%b %d'))
    style(ax)
    ax.legend(fontsize=7, frameon=False, loc='upper left')
    ax.set_title('Potato: tuber OBT at harvest / shoot TFWT at end of 1 h exposure (TLI), daytime',
                 fontsize=9, color=INK, loc='left')
    ax.set_ylabel('TLI, %', fontsize=8, color=INK2)
    ax.set_xlabel('exposure date (median and 10–90 % over 2004–2008)', fontsize=8, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig3_TLI_potato.png'), dpi=150)


def fig_chronic(chronic, bench, allowed):
    metrics = [('chronic_OBT_over_air', 'edible OBT / air HTO', [('B18', 'garden'), ('B34', 'root crops'),
                                                                 ('B35', 'grain')]),
               ('chronic_TFWT_over_air', 'edible TFWT / air HTO', [('B19', 'garden')]),
               ('chronic_leaf_OBT_over_TFWT', 'shoot OBT / shoot TFWT', [('B20', 'equilibrium')])]
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), constrained_layout=True)
    c = chronic[chronic.forcing == 'diurnal']
    for ax, (m, title, items) in zip(axes, metrics):
        groups, vals = [], []
        for crop in ('potato', 'wheat'):
            for variant in VARIANT_ORDER:
                groups.append(f'{crop[:3]}\n{variant}')
                vals.append((variant, c[(c.crop == crop) & (c.variant == variant)][m]))
        strip(ax, groups, vals)
        bands(ax, bench, allowed, items)
        ax.set_title(title, fontsize=9, color=INK, loc='left')
        ax.tick_params(axis='x', labelsize=7)
    axes[0].set_ylabel('ratio at harvest', fontsize=8, color=INK2)
    fig.supxlabel('n/a = negative or diverged (>1e6) runs, not plotted; bars = median of plotted runs', fontsize=7, color=INK2)
    fig.savefig(os.path.join(FIG, 'fig4_chronic.png'), dpi=150)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--include-test', action='store_true',
                        help='also evaluate the frozen test benchmarks (do this once, after model development)')
    args = parser.parse_args()

    os.makedirs(FIG, exist_ok=True)
    acute = pd.read_csv(os.path.join(RES, 'acute_exposures.csv'), dtype={'exposure': str})
    chronic = pd.read_csv(os.path.join(RES, 'chronic_exposure.csv'))
    bench = pd.read_csv(os.path.join(BENCH_DIR, 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index

    test_ids, test_files = load_split()
    fingerprint = test_fingerprint(bench, test_ids, test_files)
    status = 'unchanged' if fingerprint == FROZEN_TEST_FINGERPRINT else 'CHANGED since the split was frozen'
    print(f'test split: {len(test_ids)} benchmarks, {len(test_files)} datasets, fingerprint {fingerprint[:12]} ({status})')

    calib = [c for c in CHECKS if c['benchmark'] not in test_ids]
    comp = comparison(acute, chronic, bench, calib)
    comp.to_csv(os.path.join(RES, 'benchmark_comparison.csv'), index=False)
    report(comp, 'calibration benchmarks')

    allowed = set(bench.index) - test_ids
    if args.include_test:
        test = comparison(acute, chronic, bench, [c for c in CHECKS if c['benchmark'] in test_ids])
        test.to_csv(os.path.join(RES, 'benchmark_comparison_test.csv'), index=False)
        report(test, 'TEST benchmarks')
        allowed = set(bench.index)
    fig_leaf_uptake(acute, bench, allowed)
    fig_exchange(acute, bench, allowed)
    fig_tli(acute, bench, allowed)
    fig_chronic(chronic, bench, allowed)


if __name__ == '__main__':
    main()
