"""Figures 1-4 of the manuscript (Figures 5-7 are analysis/results/figures/fig_b1..b3).

Run from the repository root after validation/build_inputs.py, validation/run_validation.py and
validation/analyze.py (--include-test was run once on 2026-10-29; its CSV is committed):
    MPLBACKEND=Agg python docs/manuscript/make_figures.py
writes docs/manuscript/figures/fig1_structure.png ... fig4_tli_phase.png.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
FIG = os.path.join(HERE, 'figures')
RES = os.path.join(ROOT, 'validation', 'results')
os.makedirs(FIG, exist_ok=True)

INK, MUTED, GRID = '#222222', '#6b6b6b', '#d9d9d9'
BLUE, ORANGE, GREEN, PURPLE = '#2f6db5', '#d9822b', '#3a923a', '#8a5bb5'
plt.rcParams.update({'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.spines.top': False,
                     'axes.spines.right': False, 'savefig.dpi': 200, 'savefig.bbox': 'tight'})


def fig1_structure():
    """Compartments and the main HTO/OBT pathways; shaded boxes have rates driven by WOFOST."""
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0.4, 6.9)
    ax.axis('off')
    boxes = {
        'a2': (1.6, 6.0, 'Air moisture HTO\n(boundary condition)'),
        'bh': (1.6, 4.0, 'Leaf water\nHTO (bh)'),
        's1': (1.6, 1.8, 'Soil water HTO\nlayers 1-3 (s1-s3)'),
        'ba': (5.0, 4.0, 'Leaf assimilate\nOBT (ba)'),
        'bo': (8.4, 4.0, 'Stem/leaf\nstructural OBT (bo)'),
        'fh': (5.0, 1.8, 'Grain/tuber\nwater HTO (fh)'),
        'fo': (8.4, 1.8, 'Grain/tuber\nOBT (fo)'),
    }
    w, h = 2.2, 0.9
    for key, (x, y, text) in boxes.items():
        color = '#f4f4f4' if key in ('a2', 's1') else '#eef3fa'
        ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle='round,pad=0.04',
                                    fc=color, ec=MUTED, lw=0.8))
        ax.text(x, y, text, ha='center', va='center', color=INK, fontsize=8)

    def arrow(a, b, label, color, rad=0.0, dx=0.0, lx=0.0, ly=0.0, ha='center'):
        (x0, y0, _), (x1, y1, _) = boxes[a], boxes[b]
        ax.add_patch(FancyArrowPatch((x0 + dx, y0), (x1 + dx, y1), arrowstyle='-|>', mutation_scale=9,
                                     shrinkA=24, shrinkB=24, connectionstyle=f'arc3,rad={rad}', color=color, lw=0.9))
        if label:
            ax.text((x0 + x1) / 2 + lx, (y0 + y1) / 2 + ly, label, fontsize=7, color=color, ha=ha, va='center')

    arrow('a2', 'bh', '', BLUE, dx=-0.25)
    arrow('bh', 'a2', 'leaf-air exchange,\nsymmetric (total resistance)', BLUE, dx=0.25, lx=0.45, ly=0.3, ha='left')
    arrow('a2', 's1', 'deposition', MUTED, rad=0.55, lx=-1.25, ly=-0.3)
    arrow('s1', 'bh', 'root uptake =\ntranspiration\n+ growth', ORANGE, lx=-0.15, ha='right')
    arrow('bh', 'ba', 'gross photosynthesis', GREEN, ly=0.28)
    arrow('ba', 'bh', 'respiration', GREEN, rad=-0.45, ly=-0.6)
    arrow('ba', 'bo', 'shoot growth', GREEN, ly=0.28)
    arrow('ba', 'fo', 'grain/tuber growth', GREEN, lx=0.3, ly=0.2, ha='left')
    arrow('bh', 'fh', 'water to fruit', PURPLE, lx=0.15, ly=-0.35, ha='left')
    arrow('s1', 'fh', 'root uptake', ORANGE, ly=-0.25)
    ax.text(5.0, 0.6, 'Blue boxes: rates updated hourly from WOFOST (green-leaf water, transpiration, gross assimilation, '
                      'dry-matter partitioning).\nHourly propagator = matrix exponential of the rate matrix.',
            ha='center', fontsize=7.2, color=MUTED)
    fig.savefig(os.path.join(FIG, 'fig1_structure.png'))
    plt.close(fig)


def fig2_rates():
    """Daily mean of the growth-driven rates over one season, with the season-mean (constant) value."""
    import common as c
    from dmitte.transfer_equotions import K_INDEX
    series = {'leaf→air exchange $k_{bh,a2}$': ('kbh_a2', BLUE),
              'OBT formation $k_{bh,ba}$': ('kbh_ba', GREEN),
              'assimilate→grain/tuber $k_{ba,fo}$': ('kba_fo', ORANGE),
              'leaf→fruit water $k_{bh,fh}$': ('kbh_fh', PURPLE)}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, (crop, year) in zip(axes, (('wheat', 2006), ('potato', 2006))):
        s = c.Season(crop, year)
        idx = s.index[:s.n]
        daa = (idx - s.anthesis) / pd.Timedelta(days=1)
        for label, (name, color) in series.items():
            k = pd.Series(s.k[:, K_INDEX[name]], index=daa)
            daily = k.groupby(np.floor(k.index)).mean()
            daily = daily[daily > 0]
            ax.plot(daily.index, daily.values, color=color, lw=1.2, label=label)
            mean = s.k[:, K_INDEX[name]].mean()
            if mean > 0:
                ax.axhline(mean, color=color, lw=0.8, ls='--')
        if s.lin_start is not None:
            ax.axvspan((s.lin_start - s.anthesis) / pd.Timedelta(days=1), (s.lin_end - s.anthesis) / pd.Timedelta(days=1),
                       color=GRID, alpha=0.5, lw=0)
        ax.set_yscale('log')
        ax.set_ylim(1e-6, 3)
        ax.set_title(f'{crop} {year}', fontsize=9, color=INK)
        ax.set_xlabel('days after anthesis' if crop == 'wheat' else 'days after DVS = 1 (tuber initiation)')
        ax.grid(axis='y', color=GRID, lw=0.5)
    axes[0].set_ylabel('daily mean rate, h$^{-1}$')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, fontsize=7, frameon=False, loc='lower center', ncol=4, bbox_to_anchor=(0.5, -0.14))
    fig.text(0.5, -0.19, 'Dashed: season mean (constant-rate baseline). Grey band: linear phase of storage-organ growth.',
             ha='center', fontsize=7.5, color=MUTED)
    fig.savefig(os.path.join(FIG, 'fig2_rates.png'))
    plt.close(fig)


def signed_miss(r):
    if r.bench_low <= r.model_median <= r.bench_high:
        return 0.0
    return np.log10(r.model_median / r.bench_high) if r.model_median > r.bench_high else np.log10(r.model_median / r.bench_low)


def fig3_benchmarks():
    """Signed log10 factor by which each model median misses its published range (0 inside)."""
    cal = pd.read_csv(os.path.join(RES, 'benchmark_comparison.csv')).assign(set='calibration')
    test = pd.read_csv(os.path.join(RES, 'benchmark_comparison_test.csv')).assign(set='test')
    d = pd.concat([cal, test], ignore_index=True)
    d['miss'] = d.apply(signed_miss, axis=1)
    d['label'] = d.benchmark + ' ' + d.crop
    d = d.sort_values(['set', 'crop', 'miss']).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    x = np.arange(len(d))
    for s, marker, color in (('calibration', 'o', BLUE), ('test', '^', ORANGE)):
        m = d.set == s
        ax.scatter(x[m], d.miss[m], marker=marker, color=color, s=22, zorder=3, label=f'{s} ({m.sum()})')
    lim = np.log10(2)
    ax.axhspan(-lim, lim, color=GRID, alpha=0.6, lw=0)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(d.label, rotation=90, fontsize=6.5)
    ax.set_ylabel('log$_{10}$(model / nearest bound)')
    ax.set_ylim(-0.6, 0.6)
    ax.legend(frameon=False, fontsize=7.5, loc='upper left')
    ax.text(len(d) - 0.5, lim + 0.02, 'factor 2', ha='right', va='bottom', fontsize=7, color=MUTED)
    fig.savefig(os.path.join(FIG, 'fig3_benchmarks.png'))
    plt.close(fig)


def fig4_tli_phase():
    """Grain TLI after a 2 h exposure by grain-filling phase (test benchmarks B53-B55)."""
    a = pd.read_csv(os.path.join(RES, 'acute_exposures.csv'))
    w = a[(a.crop == 'wheat') & (a.forcing == 'diurnal') & (a.duration_h == 2)]
    ranges = {'before': (0.06, 0.2, 'B53'), 'during': (0.4, 0.9, 'B54'), 'after': (0.04, 0.2, 'B55')}
    colors = {'before': MUTED, 'during': GREEN, 'after': ORANGE}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={'width_ratios': [1.6, 1]})
    for ph, g in w.groupby('grain_phase'):
        for tod, marker in (('day', 'o'), ('night', 's')):
            gg = g[g.time == tod]
            ax1.scatter(gg.days_after_anthesis, gg.TLI_edible, s=12, marker=marker, color=colors[ph], alpha=0.8,
                        label=f'{ph}, {tod}')
    ax1.set_xlabel('exposure date, days after anthesis')
    ax1.set_ylabel('grain TLI, %')
    ax1.legend(fontsize=6.5, frameon=False, ncol=2)
    for i, ph in enumerate(('before', 'during', 'after')):
        v = w[w.grain_phase == ph].TLI_edible
        lo, hi, bid = ranges[ph]
        ax2.add_patch(plt.Rectangle((i - 0.35, lo), 0.7, hi - lo, color=GRID, alpha=0.8, lw=0))
        ax2.boxplot(v, positions=[i], widths=0.3, showfliers=False,
                    medianprops=dict(color=colors[ph], lw=1.5), boxprops=dict(color=MUTED),
                    whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
        ax2.text(i, 0.97, bid, ha='center', fontsize=7, color=MUTED)
    ax2.set_xticks(range(3))
    ax2.set_xticklabels(['before', 'during', 'after'])
    ax2.set_xlabel('linear grain-filling phase')
    ax2.set_ylim(0, 1.0)
    fig.savefig(os.path.join(FIG, 'fig4_tli_phase.png'))
    plt.close(fig)


if __name__ == '__main__':
    fig1_structure()
    fig3_benchmarks()
    fig4_tli_phase()
    fig2_rates()
