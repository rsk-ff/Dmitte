"""Supplementary analysis after the test run: the wheat chamber experiments driven by their measured hourly conditions.

The chamber benchmarks (B46-B48, B51-B55) are test data. The frozen checks took them from the generic virtual
experiment (site weather, 2 h square pulse), and the test run (validation/README.md section 4.6) found B46 low
by a factor 2.1. This script re-runs the two experiments with a complete hourly record (scenario cards
WHEAT_CH_HOURLY_DAY / _NIGHT; Diabate & Strack 1997 Tables 2 and 3) to see how much of that deviation comes from
the experimental conditions. It changes no model parameter and no benchmark verdict: it is a labelled
supplementary analysis, not an independent test.

For each experiment and each WOFOST wheat season of run_validation.py (2005-2008):

* the experiment starts at 09:00 on the day with the same number of days after anthesis (12 d day, 19 d night);
* for the 26 h in the chamber each hour uses the measured temperature, RH and PPFD (global radiation =
  PPFD / 2.2, as in scenario_field_box.py), no rain, and air-moisture HTO = hourly integrated air HTO / 3600 s /
  absolute humidity of that hour; air-to-soil deposition is switched off (soil covered with polyethylene);
* afterwards the plants are outdoors (site weather, no HTO in air) until the WOFOST harvest;
* the chamber fan speed is not reported, so the aerodynamic resistance RAM + RB in the chamber is scanned
  as in the field box: site wind, 30, 15 and 8 s/m;
* 'day_high_light' is the high-light case of B46, for which no hourly record is published: the day experiment
  with PPFD 900 instead of 120 umol m-2 s-1 whenever the chamber light is on. It is an inference from the
  conditions, not a reproduction of a measured run.

Because the air is driven in absolute units, model and measured leaf and ear TFWT (Bq/mL) and OBT (Bq/g dry
matter) are compared directly, without the normalisation to air HTO that the Table 4 ratios use.

Run from the repository root after validation/build_inputs.py:
    python validation/supplementary_chamber.py
writes validation/results/supplementary_chamber.csv, supplementary_chamber_stdout.txt and
validation/results/figures/fig6_chamber_supplementary.png.
"""
import os
import sys
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.linalg import expm  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')

import run_validation as rv  # noqa: E402
import scenario_field_box as fb  # noqa: E402
from dmitte.para_constant import COMB_WATER  # noqa: E402
from dmitte.transfer_equotions import K_INDEX, transfer_rates  # noqa: E402

DATA = fb.DATA
OUT = fb.OUT
CHAMBER_HOURS = 26
EXPOSURE_H = 2
SAMPLES = (2, 3, 4, 26)                          # sampling times in the chamber, h after the start of supply
EXPERIMENTS = {'day_low_light': ('day_low_light', 12, None),
               'day_high_light': ('day_low_light', 12, 900.0),
               'night': ('night', 19, None)}
TABLE4 = {'day_low_light': ('B47', 75, 70, 80), 'day_high_light': ('B46', 86, 84, 88), 'night': ('B48', 20, 13, 27)}


def chamber_series(name, ppfd_on):
    """Hourly chamber conditions; air_moist is the mean air-moisture HTO (Bq/mL) of the hour starting at that row."""
    d = pd.read_csv(os.path.join(DATA, 'wheat_DE_chamber_timeseries_1993.csv'))
    d = d[(d.experiment == name) & d.hours_after_start.notna()].sort_values('hours_after_start').reset_index(drop=True)
    es = 0.6108 * np.exp(17.27 * d.T_C / (d.T_C + 237.3))           # kPa
    ah = d.RH_pct / 100 * es * 1000 / (461.5 * (d.T_C + 273.15))     # kg m-3 = g/L
    conc = d.air_HTO_integrated_Bq_s_per_m3.shift(-1) / 3600        # Bq m-3 in hour [h, h+1)
    d['air_moist'] = conc / (ah * 1000)                              # Bq per g (= mL) of air moisture
    if ppfd_on is not None:
        d.loc[d.PPFD_umol_m2_s > 0, 'PPFD_umol_m2_s'] = ppfd_on
    return d.iloc[:CHAMBER_HOURS]


def observed(name):
    d = pd.read_csv(os.path.join(DATA, 'wheat_DE_chamber_timeseries_1993.csv'))
    return d[d.experiment == name]


def set_chamber(forcing, i0, ch, ra_times_u, ra):
    f = forcing.copy()
    col = {c: f.columns.get_loc(c) for c in ('Ta_Avg', 'RH_Avg', 'Pvapor_Avg', 'DR_Avg', 'Rain_Tot', 'WS10m_avg')}
    for j, r in ch.iterrows():
        i = i0 + j
        es = 0.6108 * np.exp(17.27 * r.T_C / (r.T_C + 237.3))
        f.iloc[i, col['Ta_Avg']] = r.T_C
        f.iloc[i, col['RH_Avg']] = r.RH_pct
        f.iloc[i, col['Pvapor_Avg']] = r.RH_pct / 100 * es
        f.iloc[i, col['DR_Avg']] = r.PPFD_umol_m2_s * fb.GLOBAL_PER_PPFD
        f.iloc[i, col['Rain_Tot']] = 0.0
        if ra is not None:
            f.iloc[i, col['WS10m_avg']] = ra_times_u[i] / ra
    return f


def run_chamber(k, pools, i0, air_moist):
    """Inventories from the start of supply (row 0) to harvest; air clamped to the measured series in the chamber."""
    k = k[i0:].copy()
    k[:CHAMBER_HOURS, K_INDEX['ka2_s1']] = 0.0                 # soil covered
    M = rv.generators(k)
    M[:, rv.A_A2, :] = 0                                       # air is a boundary condition
    y = np.zeros(rv.N_COMP)
    ys = [y.copy()]
    for j, Mj in enumerate(M):
        y[rv.A_A2] = pools['w_air'][i0 + j] * air_moist[j] if j < CHAMBER_HOURS else 0.0
        y = expm(Mj) @ y
        y[rv.A_A2] = 0.0
        ys.append(y.copy())
    return np.array(ys)


def metrics(ys, pools, i0, air_moist):
    n = len(pools['w_plant'])
    at = lambda j: min(i0 + j, n - 1)
    out = {}
    for j in SAMPLES:
        out[f'leaf_TWT_{j}h'] = ys[j, rv.A_BH] / pools['w_plant'][at(j)]
        out[f'ear_TWT_{j}h'] = rv.safe_div(ys[j, rv.A_FH], pools['w_fruit'][at(j)])
        out[f'leaf_OBT_{j}h'] = COMB_WATER * rv.leaf_obt_conc(ys[j], pools, at(j))
    hv = len(ys) - 1
    out['grain_OBT_harvest'] = COMB_WATER * ys[hv, rv.A_FO] / pools['comb_edible'][n - 1]
    out['leaf_TWT_harvest'] = ys[hv, rv.A_BH] / pools['w_plant'][n - 1]
    c_mean = np.mean(air_moist[:EXPOSURE_H])
    out['R0_leaf'] = 100 * out['leaf_TWT_2h'] / c_mean
    out['obt_end_over_leaf'] = 100 * rv.leaf_obt_conc(ys[EXPOSURE_H], pools, at(EXPOSURE_H)) / (ys[EXPOSURE_H, rv.A_BH] / pools['w_plant'][at(EXPOSURE_H)])
    out['TLI'] = 100 * ys[hv, rv.A_FO] / pools['comb_edible'][n - 1] / out['leaf_TWT_2h']
    return out


def run(cases=fb.RA_CASES):
    rows = []
    for year in rv.CROPS['wheat']['years']:
        pr, forcing, start, end, anthesis = fb.season(year)
        idx = forcing.index
        air0 = fb.model_state(pr, forcing, start, end)[0]
        ra_times_u = (air0.RAM + air0.RB) * air0.U10
        for exp, (name, daa, ppfd_on) in EXPERIMENTS.items():
            ch = chamber_series(name, ppfd_on)
            i0 = idx.get_loc(anthesis.normalize() + pd.Timedelta(days=daa, hours=9))
            for case, ra in cases.items():
                air, soil, plant = fb.model_state(pr, set_chamber(forcing, i0, ch, ra_times_u, ra), start, end)
                k = transfer_rates('HTO', 'CEREAL', air, soil, plant)
                pools = rv.season_arrays(air, soil, plant)
                ys = run_chamber(k, pools, i0, ch.air_moist.values)
                rows.append(dict(experiment=exp, year=year, case=case,
                                 ra_exposure=float(np.mean(air.RAM[i0:i0 + 2] + air.RB[i0:i0 + 2])),
                                 rc_exposure=float(np.mean(plant.RC1[i0:i0 + 2]) * 100),
                                 **metrics(ys, pools, i0, ch.air_moist.values)))
        print('wheat', year, 'done', flush=True)
    return pd.DataFrame(rows)


def report(res):
    med = res.groupby(['experiment', 'case'], sort=False).median(numeric_only=True)
    print('\nmodel (median over 4 seasons) / measured, absolute concentrations (Bq/mL water, Bq/g dm):')
    for exp, (name, _, _) in EXPERIMENTS.items():
        o = observed(name).set_index('hours_after_start')
        harvest = observed(name)[observed(name).hours_after_start.isna()].iloc[0]
        ch = chamber_series(name, None)
        c_mean = ch.air_moist.iloc[:EXPOSURE_H].mean()
        rows = []
        for case in fb.RA_CASES:
            m = med.loc[(exp, case)]
            r = dict(case=case, RAMRB=round(m.ra_exposure), RC=round(m.rc_exposure))
            for j in SAMPLES:
                r[f'leaf {j}h'] = m[f'leaf_TWT_{j}h'] / o.loc[j, 'TWT_leaves_Bq_per_mL']
            for j in SAMPLES:
                r[f'ear {j}h'] = m[f'ear_TWT_{j}h'] / o.loc[j, 'TWT_ears_Bq_per_mL']
            for j in SAMPLES:
                r[f'leafOBT {j}h'] = m[f'leaf_OBT_{j}h'] / o.loc[j, 'OBT_leaves_Bq_per_g_dm']
            r['grain OBT hv'] = m.grain_OBT_harvest / harvest.OBT_grains_Bq_per_g_dm
            rows.append(r)
        t = pd.DataFrame(rows)
        note = ' (measured: day_low_light record; only the model light changes)' if exp == 'day_high_light' else ''
        print(f'--- {exp}{note}; measured mean air moisture over the 2 h supply {c_mean:.3g} Bq/mL')
        print(t.round(2).to_string(index=False))
    print('\nbenchmark-style ratios (model median over seasons; 5-95 % over seasons and cases in brackets):')
    print('R0_leaf = leaf TFWT at 2 h / mean air-moisture HTO over the 2 h supply, both computed from the hourly record')
    out = []
    for exp, (name, _, _) in EXPERIMENTS.items():
        o = observed(name).set_index('hours_after_start')
        harvest = observed(name)[observed(name).hours_after_start.isna()].iloc[0]
        c_mean = chamber_series(name, None).air_moist.iloc[:EXPOSURE_H].mean()
        bid, t4, lo, hi = TABLE4[exp]
        g = res[res.experiment == exp]
        q = lambda c: f'{g[c].median():.3g} ({g[c].quantile(.05):.3g}-{g[c].quantile(.95):.3g})'
        out.append(dict(experiment=exp, benchmark=bid, Table4_R0=f'{t4} ({lo}-{hi})',
                        record_R0=('-' if exp == 'day_high_light' else round(100 * o.loc[2, 'TWT_leaves_Bq_per_mL'] / c_mean)),
                        model_R0=q('R0_leaf'),
                        record_OBT_over_TWT_2h=('-' if exp == 'day_high_light' else
                                                round(100 * o.loc[2, 'OBT_leaves_Bq_per_g_dm'] / COMB_WATER / o.loc[2, 'TWT_leaves_Bq_per_mL'], 2)),
                        model_OBT_over_TWT_2h=q('obt_end_over_leaf'),
                        record_TLI=('-' if exp == 'day_high_light' else
                                    round(100 * harvest.OBT_grains_Bq_per_g_dm / COMB_WATER / o.loc[2, 'TWT_leaves_Bq_per_mL'], 2)),
                        model_TLI=q('TLI')))
    print(pd.DataFrame(out).to_string(index=False))
    print('\nmodel R0_leaf by box ventilation case (median over seasons):')
    print(res.groupby(['case', 'experiment'], sort=False).R0_leaf.median().unstack().round(1).to_string())
    hl = res[res.experiment == 'day_high_light'].set_index(['year', 'case']).R0_leaf
    ll = res[res.experiment == 'day_low_light'].set_index(['year', 'case']).R0_leaf
    print(f'\nhigh / low light R0_leaf: model median {(hl / ll).median():.2f} ({(hl / ll).min():.2f}-{(hl / ll).max():.2f}); '
          f'Table 4: {86 / 75:.2f}')


def figure(res):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
    colors = dict(zip(fb.RA_CASES, ['#52514e', '#2a78d6', '#eb6834', '#1baf7a']))
    for ax, (exp, (name, _, _)) in zip(axes, EXPERIMENTS.items()):
        o = observed(name).set_index('hours_after_start')
        med = res[res.experiment == exp].groupby('case', sort=False).median(numeric_only=True)
        for case, m in med.iterrows():
            ax.plot(SAMPLES, [m[f'leaf_TWT_{j}h'] for j in SAMPLES], '-o', ms=3, color=colors[case], label=case, linewidth=1.2)
        if exp != 'day_high_light':
            ax.errorbar(SAMPLES, [o.loc[j, 'TWT_leaves_Bq_per_mL'] for j in SAMPLES],
                        yerr=[o.loc[j, 'TWT_leaves_sd'] for j in SAMPLES], fmt='s', color='#0b0b0b', ms=5,
                        label='measured', zorder=5)
        else:
            ax.text(0.03, 0.04, 'no measured series (inferred case)', transform=ax.transAxes, fontsize=8, color='#52514e')
        ax.set_yscale('log')
        ax.set_xlabel('hours after start of the 2 h supply')
        ax.set_ylabel('leaf TFWT, Bq/mL')
        ax.set_title(exp.replace('_', ' '), fontsize=10)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=8, frameon=False, title='chamber ventilation', title_fontsize=8)
    fig.suptitle('Wheat chamber (Diabate & Strack 1997) driven by the measured hourly record: supplementary analysis after the test run',
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'figures', 'fig6_chamber_supplementary.png'), dpi=150)


def main():
    res = run()
    res.to_csv(os.path.join(OUT, 'supplementary_chamber.csv'), index=False)
    report(res)
    figure(res)


if __name__ == '__main__':
    main()
