"""Drive the wheat field-box experiments (scenario card WHEAT_FB) with their measured conditions.

TECDOC-1738 Table 22 / Strack et al. (2005): 14 one-hour HTO exposures of field-grown winter wheat in a
ventilated acrylic box (30 x 30 x 100 cm), 07:00-23:00, 12-28 days after anthesis. For each exposure and
each WOFOST wheat season of run_validation.py (2005-2008):

* the exposure starts at the experiment's hour on the day with the same number of days after anthesis
  (anthesis = WOFOST DVS 1);
* the exposure hour uses the measured box temperature, RH and PPFD (global radiation = PPFD / 2.2,
  i.e. PAR = 0.48 of global radiation and 4.57 umol/J) and no rain; the hour after uses the same values
  (the plants stay in the field, outdoor conditions are not reported) with the site wind;
* air HTO follows the measured time course in the box (Fig. 26, read every 5 min; experiments without a
  curve use the mean shape of the five curves), normalised to its 1 h mean, and the exposure hour is
  integrated in 5 min steps;
* the box air speed is not reported, so the aerodynamic resistance RAM + RB during the exposure hour is
  scanned: the site wind, and 30, 15 and 8 s/m.

Model output is compared exposure by exposure with leaf TFWT / mean air HTO at the end of exposure
(Fig. 25), leaf OBT relative to that TFWT 1 h, 2 h, 4 h and 1 d after the start (Fig. 27), grain TLI at
harvest, and by time of day with the leaf TFWT half-lives within the hour after exposure (Table 23).
All of these are calibration data (benchmarks/split.csv).

Run from the repository root after validation/build_inputs.py:
    python validation/scenario_field_box.py
writes validation/results/scenario_field_box.csv and validation/results/figures/fig5_field_box.png.
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
from dmitte import calc_para  # noqa: E402
from dmitte.calc_para import Calc_air, Calc_plant, Calc_soil  # noqa: E402
from dmitte.transfer_equotions import transfer_rates  # noqa: E402

DATA = os.path.join(HERE, 'benchmarks', 'datasets')
OUT = os.path.join(HERE, 'results')
GLOBAL_PER_PPFD = 1 / 2.2          # W m-2 global radiation per umol m-2 s-1 PPFD
RA_CASES = {'site wind': None, 'RAM+RB 30 s/m': 30.0, 'RAM+RB 15 s/m': 15.0, 'RAM+RB 8 s/m': 8.0}
SUBSTEPS = 12                      # 5 min steps in the exposure hour
# time-of-day classes of TECDOC-1738 Table 23 (leaf TFWT half-life within the hour after exposure, min)
HALF_LIFE = {'dawn': (40, 60), 'day': (25, 49), 'dusk': (230, 660), 'night': (110, 170)}


def time_class(hour):
    return 'dawn' if hour <= 8 else 'day' if hour <= 15 else 'dusk' if hour <= 20 else 'night'


def air_shapes(exps):
    """Air HTO in each 5 min step relative to the 1 h mean (Fig. 26; mean shape where no curve exists)."""
    a = pd.read_csv(os.path.join(DATA, 'wheat_DE_field_box_air_HTO.csv'))
    shapes = {}
    for e, g in a.groupby('experiment'):
        v = g.sort_values('minute').air_HTO_kBq_per_mL.values
        step = (v[:-1] + v[1:]) / 2
        shapes[e] = step / step.mean()
    mean_shape = np.mean(list(shapes.values()), axis=0)
    return {e: shapes.get(e, mean_shape) for e in exps}, set(shapes)


def season(year):
    cfg = rv.CROPS['wheat']
    start, end = (f'{year}-{md}' for md in cfg['window'])
    plantmodel = [*cfg['crop'], 'ec2.soil', cfg['agro'].format(y=year)]
    pr = rv.run_wofost.plantModel(cfg['plant_type'], plantmodel, [start, end])
    forcing = rv.diurnal_meteo(start, end, pr['SM'])
    anthesis, _, _ = rv.grain_phases(pr)
    return pr, forcing, start, end, anthesis


def model_state(pr, forcing, start, end):
    stock = calc_para.meteodata
    calc_para.meteodata = lambda s, e: forcing
    try:
        args = (pr, start, end, 'CEREAL', rv.STABILITY, rv.HEG, rv.X)
        return Calc_air(*args), Calc_soil(*args), Calc_plant(*args, rv.CROPS['wheat']['drymatter'])
    finally:
        calc_para.meteodata = stock


def set_box_conditions(forcing, i0, exp, u10_box):
    f = forcing.copy()
    es = 0.6108 * np.exp(17.27 * exp.temperature_C / (exp.temperature_C + 237.3))
    for i in (i0, i0 + 1):
        f.iloc[i, f.columns.get_loc('Ta_Avg')] = exp.temperature_C
        f.iloc[i, f.columns.get_loc('RH_Avg')] = exp.RH_pct
        f.iloc[i, f.columns.get_loc('Pvapor_Avg')] = exp.RH_pct / 100 * es
        f.iloc[i, f.columns.get_loc('DR_Avg')] = exp.PPFD_umol_m2_s * GLOBAL_PER_PPFD
        f.iloc[i, f.columns.get_loc('Rain_Tot')] = 0.0
    if u10_box is not None:
        f.iloc[i0, f.columns.get_loc('WS10m_avg')] = u10_box
    return f


def run_exposure(k, pools, i0, shape):
    """Inventories (Bq m-2) from the start of exposure (row 0) to harvest, air at concentration shape(t)."""
    M = rv.generators(k[i0:])
    M[:, rv.A_A2, :] = 0                       # air is a boundary condition
    P_sub = expm(M[0] / SUBSTEPS)
    P = expm(M[1:])
    y = np.zeros(10)
    ys = [y.copy()]
    for c in shape:
        y[rv.A_A2] = pools['w_air'][i0] * c
        y = P_sub @ y
    y[rv.A_A2] = 0.0
    ys.append(y.copy())
    for Pi in P:
        y = Pi @ y
        y[rv.A_A2] = 0.0
        ys.append(y.copy())
    return np.array(ys)


def metrics(ys, pools, i0):
    n = len(pools['w_plant'])
    leaf = lambda j: ys[j, rv.A_BH] / pools['w_plant'][min(i0 + j, n - 1)]
    c_end = leaf(1)
    obt = lambda j: 100 * ys[j, rv.A_BO] / pools['comb_body'][min(i0 + j, n - 1)] / c_end
    hv = len(ys) - 1
    return dict(R0_leaf=100 * c_end,
                k_after=np.log(leaf(1) / leaf(2)),
                obt_1h=obt(1), obt_2h=obt(2), obt_4h=obt(4), obt_1d=obt(24),
                TLI=100 * ys[hv, rv.A_FO] / pools['comb_edible'][n - 1] / c_end)


def main():
    exps = pd.read_csv(os.path.join(DATA, 'wheat_DE_field_box_1h.csv'))
    shapes, measured_curve = air_shapes(exps.experiment)
    rows = []
    for year in rv.CROPS['wheat']['years']:
        pr, forcing, start, end, anthesis = season(year)
        idx = forcing.index
        air0 = model_state(pr, forcing, start, end)[0]
        ra_times_u = (air0.RAM + air0.RB) * air0.U10      # RAM + RB is proportional to 1 / U10
        for _, e in exps.iterrows():
            t0 = anthesis.normalize() + pd.Timedelta(days=int(e.days_after_anthesis), hours=int(e.start_hour))
            i0 = idx.get_loc(t0)
            for case, ra in RA_CASES.items():
                u10 = None if ra is None else ra_times_u[i0] / ra
                air, soil, plant = model_state(pr, set_box_conditions(forcing, i0, e, u10), start, end)
                k = transfer_rates('HTO', 'CEREAL', air, soil, plant)
                pools = rv.season_arrays(air, soil, plant)
                ys = run_exposure(k, pools, i0, shapes[e.experiment])
                rows.append(dict(experiment=e.experiment, year=year, case=case, time_class=time_class(e.start_hour),
                                 date=t0, ra_exposure=float(air.RAM[i0] + air.RB[i0]),
                                 rc_exposure=float(plant.RC1[i0] * 100), measured_air_curve=e.experiment in measured_curve,
                                 **metrics(ys, pools, i0)))
        print('wheat', year, 'done', flush=True)
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(OUT, 'scenario_field_box.csv'), index=False)
    report(res, exps)
    figure(res, exps)


OBS = {'R0_leaf': 'leaf_TFWT_rel_end_pct', 'obt_1h': 'leaf_OBT_rel_1h_pct', 'obt_2h': 'leaf_OBT_rel_2h_pct',
       'obt_4h': 'leaf_OBT_rel_4h_pct', 'obt_1d': 'leaf_OBT_rel_1d_pct', 'TLI': 'grain_TLI_harvest_pct'}


def paired(res, exps):
    med = res.groupby(['case', 'experiment'], sort=False).median(numeric_only=True).reset_index()
    return med.merge(exps, on='experiment')


def report(res, exps):
    m = paired(res, exps)
    print('\nmodel (median over 4 seasons) / observed, by exposure:')
    for case, g in m.groupby('case', sort=False):
        t = pd.DataFrame({'exp': g.experiment, 'hour': g.start_hour})
        for metric, col in OBS.items():
            t[metric] = (g[metric] / g[col]).round(2).values
        print(f'--- {case} (RAM+RB at exposure: median {g.ra_exposure.median():.0f} s/m)')
        print(t.to_string(index=False))
    print('\nsummary: geometric mean of model/observed, share within a factor 2, Spearman r across exposures')
    out = []
    for case, g in m.groupby('case', sort=False):
        for metric, col in OBS.items():
            ok = g[col].notna()
            r = g.loc[ok, metric] / g.loc[ok, col]
            out.append(dict(case=case, metric=metric, n=int(ok.sum()), gmean_ratio=float(np.exp(np.log(r).mean())),
                            within_2x=float(((r > 0.5) & (r < 2)).mean()),
                            spearman=float(g.loc[ok, [metric, col]].corr(method='spearman').iloc[0, 1])))
    print(pd.DataFrame(out).round(2).to_string(index=False))
    print('\nleaf TFWT loss rate within 1 h after exposure (1/h), median by time of day; Table 23 range:')
    k = res.groupby(['case', 'time_class'], sort=False).k_after.median().unstack()
    lit = {c: f'{np.log(2) / (hi / 60):.2f}-{np.log(2) / (lo / 60):.2f}' for c, (lo, hi) in HALF_LIFE.items()}
    k.loc['Table 23'] = pd.Series(lit)
    print(k[list(HALF_LIFE)].to_string())
    day = m.start_hour.between(9, 15)
    night = m.start_hour == 23
    print('\nB49 (day 09-15 h, 54-92 %) and B50 (night 23 h, 17-19 %), mean of the exposure medians:')
    print(m.assign(group=np.where(day, 'B49 day', np.where(night, 'B50 night', '')))
          .query('group != ""').groupby(['case', 'group'], sort=False).R0_leaf.mean().round(1).unstack().to_string())


def figure(res, exps):
    m = paired(res, exps)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    colors = dict(zip(RA_CASES, ['#52514e', '#2a78d6', '#eb6834', '#1baf7a']))
    for ax, (metric, label) in zip(axes, [('R0_leaf', 'leaf TFWT / mean air HTO at end, %'),
                                          ('obt_1h', 'leaf OBT / leaf TFWT at end, %'),
                                          ('TLI', 'grain TLI at harvest, %')]):
        col = OBS[metric]
        for case, g in m.groupby('case', sort=False):
            ax.scatter(g[col], g[metric], s=22, color=colors[case], label=case, alpha=0.85, linewidths=0)
        lo = min(m[col].min(), m[metric].min()) * 0.7
        hi = max(m[col].max(), m[metric].max()) * 1.4
        x = np.array([lo, hi])
        ax.plot(x, x, color='#0b0b0b', linewidth=1)
        ax.fill_between(x, x / 2, x * 2, color='#c9c8c2', alpha=0.4, linewidth=0)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_xlabel('observed')
        ax.set_ylabel('model')
        ax.set_title(label, fontsize=10)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    axes[0].legend(fontsize=8, frameon=False, title='box ventilation', title_fontsize=8)
    fig.suptitle('Wheat field box, 14 one-hour exposures driven by the measured conditions (grey: factor 2)', fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'figures', 'fig5_field_box.png'), dpi=150)


if __name__ == '__main__':
    main()
