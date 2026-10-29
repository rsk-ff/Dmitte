"""Hypothetical accidental release (EMRAS TC14 chapter 9 style), dynamic versus constant rates (thesis analysis B2).

Scenario settings are those of the repository's case scripts (iaea-case1-HTO.py, cereals_HTO.py,
potato_HTO.py): 3.7e15 Bq HTO released over 1 h at an effective height of 20 m, receptors 1, 3, 10 and 30 km
downwind on the plume axis, and three weather cases

    Case 1: stability A (1), wind 2 m/s, release at 12:00
    Case 2: stability D (4), wind 5 m/s, release at 12:00
    Case 3: stability F (6), wind 2 m/s, release at 00:00

(A occurs only by day and F only by night; D is taken by day). Air HTO at 1 m height follows
dmitte.guassian.gaussian_puff_model, is integrated per hour, and is converted to air-moisture HTO with the
absolute humidity of the site weather of that hour. The release is placed on every day between emergence and
maturity of each WOFOST season, so the scenario result is a distribution over release dates; the crop
response uses the harvest endpoints of common.Season (superposition of the one-hour responses, exact because
the model is linear in the air concentration), for the variants of common.ALL_VARIANTS.

Coupling of dispersion and compartment model (design review items 13 and 14): the plume concentration is a
forcing of the air compartment (boundary condition, hour by hour over the whole passage), not an initial
inventory, so no mixing-layer height, integration height or forced 1 h cut-off enters. The stability class of a
case is used only for the dispersion; the compartment model's own air parameters (RAM, deposition) come from the
hourly site weather with neutral stability (run_validation.STABILITY = 4) for the whole season, including the
plume hours. Limitations: the dispersion uses the case wind speed while leaf exchange uses the site wind of that
hour, and HTO re-emitted from soil and plants is not returned to the air (no local air compartment).

Run from the repository root after validation/build_inputs.py:
    python analysis/emras_scenarios.py
writes analysis/results/emras_scenarios.csv.gz, emras_summary.csv and figures/fig_b2_emras.png.
"""
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import common as c  # noqa: E402
from dmitte.guassian import gaussian_puff_model  # noqa: E402

Q_TOTAL = 3.7e15            # Bq
T_RELEASE = 3600            # s
HEG = 20                    # m
Z_RECEPTOR = 1.0            # m
DISTANCES = (1e3, 3e3, 10e3, 30e3)
CASES = {'Case 1 (A, 2 m/s, 12:00)': (1, 2.0, 12),
         'Case 2 (D, 5 m/s, 12:00)': (4, 5.0, 12),
         'Case 3 (F, 2 m/s, 00:00)': (6, 2.0, 0)}
CONSUMPTION = {'potato': 100.0, 'wheat': 100.0}   # kg per year eaten, for the illustrative dose only


def hourly_air(stability, u, x, hours):
    """Hourly mean air HTO (Bq m-3) at the receptor, hour 0 = the release hour."""
    t = np.arange(0, hours * 3600, 60.0) + 30.0
    con = gaussian_puff_model(np.full_like(t, x), 0.0, Z_RECEPTOR, t, Q_total=Q_TOTAL, UU=u, HEG=HEG,
                              stability=stability, t_release=T_RELEASE, puff_num=120)
    return con.reshape(hours, 60).mean(axis=1)


def absolute_humidity(s):
    """kg water per m3 air, hourly, from the site forcing of the season."""
    return np.asarray(s.air.a_h, dtype=float)[:s.n] / 1000         # calc_para a_h is in g m-3


def run():
    plumes = {}
    for case, (stab, u, _) in CASES.items():
        for x in DISTANCES:
            hours = int(np.ceil((x / u + T_RELEASE) / 3600)) + 3
            plumes[case, x] = hourly_air(stab, u, x, hours)
    rows = []
    for s in c.seasons():
        dvs = s.pr['DVS'].reindex(s.index[:s.n]).values
        ah = absolute_humidity(s)
        for v in c.VARIANTS:
            names, R = s.release_response(v)
            for case, (stab, u, hour) in CASES.items():
                starts = [j for j in range(s.n) if s.index[j].hour == hour and np.isfinite(dvs[j]) and 0 <= dvs[j] < 2]
                for x in DISTANCES:
                    air = plumes[case, x]
                    for j in starts:
                        h = np.arange(j, min(j + len(air), s.n))
                        moist = air[:len(h)] / ah[h]                     # Bq per kg air moisture
                        out = moist @ R[h]
                        rows.append(dict(crop=s.crop, year=s.year, variant=v, case=case, distance_km=x / 1e3,
                                         release=s.index[j], days_after_anthesis=s.days_after_anthesis(j),
                                         TIC_air_Bq_h_m3=air.sum(), TIC_moist_Bq_h_kg=moist.sum(),
                                         **dict(zip(names, out))))
        for case, (stab, u, hour) in CASES.items():
            starts = [j for j in range(s.n) if s.index[j].hour == hour and np.isfinite(dvs[j]) and 0 <= dvs[j] < 2]
            for j in starts:
                longest = max(len(plumes[case, x]) for x in DISTANCES)
                names, Rf = s.release_frozen_response(j, longest)
                for x in DISTANCES:
                    air = plumes[case, x]
                    h = np.arange(j, min(j + len(air), s.n))
                    moist = air[:len(h)] / ah[h]
                    rows.append(dict(crop=s.crop, year=s.year, variant='release_frozen', case=case, distance_km=x / 1e3,
                                     release=s.index[j], days_after_anthesis=s.days_after_anthesis(j),
                                     TIC_air_Bq_h_m3=air.sum(), TIC_moist_Bq_h_kg=moist.sum(),
                                     **dict(zip(names, moist @ Rf[:len(h)]))))
        print(s.crop, s.year, 'done', flush=True)
    d = pd.DataFrame(rows)
    d['annual_dose_mSv'] = d.dose_per_kg * d.crop.map(CONSUMPTION)
    return d


def summary(d):
    g = d.groupby(['crop', 'case', 'distance_km', 'variant'], sort=False)
    s = g.agg(TIC_moist=('TIC_moist_Bq_h_kg', 'median'),
              dose_p5=('annual_dose_mSv', lambda x: x.quantile(.05)),
              dose_median=('annual_dose_mSv', 'median'),
              dose_p95=('annual_dose_mSv', lambda x: x.quantile(.95)),
              dose_max=('annual_dose_mSv', 'max'),
              edible_fresh_median=('edible_fresh', 'median')).reset_index()
    dyn = s[s.variant == 'dynamic'].set_index(['crop', 'case', 'distance_km'])
    key = lambda r: (r.crop, r.case, r.distance_km)
    s['median_vs_dynamic'] = s.apply(lambda r: r.dose_median / dyn.loc[key(r), 'dose_median'], axis=1)
    s['max_vs_dynamic'] = s.apply(lambda r: r.dose_max / dyn.loc[key(r), 'dose_max'], axis=1)
    return s


def figure(d):
    sub = d[d.distance_km == 1]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4))
    colors = {'dynamic': '#2a78d6', 'diurnal_only': '#eb6834', 'constant': '#52514e', 'release_frozen': '#1baf7a'}
    for ax, crop in zip(axes, c.CROPS):
        g = sub[sub.crop == crop]
        pos = 0
        ticks, labels = [], []
        for case in CASES:
            for k, v in enumerate(c.ALL_VARIANTS):
                vals = g[(g.case == case) & (g.variant == v)].annual_dose_mSv.values
                b = ax.boxplot(vals, positions=[pos + k * 0.8], widths=0.65, patch_artist=True, showfliers=False)
                for p in b['boxes']:
                    p.set_facecolor(colors[v])
                    p.set_alpha(0.55)
                for m in b['medians']:
                    m.set_color('#0b0b0b')
            ticks.append(pos + 1.2)
            labels.append(case.replace(' (', '\n('))
            pos += 4.0
        ax.set_xticks(ticks)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_yscale('log')
        ax.set_ylabel(f'ingestion dose, mSv ({CONSUMPTION[crop]:.0f} kg eaten)')
        ax.set_title(f'{crop} at 1 km: distribution over release dates and seasons', fontsize=10)
        for s in ('top', 'right'):
            ax.spines[s].set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[v], alpha=0.55) for v in c.ALL_VARIANTS]
    axes[0].legend(handles, [c.VARIANT_LABEL[v] for v in c.ALL_VARIANTS], fontsize=7.5, frameon=False, loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(c.FIG, 'fig_b2_emras.png'), dpi=150)


def main():
    os.makedirs(c.FIG, exist_ok=True)
    d = run()
    d.to_csv(os.path.join(c.OUT, 'emras_scenarios.csv.gz'), index=False)
    s = summary(d)
    s.to_csv(os.path.join(c.OUT, 'emras_summary.csv'), index=False)
    pd.set_option('display.width', 220)
    print('\nannual ingestion dose (mSv, 100 kg eaten) over release dates and seasons; TIC in Bq h per kg air moisture:')
    print(s.to_string(index=False, float_format=lambda x: f'{x:.3g}'))
    figure(d)


if __name__ == '__main__':
    main()
