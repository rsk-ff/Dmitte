"""Shared machinery of the thesis analyses (analysis/*.py): seasons, frozen-rate baselines, release responses.

All analyses use the validated set-up of validation/run_validation.py (Wageningen weather with a diurnal
cycle, WOFOST potato 2004-2008 and winter wheat 2005-2008) and treat air-moisture HTO as a boundary
condition, as the validation experiments do.

Release response. For a one-hour exposure starting at hour j with air-moisture HTO 1 Bq/kg (time-integrated
concentration 1 Bq h kg-1), every harvest endpoint c^T y(harvest) is linear in the exposure, so the responses
for all release hours of a season follow from one backward (adjoint) sweep:
    w_n = c,   w_j = Z P_j^T Z w_{j+1},   response(j) = w_{j+1} . (Z P_j e_air) * W_air(j)
where P_j = expm(M_j) is the hourly propagator with the air row zeroed (air constant within the hour) and Z
sets the air compartment to zero (no HTO in air outside the exposure hour).
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd
from scipy.linalg import expm

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'validation'))
warnings.filterwarnings('ignore')

import run_validation as rv  # noqa: E402
from dmitte.para_constant import DOSE_COEF_HTO, DOSE_COEF_OBT  # noqa: E402
from dmitte.transfer_equotions import transfer_rates  # noqa: E402

OUT = os.path.join(ROOT, 'analysis', 'results')
FIG = os.path.join(OUT, 'figures')
A_A2, A_BH, A_FH, A_FO = rv.A_A2, rv.A_BH, rv.A_FH, rv.A_FO
CROPS = rv.CROPS

# rate variants: the dynamic model and baselines with its time variation frozen. 'release_frozen' (rates of
# the release hour held for the rest of the season) differs per release, so it has no adjoint sweep and is
# computed per release with Season.release_frozen_response.
VARIANTS = ('dynamic', 'diurnal_only', 'constant')
ALL_VARIANTS = VARIANTS + ('release_frozen',)
VARIANT_LABEL = {'dynamic': 'dynamic (hourly rates from the crop model)',
                 'diurnal_only': 'growth frozen (season-mean diurnal cycle of rates)',
                 'constant': 'constant (season-mean rates)',
                 'release_frozen': 'frozen at release (rates of the release hour all season)'}


class Season:
    """One crop season: WOFOST run, calc_para state, rates, pools and propagators of each variant."""

    def __init__(self, crop, year, rates=None, state=None):
        cfg = CROPS[crop]
        self.crop, self.year, self.cfg = crop, year, cfg
        if state is None:
            state = rv.model_state(crop, year, 'diurnal')
        self.pr, self.air, self.soil, self.plant = state
        self.index = pd.date_range(f'{year}-{cfg["window"][0]}', f'{year}-{cfg["window"][1]}', freq='h')
        self.pools = rv.season_arrays(self.air, self.soil, self.plant)
        k = transfer_rates('HTO', cfg['plant_type'], self.air, self.soil, self.plant, **(rates or {}))
        self.k = rv.clean_rates(k)[0]
        self.n = len(self.k)
        self.anthesis, self.lin_start, self.lin_end = rv.grain_phases(self.pr)
        self._P = {}

    def rates(self, variant):
        k = self.k
        if variant == 'dynamic':
            return k
        if variant == 'constant':
            return np.repeat(k.mean(axis=0, keepdims=True), self.n, axis=0)
        if variant == 'diurnal_only':
            hod = self.index[:self.n].hour
            mean = pd.DataFrame(k).groupby(hod).mean().values
            return mean[hod]
        raise ValueError(variant)

    def propagators(self, variant):
        if variant not in self._P:
            M = rv.generators(self.rates(variant))
            M[:, A_A2, :] = 0
            self._P[variant] = expm(M)
        return self._P[variant]

    def fresh_edible(self):
        """Fresh mass of the edible part at harvest, kg m-2."""
        return self.pools['comb_edible'][-1] / rv.COMB_WATER / self.cfg['drymatter']

    def endpoints(self):
        """Harvest endpoints as linear functionals of the state (rows) and their names."""
        n = len(self.pools['w_plant']) - 1
        c = {}
        e = np.zeros(rv.N_COMP)
        e[A_FO] = 1 / self.pools['comb_edible'][n]
        c['edible_OBT'] = e                                   # Bq per kg combustion water
        e = np.zeros(rv.N_COMP)
        e[A_FH] = 1 / self.pools['w_fruit'][n] if self.pools['w_fruit'][n] > 0 else 0
        c['edible_TFWT'] = e                                  # Bq per kg tissue water
        fresh = self.fresh_edible()
        e = np.zeros(rv.N_COMP)
        e[A_FO], e[A_FH] = 1 / fresh, 1 / fresh
        c['edible_fresh'] = e                                 # Bq per kg fresh weight
        e = np.zeros(rv.N_COMP)
        e[A_FO], e[A_FH] = DOSE_COEF_OBT / fresh, DOSE_COEF_HTO / fresh
        c['dose_per_kg'] = e                                  # mSv per kg eaten
        return list(c), np.array(list(c.values()))

    def release_response(self, variant):
        """Harvest endpoints for a 1 h exposure (air moisture 1 Bq/kg) starting at every hour of the season.

        Returns an array (hours, endpoints); response per unit time-integrated air-moisture concentration.
        """
        P = self.propagators(variant)
        names, C = self.endpoints()
        W = C.T.copy()                                        # (N_COMP, n_endpoints), w_n
        out = np.zeros((self.n, len(names)))
        for j in range(self.n - 1, -1, -1):
            col = P[j][:, A_A2] * self.pools['w_air'][j]
            col[A_A2] = 0.0
            out[j] = col @ W
            W[A_A2] = 0.0
            W = P[j].T @ W
            W[A_A2] = 0.0
        return names, out

    def release_frozen_response(self, j, hours=1):
        """Harvest endpoints for 1 h exposures in hours j, j+1, ..., j+hours-1 with the rates of hour j held
        constant to harvest (a quasi-static model evaluated with the conditions at release). Returns (hours, endpoints)."""
        M = rv.generators(self.k[j:j + 1])[0]
        M[A_A2, :] = 0
        names, C = self.endpoints()
        P1 = expm(M)
        out = np.zeros((hours, len(names)))
        for t in range(min(hours, self.n - j)):
            v = P1[:, A_A2] * self.pools['w_air'][j + t]
            v[A_A2] = 0.0
            out[t] = C @ (expm(M * (self.n - j - t - 1)) @ v)
        return names, out

    def forward(self, variant, air_moist, i0):
        """State from hour i0 to harvest with air-moisture HTO air_moist[t] (Bq/kg) in hours i0, i0+1, ..."""
        P = self.propagators(variant)
        y = np.zeros(rv.N_COMP)
        ys = [y.copy()]
        for j in range(i0, self.n):
            t = j - i0
            y[A_A2] = self.pools['w_air'][j] * air_moist[t] if t < len(air_moist) else 0.0
            y = P[j] @ y
            y[A_A2] = 0.0
            ys.append(y.copy())
        return np.array(ys)

    def harvest(self, y):
        names, C = self.endpoints()
        return dict(zip(names, C @ y))

    def days_after_anthesis(self, j):
        return (self.index[j] - self.anthesis) / pd.Timedelta(days=1) if self.anthesis is not None else np.nan


def seasons(crops=CROPS, rates=None):
    for crop in crops:
        for year in CROPS[crop]['years']:
            yield Season(crop, year, rates)
