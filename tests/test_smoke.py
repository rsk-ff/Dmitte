"""Smoke tests that run without the gitignored ./data directory."""
import os

import numpy as np
import pytest

from dmitte import para_constant
from dmitte.guassian import calc_sigmaxyz, gaussian_puff_model
from dmitte.transfer_equotions import solve_con

DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'data')


def test_package_imports():
    # run_wofost and calc_para pull in pcse and pandas; importing them catches version breaks
    from dmitte import calc_para, run_wofost  # noqa: F401


@pytest.mark.parametrize('stability', range(1, 7))
def test_sigmaxyz_positive(stability):
    sx, sy, sz = calc_sigmaxyz(stability, np.array([100.0, 1000.0]))
    assert np.all(sx > 0) and np.all(sy > 0) and np.all(sz > 0)
    assert np.all(sx == sy)


def test_sigmaxyz_invalid_stability():
    with pytest.raises(ValueError):
        calc_sigmaxyz(7, 100.0)


def test_gaussian_puff_centerline_peak():
    y = np.array([0.0, 50.0, 200.0])
    x = np.full_like(y, 500.0)
    z = np.full_like(y, 1.0)
    con = gaussian_puff_model(x, y, z, t=np.full_like(y, 600.0), Q_total=1e12, UU=2.0, HEG=50,
                              stability=4, t_release=600, puff_num=50)
    assert np.all(np.isfinite(con))
    assert np.all(con >= 0)
    assert con[0] > con[1] > con[2]


def test_solve_con_conserves_activity():
    # Only a2 -> s1 transfer (index 8); total activity changes only by decay
    n_steps = 24
    k_array = np.zeros((n_steps, 29))
    k_array[:, 8] = 0.1
    As = solve_con('HTO', 1.0, k_array, para_constant.LAMBDA_T_H)
    assert As.shape == (n_steps + 1, 10)
    totals = As[1:].sum(axis=1)
    decay = np.exp(-para_constant.LAMBDA_T_H * np.arange(n_steps))
    np.testing.assert_allclose(totals, totals[0] * decay, rtol=1e-6)
    assert As[1, 3] > 0


def test_solve_con_rejects_unknown_char():
    with pytest.raises(ValueError):
        solve_con('XX', 1.0, np.zeros((2, 29)), para_constant.LAMBDA_T_H)


@pytest.mark.skipif(not os.path.isdir(DATA_DIR), reason='./data not present')
def test_meteodata_loads(monkeypatch):
    from dmitte.calc_para import meteodata
    monkeypatch.chdir(os.path.join(DATA_DIR, '..'))
    df = meteodata('2021-04-20', '2021-04-22')
    assert {'Ta_Avg', 'WS10m_avg', 'RH_Avg'} <= set(df.columns)
