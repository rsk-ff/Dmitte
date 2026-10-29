"""Smoke tests that run without the gitignored ./data directory."""
import os

import numpy as np
import pytest

from dmitte import para_constant
from dmitte.guassian import calc_sigmaxyz, gaussian_puff_model
from dmitte.transfer_equotions import A_NAMES, K_NAMES, solve_con

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
    k_array = np.zeros((n_steps, len(K_NAMES)))
    k_array[:, 8] = 0.1
    As = solve_con('HTO', 1.0, k_array, para_constant.LAMBDA_T_H)
    assert As.shape == (n_steps + 1, len(A_NAMES))
    totals = As[1:].sum(axis=1)
    decay = np.exp(-para_constant.LAMBDA_T_H * np.arange(n_steps))
    np.testing.assert_allclose(totals, totals[0] * decay, rtol=1e-6)
    assert As[1, 3] > 0


def test_solve_con_rejects_unknown_char():
    with pytest.raises(ValueError):
        solve_con('XX', 1.0, np.zeros((2, len(K_NAMES))), para_constant.LAMBDA_T_H)


@pytest.mark.skipif(not os.path.isdir(DATA_DIR), reason='./data not present')
def test_meteodata_loads(monkeypatch):
    from dmitte.calc_para import meteodata
    monkeypatch.chdir(os.path.join(DATA_DIR, '..'))
    df = meteodata('2021-04-20', '2021-04-22')
    assert {'Ta_Avg', 'WS10m_avg', 'RH_Avg'} <= set(df.columns)


def test_wofost_conf_outputs_what_calc_plant_reads():
    # Calc_plant reads PGASS and PMRES, which the stock PCSE config does not output
    from pcse.engine import Engine  # noqa: F401
    from dmitte import run_wofost
    conf = {}
    with open(run_wofost.WOFOST_CONF) as f:
        exec(compile(f.read(), run_wofost.WOFOST_CONF, 'exec'), conf)
    assert {'PGASS', 'PMRES', 'LAI', 'TAGP', 'TWSO', 'TWLV', 'SM', 'WWLOW'} <= set(conf['OUTPUT_VARS'])
    assert run_wofost.WOFOST72_CROP_REPO.endswith('/wofost72/')


def test_benchmark_test_split_is_frozen():
    # The test benchmarks are evaluated once after model development; editing them changes the fingerprint
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'validation'))
    import pandas as pd
    import analyze
    bench = pd.read_csv(os.path.join(analyze.BENCH_DIR, 'literature_benchmarks.csv')).set_index('id')
    bench['id'] = bench.index
    test_ids, test_files = analyze.load_split()
    assert len(test_ids) == 15 and len(test_files) == 4
    assert analyze.test_fingerprint(bench, test_ids, test_files) == analyze.FROZEN_TEST_FINGERPRINT
    # every check refers to an existing benchmark, and no test benchmark is missing a split entry
    assert {c['benchmark'] for c in analyze.CHECKS} <= set(bench.index)
    split = pd.read_csv(os.path.join(analyze.BENCH_DIR, 'split.csv'))
    assert set(split.id[~split.id.str.startswith('datasets/')]) == set(bench.index)


def test_solve_con_is_exact_for_stiff_rates():
    # a2 <-> s1 exchange much faster than the hourly step: equilibrium within the first hour
    k_array = np.zeros((3, len(K_NAMES)))
    k_array[:, 8] = 5000.0                      # ka2_s1
    k_array[:, 7] = 1000.0                      # ks1_a2
    As = solve_con('HTO', 1.0, k_array, 0.0)
    np.testing.assert_allclose(As[1, [2, 3]], [0.0, 5 / 6], atol=1e-9)   # air reset to 0 after the first hour
    np.testing.assert_allclose(As[2, [2, 3]], [5 / 6 / 6, 5 / 6 * 5 / 6], rtol=1e-9)
