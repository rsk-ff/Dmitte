"""The adjoint release response of analysis/common.py equals forward runs (needs ./data, built by validation/build_inputs.py)."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.join(os.path.dirname(__file__), '..')
DATA_DIR = os.path.join(ROOT, 'data')
pytestmark = pytest.mark.skipif(not os.path.isdir(os.path.join(DATA_DIR, 'agro')), reason='./data not present')


@pytest.fixture(scope='module')
def season():
    sys.path.insert(0, os.path.join(ROOT, 'analysis'))
    import common
    return common, common.Season('wheat', 2006)


@pytest.mark.parametrize('variant', ['dynamic', 'diurnal_only', 'constant'])
def test_release_response_matches_forward(season, variant):
    common, s = season
    names, R = s.release_response(variant)
    for j in (400, 1500, 2100):
        h = s.harvest(s.forward(variant, [1.0], j)[-1])
        np.testing.assert_allclose(R[j], [h[n] for n in names], rtol=1e-9)


def test_superposition_of_hourly_responses(season):
    common, s = season
    names, R = s.release_response('dynamic')
    air = np.array([0.3, 1.0, 0.5])
    h = s.harvest(s.forward('dynamic', air, 1500)[-1])
    np.testing.assert_allclose(air @ R[1500:1503], [h[n] for n in names], rtol=1e-9)


def test_frozen_variants_have_no_seasonal_variation(season):
    common, s = season
    k = s.rates('constant')
    assert np.allclose(k, k[0])
    d = s.rates('diurnal_only')
    assert np.allclose(d[24:48], d[48:72])


def test_release_frozen_matches_forward_with_constant_rates(season):
    common, s = season
    from scipy.linalg import expm
    j = 1500
    M = common.rv.generators(s.k[j:j + 1])[0]
    M[common.A_A2, :] = 0
    s._P['frozen_test'] = np.repeat(expm(M)[None], s.n, axis=0)
    names, R = s.release_frozen_response(j, hours=2)
    h = s.harvest(s.forward('frozen_test', [1.0, 1.0], j)[-1])
    np.testing.assert_allclose(R.sum(axis=0), [h[n] for n in names], rtol=1e-6)
