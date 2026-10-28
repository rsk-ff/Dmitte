"""Conservation, rate-sign, leaf-exchange and plant-process tests (review step 2, item 8 and step 3)."""
from types import SimpleNamespace

import numpy as np
import pytest

from dmitte import para_constant
from dmitte.calc_para import Calc_plant
from dmitte.transfer_equotions import (K_INDEX, K_NAMES, ROOT_LAYER_SHARE, _rate_of_change, _split_darcy,
                                       check_rates, leaf_exchange_rate, net_transpiration, solve_con,
                                       transfer_equotions, transfer_rates, transfer_rates_baomi,
                                       transfer_rates_Korea, transfer_rates_morris, water_balance_residual)

SYMMETRIC = [transfer_rates, transfer_rates_Korea, transfer_rates_baomi, transfer_rates_morris]
EXTERNAL_SINKS = ['ka1_a1', 'ka2_a2', 'ks3_s3']    # HT/HTO dispersion out of the column, deep percolation


def stub_state(n=4, friut_w=1.0, plant_w=2.0, rho_v=12.0):
    ones = np.ones(n)
    ml = para_constant.mixLayerH[3]
    air = SimpleNamespace(ML=ml, rainfall=0 * ones, atm_h=ml * rho_v / 1000 * 0.11 * ones,
                          U10=2.0 * ones, sta=4, ICOMP=para_constant.Plant_dict['ROOT_VEG']['ICOMP'])
    soil = SimpleNamespace(VDSO=4e-3 * ones, ESOIL=0.1 * ones, Va_b1=0.01 * ones, Va_b2=0.01 * ones,
                           soil1w=15.0 * ones, soil2w=30.0 * ones, soil3w=45.0 * ones)
    plant = SimpleNamespace(VDPF=1e-2 * ones, ROW=17.3 * ones, ETRM=0.2 * ones,
                            plant_w=plant_w * ones, plant_wh=0.11 * plant_w * ones,
                            friut_w=friut_w * ones, friut_wh=0.11 * friut_w * ones, TROBT=1e-4 * ones,
                            TROBT_FRUIT=2e-4 * ones)
    return air, soil, plant


SOIL_W = (15.0, 30.0, 45.0)


def root_uptake(k):
    '''Water flux from the three soil layers into plant and fruit water, kg m-2 h-1'''
    to_bh = sum(k[:, K_INDEX[f'ks{layer}_bh']] * w for layer, w in zip('123', SOIL_W))
    to_fh = sum(k[:, K_INDEX[f'ks{layer}_fh']] * w for layer, w in zip('123', SOIL_W))
    return to_bh, to_fh


def test_ode_loses_activity_only_through_external_sinks():
    rng = np.random.default_rng(0)
    rates = rng.random(len(K_NAMES))
    for name in EXTERNAL_SINKS:
        rates[K_INDEX[name]] = 0
    for _ in range(5):
        y = rng.random(10)
        dy = transfer_equotions(0, y, rates, 0.0)
        assert sum(dy) == pytest.approx(0.0, abs=1e-12)


def test_plant_to_soil_transfer_lands_in_soil_layer_1():
    rates = np.zeros(len(K_NAMES))
    rates[K_INDEX['kbh_so']] = 0.1
    y = np.zeros(10)
    y[6] = 1.0                                         # A_bh
    dy = transfer_equotions(0, y, rates, 0.0)
    assert dy[6] == pytest.approx(-0.1)
    assert dy[3] == pytest.approx(0.1)                 # A_s1


def test_check_rates_rejects_negative_and_nonfinite():
    k = np.zeros((3, len(K_NAMES)))
    check_rates(k)
    k[1, K_INDEX['kbh_bo']] = -1e-5
    k[2, K_INDEX['kfh_bh']] = np.inf
    with pytest.raises(ValueError, match='kbh_bo: 1 negative.*kfh_bh: 1 non-finite'):
        check_rates(k)
    with pytest.raises(ValueError):
        solve_con('HTO', 1.0, k, para_constant.LAMBDA_T_H)


def test_split_darcy_turns_upward_flow_into_a_rate_from_the_lower_layer():
    w1, w2 = np.array([15.0, 15.0]), np.array([30.0, 30.0])
    k_down, k_up = _split_darcy(np.array([0.02, -0.02]), 1e-3, w1, w2)
    np.testing.assert_allclose(k_down, [0.02, 0.0])
    np.testing.assert_allclose(k_up, [1e-3, 1e-3 + 0.02 * 15 / 30])


@pytest.mark.parametrize('rate_fn', SYMMETRIC, ids=lambda f: f.__name__)
def test_rates_are_valid_before_the_fruit_and_the_plant_exist(rate_fn):
    air, soil, plant = stub_state(n=6)
    plant.friut_w = plant.friut_wh = np.array([0.0, 0.0, 1.0, 1.0, 1.0, 1.0])
    plant.plant_w = plant.plant_wh = np.array([0.0, 2.0, 2.0, 2.0, 2.0, 2.0])
    k_ref = rate_fn('HTO', 'ROOT_VEG', *stub_state(n=6))
    k = rate_fn('HTO', 'ROOT_VEG', air, soil, plant)
    check_rates(k)
    # no plant yet: nothing goes into or out of the plant compartments
    plant_cols = [K_INDEX[n] for n in K_NAMES if n.endswith(('_bh', '_fh')) or n.startswith(('kbh', 'kfh', 'kbo'))]
    assert np.all(k[0, plant_cols] == 0)
    # plant but no fruit yet: root uptake all goes to the plant water, including what the fruit would take
    for layer in '123':
        assert k[1, K_INDEX[f'ks{layer}_fh']] == 0
    to_bh, _ = root_uptake(k)
    growth = _rate_of_change(plant.plant_w) + _rate_of_change(plant.friut_w)
    expected = net_transpiration(air, plant) + k[:, K_INDEX['kbh_so']] * plant.plant_w + growth
    assert to_bh[1] == pytest.approx(expected[1])
    # once both exist and stop growing, the rates are those of a plant that was always there
    np.testing.assert_allclose(k[4:], k_ref[4:])


@pytest.mark.parametrize('rate_fn', SYMMETRIC, ids=lambda f: f.__name__)
@pytest.mark.parametrize('plant_type', ['CEREAL', 'ROOT_VEG', 'LEAFY_PLANT'])
def test_plant_and_fruit_water_balance_closes(rate_fn, plant_type):
    # growing plant and fruit, transpiring (rho_v 9 g/m3 < rho_sat 17.3 g/m3)
    n = 48
    air, soil, plant = stub_state(n=n, rho_v=9.0)
    plant.plant_w = 1.0 + 0.01 * np.arange(n)
    plant.friut_w = 0.5 + 0.02 * np.arange(n)
    plant.plant_wh, plant.friut_wh = 0.11 * plant.plant_w, 0.11 * plant.friut_w
    k = rate_fn('HTO', plant_type, air, soil, plant)
    check_rates(k)
    res = water_balance_residual(k, air, soil, plant)
    np.testing.assert_allclose(res['bh']['residual'], _rate_of_change(plant.plant_w) / plant.plant_w, atol=1e-12)
    has_fruit = plant_type != 'LEAFY_PLANT' or rate_fn is transfer_rates_morris   # morris ignores plant_type
    if has_fruit:
        np.testing.assert_allclose(res['fh']['residual'], _rate_of_change(plant.friut_w) / plant.friut_w,
                                   atol=1e-12)
    # root uptake supplies net transpiration, the plant -> soil flow and the growth of plant and fruit water
    to_bh, to_fh = root_uptake(k)
    fruit = _rate_of_change(plant.friut_w) if has_fruit else 0
    expected = net_transpiration(air, plant) + k[:, K_INDEX['kbh_so']] * plant.plant_w + \
        _rate_of_change(plant.plant_w) + fruit
    np.testing.assert_allclose(to_bh + to_fh, expected)
    # each soil layer supplies its share
    for layer, share, w in zip('123', ROOT_LAYER_SHARE, SOIL_W):
        flux = (k[:, K_INDEX[f'ks{layer}_bh']] + k[:, K_INDEX[f'ks{layer}_fh']]) * w
        np.testing.assert_allclose(flux, share * expected)


def test_rain_is_not_intercepted_by_the_plant():
    air, soil, plant = stub_state()
    k_dry = transfer_rates_baomi('HTO', 'ROOT_VEG', air, soil, plant)
    air.rainfall = np.full(4, 2.0)                      # mm/h
    k_wet = transfer_rates_baomi('HTO', 'ROOT_VEG', air, soil, plant)
    np.testing.assert_array_equal(k_wet[:, K_INDEX['ka2_bh']], k_dry[:, K_INDEX['ka2_bh']])
    assert np.all(k_wet[:, K_INDEX['ka2_s1']] > k_dry[:, K_INDEX['ka2_s1']])


def growing_plant(days=3, dark=para_constant.OBT_DARK):
    '''Calc_plant with a diurnal radiation cycle; shoot grows on day 1, fruit on days 1-2, shoot senesces on day 3.'''
    import pandas as pd
    idx = pd.date_range('2021-07-01', periods=24 * days, freq='h')
    hour = idx.hour.values
    plant = Calc_plant.__new__(Calc_plant)
    plant.meteo = pd.DataFrame(index=idx)
    plant.PAR = np.where((hour >= 6) & (hour < 20), 500 * np.sin(np.pi * (hour - 6 + 0.5) / 14), 0.001)
    t = np.arange(len(idx)) / 24
    plant.twlv = 2000 + 500 * np.clip(t, 0, 1) - 300 * np.clip(t - 2, 0, 1)
    plant.twst = np.full(len(idx), 1000.0)
    plant.twso = 1000 * np.clip(t, 0, 2)
    plant.plant_w = np.full(len(idx), 1.0)
    return plant, idx


def test_obt_formation_follows_dry_matter_growth():
    plant, idx = growing_plant()
    plant.calc_transOBT()
    assert np.all(plant.TROBT >= 0) and np.all(plant.TROBT_FRUIT >= 0)
    day = idx.normalize()
    coef = para_constant.ID_OBT * para_constant.COMB_WATER
    # daily OBT formed in the fruit = ID * combustion water * daily fruit growth (per unit leaf HTO concentration)
    formed = pd_sum(plant.TROBT_FRUIT * plant.plant_w, day)
    one_hour = coef * 0.1 / 24                          # the hourly growth rate is smeared over the hour of change
    np.testing.assert_allclose(formed, coef * np.array([0.1, 0.1, 0.0]), rtol=0.05, atol=one_hour)
    # shoot: growth plus replacement of the turned-over OBT; senescence does not destroy OBT
    m_veg = (plant.twlv + plant.twst) / 1e4
    turnover = pd_sum(para_constant.KBO_BH * m_veg, day)
    np.testing.assert_allclose(pd_sum(plant.TROBT * plant.plant_w, day),
                               coef * (np.array([0.05, 0.0, 0.0]) + turnover), rtol=0.05, atol=one_hour)
    # night formation is OBT_DARK times the mean daylight rate
    is_day = plant.PAR > 2
    d1 = day == day[0]
    ratio = plant.TROBT_FRUIT[d1 & ~is_day].mean() / plant.TROBT_FRUIT[d1 & is_day].mean()
    assert ratio == pytest.approx(para_constant.OBT_DARK, rel=1e-6)


def pd_sum(values, day):
    import pandas as pd
    return pd.Series(values).groupby(day).sum().values


def test_chronic_shoot_obt_reaches_discrimination_factor():
    # leaf water held at a constant concentration: shoot OBT / combustion water -> ID_OBT * C_leaf
    plant, idx = growing_plant(days=60)
    plant.twlv = np.full(len(idx), 2500.0)
    plant.calc_transOBT()
    c_leaf, a_bo = 1.0, 0.0
    m_veg = (plant.twlv + plant.twst) / 1e4
    for i in range(len(idx)):
        a_bo += plant.TROBT[i] * plant.plant_w[i] * c_leaf - para_constant.KBO_BH * a_bo
    assert a_bo / (para_constant.COMB_WATER * m_veg[-1]) == pytest.approx(para_constant.ID_OBT, rel=0.05)


@pytest.mark.parametrize('rate_fn', SYMMETRIC, ids=lambda f: f.__name__)
def test_leaf_exchange_reaches_isotopic_equilibrium_with_air(rate_fn):
    rho_v = 12.0
    air, soil, plant = stub_state(rho_v=rho_v)
    k = rate_fn('HTO', 'ROOT_VEG', air, soil, plant)
    w_air = air.ML * rho_v / 1000                      # kg water per m2 in the mixing layer
    # steady state of air -> leaf -> air alone: C_leaf / C_air = ka2_bh * W_air / (kbh_a2 * W_leaf)
    ratio = k[:, K_INDEX['ka2_bh']] * w_air / (k[:, K_INDEX['kbh_a2']] * plant.plant_w)
    np.testing.assert_allclose(ratio, para_constant.ALF_HTO * rho_v / plant.ROW)


@pytest.mark.parametrize('rate_fn', SYMMETRIC, ids=lambda f: f.__name__)
def test_leaf_releases_hto_without_transpiration(rate_fn):
    air, soil, plant = stub_state()
    plant.ETRM = np.zeros(4)                           # night: no net transpiration
    k = rate_fn('HTO', 'ROOT_VEG', air, soil, plant)
    np.testing.assert_allclose(k[:, K_INDEX['kbh_a2']], leaf_exchange_rate(plant))
    assert np.all(k[:, K_INDEX['kbh_a2']] > 0.01)


def test_water_balance_residual_is_zero_for_balanced_fluxes():
    air, soil, plant = stub_state(n=3)
    k = np.zeros((3, len(K_NAMES)))
    k[:, K_INDEX['ks1_s2']] = 0.02
    k[:, K_INDEX['ks2_s1']] = 0.02 * 15 / 30           # same water flux back up
    res = water_balance_residual(k, air, soil, plant, alpha=1.0)
    np.testing.assert_allclose(res['s1']['residual'], 0, atol=1e-15)
    np.testing.assert_allclose(res['s2']['residual'], 0, atol=1e-15)
    np.testing.assert_allclose(res['s1']['k_out'], 0.02)

    k[:, K_INDEX['ks2_s1']] = 0                         # one-way drainage: s1 loses, s2 gains
    res = water_balance_residual(k, air, soil, plant, alpha=1.0)
    np.testing.assert_allclose(res['s1']['residual'], -0.02)
    np.testing.assert_allclose(res['s2']['residual'], 0.02 * 15 / 30)
