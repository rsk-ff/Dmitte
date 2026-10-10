"""Regression tests for the unit and indexing fixes in the transfer rates (review step 1)."""
from types import SimpleNamespace

import numpy as np
import pytest

from dmitte import para_constant
from dmitte.calc_para import Calc_plant, Calc_soil
from dmitte.transfer_equotions import (transfer_rates, transfer_rates_baomi, transfer_rates_Korea,
                                       transfer_rates_morris, transfer_rates_UFORTI)

K_NAMES = ['ka1_a1', 'ks0_a1', 'ka1_s0', 'ka1_a2', 'ks0_s1', 'ks0_s2', 'ks0_s3', 'ks1_a2', 'ka2_s1', 'ks2_s1',
           'ks1_s2', 'ks3_s2', 'ks2_s3', 'ks3_s3', 'ks1_bh', 'ks2_bh', 'ks3_bh', 'ka2_a2', 'kbh_a2', 'ka2_bh',
           'kbh_so', 'kbh_bo', 'kbo_bh', 'kfh_bh', 'kbh_fh', 'kbh_fo', 'ks1_fh', 'ks2_fh', 'ks3_fh']
K = {name: i for i, name in enumerate(K_NAMES)}

RATE_FUNCTIONS = [transfer_rates, transfer_rates_UFORTI, transfer_rates_Korea,
                  transfer_rates_baomi, transfer_rates_morris]
PLANT_TYPES = ['LEAFY_PLANT', 'CEREAL', 'ROOT_VEG']


def penman_monteith(temp, rh, q, ra, rs):
    """Reference Penman-Monteith evaporation, kg m-2 h-1 (FAO-56 units: kPa, J/(kg K), s/m)."""
    es = 0.6108 * np.exp(17.27 * temp / (temp + 237.3))
    vpd = es * (1 - rh)
    delta = 4098 * es / (temp + 237.3) ** 2
    gamma = 0.067
    rho, cp, lv = 1.2923, 1013.0, 2.45e6
    return (delta * q + rho * cp * vpd / ra) / (delta + gamma * (1 + rs / ra)) / lv * 3600


def stub_state(plant_type, sta, n=3):
    """Minimal air/soil/plant state; soil layers deliberately hold different amounts of water."""
    ones = np.ones(n)
    air = SimpleNamespace(ML=para_constant.mixLayerH[sta - 1], rainfall=0.1 * ones, atm_h=2.0 * ones,
                          atm_w=1e9 * ones, U10=2.0 * ones, sta=sta,
                          ICOMP=para_constant.Plant_dict[plant_type]['ICOMP'])
    soil = SimpleNamespace(VDSO=4e-3 * ones, ESOIL=0.1 * ones, Va_b1=0.01 * ones, Va_b2=0.01 * ones,
                           BODW1=0.3 * ones, soil1w=15.0 * ones, soil2w=30.0 * ones, soil3w=45.0 * ones)
    plant = SimpleNamespace(VDPF=1e-2 * ones, ETRM=0.2 * ones, plant_w=2.0 * ones, plant_wh=0.22 * ones,
                            friut_w=1.0 * ones, friut_wh=0.11 * ones, TROBT=1e-4 * ones)
    return air, soil, plant


def test_ingestion_dose_coefficients_match_icrp72():
    sv_to_msv = 1e3
    assert para_constant.DOSE_COEF_HTO == pytest.approx(1.8e-11 * sv_to_msv)
    assert para_constant.DOSE_COEF_OBT == pytest.approx(4.2e-11 * sv_to_msv)


def test_esoil_matches_penman_monteith():
    temp, rh, q, ra, rs = 20.0, 0.6, 400.0, 50.0, 100.0
    es = 0.6108 * np.exp(17.27 * temp / (temp + 237.3))
    soil = Calc_soil.__new__(Calc_soil)
    soil.RAM, soil.RB, soil.RSOIL = np.array([ra / 2]), np.array([ra / 2]), np.array([rs])
    soil.TEMP, soil.ea = np.array([temp]), np.array([rh * es])
    soil.PAR, soil.LAI = np.array([q]), np.array([0.0])
    soil.calc_ESOIL()
    assert soil.ESOIL[0] == pytest.approx(penman_monteith(temp, rh, q, ra, rs), rel=1e-6)
    assert soil.ESOIL[0] == pytest.approx(0.35, rel=0.02)


def test_etrm_keeps_vapour_pressure_deficit_term_at_night():
    temp, rh, ra, rc1 = 15.0, 0.7, 60.0, 1.0      # RC1 is scaled by 100 inside calc_ETRM
    es = 0.6108 * np.exp(17.27 * temp / (temp + 237.3))
    plant = Calc_plant.__new__(Calc_plant)
    plant.RAM, plant.RB, plant.RC1 = np.array([ra / 2]), np.array([ra / 2]), np.array([rc1])
    plant.TEMP, plant.ea = np.array([temp]), np.array([rh * es])
    plant.PAR, plant.LAI = np.array([0.001]), np.array([3.0])
    plant.calc_ETRM()
    expected = penman_monteith(temp, rh, plant.QSTR[0], ra, rc1 * 100)
    assert plant.ETRM[0] == pytest.approx(expected, rel=1e-6)
    assert plant.ETRM[0] > 0.01


@pytest.mark.parametrize('rate_fn', RATE_FUNCTIONS, ids=lambda f: f.__name__)
@pytest.mark.parametrize('plant_type', PLANT_TYPES)
def test_rates_are_finite_nonnegative_and_independent_of_stability(rate_fn, plant_type):
    rates = []
    for sta in range(1, 7):
        air, soil, plant = stub_state(plant_type, sta)
        air.ML = para_constant.mixLayerH[3]     # isolate the explicit use of sta
        k = rate_fn('HTO', plant_type, air, soil, plant)
        assert np.all(np.isfinite(k))
        rates.append(k)
    for k in rates[1:]:
        np.testing.assert_array_equal(k, rates[0])
    # Va_b1/Va_b2 (calc_soil_ab) can still be negative; everything else must not be
    k = rates[0]
    signed = [K['ks1_s2'], K['ks2_s3']]
    assert np.all(np.delete(k, signed, axis=1) >= 0)


@pytest.mark.parametrize('rate_fn', RATE_FUNCTIONS, ids=lambda f: f.__name__)
@pytest.mark.parametrize('plant_type', ['CEREAL', 'ROOT_VEG'])
def test_layer3_root_uptake_uses_layer3_water(rate_fn, plant_type):
    air, soil, plant = stub_state(plant_type, sta=4)
    k = rate_fn('HTO', plant_type, air, soil, plant)
    if rate_fn is transfer_rates_UFORTI and plant_type == 'ROOT_VEG':
        pytest.skip('UFOTRI ROOT_VEG uses constant rates')
    # flux out of layer 3 relative to flux out of layer 2, both per unit of plant uptake
    flux3 = (k[:, K['ks3_bh']] + k[:, K['ks3_fh']]) * soil.soil3w
    flux2 = (k[:, K['ks2_bh']] + k[:, K['ks2_fh']]) * soil.soil2w
    if np.all(flux2 > 0) and np.all(flux3 > 0):
        np.testing.assert_allclose(flux3, flux2)       # both layers supply 40 % of uptake


@pytest.mark.parametrize('rate_fn', [transfer_rates, transfer_rates_Korea, transfer_rates_baomi,
                                     transfer_rates_morris], ids=lambda f: f.__name__)
def test_soil_evaporation_rate_divides_by_layer1_water(rate_fn):
    air, soil, plant = stub_state('ROOT_VEG', sta=4)
    k = rate_fn('HTO', 'ROOT_VEG', air, soil, plant)
    np.testing.assert_allclose(k[:, K['ks1_a2']], soil.ESOIL / soil.soil1w * 1.1)


def test_cereal_fruit_formation_rate_uses_crop_half_life():
    air, soil, plant = stub_state('CEREAL', sta=6)
    for rate_fn in (transfer_rates, transfer_rates_UFORTI, transfer_rates_Korea):
        k = rate_fn('HTO', 'CEREAL', air, soil, plant)
        hwz = para_constant.HWZ[para_constant.Plant_dict['CEREAL']['ICOMP'] - 1]
        expected = np.log(2) / (hwz / 2) * (plant.friut_wh / plant.plant_wh) / 24
        np.testing.assert_allclose(k[:, K['kbh_fo']], expected)
        np.testing.assert_allclose(k[:, K['kbh_fh']], np.log(2) / 2)
