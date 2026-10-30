"""Physics and maths checks: apparent resistivity, phase, phase tensors, rotation, 1-D modelling.

These test against known results (a uniform half-space, an independent 1-D
recursion, rotation invariants) rather than against stored outputs.
"""
from copy import deepcopy

import numpy as np
import pytest

import gateway_mt.utils as utils
from gateway_mt.modelling_1D import wait_1D

from helpers import make_site, one_d_impedance, one_d_tensor, rotate_axes

PERIODS = np.logspace(-3, 3, 13)


# --------------------------------------------------------------------------- #
# Apparent resistivity and phase
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize('comp', ['xy', 'yx', 'det', 'ssq', 'gav', 'aav'])
def test_half_space_gives_its_resistivity_and_45_degrees(comp):
    """Over a uniform 100 ohm-m half-space every definition of rho_a is 100 and every phase 45 degrees."""
    site = make_site(PERIODS, one_d_tensor(one_d_impedance(PERIODS, [100.0], [])))
    rho, _, _ = utils.compute_rho(site, calc_comp=comp)
    phase, _ = utils.compute_phase(site, calc_comp=comp, wrap=1)
    np.testing.assert_allclose(np.real(rho), 100.0, rtol=1e-4)   # periods are stored to 5 significant figures
    np.testing.assert_allclose(phase, 45.0, atol=0.01)


def test_layered_earth_apparent_resistivity_limits():
    """A resistive layer over a conductor: rho_a tends to the top layer at short periods and falls at long periods."""
    Z = one_d_impedance(PERIODS, [1000.0, 10.0], [2000.0])
    site = make_site(PERIODS, one_d_tensor(Z))
    rho, _, _ = utils.compute_rho(site, calc_comp='xy')
    phase, _ = utils.compute_phase(site, calc_comp='xy')
    assert np.real(rho[0]) == pytest.approx(1000.0, rel=0.01)
    assert np.real(rho[-1]) < 20.0
    assert phase.max() > 60.0          # phase > 45 degrees over the conductor


def test_rho_errors_are_positive_and_flagged_points_stay_flagged():
    site = make_site(PERIODS, one_d_tensor(one_d_impedance(PERIODS, [100.0], [])))
    site.used_error['ZXYR'][2] = site.REMOVE_FLAG
    _, rho_err, log_err = utils.compute_rho(site, calc_comp='xy', errtype='used_error')
    assert rho_err[2] == site.REMOVE_FLAG
    assert np.all(np.delete(rho_err, 2) > 0) and np.all(np.delete(log_err, 2) > 0)


def test_skin_depth():
    assert utils.skin_depth(100, 1) == pytest.approx(5000.0)
    assert utils.skin_depth(1, 100) == pytest.approx(5000.0)


# --------------------------------------------------------------------------- #
# 1-D forward modelling
# --------------------------------------------------------------------------- #
def test_wait_1d_matches_independent_recursion():
    """wait_1D (layer thicknesses as ``dz``) agrees with a textbook implementation for a 5-layer earth.

    This is also how the synthetic test data is generated.
    """
    rho = [300.0, 3000.0, 10.0, 2000.0, 100.0]
    thicknesses = [1000.0, 12000.0, 5000.0, 100000.0]
    np.testing.assert_allclose(wait_1D(PERIODS, rho, np.array(thicknesses)),
                               one_d_impedance(PERIODS, rho, thicknesses), rtol=1e-10)


def test_wait_1d_half_space():
    np.testing.assert_allclose(wait_1D(PERIODS, [100.0], np.array([])),
                               one_d_impedance(PERIODS, [100.0], []), rtol=1e-10)


def test_wait_1d_two_layers():
    """dz has one entry per layer, as the 1-D modelling window passes it; the half-space's entry is unused."""
    np.testing.assert_allclose(wait_1D(PERIODS, [1000.0, 10.0], np.array([2000.0, 1e6])),
                               one_d_impedance(PERIODS, [1000.0, 10.0], [2000.0]), rtol=1e-10)


# --------------------------------------------------------------------------- #
# Phase tensors
# --------------------------------------------------------------------------- #
def test_phase_tensor_of_a_1d_earth_is_a_circle():
    """Over a 1-D earth the phase tensor is isotropic: phi_max = phi_min, no skew (beta = 0)."""
    Z = one_d_impedance(PERIODS, [300.0, 10.0, 1000.0], [1000.0, 5000.0])
    site = make_site(PERIODS, one_d_tensor(Z))
    for pt in site.phase_tensors:
        assert pt.phi_max == pytest.approx(pt.phi_min, rel=1e-9)
        assert pt.beta == pytest.approx(0.0, abs=1e-9)


def test_phase_tensor_of_a_half_space_is_45_degrees():
    site = make_site(PERIODS, one_d_tensor(one_d_impedance(PERIODS, [100.0], [])))
    for pt in site.phase_tensors:
        # phi_max / phi_min are stored as tangents of the phase
        assert np.degrees(np.arctan(pt.phi_max)) == pytest.approx(45.0, abs=1e-3)


# --------------------------------------------------------------------------- #
# Rotation
# --------------------------------------------------------------------------- #
def _two_d_site(strike_deg):
    """An ideal 2-D response (strike along x'), observed in axes rotated by ``strike_deg``."""
    Zte = one_d_impedance(PERIODS, [300.0, 10.0, 1000.0], [1000.0, 5000.0])
    Ztm = one_d_impedance(PERIODS, [300.0, 80.0, 1000.0], [1000.0, 5000.0])
    Z = np.zeros((len(PERIODS), 2, 2), dtype=complex)
    Z[:, 0, 1], Z[:, 1, 0] = Zte, -Ztm
    T = np.zeros((len(PERIODS), 2), dtype=complex)
    T[:, 1] = 0.2 + 0.1j                       # 2-D: Hz depends only on the across-strike field
    Zr, Tr = rotate_axes(Z, T, -strike_deg)
    return make_site(PERIODS, Zr, Tr)


def test_rotation_is_reversible():
    site = _two_d_site(25)
    original = deepcopy(site.data)
    site.rotate_data(azi=40)
    assert site.azimuth == 40
    site.rotate_data(azi=0)
    for comp in original:
        np.testing.assert_allclose(site.data[comp], original[comp], atol=1e-12)


@pytest.mark.parametrize('strike', [20, 45, 70])
def test_rotating_onto_strike_removes_diagonals_and_along_strike_tipper(strike):
    """Impedance and tipper rotate consistently: on strike, Zxx = Zyy = 0 and Tzx = 0 together."""
    site = _two_d_site(strike)
    assert np.abs(site.data['ZXXR']).max() > 1e-3          # not aligned before rotating
    site.rotate_data(azi=-strike)
    for comp in ('ZXXR', 'ZXXI', 'ZYYR', 'ZYYI', 'TZXR', 'TZXI'):
        np.testing.assert_allclose(site.data[comp], 0, atol=1e-12)
    np.testing.assert_allclose(site.data['TZYR'], 0.2, rtol=1e-9)


def test_rotation_invariants():
    """The impedance determinant and the phase-tensor invariants do not depend on the axes."""
    a, b = _two_d_site(30), _two_d_site(30)
    b.rotate_data(azi=55)
    det_a, _ = utils.compute_MT_determinant(a)
    det_b, _ = utils.compute_MT_determinant(b)
    np.testing.assert_allclose(det_a, det_b, rtol=1e-10)
    for pa, pb in zip(a.phase_tensors, b.phase_tensors):
        assert pa.phi_max == pytest.approx(pb.phi_max, rel=1e-9)
        assert pa.phi_min == pytest.approx(pb.phi_min, rel=1e-9)
        assert pa.beta == pytest.approx(pb.beta, abs=1e-9)


# --------------------------------------------------------------------------- #
# Coordinates
# --------------------------------------------------------------------------- #
def _as_dms_text(dd):
    d, m, s = utils.dd2dms(dd)
    return '{}{}:{}:{}'.format('-' if dd < 0 else '', abs(int(d)), int(m), s)


@pytest.mark.parametrize('dd', [49.123456, -84.5123, 1.5, -1.25])
def test_degrees_minutes_seconds_round_trip(dd):
    assert utils.dms2dd(_as_dms_text(dd)) == pytest.approx(dd, abs=1e-9)


@pytest.mark.parametrize('dd', [0.5, -0.25])
def test_degrees_minutes_seconds_round_trip_below_one_degree(dd):
    """Within a degree of the equator / Greenwich the degrees field is (-)0, so the sign comes from the text."""
    assert utils.dms2dd(_as_dms_text(dd)) == pytest.approx(dd, abs=1e-9)


def test_utm_projection_matches_pyproj():
    pyproj = pytest.importorskip('pyproj')
    zone, letter, easting, northing = utils.project((-84.5, 49.0))
    assert (zone, letter) == (16, 'U')
    expected = pyproj.Transformer.from_crs('epsg:4326', 'epsg:32616', always_xy=True).transform(-84.5, 49.0)
    assert (easting, northing) == pytest.approx(expected, abs=0.01)
    # geo2utm is the faster built-in transform used when reading EDIs; it returns (easting, northing, origin)
    easting, northing, _ = utils.geo2utm(49.0, -84.5)
    assert (easting, northing) == pytest.approx(expected, abs=0.1)
