"""Recombination checks (``03_numerics_spec.md`` §5.1-5.2)."""

import numpy as np
import pytest

from cmbcore import Params, BackgroundCosmology, Recombination
from cmbcore import constants as const


@pytest.fixture(scope="module")
def rec():
    p = Params.fiducial()
    bg = BackgroundCosmology(p)
    return Recombination(bg, p), bg, p


def test_visibility_normalized(rec):
    r, _, _ = rec
    assert abs(r.visibility_norm() - 1.0) < 1e-3


def test_z_star(rec):
    r, _, _ = rec
    assert abs(r.z_star() - 1090) < 15


def test_sound_horizon(rec):
    r, _, _ = rec
    rs = r.r_s(r.x_star()) / const.Mpc
    assert abs(rs - 144.5) < 4.0


def test_comoving_distance_to_star(rec):
    r, bg, _ = rec
    chi = bg.comoving_distance(r.x_star()) / const.Gpc
    assert abs(chi - 13.9) < 0.3


def test_theta_star(rec):
    r, bg, _ = rec
    rs = r.r_s(r.x_star())
    chi = bg.comoving_distance(r.x_star())
    theta100 = 100.0 * rs / chi
    assert abs(theta100 - 1.041) < 0.02


def test_Xe_monotone_after_recomb(rec):
    r, _, _ = rec
    x = np.linspace(-7.0, -4.0, 50)
    Xe = r.Xe(x)
    assert np.all(np.diff(Xe) <= 1e-6)  # non-increasing through recombination


def test_dtau_analytic_matches_spline(rec):
    r, _, _ = rec
    x = -7.0
    dx = 1e-4
    num = (r.tau(x + dx) - r.tau(x - dx)) / (2 * dx)
    assert abs(r.dtau(x) - num) / abs(num) < 1e-2
