"""Background-cosmology checks (``03_numerics_spec.md`` §5.1-5.2)."""

import numpy as np
import pytest

from cmbcore import Params, BackgroundCosmology
from cmbcore import constants as const


@pytest.fixture(scope="module")
def bg():
    return BackgroundCosmology(Params.fiducial())


def test_flat_closure():
    p = Params.fiducial()
    total = (p.Omega_b + p.Omega_c + p.Omega_gamma + p.Omega_nu
             + p.Omega_Lambda + p.Omega_k)
    assert abs(total - 1.0) < 1e-12


def test_z_eq(bg):
    assert abs(bg.z_eq() - 3400) < 100


def test_age_of_universe(bg):
    age_gyr = bg.t0 / const.Gyr
    assert 13.0 < age_gyr < 14.5


def test_eta_monotonic(bg):
    x = np.linspace(-15, 0, 100)
    eta = bg.eta(x)
    assert np.all(np.diff(eta) > 0)


def test_Hp_derivative_matches_numeric(bg):
    """Analytic Hp' agrees with a finite difference (sign/coefficient guard)."""
    x = -5.0
    dx = 1e-5
    num = (bg.Hp(x + dx) - bg.Hp(x - dx)) / (2 * dx)
    assert abs(bg.dHp(x) - num) / abs(num) < 1e-4


def test_dimensionless_invariance(bg):
    """The reduced expansion rate Hp(x)/H0 is the dimensionless sqrt(E(x)).

    This is the §5.1 dimensional-consistency check: Hp carries the only
    dimensionful scale H0, so Hp/H0 must equal the dimensionless density sum.
    """
    p = bg.p
    x = np.array([-12.0, -8.0, -4.0, 0.0])
    E = ((p.Omega_b + p.Omega_c) * np.exp(-x)
         + (p.Omega_gamma + p.Omega_nu) * np.exp(-2 * x)
         + p.Omega_k + p.Omega_Lambda * np.exp(2 * x))
    assert np.allclose(bg.Hp(x) / p.H0, np.sqrt(E), rtol=1e-12)
