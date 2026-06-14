"""Spectrum-level checks (``03_numerics_spec.md`` §5.2).

These are *fast* sanity checks on a coarse grid; the full CLASS comparison
(§5.3) lives in the validation pipeline, not the unit-test suite.
"""

import numpy as np
import pytest

from cmbcore import Params, PowerSpectrum
from cmbcore import constants as const


@pytest.fixture(scope="module")
def coarse_spectrum():
    from cmbcore.spectrum import k_grid
    p = Params.fiducial()
    ps = PowerSpectrum(p)
    ells = np.array([2, 20, 100, 180, 200, 220, 240, 300, 400, 500])
    ks = k_grid(nk=120) / const.Mpc  # quadratic grid (low-k dense)
    ells, Dl = ps.dls(ells=ells, ks_si=ks)
    return ells, Dl


def test_positive_spectrum(coarse_spectrum):
    _, Dl = coarse_spectrum
    assert np.all(Dl > 0)


def test_first_peak_near_220(coarse_spectrum):
    ells, Dl = coarse_spectrum
    peak_ell = ells[np.argmax(Dl)]
    assert 180 <= peak_ell <= 260


def test_acoustic_scale_sane():
    from cmbcore.background import BackgroundCosmology
    from cmbcore.recombination import Recombination
    from cmbcore import analytic
    p = Params.fiducial()
    bg = BackgroundCosmology(p)
    rec = Recombination(bg, p)
    lA = analytic.acoustic_scale(bg, rec)
    # l_A ~ pi chi / r_s ~ 300; first peak ~ 0.75 l_A ~ 220.
    assert 280 < lA < 320


def test_polarization_TTEETE():
    """EE is positive and out of phase with TT; TE oscillates in sign."""
    from cmbcore.spectrum import k_grid
    ps = PowerSpectrum(Params.fiducial())
    ells = np.array([2, 100, 200, 300, 400, 500, 700])
    ks = k_grid(nk=100) / const.Mpc
    ells, D = ps.cls_all(ells=ells, ks_si=ks)
    assert np.all(D["EE"] >= 0)                 # EE auto-spectrum non-negative
    assert D["EE"][ells == 400] > D["EE"][ells == 200]  # EE rises (TT trough ~ EE peak)
    assert np.any(D["TE"] < 0) and np.any(D["TE"] > 0)  # TE changes sign
    assert np.all(D["TT"] > 0)
