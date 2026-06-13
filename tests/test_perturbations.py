"""Perturbation / Boltzmann-hierarchy checks (``03_numerics_spec.md`` §5.1)."""

import numpy as np
import pytest

from cmbcore import (Params, BackgroundCosmology, Recombination,
                     PerturbationSolver)
from cmbcore import constants as const


@pytest.fixture(scope="module")
def setup():
    p = Params.fiducial()
    bg = BackgroundCosmology(p)
    rec = Recombination(bg, p)
    solver = PerturbationSolver(bg, rec, p)
    return solver, bg, rec, p


def test_initial_adiabatic_relations(setup):
    """f_nu->0 gives Phi=-Psi initially (no anisotropic stress, §3.2)."""
    p = Params(N_eff=1e-6)  # nearly massless-free -> f_nu ~ 0
    bg = BackgroundCosmology(p)
    rec = Recombination(bg, p)
    solver = PerturbationSolver(bg, rec, p)
    k = 0.05 / const.Mpc
    y0 = solver._initial_conditions(solver.x_start, k)
    i = solver._idx()
    Phi = y0[i["Phi"]]
    # Analytic super-horizon Psi (the IC is built from this); with f_nu -> 0 the
    # adiabatic relation reduces to Phi = -Psi exactly.
    Psi = -1.0 / (1.5 + 0.4 * p.f_nu)
    assert abs(Phi + Psi) / abs(Phi) < 1e-3
    # And the photon/neutrino monopoles satisfy Theta0 = N0 = -Psi/2.
    assert abs(y0[i["Theta"]] - (-0.5 * Psi)) < 1e-6


def test_matter_growth(setup):
    """delta_c grows by orders of magnitude from x_start to today."""
    solver, _, _, _ = setup
    k = 0.05 / const.Mpc
    res = solver.solve(k)
    assert abs(res["delta_c"](0.0)) > 100 * abs(res["delta_c"](-15))


def test_tc_continuity(setup):
    """Theta0, Theta1 are continuous across the TC handoff."""
    solver, _, _, _ = setup
    k = 0.1 / const.Mpc
    res = solver.solve(k)
    x_tc = solver._x_tc_end(k)
    eps = 1e-3
    # The handoff copies the TC state into the full system, so the variables are
    # continuous; allow either a small absolute or a small relative mismatch
    # (Theta1 crosses zero near the handoff, so a pure relative test is unfair).
    scale = max(abs(res["Phi"](x_tc)), abs(res["Theta0"](x_tc)))
    # No gross discontinuity (the handoff copies state exactly; the few-percent
    # tolerance absorbs the TC approximation's higher-derivative kink across the
    # ~eps window straddling the switch).
    for nm in ["Theta0", "Theta1", "Phi"]:
        left = res[nm](x_tc - eps)
        right = res[nm](x_tc + eps)
        assert abs(left - right) < 5e-2 * scale + 1e-8


def test_superhorizon_psi_constant(setup):
    """On super-horizon, Psi is roughly constant (drift small) early on."""
    solver, bg, _, _ = setup
    k = 1e-4 / const.Mpc  # very large scale
    res = solver.solve(k)
    psi_early = res["Psi"](-15)
    psi_later = res["Psi"](-10)
    assert abs(psi_later - psi_early) / abs(psi_early) < 0.1
