"""Analytic estimators (``03_numerics_spec.md`` §6, NB7).

Closed-form approximations used for cross-checks and pedagogy: the acoustic
scale and peak positions, the Silk damping scale, and the sound horizon.
"""

from __future__ import annotations

import numpy as np

from . import constants as const
from .background import BackgroundCosmology
from .recombination import Recombination


def acoustic_scale(bg: BackgroundCosmology, rec: Recombination) -> float:
    """Acoustic angular scale :math:`\\ell_A=\\pi\\,\\chi(x_*)/r_s(x_*)`."""
    x_star = rec.x_star()
    rs = float(rec.r_s(x_star))
    chi = float(bg.comoving_distance(x_star))
    return np.pi * chi / rs


def peak_positions(bg: BackgroundCosmology, rec: Recombination,
                   n: int = 5) -> np.ndarray:
    """Approximate acoustic peak multipoles :math:`\\ell_m\\approx \\ell_A(m-\\phi)`.

    Uses a phase shift :math:`\\phi\\approx0.27` typical for the first peak.
    """
    lA = acoustic_scale(bg, rec)
    phi = 0.27
    m = np.arange(1, n + 1)
    return lA * (m - phi)


def silk_scale(bg: BackgroundCosmology, rec: Recombination) -> float:
    """Silk damping wavenumber :math:`k_D` [1/m] (diffusion integral).

    :math:`k_D^{-2}=\\int \\frac{dx}{6(1+R_s)\\,n_e\\sigma_T a}
    \\left(\\frac{R_s^2}{1+R_s}+\\frac{16}{15}\\right)` evaluated to :math:`x_*`.
    """
    x_star = rec.x_star()
    xs = np.linspace(bg.x_start + 1.0, x_star, 2000)
    a = np.exp(xs)
    Rs = 3.0 * rec.p.Omega_b / (4.0 * rec.p.Omega_gamma) * a
    Hp = bg.Hp(xs)
    dtau = np.abs(rec.dtau(xs))  # |tau'| = c n_e sigma_T / H
    # 1/k_D^2 = int deta /(6(1+R)) [R^2/(1+R)+16/15] / (n_e sigma_T a),
    # with deta = c dx / Hp (length convention) and
    # 1/(n_e sigma_T a) = c / (|tau'| Hp).  Hence the c^2/Hp^2 factor.
    integrand = (const.c ** 2 / (Hp ** 2 * dtau)) \
        / (6.0 * (1.0 + Rs)) * (Rs ** 2 / (1.0 + Rs) + 16.0 / 15.0)
    kD2_inv = np.trapezoid(integrand, xs)
    return 1.0 / np.sqrt(kD2_inv)


def theta_star(bg: BackgroundCosmology, rec: Recombination) -> float:
    """Sound-horizon angular scale :math:`100\\,\\theta_*=100\\,r_s/\\chi`."""
    x_star = rec.x_star()
    rs = float(rec.r_s(x_star))
    chi = float(bg.comoving_distance(x_star))
    return 100.0 * rs / chi
