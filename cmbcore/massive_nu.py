"""Massive-neutrino background (Sigma m_nu; ch.13.8, 04 §7).

A massive neutrino species has comoving energy
:math:`\\varepsilon=\\sqrt{q^2+a^2m^2}`, so its energy density transitions from
radiation-like (early) to matter-like (late). This module provides the
phase-space integrals and the today-normalized density
:math:`f_{\\rm ncdm}(a)=\\rho_{\\rm ncdm}(a)/\\rho_{\\rm crit,0}` used by
:mod:`cmbcore.background`. The perturbation-level free-streaming (the q-grid
hierarchy) is the documented next step; the dominant CMB-TT effect is this
background geometry change.
"""

from __future__ import annotations

import numpy as np
from scipy import integrate

from . import constants as const

# Fermi-Dirac phase-space grid in x = q c / (k_B T_nu) (comoving momentum in
# temperature units). f0(x) = 1/(e^x + 1).
_X = np.linspace(1e-3, 30.0, 600)
_F0 = 1.0 / (np.exp(_X) + 1.0)
_I_RHO_0 = float(np.trapezoid(_X ** 3 * _F0, _X))  # = 7 pi^4 / 120


def _y(a: float, m_eV: float, T_ncdm_eV: float):
    """Dimensionless mass parameter y = a m c^2 / (k_B T_nu)."""
    return a * m_eV / T_ncdm_eV


def _I_rho(y):
    """:math:`\\int x^2\\sqrt{x^2+y^2} f_0\\,dx` (vectorized over y)."""
    y = np.atleast_1d(np.asarray(y, dtype=float))
    integ = (_X[:, None] ** 2 * np.sqrt(_X[:, None] ** 2 + y[None, :] ** 2)
             * _F0[:, None])
    return np.trapezoid(integ, _X, axis=0)


def T_ncdm_eV(T_cmb_K: float) -> float:
    """Neutrino temperature today in eV: (4/11)^{1/3} k_B T_CMB."""
    kT_cmb = const.k_B * T_cmb_K / const.eV
    return (4.0 / 11.0) ** (1.0 / 3.0) * kT_cmb


class MassiveNu:
    """Today-normalized massive-neutrino density f_ncdm(a) = rho/rho_crit0.

    ``Omega_gamma`` is today's photon density parameter; one massive species in
    the relativistic limit matches one massless species,
    :math:`(7/8)(4/11)^{4/3}\\Omega_\\gamma`.
    """

    def __init__(self, m_eV: float, N_ncdm: int, T_cmb_K: float,
                 Omega_gamma: float):
        self.m_eV = m_eV
        self.N_ncdm = N_ncdm
        self.T_nu = T_ncdm_eV(T_cmb_K)
        # Normalization so that y->0 reproduces N_ncdm massless species.
        omega_1nu_rel = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * Omega_gamma
        self._A = N_ncdm * omega_1nu_rel / _I_RHO_0

    def f_ncdm(self, a):
        """:math:`\\rho_{\\rm ncdm}(a)/\\rho_{\\rm crit,0}` (today-normalized)."""
        a = np.asarray(a, dtype=float)
        y = self.m_eV * a / self.T_nu
        return self._A * a ** -4 * _I_rho(y).reshape(a.shape)

    def Omega_ncdm0(self) -> float:
        """Today's massive-neutrino density parameter."""
        return float(self.f_ncdm(np.array([1.0]))[0])
