"""Recombination history (``03_numerics_spec.md`` §2).

Computes the free-electron fraction :math:`X_e(x)` via Saha then Peebles, the
optical depth :math:`\\tau(x)`, the visibility function :math:`\\tilde g(x)` and
the sound horizon :math:`r_s(x)`.
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from . import constants as const
from .params import Params
from .background import BackgroundCosmology


class Recombination:
    """Solve and spline recombination and reionization history."""

    def __init__(self, bg: BackgroundCosmology, params: Params,
                 x_start: float = -12.0, x_end: float = 0.0,
                 n: int = 4000) -> None:
        self.bg = bg
        self.p = params
        self.x = np.linspace(x_start, x_end, n)
        self._solve_Xe()
        self._solve_tau()
        self._solve_sound_horizon()

    # --- helpers ---------------------------------------------------------------
    def _n_H(self, x):
        """Hydrogen number density [1/m^3], :math:`Y_p=0` default."""
        p = self.p
        a = np.exp(x)
        rho_b = p.Omega_b * 3.0 * p.H0 ** 2 / (8.0 * np.pi * const.G) * a ** -3
        return (1.0 - p.Yp) * rho_b / const.m_H

    def _T_b(self, x):
        """Baryon temperature, :math:`T_b=T_\\gamma=T_{CMB}/a` [K]."""
        Tb = self.p.T_CMB * np.exp(-x)
        return Tb * (1.0 + self.p.recomb_shift)

    # --- Saha ------------------------------------------------------------------
    def _saha_Xe(self, x):
        """Saha equilibrium :math:`X_e` (stable quadratic root)."""
        Tb = self._T_b(x)
        nH = self._n_H(x)
        # RHS = (1/nb)(m_e k_B T / 2 pi hbar^2)^{3/2} exp(-eps0/kT).
        pref = (const.m_e * const.k_B * Tb /
                (2.0 * np.pi * const.hbar ** 2)) ** 1.5
        y = pref * np.exp(-const.eps_0 / (const.k_B * Tb)) / nH
        # X_e^2/(1-X_e) = y  ->  X_e = 2/(1+sqrt(1+4/y)) (overflow safe).
        with np.errstate(over="ignore", divide="ignore"):
            Xe = 2.0 / (1.0 + np.sqrt(1.0 + 4.0 / y))
        return np.where(np.isfinite(Xe), Xe, 1.0)

    # --- Peebles ---------------------------------------------------------------
    def _peebles_rhs(self, x, Xe):
        """RHS of the Peebles ODE :math:`dX_e/dx`."""
        Xe = float(Xe[0])
        Xe = min(max(Xe, 1e-30), 1.0)
        p = self.p
        Tb = self._T_b(x)
        nH = self._n_H(x)
        H = float(self.bg.H(x))
        kT = const.k_B * Tb

        phi2 = 0.448 * np.log(const.eps_0 / kT)
        alpha2 = (8.0 / np.sqrt(3.0 * np.pi)) * const.sigma_T * const.c \
            * np.sqrt(const.eps_0 / kT) * phi2
        beta = alpha2 * (const.m_e * kT /
                         (2.0 * np.pi * const.hbar ** 2)) ** 1.5 \
            * np.exp(-const.eps_0 / kT)
        # beta^(2) = beta * exp(3 eps0 / 4 kT), built as analytic composite to
        # avoid the e^{-eps0/kT} underflow cancelling e^{3eps0/4kT}.
        beta2 = alpha2 * (const.m_e * kT /
                          (2.0 * np.pi * const.hbar ** 2)) ** 1.5 \
            * np.exp(-const.eps_0 / (4.0 * kT))

        n1s = (1.0 - Xe) * nH
        Lambda_alpha = H * (3.0 * const.eps_0) ** 3 \
            / ((8.0 * np.pi) ** 2 * (const.hbar * const.c) ** 3 * n1s)
        Cr = (const.Lambda_2s1s + Lambda_alpha) \
            / (const.Lambda_2s1s + Lambda_alpha + beta2)

        dXedx = (Cr / H) * (beta * (1.0 - Xe) - nH * alpha2 * Xe ** 2)
        return [dXedx]

    def _solve_Xe(self) -> None:
        x = self.x
        Xe = np.empty_like(x)
        # Phase 1: Saha until Xe drops below 0.99.
        saha = self._saha_Xe(x)
        switch = np.argmax(saha < 0.99)
        if saha[switch] >= 0.99:  # never drops (shouldn't happen) -> all Saha
            switch = len(x) - 1
        Xe[:switch] = saha[:switch]

        # Phase 2: Peebles ODE from the switch point onward.
        x_p = x[switch:]
        sol = solve_ivp(
            self._peebles_rhs, (x_p[0], x_p[-1]), [saha[switch]],
            t_eval=x_p, method="LSODA", rtol=1e-8, atol=1e-12,
        )
        Xe[switch:] = sol.y[0]

        # Reionization (tanh model) added on top, if enabled (tau_reio != 0
        # is used as the on/off flag; chapter 13 / NB6).
        self._Xe_rec = np.clip(Xe, 1e-30, 1.0 + 1e-12)
        if self.p.tau_reio != 0.0:
            self._Xe_rec = self._Xe_rec + self._reionization_Xe(x)
        self.Xe_spline = CubicSpline(x, np.log(np.clip(self._Xe_rec, 1e-30, None)))
        self.x_switch = x[switch]

    def _reionization_Xe(self, x):
        """tanh reionization contribution (``03`` §2.4)."""
        p = self.p
        z = np.exp(-x) - 1.0
        f_He = 0.0  # Yp=0 default
        dy = ((1.0 + p.z_reio) ** 1.5) - ((1.0 + p.z_reio - p.delta_z_reio) ** 1.5)
        dy = max(dy, 1e-6)
        arg = (((1.0 + p.z_reio) ** 1.5) - ((1.0 + z) ** 1.5)) / dy
        return (1.0 + f_He) / 2.0 * (1.0 + np.tanh(arg))

    def Xe(self, x):
        """Free-electron fraction :math:`X_e(x)`."""
        return np.exp(self.Xe_spline(np.asarray(x, dtype=float)))

    def n_e(self, x):
        """Free-electron number density [1/m^3]."""
        return self.Xe(x) * self._n_H(x)

    # --- optical depth & visibility -------------------------------------------
    def _solve_tau(self) -> None:
        x = self.x
        # dtau/dx = -c n_e sigma_T / H ; integrate backward from x=0 (tau=0).
        def rhs(xx, y):
            ne = float(self.n_e(xx))
            H = float(self.bg.H(xx))
            return [-const.c * ne * const.sigma_T / H]

        sol = solve_ivp(
            rhs, (x[-1], x[0]), [0.0], t_eval=x[::-1],
            method="LSODA", rtol=1e-8, atol=1e-12,
        )
        tau = sol.y[0][::-1]
        self.tau_spline = CubicSpline(x, tau)
        self._tau_arr = tau
        # Visibility g~ = -tau' e^{-tau}.
        taup = self.tau_spline(x, 1)
        g = -taup * np.exp(-tau)
        self.g_spline = CubicSpline(x, g)
        self._g_arr = g

    def tau(self, x):
        """Optical depth :math:`\\tau(x)`."""
        return self.tau_spline(np.asarray(x, dtype=float))

    def dtau(self, x):
        """:math:`\\tau'=d\\tau/dx = -c\\,n_e\\sigma_T/H` (analytic)."""
        x = np.asarray(x, dtype=float)
        return -const.c * self.n_e(x) * const.sigma_T / self.bg.H(x)

    def ddtau(self, x):
        """:math:`\\tau''` from the spline of :math:`\\tau`."""
        return self.tau_spline(np.asarray(x, dtype=float), 2)

    def g_tilde(self, x):
        """Visibility function :math:`\\tilde g(x)=-\\tau'e^{-\\tau}`."""
        return self.g_spline(np.asarray(x, dtype=float))

    def dg_tilde(self, x):
        return self.g_spline(np.asarray(x, dtype=float), 1)

    def ddg_tilde(self, x):
        return self.g_spline(np.asarray(x, dtype=float), 2)

    def visibility_norm(self) -> float:
        """:math:`\\int\\tilde g\\,dx`; should be 1 (``03`` §2.3 check)."""
        return float(np.trapezoid(self._g_arr, self.x))

    # --- decoupling & sound horizon -------------------------------------------
    def x_star(self) -> float:
        """:math:`x_*` where :math:`\\tau(x_*)=1`."""
        f = lambda xx: float(self.tau(xx)) - 1.0
        # tau decreases from large (early) to 0 (today); bracket it.
        return brentq(f, self.x[0], -0.1)

    def z_star(self) -> float:
        return float(np.exp(-self.x_star()) - 1.0)

    def _Rs(self, x):
        """Baryon-photon ratio :math:`R_s=3\\Omega_{b0}/(4\\Omega_{\\gamma0})\\,a`."""
        return 3.0 * self.p.Omega_b / (4.0 * self.p.Omega_gamma) * np.exp(x)

    def _solve_sound_horizon(self) -> None:
        # Comoving sound horizon r_s = int c_s d(eta_conf).  With the length
        # convention eta carries a factor c (deta_len/dx = c/Hp), so
        # d(eta_conf) = deta_len/c and the integrand becomes
        #   (c_s/c) * (c/Hp) = sqrt(1/(3(1+R_s))) * c/Hp   [metres].
        x = self.x

        def rhs(xx, y):
            Rs = self._Rs(xx)
            cs_over_c = np.sqrt(1.0 / (3.0 * (1.0 + Rs)))
            return [cs_over_c * const.c / float(self.bg.Hp(xx))]

        sol = solve_ivp(
            rhs, (x[0], x[-1]), [0.0], t_eval=x,
            method="LSODA", rtol=1e-8, atol=1e-12,
        )
        self.rs_spline = CubicSpline(x, sol.y[0])

    def r_s(self, x):
        """Comoving sound horizon :math:`r_s(x)` [m]."""
        return self.rs_spline(np.asarray(x, dtype=float))
