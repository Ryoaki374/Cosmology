"""Background FLRW cosmology (``03_numerics_spec.md`` §1).

Independent variable is :math:`x=\\ln a` with :math:`a_0=1`. Primes denote
``d/dx``. The conformal Hubble rate is :math:`\\mathcal{H}=aH` (code name ``Hp``).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from . import constants as const
from .params import Params


class BackgroundCosmology:
    """Solve and spline the homogeneous background quantities.

    Provides :math:`\\mathcal{H}(x)` and derivatives, conformal time
    :math:`\\eta(x)`, cosmic time :math:`t(x)`, density fractions and
    comoving distance over :math:`x\\in[-20,5]`.
    """

    def __init__(self, params: Params, x_start: float = -20.0,
                 x_end: float = 5.0, n: int = 4000) -> None:
        self.p = params
        self.x_start = x_start
        self.x_end = x_end
        self.x = np.linspace(x_start, x_end, n)
        # Massive-neutrino contribution to E = (Hp/H0)^2 is a^2 f_ncdm(a)
        # = e^{2x} rho_ncdm/rho_crit0. Spline it (and use spline derivatives)
        # only when present; massless runs keep the analytic path bit-for-bit
        # (Omega_ur == Omega_nu when Sigma_mnu = 0).
        self._mnu_g = None
        if getattr(params, "N_ncdm", 0) > 0 and params._massive is not None:
            xg = np.linspace(x_start, x_end, n)
            g = np.exp(2.0 * xg) * params._massive.f_ncdm(np.exp(xg))
            self._mnu_g = CubicSpline(xg, g)
        self._solve_eta_t()

    # --- Hubble rate -----------------------------------------------------------
    def Hp(self, x):
        """Conformal Hubble rate :math:`\\mathcal{H}=aH` [1/s]."""
        E, _, _ = self._e_and_deriv(x)
        return self.p.H0 * np.sqrt(E)

    def _e_and_deriv(self, x):
        """Return (E, dE/dx, d2E/dx2) where Hp = H0 sqrt(E)."""
        p = self.p
        x = np.asarray(x, dtype=float)
        mat = (p.Omega_b + p.Omega_c) * np.exp(-x)
        rad = (p.Omega_gamma + p.Omega_ur) * np.exp(-2.0 * x)
        lam = p.Omega_Lambda * np.exp(2.0 * x)
        E = mat + rad + p.Omega_k + lam
        dE = -mat - 2.0 * rad + 2.0 * lam
        d2E = mat + 4.0 * rad + 4.0 * lam
        if self._mnu_g is not None:
            E = E + self._mnu_g(x)
            dE = dE + self._mnu_g(x, 1)
            d2E = d2E + self._mnu_g(x, 2)
        return E, dE, d2E

    def dHp(self, x):
        """:math:`\\mathcal{H}'=d\\mathcal{H}/dx` (analytic)."""
        E, dE, _ = self._e_and_deriv(x)
        return self.p.H0 * dE / (2.0 * np.sqrt(E))

    def ddHp(self, x):
        """:math:`\\mathcal{H}''` (analytic)."""
        E, dE, d2E = self._e_and_deriv(x)
        return self.p.H0 * (2.0 * E * d2E - dE ** 2) / (4.0 * E ** 1.5)

    def H(self, x):
        """Hubble rate :math:`H=\\mathcal{H}/a` [1/s]."""
        x = np.asarray(x, dtype=float)
        return self.Hp(x) * np.exp(-x)

    # --- density fractions -----------------------------------------------------
    def Omega(self, x):
        """Return a dict of time-dependent density fractions :math:`\\Omega_i(x)`."""
        p = self.p
        x = np.asarray(x, dtype=float)
        # Omega_i(x) = rho_i(x)/rho_crit(x) = (rho_i/rho_crit0)/(H/H0)^2.
        a = np.exp(x)
        fac = p.H0 ** 2 / (self.H(x)) ** 2
        # Neutrinos: massless part scales as a^-4; massive part via f_ncdm(a).
        nu = p.Omega_ur * a ** -4
        if self._mnu_g is not None:
            nu = nu + self._mnu_g(x) * a ** -2   # a^2 f_ncdm -> f_ncdm = (a^2 f)/a^2
        return {
            "b": p.Omega_b * a ** -3 * fac,
            "c": p.Omega_c * a ** -3 * fac,
            "gamma": p.Omega_gamma * a ** -4 * fac,
            "nu": nu * fac,
            "k": p.Omega_k * a ** -2 * fac,
            "Lambda": p.Omega_Lambda * fac,
        }

    # --- conformal & cosmic time ----------------------------------------------
    def _solve_eta_t(self) -> None:
        # deta/dx = c/Hp ; dt/dx = 1/H.  Radiation-dominated initial values.
        x0 = self.x_start
        eta0 = const.c / self.Hp(x0)
        t0 = 1.0 / (2.0 * self.H(x0))

        def rhs(x, y):
            return [const.c / self.Hp(x), 1.0 / self.H(x)]

        sol = solve_ivp(
            rhs, (x0, self.x_end), [eta0, t0], t_eval=self.x,
            method="RK45", rtol=1e-10, atol=1e-12, dense_output=True,
        )
        self._eta = CubicSpline(self.x, sol.y[0])
        self._t = CubicSpline(self.x, sol.y[1])
        self.eta0 = float(self._eta(0.0))  # conformal time today
        self.t0 = float(self._t(0.0))

    def eta(self, x):
        """Conformal time :math:`\\eta(x)` [s]."""
        return self._eta(np.asarray(x, dtype=float))

    def t(self, x):
        """Cosmic time :math:`t(x)` [s]."""
        return self._t(np.asarray(x, dtype=float))

    # --- distances -------------------------------------------------------------
    def comoving_distance(self, x):
        """Comoving distance :math:`\\chi(x)=\\eta_0-\\eta(x)` [m].

        Here :math:`\\eta` already carries a factor of :math:`c`
        (``deta/dx = c/Hp``), so it is a comoving length, not a time.
        """
        return self.eta0 - self.eta(x)

    def angular_distance_factor(self, chi):
        """Curvature-corrected radial function :math:`S_k(\\chi)` [m].

        For flat space returns ``chi``; otherwise the ``sinh``/``sin`` form
        (chapter 13). ``chi`` is a *comoving distance* in metres.
        """
        p = self.p
        if abs(p.Omega_k) < 1e-12:
            return chi
        # Curvature scale: H0 sqrt(|Omega_k|)/c.
        kappa = p.H0 * np.sqrt(abs(p.Omega_k)) / const.c
        if p.Omega_k > 0:  # open
            return np.sinh(kappa * chi) / kappa
        return np.sin(kappa * chi) / kappa

    # --- equality redshift -----------------------------------------------------
    def z_eq(self) -> float:
        """Matter-radiation equality redshift."""
        a_eq = self.p.Omega_r / self.p.Omega_m
        return 1.0 / a_eq - 1.0
