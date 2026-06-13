"""Linear perturbations / Boltzmann hierarchy (``03_numerics_spec.md`` §3).

Scalar perturbations, conformal-Newtonian gauge, temperature only (no
polarization). For each wavenumber ``k`` the system is integrated in a
tight-coupling phase followed by the full hierarchy, from ``x_start`` to ``x=0``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from . import constants as const
from .params import Params
from .background import BackgroundCosmology
from .recombination import Recombination


@dataclass
class PerturbationResult:
    """Splined perturbation variables for a single ``k`` (SI ``k`` in 1/m)."""

    k: float
    x: np.ndarray
    splines: dict  # name -> CubicSpline

    def __getitem__(self, name: str):
        return self.splines[name]


class PerturbationSolver:
    """Integrate the Boltzmann hierarchy for individual wavenumbers."""

    def __init__(self, bg: BackgroundCosmology, rec: Recombination,
                 params: Params, lmax: int = 8, lmax_nu: int = 10,
                 x_start: float = -18.0, n_x: int = 2000) -> None:
        self.bg = bg
        self.rec = rec
        self.p = params
        self.lmax = lmax
        self.lmax_nu = lmax_nu
        self.x_start = x_start
        self.n_x = n_x

    # --- variable layout -------------------------------------------------------
    # Full system state vector:
    #   [delta_c, v_c, delta_b, v_b, Phi,
    #    Theta_0..Theta_lmax,
    #    N_0..N_lmax_nu]
    def _idx(self):
        i = {}
        i["delta_c"], i["v_c"], i["delta_b"], i["v_b"], i["Phi"] = 0, 1, 2, 3, 4
        base = 5
        i["Theta"] = base
        base += self.lmax + 1
        i["N"] = base
        base += self.lmax_nu + 1
        i["n_full"] = base
        return i

    # --- algebraic closure for Psi --------------------------------------------
    def _Psi(self, x, Phi, Theta2, N2, k):
        """(G2): :math:`\\Psi=-\\Phi-\\frac{12H_0^2}{c^2k^2a^2}(\\Omega_\\gamma\\Theta_2+\\Omega_\\nu\\mathcal N_2)`."""
        p = self.p
        a = np.exp(x)
        return -Phi - 12.0 * p.H0 ** 2 / (const.c ** 2 * k ** 2 * a ** 2) \
            * (p.Omega_gamma * Theta2 + p.Omega_nu * N2)

    # --- initial conditions (§3.2) --------------------------------------------
    def _initial_conditions(self, x0, k):
        p = self.p
        Hp = float(self.bg.Hp(x0))
        dtau = float(self.rec.dtau(x0))
        ckH = const.c * k / Hp
        f_nu = p.f_nu

        Psi = -1.0 / (1.5 + 0.4 * f_nu)
        Phi = -(1.0 + 0.4 * f_nu) * Psi
        delta_c = delta_b = -1.5 * Psi
        v_c = v_b = -0.5 * ckH * Psi

        Theta = np.zeros(self.lmax + 1)
        Theta[0] = -0.5 * Psi
        Theta[1] = ckH / 6.0 * Psi
        Theta[2] = -20.0 * ckH / (45.0 * dtau) * Theta[1]
        for l in range(3, self.lmax + 1):
            Theta[l] = -l / (2.0 * l + 1.0) * ckH / dtau * Theta[l - 1]

        N = np.zeros(self.lmax_nu + 1)
        N[0] = -0.5 * Psi
        N[1] = ckH / 6.0 * Psi
        a = np.exp(x0)
        N[2] = -(const.c ** 2 * k ** 2 * a ** 2 * (Phi + Psi)) \
            / (12.0 * p.H0 ** 2 * p.Omega_nu)
        for l in range(3, self.lmax_nu + 1):
            N[l] = ckH / ((2.0 * l + 1.0)) * N[l - 1]

        y = np.zeros(self._idx()["n_full"])
        i = self._idx()
        y[i["delta_c"]] = delta_c
        y[i["v_c"]] = v_c
        y[i["delta_b"]] = delta_b
        y[i["v_b"]] = v_b
        y[i["Phi"]] = Phi
        y[i["Theta"]:i["Theta"] + self.lmax + 1] = Theta
        y[i["N"]:i["N"] + self.lmax_nu + 1] = N
        return y

    # --- common gravity & RHS pieces ------------------------------------------
    def _R_coupling(self, x):
        """Coupling coefficient :math:`R=4\\Omega_{\\gamma0}/(3\\Omega_{b0})e^{-x}` (= 1/R_s)."""
        return 4.0 * self.p.Omega_gamma / (3.0 * self.p.Omega_b) * np.exp(-x)

    def _Phi_prime(self, x, k, Phi, Psi, delta_c, delta_b, Theta0, N0):
        """(G1)."""
        p = self.p
        Hp = float(self.bg.Hp(x))
        a = np.exp(x)
        ck2H2 = (const.c * k / Hp) ** 2
        src = (
            p.Omega_c * np.exp(-x) * delta_c
            + p.Omega_b * np.exp(-x) * delta_b
            + 4.0 * p.Omega_gamma * np.exp(-2.0 * x) * Theta0
            + 4.0 * p.Omega_nu * np.exp(-2.0 * x) * N0
        )
        return Psi - ck2H2 / 3.0 * Phi + (p.H0 ** 2) / (2.0 * Hp ** 2) * src

    # --- full hierarchy RHS ----------------------------------------------------
    def _rhs_full(self, x, y, k):
        p = self.p
        i = self._idx()
        Hp = float(self.bg.Hp(x))
        dtau = float(self.rec.dtau(x))
        ckH = const.c * k / Hp

        delta_c = y[i["delta_c"]]
        v_c = y[i["v_c"]]
        delta_b = y[i["delta_b"]]
        v_b = y[i["v_b"]]
        Phi = y[i["Phi"]]
        Th = y[i["Theta"]:i["Theta"] + self.lmax + 1]
        N = y[i["N"]:i["N"] + self.lmax_nu + 1]

        Psi = self._Psi(x, Phi, Th[2], N[2], k)
        Phip = self._Phi_prime(x, k, Phi, Psi, delta_c, delta_b, Th[0], N[0])

        dy = np.zeros_like(y)
        # Matter (M1-M4).
        dy[i["delta_c"]] = ckH * v_c - 3.0 * Phip
        dy[i["v_c"]] = -v_c - ckH * Psi
        dy[i["delta_b"]] = ckH * v_b - 3.0 * Phip
        R = self._R_coupling(x)
        dy[i["v_b"]] = -v_b - ckH * Psi + dtau * R * (3.0 * Th[1] + v_b)

        # Photons (P1-P4).
        dTh = np.zeros(self.lmax + 1)
        dTh[0] = -ckH * Th[1] - Phip
        dTh[1] = (ckH / 3.0) * Th[0] - (2.0 * ckH / 3.0) * Th[2] \
            + (ckH / 3.0) * Psi + dtau * (Th[1] + v_b / 3.0)
        for l in range(2, self.lmax):
            src = dtau * (Th[l] - 0.1 * Th[2] * (1 if l == 2 else 0))
            dTh[l] = (l * ckH / (2.0 * l + 1.0)) * Th[l - 1] \
                - ((l + 1.0) * ckH / (2.0 * l + 1.0)) * Th[l + 1] + src
        L = self.lmax
        eta = float(self.bg.eta(x))
        dTh[L] = ckH * Th[L - 1] \
            - const.c * (L + 1.0) / (Hp * eta) * Th[L] + dtau * Th[L]
        dy[i["Theta"]:i["Theta"] + self.lmax + 1] = dTh

        # Neutrinos (N1-N4): photon eqns with tau'->0, Theta->N.
        dN = np.zeros(self.lmax_nu + 1)
        dN[0] = -ckH * N[1] - Phip
        dN[1] = (ckH / 3.0) * N[0] - (2.0 * ckH / 3.0) * N[2] + (ckH / 3.0) * Psi
        for l in range(2, self.lmax_nu):
            dN[l] = (l * ckH / (2.0 * l + 1.0)) * N[l - 1] \
                - ((l + 1.0) * ckH / (2.0 * l + 1.0)) * N[l + 1]
        Ln = self.lmax_nu
        dN[Ln] = ckH * N[Ln - 1] - const.c * (Ln + 1.0) / (Hp * eta) * N[Ln]
        dy[i["N"]:i["N"] + self.lmax_nu + 1] = dN

        dy[i["Phi"]] = Phip
        return dy

    # --- tight-coupling RHS ----------------------------------------------------
    # TC state: [delta_c, v_c, delta_b, v_b, Phi, Theta0, Theta1, N_0..N_lmax_nu]
    def _tc_idx(self):
        i = {}
        i["delta_c"], i["v_c"], i["delta_b"], i["v_b"], i["Phi"] = 0, 1, 2, 3, 4
        i["Theta0"], i["Theta1"] = 5, 6
        i["N"] = 7
        i["n_tc"] = 7 + self.lmax_nu + 1
        return i

    def _tc_theta2(self, x, k, Theta1):
        Hp = float(self.bg.Hp(x))
        dtau = float(self.rec.dtau(x))
        ckH = const.c * k / Hp
        return -20.0 * ckH / (45.0 * dtau) * Theta1

    def _rhs_tc(self, x, y, k):
        p = self.p
        i = self._tc_idx()
        bg = self.bg
        Hp = float(bg.Hp(x))
        dHp = float(bg.dHp(x))
        dtau = float(self.rec.dtau(x))
        ddtau = float(self.rec.ddtau(x))
        ckH = const.c * k / Hp

        delta_c = y[i["delta_c"]]
        v_c = y[i["v_c"]]
        delta_b = y[i["delta_b"]]
        v_b = y[i["v_b"]]
        Phi = y[i["Phi"]]
        Th0 = y[i["Theta0"]]
        Th1 = y[i["Theta1"]]
        N = y[i["N"]:i["N"] + self.lmax_nu + 1]
        Th2 = self._tc_theta2(x, k, Th1)

        Psi = self._Psi(x, Phi, Th2, N[2], k)
        Phip = self._Phi_prime(x, k, Phi, Psi, delta_c, delta_b, Th0, N[0])

        dy = np.zeros_like(y)
        R = self._R_coupling(x)

        # Theta0' from (P1).
        dTh0 = -ckH * Th1 - Phip

        # q and v_b', Theta1' from TC equations (Callin 70-72).
        num = (
            -((1.0 - R) * dtau + (1.0 + R) * ddtau) * (3.0 * Th1 + v_b)
            - ckH * Psi
            + (1.0 - dHp / Hp) * ckH * (-Th0 + 2.0 * Th2)
            - ckH * dTh0
        )
        den = (1.0 + R) * dtau + dHp / Hp - 1.0
        q = num / den

        dv_b = (1.0 / (1.0 + R)) * (
            -v_b - ckH * Psi
            + R * (q + ckH * (-Th0 + 2.0 * Th2) - ckH * Psi)
        )
        dTh1 = (q - dv_b) / 3.0

        dy[i["delta_c"]] = ckH * v_c - 3.0 * Phip
        dy[i["v_c"]] = -v_c - ckH * Psi
        dy[i["delta_b"]] = ckH * v_b - 3.0 * Phip
        dy[i["v_b"]] = dv_b
        dy[i["Phi"]] = Phip
        dy[i["Theta0"]] = dTh0
        dy[i["Theta1"]] = dTh1

        # Neutrinos evolve fully even in TC.
        dN = np.zeros(self.lmax_nu + 1)
        eta = float(bg.eta(x))
        dN[0] = -ckH * N[1] - Phip
        dN[1] = (ckH / 3.0) * N[0] - (2.0 * ckH / 3.0) * N[2] + (ckH / 3.0) * Psi
        for l in range(2, self.lmax_nu):
            dN[l] = (l * ckH / (2.0 * l + 1.0)) * N[l - 1] \
                - ((l + 1.0) * ckH / (2.0 * l + 1.0)) * N[l + 1]
        Ln = self.lmax_nu
        dN[Ln] = ckH * N[Ln - 1] - const.c * (Ln + 1.0) / (Hp * eta) * N[Ln]
        dy[i["N"]:i["N"] + self.lmax_nu + 1] = dN
        return dy

    # --- tight coupling end ----------------------------------------------------
    def _x_tc_end(self, k):
        """Largest x where TC still holds (§3.3 conditions)."""
        xs = np.linspace(self.x_start, 0.0, 4000)
        Hp = self.bg.Hp(xs)
        dtau = np.abs(self.rec.dtau(xs))
        ckH = const.c * k / Hp
        Xe = self.rec.Xe(xs)
        cond = (dtau > 10.0) & (dtau > 10.0 * ckH) & (Xe > 0.99)
        if not cond.any():
            return self.x_start  # TC never valid
        # First x where condition fails after start.
        idx = np.argmax(~cond)
        if idx == 0:
            return self.x_start
        return xs[idx]

    # --- solve a single k ------------------------------------------------------
    def solve(self, k: float) -> PerturbationResult:
        """Integrate the hierarchy for SI wavenumber ``k`` [1/m].

        Two phases (tight coupling, then full hierarchy) are integrated and the
        recorded variables are reconstructed with vectorized array operations
        (no per-point Python loop) for speed.
        """
        x_grid = np.linspace(self.x_start, 0.0, self.n_x)
        x_tc = self._x_tc_end(k)

        i_tc = self._tc_idx()
        i_full = self._idx()
        nnu = self.lmax_nu + 1

        # Initial conditions (full layout), then project to TC layout.
        y0_full = self._initial_conditions(self.x_start, k)
        y0_tc = np.zeros(i_tc["n_tc"])
        y0_tc[i_tc["delta_c"]] = y0_full[i_full["delta_c"]]
        y0_tc[i_tc["v_c"]] = y0_full[i_full["v_c"]]
        y0_tc[i_tc["delta_b"]] = y0_full[i_full["delta_b"]]
        y0_tc[i_tc["v_b"]] = y0_full[i_full["v_b"]]
        y0_tc[i_tc["Phi"]] = y0_full[i_full["Phi"]]
        y0_tc[i_tc["Theta0"]] = y0_full[i_full["Theta"] + 0]
        y0_tc[i_tc["Theta1"]] = y0_full[i_full["Theta"] + 1]
        y0_tc[i_tc["N"]:i_tc["N"] + nnu] = \
            y0_full[i_full["N"]:i_full["N"] + nnu]

        names = (["delta_c", "v_c", "delta_b", "v_b", "Phi"]
                 + [f"Theta{l}" for l in range(self.lmax + 1)]
                 + [f"N{l}" for l in range(self.lmax_nu + 1)]
                 + ["Psi"])
        cols = {nm: [] for nm in names}
        x_chunks = []

        mask_tc = x_grid <= x_tc
        have_tc = mask_tc.sum() >= 2

        # Phase 1: tight coupling.
        if have_tc:
            xt = x_grid[mask_tc]
            sol = solve_ivp(
                lambda xx, yy: self._rhs_tc(xx, yy, k),
                (xt[0], xt[-1]), y0_tc, t_eval=xt,
                method="LSODA", rtol=1e-8, atol=1e-10,
            )
            self._collect_tc(cols, x_chunks, sol.t, sol.y, k)
            y0_full = self._tc_to_full(xt[-1], sol.y[:, -1], k)
            x_full_start = xt[-1]
        else:
            x_full_start = self.x_start

        # Phase 2: full hierarchy.
        xf = x_grid[x_grid > x_full_start]
        if len(xf) >= 1:
            xf = np.concatenate(([x_full_start], xf))
            sol = solve_ivp(
                lambda xx, yy: self._rhs_full(xx, yy, k),
                (xf[0], xf[-1]), y0_full, t_eval=xf,
                method="LSODA", rtol=1e-8, atol=1e-10,
            )
            sl = slice(1, None) if have_tc else slice(None)
            self._collect_full(cols, x_chunks, sol.t[sl], sol.y[:, sl], k)

        x_out = np.concatenate(x_chunks)
        splines = {nm: CubicSpline(x_out, np.concatenate(cols[nm]))
                   for nm in names}
        return PerturbationResult(k=k, x=x_out, splines=splines)

    def _tc_to_full(self, x, y_tc, k):
        i_tc = self._tc_idx()
        i_full = self._idx()
        y = np.zeros(i_full["n_full"])
        y[i_full["delta_c"]] = y_tc[i_tc["delta_c"]]
        y[i_full["v_c"]] = y_tc[i_tc["v_c"]]
        y[i_full["delta_b"]] = y_tc[i_tc["delta_b"]]
        y[i_full["v_b"]] = y_tc[i_tc["v_b"]]
        y[i_full["Phi"]] = y_tc[i_tc["Phi"]]
        Th0 = y_tc[i_tc["Theta0"]]
        Th1 = y_tc[i_tc["Theta1"]]
        y[i_full["Theta"] + 0] = Th0
        y[i_full["Theta"] + 1] = Th1
        Hp = float(self.bg.Hp(x))
        dtau = float(self.rec.dtau(x))
        ckH = const.c * k / Hp
        Th2 = -20.0 * ckH / (45.0 * dtau) * Th1
        y[i_full["Theta"] + 2] = Th2
        for l in range(3, self.lmax + 1):
            y[i_full["Theta"] + l] = \
                -l / (2.0 * l + 1.0) * ckH / dtau * y[i_full["Theta"] + l - 1]
        y[i_full["N"]:i_full["N"] + self.lmax_nu + 1] = \
            y_tc[i_tc["N"]:i_tc["N"] + self.lmax_nu + 1]
        return y

    def _collect_tc(self, cols, x_chunks, xs, Y, k):
        """Vectorized recording of a tight-coupling solution block."""
        i = self._tc_idx()
        nnu = self.lmax_nu + 1
        Hp = self.bg.Hp(xs)
        dtau = self.rec.dtau(xs)
        ckH = const.c * k / Hp

        Th1 = Y[i["Theta1"]]
        Th2 = -20.0 * ckH / (45.0 * dtau) * Th1
        N = Y[i["N"]:i["N"] + nnu]
        Phi = Y[i["Phi"]]
        a = np.exp(xs)
        Psi = -Phi - 12.0 * self.p.H0 ** 2 / (const.c ** 2 * k ** 2 * a ** 2) \
            * (self.p.Omega_gamma * Th2 + self.p.Omega_nu * N[2])

        x_chunks.append(xs)
        cols["delta_c"].append(Y[i["delta_c"]])
        cols["v_c"].append(Y[i["v_c"]])
        cols["delta_b"].append(Y[i["delta_b"]])
        cols["v_b"].append(Y[i["v_b"]])
        cols["Phi"].append(Phi)
        cols["Theta0"].append(Y[i["Theta0"]])
        cols["Theta1"].append(Th1)
        cols["Theta2"].append(Th2)
        prev = Th2
        for l in range(3, self.lmax + 1):
            prev = -l / (2.0 * l + 1.0) * ckH / dtau * prev
            cols[f"Theta{l}"].append(prev)
        for l in range(nnu):
            cols[f"N{l}"].append(N[l])
        cols["Psi"].append(Psi)

    def _collect_full(self, cols, x_chunks, xs, Y, k):
        """Vectorized recording of a full-hierarchy solution block."""
        i = self._idx()
        nnu = self.lmax_nu + 1
        Th = Y[i["Theta"]:i["Theta"] + self.lmax + 1]
        N = Y[i["N"]:i["N"] + nnu]
        Phi = Y[i["Phi"]]
        a = np.exp(xs)
        Psi = -Phi - 12.0 * self.p.H0 ** 2 / (const.c ** 2 * k ** 2 * a ** 2) \
            * (self.p.Omega_gamma * Th[2] + self.p.Omega_nu * N[2])

        x_chunks.append(xs)
        cols["delta_c"].append(Y[i["delta_c"]])
        cols["v_c"].append(Y[i["v_c"]])
        cols["delta_b"].append(Y[i["delta_b"]])
        cols["v_b"].append(Y[i["v_b"]])
        cols["Phi"].append(Phi)
        for l in range(self.lmax + 1):
            cols[f"Theta{l}"].append(Th[l])
        for l in range(nnu):
            cols[f"N{l}"].append(N[l])
        cols["Psi"].append(Psi)
