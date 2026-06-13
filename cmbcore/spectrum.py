"""CMB temperature power spectrum (``03_numerics_spec.md`` §4).

Builds the line-of-sight source function, performs the Bessel integral to get
the transfer functions :math:`\\Theta_\\ell(k)`, and assembles
:math:`C_\\ell^{TT}`.
"""

from __future__ import annotations

import os
from functools import partial
from multiprocessing import Pool

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn

from . import constants as const
from .params import Params
from .background import BackgroundCosmology
from .recombination import Recombination
from .perturbations import PerturbationSolver


class BesselBank:
    """Pre-splined spherical Bessel functions :math:`j_\\ell(z)` (§4.2).

    Evaluating ``spherical_jn`` directly for high orders over many points is the
    dominant cost, so we tabulate each :math:`j_\\ell` on a shared ``z`` grid
    (:math:`\\Delta z\\le 2\\pi/16`) once and spline-interpolate thereafter.
    Below ``z ~ l`` the function is exponentially small and set to zero.
    """

    def __init__(self, ells, z_max, dz=2.0 * np.pi / 16.0):
        self.ells = np.asarray(ells)
        z = np.arange(0.0, z_max + dz, dz)
        self._z = z
        self._splines = {}
        for l in self.ells:
            jl = spherical_jn(int(l), z)
            self._splines[int(l)] = CubicSpline(z, jl)

    def eval(self, l, zq):
        """Evaluate :math:`j_\\ell` at query points ``zq`` (clipped to table)."""
        sp = self._splines[int(l)]
        zq = np.clip(zq, 0.0, self._z[-1])
        return sp(zq)


# --- worker for k-parallel transfer (module level so it is picklable) --------
_WORKER = {}


def _init_worker(params, lmax, lmax_nu, x_los, eta0, eta_x, ells, z_max):
    bg = BackgroundCosmology(params)
    rec = Recombination(bg, params)
    solver = PerturbationSolver(bg, rec, params, lmax=lmax, lmax_nu=lmax_nu)
    _WORKER.update(dict(params=params, bg=bg, rec=rec, solver=solver,
                        x_los=x_los, eta0=eta0, eta_x=eta_x,
                        ells=np.asarray(ells),
                        bessel=BesselBank(ells, z_max)))


def _transfer_one_k(k):
    """Return Theta_ell(k) for all ells (R-normalized) for one wavenumber."""
    ps = PowerSpectrum.__new__(PowerSpectrum)
    ps.p = _WORKER["params"]
    ps.bg = _WORKER["bg"]
    ps.rec = _WORKER["rec"]
    ps.solver = _WORKER["solver"]
    ps.lmax = _WORKER["solver"].lmax
    x_los = _WORKER["x_los"]
    eta0 = _WORKER["eta0"]
    eta_x = _WORKER["eta_x"]
    ells = _WORKER["ells"]
    bessel = _WORKER["bessel"]

    res = ps.solver.solve(k)
    S = ps._source(res, k, x_los) / ps._R_ini(k)
    kchi = k * (eta0 - eta_x)
    out = np.empty(len(ells))
    for il, l in enumerate(ells):
        out[il] = np.trapezoid(S * bessel.eval(l, kchi), x_los)
    return out


# Default multipole sampling (§4.3): dense at low ell, coarser at high ell.
def default_ells() -> np.ndarray:
    ells = list(range(2, 11))
    ells += list(range(12, 31, 2))
    ells += list(range(35, 101, 5))
    ells += list(range(110, 301, 10))
    ells += list(range(325, 1501, 25))
    return np.array(sorted(set(ells)))


def k_grid(kmin=5e-5, kmax=0.35, nk=250):
    """Quadratic ``k`` grid in 1/Mpc (low-k dense), returned in 1/Mpc (§4.3).

    Low-k density captures the SW/ISW contributions; the acoustic range is
    resolved by raising ``nk`` (convergence study, 04 §3).
    """
    i = np.arange(nk + 1)
    return kmin + (kmax - kmin) * (i / nk) ** 2


class PowerSpectrum:
    """End-to-end :math:`C_\\ell^{TT}` calculator."""

    def __init__(self, params: Params,
                 bg: BackgroundCosmology | None = None,
                 rec: Recombination | None = None,
                 lmax: int = 8, lmax_nu: int = 10) -> None:
        self.p = params
        self.bg = bg or BackgroundCosmology(params)
        self.rec = rec or Recombination(self.bg, params)
        self.solver = PerturbationSolver(self.bg, self.rec, params,
                                         lmax=lmax, lmax_nu=lmax_nu)
        self.lmax = lmax

    # --- primordial spectrum ---------------------------------------------------
    def P_R(self, k_si):
        """Dimensionless primordial spectrum :math:`\\mathcal P_\\mathcal R(k)`.

        ``k_si`` in 1/m; pivot is in 1/Mpc.
        """
        k_mpc = k_si * const.Mpc
        return self.p.A_s * (k_mpc / self.p.k_pivot) ** (self.p.n_s - 1.0)

    # --- comoving curvature at the initial slice (§4.3 normalization) ---------
    def _R_ini(self, k_si) -> float:
        """Comoving curvature :math:`\\mathcal R_{\\rm ini}` for the standard ICs.

        The §3.2 initial conditions set
        :math:`\\Psi_{\\rm ini}=-1/(\\tfrac32+\\tfrac25 f_\\nu)
        =-\\tfrac23\\mathcal R/(1+\\tfrac{4}{15}f_\\nu)`, i.e. the standard
        super-horizon adiabatic relation with :math:`\\mathcal R=1`. The modes
        are therefore already normalized to unit comoving curvature, so the
        conversion factor is exactly 1. (Confirmed against the analytic
        super-horizon value; the final amplitude defence is the CLASS
        comparison, ``03`` §5.3.)
        """
        return 1.0

    # --- source function (§4.1) ------------------------------------------------
    def _source(self, res, k_si, x):
        """Line-of-sight source :math:`\\tilde S(k,x)` on grid ``x``.

        Returns the total source (sum of the four terms).
        """
        bg, rec = self.bg, self.rec
        Hp = bg.Hp(x)
        g = rec.g_tilde(x)
        etau = np.exp(-rec.tau(x))

        Th0 = res["Theta0"](x)
        Th2 = res["Theta2"](x)
        Psi = res["Psi"](x)
        Phi = res["Phi"](x)
        Psip = res["Psi"](x, 1)
        Phip = res["Phi"](x, 1)
        v_b = res["v_b"](x)

        ck = const.c * k_si

        # Term 1: g (Theta0 + Psi + Theta2/4).
        T1 = g * (Th0 + Psi + 0.25 * Th2)
        # Term 2: e^{-tau} (Psi' - Phi').
        T2 = etau * (Psip - Phip)
        # Term 3: -(1/ck) d/dx (Hp g v_b).
        Hgv = Hp * g * v_b
        Hgv_sp = CubicSpline(x, Hgv)
        T3 = -(1.0 / ck) * Hgv_sp(x, 1)
        # Term 4: (3/4 c^2 k^2) d/dx [ Hp d/dx (Hp g Theta2) ].
        HgT2 = Hp * g * Th2
        HgT2_sp = CubicSpline(x, HgT2)
        inner = Hp * HgT2_sp(x, 1)
        inner_sp = CubicSpline(x, inner)
        T4 = (3.0 / (4.0 * const.c ** 2 * k_si ** 2)) * inner_sp(x, 1)

        return T1 + T2 + T3 + T4

    # --- transfer functions (§4.2) --------------------------------------------
    def _default_x_los(self, n=1600):
        """Single monotonic line-of-sight x-grid, clustered near recombination.

        Avoids the abrupt density jumps of a concatenated grid (which inject
        spline-derivative artifacts into the source). The recombination window
        gets Δx ~ 10⁻³ (spec §4.2) while the late-time ISW tail stays sampled
        at Δx ≲ 0.05.
        """
        x0 = self.solver.x_start
        # ~half the points densely across the visibility window, the rest
        # spread over the early and late epochs.
        x_rec = np.linspace(-7.8, -6.2, n // 2)            # Δx ~ 1.6/1500 ~1e-3
        x_early = np.linspace(x0, -7.8, n // 6, endpoint=False)
        x_late = np.linspace(-6.2, 0.0, n // 3)[1:]
        return np.unique(np.concatenate([x_early, x_rec, x_late]))

    def transfer(self, ells, ks_si=None, x_los=None, verbose=False,
                 nproc=None):
        """Compute :math:`\\Theta_\\ell(k)` for all ``ells`` and ``ks_si``.

        Returns ``(ks_si, Theta)`` with ``Theta`` of shape ``(len(ells), nk)``,
        already normalized to :math:`\\mathcal R=1`. ``nproc`` parallelizes over
        ``k`` (default: CPU count; set 1 to disable).
        """
        if ks_si is None:
            ks_si = k_grid() / const.Mpc
        if x_los is None:
            x_los = self._default_x_los()

        eta0 = self.bg.eta0
        eta_x = self.bg.eta(x_los)
        ells = np.asarray(ells)
        nk = len(ks_si)
        z_max = float(ks_si.max() * (eta0 - eta_x.min())) * 1.01

        if nproc is None:
            nproc = max(1, os.cpu_count() or 1)

        if nproc == 1:
            bessel = BesselBank(ells, z_max)
            Theta = np.zeros((len(ells), nk))
            for ik, k in enumerate(ks_si):
                res = self.solver.solve(k)
                S = self._source(res, k, x_los) / self._R_ini(k)
                kchi = k * (eta0 - eta_x)
                for il, l in enumerate(ells):
                    Theta[il, ik] = np.trapezoid(
                        S * bessel.eval(l, kchi), x_los)
                if verbose and ik % 25 == 0:
                    print(f"  transfer k {ik+1}/{nk}")
            return ks_si, Theta

        init_args = (self.p, self.solver.lmax, self.solver.lmax_nu,
                     x_los, eta0, eta_x, ells, z_max)
        with Pool(processes=nproc, initializer=_init_worker,
                  initargs=init_args) as pool:
            results = pool.map(_transfer_one_k, list(ks_si),
                               chunksize=max(1, nk // (nproc * 4)))
        Theta = np.array(results).T  # (n_ell, n_k)
        return ks_si, Theta

    # --- C_ell assembly (§4.3) -------------------------------------------------
    def cls(self, ells=None, ks_si=None, verbose=False):
        """Return ``(ells, C_ell)`` in dimensionless units (Theta^2)."""
        if ells is None:
            ells = default_ells()
        ks_si, Theta = self.transfer(ells, ks_si, verbose=verbose)
        # C_ell = 4 pi int P_R(k) |Theta_ell|^2 dk/k = 4 pi int P_R |Theta|^2 dlnk.
        # A plain trapezoid is used: for the (necessarily) under-resolved
        # oscillatory integrand |Theta_ell(k)|^2, spline/Simpson schemes overshoot
        # between samples and amplify l-to-l jitter, whereas the trapezoid is
        # stable. Residual mid-/high-l jitter is reduced by raising n_k (the
        # convergence study, 04 §3).
        P = self.P_R(ks_si)
        lnk = np.log(ks_si)
        Cl = np.zeros(len(ells))
        for il in range(len(ells)):
            Cl[il] = 4.0 * np.pi * np.trapezoid(P * Theta[il] ** 2, lnk)
        return np.asarray(ells), Cl

    def dls(self, ells=None, ks_si=None, verbose=False):
        """Return ``(ells, D_ell)`` with :math:`D_\\ell=\\ell(\\ell+1)C_\\ell T_{CMB}^2/2\\pi` [μK²]."""
        ells, Cl = self.cls(ells, ks_si, verbose=verbose)
        T2 = (self.p.T_CMB * 1e6) ** 2  # K -> uK, squared
        Dl = ells * (ells + 1.0) / (2.0 * np.pi) * Cl * T2
        return ells, Dl
