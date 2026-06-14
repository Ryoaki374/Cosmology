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
    (:math:`\\Delta z=2\\pi/20`) once and spline-interpolate thereafter.
    (Tested: refining to 2pi/40 leaves C_ell unchanged, so the table is not the
    jitter source.)
    """

    def __init__(self, ells, z_max, dz=2.0 * np.pi / 20.0):
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


def _init_worker(params, lmax, lmax_nu, x_los):
    bg = BackgroundCosmology(params)
    rec = Recombination(bg, params)
    solver = PerturbationSolver(bg, rec, params, lmax=lmax, lmax_nu=lmax_nu)
    _WORKER.update(dict(params=params, bg=bg, rec=rec, solver=solver,
                        x_los=x_los))


def _worker_ps():
    ps = PowerSpectrum.__new__(PowerSpectrum)
    ps.p = _WORKER["params"]
    ps.bg = _WORKER["bg"]
    ps.rec = _WORKER["rec"]
    ps.solver = _WORKER["solver"]
    ps.lmax = _WORKER["solver"].lmax
    return ps


def _sources_TE_one_k(k):
    """Return the temperature and E-mode sources (2, n_x) for one wavenumber."""
    ps = _worker_ps()
    x_los = _WORKER["x_los"]
    res = ps.solver.solve(k)
    norm = ps._R_ini(k)
    S_T = ps._source(res, k, x_los) / norm
    S_E = ps._source_E(res, k, x_los) / norm
    return np.stack([S_T, S_E])


def _source_one_k(k):
    """Return the (R-normalized) temperature source S(x) for one wavenumber.

    The source is smooth in k, so it is computed on the (coarse) solve grid and
    later spline-interpolated onto a fine k-grid for the Bessel integral.
    """
    ps = _worker_ps()
    x_los = _WORKER["x_los"]
    res = ps.solver.solve(k)
    return ps._source(res, k, x_los) / ps._R_ini(k)


# --- worker for ell-parallel line-of-sight integral --------------------------
_LOS = {}


def _los_init(S_fine, x_los, k_fine, chi, z_max, dz):
    """Store the (shared, constant) LOS inputs in each worker once."""
    _LOS.update(dict(S_fine=S_fine, x_los=x_los, z_max=z_max, dz=dz,
                     arg=k_fine[:, None] * chi[None, :]))


def _los_one_ell(l):
    """Theta_ell(k_fine) for one ell: build j_ell table then integrate over x."""
    dz = _LOS["dz"]
    z = np.arange(0.0, _LOS["z_max"] + dz, dz)
    sp = CubicSpline(z, spherical_jn(int(l), z))
    jl = sp(np.clip(_LOS["arg"], 0.0, z[-1]))
    return np.trapezoid(_LOS["S_fine"] * jl, _LOS["x_los"], axis=1)


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
        Pi = res["Pi"](x)            # = Theta_2 + ThetaP_0 + ThetaP_2
        Psi = res["Psi"](x)
        Phi = res["Phi"](x)
        Psip = res["Psi"](x, 1)
        Phip = res["Phi"](x, 1)
        v_b = res["v_b"](x)

        ck = const.c * k_si

        # Term 1: g (Theta0 + Psi + Pi/4)  [polarization enters via Pi].
        T1 = g * (Th0 + Psi + 0.25 * Pi)
        # Term 2: e^{-tau} (Psi' - Phi').
        T2 = etau * (Psip - Phip)
        # Term 3: -(1/ck) d/dx (Hp g v_b).
        Hgv = Hp * g * v_b
        Hgv_sp = CubicSpline(x, Hgv)
        T3 = -(1.0 / ck) * Hgv_sp(x, 1)
        # Term 4: (3/4 c^2 k^2) d/dx [ Hp d/dx (Hp g Pi) ].
        HgT2 = Hp * g * Pi
        HgT2_sp = CubicSpline(x, HgT2)
        inner = Hp * HgT2_sp(x, 1)
        inner_sp = CubicSpline(x, inner)
        T4 = (3.0 / (4.0 * const.c ** 2 * k_si ** 2)) * inner_sp(x, 1)

        return T1 + T2 + T3 + T4

    def _source_E(self, res, k_si, x):
        """E-mode polarization source (§ appendix E / Callin).

        :math:`\\tilde S_E(k,x)=\\dfrac{3\\,\\tilde g(x)\\,\\Pi(x)}{4\\,(k\\chi)^2}`,
        with :math:`\\chi=\\eta_0-\\eta`. The geometric factor
        :math:`\\sqrt{(\\ell+2)!/(\\ell-2)!}` is applied per-ell in the transfer.
        """
        g = self.rec.g_tilde(x)
        Pi = res["Pi"](x)
        chi = self.bg.eta0 - self.bg.eta(x)
        kchi = k_si * chi
        out = np.zeros_like(x, dtype=float)
        m = kchi > 1e-8
        out[m] = 0.75 * g[m] * Pi[m] / kchi[m] ** 2
        return out

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

    def sources(self, ks_si, x_los=None, verbose=False, nproc=None):
        """Compute the line-of-sight source ``S(k,x)`` on the solve grid.

        Returns ``(ks_si, x_los, S)`` with ``S`` of shape ``(len(ks_si),
        len(x_los))``, R-normalized. This is the expensive step (one ODE solve
        per k) and is k-parallelized.
        """
        if x_los is None:
            x_los = self._default_x_los()
        ks_si = np.asarray(ks_si)
        nk = len(ks_si)
        if nproc is None:
            nproc = max(1, os.cpu_count() or 1)

        if nproc == 1:
            S = np.empty((nk, len(x_los)))
            for ik, k in enumerate(ks_si):
                res = self.solver.solve(k)
                S[ik] = self._source(res, k, x_los) / self._R_ini(k)
                if verbose and ik % 25 == 0:
                    print(f"  source k {ik+1}/{nk}")
            return ks_si, x_los, S

        init_args = (self.p, self.solver.lmax, self.solver.lmax_nu, x_los)
        with Pool(processes=nproc, initializer=_init_worker,
                  initargs=init_args) as pool:
            results = pool.map(_source_one_k, list(ks_si),
                               chunksize=max(1, nk // (nproc * 4)))
        return ks_si, x_los, np.array(results)

    def sources_TE(self, ks_si, x_los=None, verbose=False, nproc=None):
        """Temperature and E-mode sources from one solve per k.

        Returns ``(ks_si, x_los, S_T, S_E)``; both have shape
        ``(len(ks_si), len(x_los))``, R-normalized.
        """
        if x_los is None:
            x_los = self._default_x_los()
        ks_si = np.asarray(ks_si)
        nk = len(ks_si)
        if nproc is None:
            nproc = max(1, os.cpu_count() or 1)

        if nproc == 1:
            S_T = np.empty((nk, len(x_los)))
            S_E = np.empty((nk, len(x_los)))
            for ik, k in enumerate(ks_si):
                res = self.solver.solve(k)
                norm = self._R_ini(k)
                S_T[ik] = self._source(res, k, x_los) / norm
                S_E[ik] = self._source_E(res, k, x_los) / norm
            return ks_si, x_los, S_T, S_E

        init_args = (self.p, self.solver.lmax, self.solver.lmax_nu, x_los)
        with Pool(processes=nproc, initializer=_init_worker,
                  initargs=init_args) as pool:
            res = pool.map(_sources_TE_one_k, list(ks_si),
                           chunksize=max(1, nk // (nproc * 4)))
        res = np.array(res)  # (nk, 2, n_x)
        return ks_si, x_los, res[:, 0, :], res[:, 1, :]

    def los_integral(self, ells, ks_si, x_los, S, k_fine=None, nproc=None,
                     dz=2.0 * np.pi / 20.0):
        """Line-of-sight Bessel integral on a fine k-grid (§4.2).

        :math:`\\Theta_\\ell(k)=\\int S(k,x)\\,j_\\ell[k(\\eta_0-\\eta)]\\,dx`.

        The source ``S`` (smooth in k) is spline-interpolated from the coarse
        solve grid ``ks_si`` onto ``k_fine`` before the integral, because
        :math:`\\Theta_\\ell` itself oscillates rapidly in k (period
        :math:`\\sim2\\pi/\\chi_*`) and must be sampled finely for a clean C_ell.
        The per-ell work (the ``spherical_jn`` table build, which dominates at
        high ell, plus the x-integral) is parallelized over ``ells``.
        Returns ``(k_fine, Theta)`` with ``Theta`` shape ``(len(ells), len(k_fine))``.
        """
        ells = np.asarray(ells)
        eta0 = self.bg.eta0
        chi = eta0 - self.bg.eta(x_los)          # comoving distance [m]
        if k_fine is None:
            k_fine = self._fine_k_grid(ks_si)

        # Spline the (smooth-in-k) source onto the fine grid, in ln k.
        lnk = np.log(ks_si)
        S_fine = CubicSpline(lnk, S, axis=0)(np.log(k_fine))  # (n_fine, n_x)
        z_max = float(k_fine.max() * chi.max()) * 1.01

        if nproc is None:
            nproc = max(1, os.cpu_count() or 1)

        if nproc == 1 or len(ells) < 8:
            bank = BesselBank(ells, z_max, dz=dz)
            arg = k_fine[:, None] * chi[None, :]
            Theta = np.empty((len(ells), len(k_fine)))
            for il, l in enumerate(ells):
                Theta[il] = np.trapezoid(S_fine * bank.eval(l, arg), x_los, axis=1)
            return k_fine, Theta

        init = (S_fine, x_los, k_fine, chi, z_max, dz)
        with Pool(processes=nproc, initializer=_los_init, initargs=init) as pool:
            rows = pool.map(_los_one_ell, [int(l) for l in ells],
                            chunksize=max(1, len(ells) // (nproc * 4)))
        return k_fine, np.array(rows)

    def _fine_k_grid(self, ks_si):
        """Uniform fine k-grid resolving the j_l(k chi) oscillation (~2pi/chi)."""
        chi_max = self.bg.eta0 - float(self.bg.eta(self.solver.x_start))
        dk = (2.0 * np.pi / chi_max) / 8.0        # ~8 samples per oscillation
        kmin, kmax = float(ks_si.min()), float(ks_si.max())
        n_fine = int(np.ceil((kmax - kmin) / dk)) + 1
        return np.linspace(kmin, kmax, max(n_fine, len(ks_si)))

    def transfer(self, ells, ks_si=None, x_los=None, verbose=False,
                 nproc=None, k_fine=None):
        """Compute :math:`\\Theta_\\ell(k)` (R-normalized) via source-spline + LOS.

        Returns ``(k_fine, Theta)`` with ``Theta`` of shape
        ``(len(ells), len(k_fine))``. ``ks_si`` is the (coarse) solve grid;
        the C_ell-resolution fine grid is built internally (or passed).
        """
        if ks_si is None:
            ks_si = k_grid() / const.Mpc
        ells = np.asarray(ells)
        ks_si, x_los, S = self.sources(ks_si, x_los, verbose=verbose,
                                       nproc=nproc)
        return self.los_integral(ells, ks_si, x_los, S, k_fine=k_fine,
                                 nproc=nproc)

    # --- C_ell assembly (§4.3) -------------------------------------------------
    def integrate_cls(self, ells, ks_si, Theta):
        """C_ell from a transfer evaluated on a (fine) k-grid.

        C_ell = 4 pi int P_R(k) |Theta_ell|^2 dk/k = 4 pi int P_R |Theta|^2 dlnk.
        ``Theta`` here is sampled on the fine grid returned by :meth:`transfer`,
        which resolves the rapid j_l(k chi) oscillation, so a plain trapezoid is
        accurate. ``ks_si`` must be that same fine grid.
        """
        P = self.P_R(ks_si)
        lnk = np.log(ks_si)
        Cl = np.array([4.0 * np.pi * np.trapezoid(P * Theta[il] ** 2, lnk)
                       for il in range(len(ells))])
        return np.asarray(ells), Cl

    def cls(self, ells=None, ks_si=None, verbose=False, return_transfer=False):
        """Return ``(ells, C_ell)`` in dimensionless units (Theta^2).

        With ``return_transfer=True`` also returns ``(ks_si, Theta)`` so the
        transfer can be cached and re-integrated cheaply.
        """
        if ells is None:
            ells = default_ells()
        ks_si, Theta = self.transfer(ells, ks_si, verbose=verbose)
        ells, Cl = self.integrate_cls(ells, ks_si, Theta)
        if return_transfer:
            return ells, Cl, ks_si, Theta
        return ells, Cl

    def dls(self, ells=None, ks_si=None, verbose=False):
        """Return ``(ells, D_ell)`` with :math:`D_\\ell=\\ell(\\ell+1)C_\\ell T_{CMB}^2/2\\pi` [μK²]."""
        ells, Cl = self.cls(ells, ks_si, verbose=verbose)
        T2 = (self.p.T_CMB * 1e6) ** 2  # K -> uK, squared
        Dl = ells * (ells + 1.0) / (2.0 * np.pi) * Cl * T2
        return ells, Dl

    # --- temperature + polarization (TT, EE, TE) -------------------------------
    def cls_all(self, ells=None, ks_si=None, verbose=False, nproc=None):
        """Return ``(ells, {'TT','EE','TE'})`` (dimensionless) including E-mode.

        One ODE solve per k feeds both the temperature and E-mode sources; the
        E transfer carries the geometric factor :math:`\\sqrt{(\\ell+2)!/(\\ell-2)!}`.
        """
        if ells is None:
            ells = default_ells()
        ells = np.asarray(ells)
        if ks_si is None:
            ks_si = k_grid() / const.Mpc
        ks, x_los, S_T, S_E = self.sources_TE(ks_si, verbose=verbose, nproc=nproc)
        k_fine, Th_T = self.los_integral(ells, ks, x_los, S_T, nproc=nproc)
        _, Th_E = self.los_integral(ells, ks, x_los, S_E, k_fine=k_fine,
                                    nproc=nproc)
        pref = np.sqrt(np.clip((ells + 2.0) * (ells + 1.0) * ells * (ells - 1.0),
                               0.0, None))
        Th_E = pref[:, None] * Th_E

        P = self.P_R(k_fine)
        lnk = np.log(k_fine)

        def integ(A, B):
            return np.array([4.0 * np.pi * np.trapezoid(P * A[i] * B[i], lnk)
                             for i in range(len(ells))])
        return ells, {"TT": integ(Th_T, Th_T), "EE": integ(Th_E, Th_E),
                      "TE": integ(Th_T, Th_E)}

    def dls_all(self, ells=None, ks_si=None, verbose=False, nproc=None):
        """Like :meth:`cls_all` but D_ell = l(l+1)C_l T_CMB^2/2pi [μK²] (TE keeps sign)."""
        ells, C = self.cls_all(ells, ks_si, verbose=verbose, nproc=nproc)
        T2 = (self.p.T_CMB * 1e6) ** 2
        fac = ells * (ells + 1.0) / (2.0 * np.pi) * T2
        return ells, {k: fac * v for k, v in C.items()}
