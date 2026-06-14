#!/usr/bin/env python3
"""Massive-neutrino TT response vs CLASS (Sigma m_nu; ch.13.8, NB6).

Validates the native massive-neutrino *background* implementation by comparing,
against CLASS (full ncdm), the TT spectrum and especially the response ratio
C_l(Sigma m_nu)/C_l(massless) — the latter isolates the neutrino-mass effect
from the absolute calibration. h is held fixed (Omega_Lambda absorbs the change).
"""
from __future__ import annotations
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from cmbcore import Params, PowerSpectrum
from cmbcore.spectrum import k_grid, default_ells
from cmbcore import constants as const

OUT = Path(__file__).resolve().parent.parent / "figures"
MASSES = [0.0, 0.12]
NK = 150
LMAX = 1500


def cmbcore_TT(Sigma, ells):
    ps = PowerSpectrum(Params(Sigma_mnu=Sigma))
    _, Dl = ps.dls(ells=ells, ks_si=k_grid(nk=NK) / const.Mpc)
    return Dl


def class_TT(Sigma, ells):
    from classy import Class
    p = Params(Sigma_mnu=Sigma)
    base = dict(output="tCl", lensing="no", modes="s",
                h=p.h, omega_b=p.Omega_b * p.h ** 2,
                omega_cdm=p.Omega_c * p.h ** 2, T_cmb=p.T_CMB,
                n_s=p.n_s, A_s=p.A_s, k_pivot=p.k_pivot,
                reio_parametrization="reio_none", l_max_scalars=LMAX + 200,
                YHe=0.01)
    if Sigma > 0:
        base.update(N_ur=p.N_ur, N_ncdm=1, m_ncdm=Sigma, T_ncdm=0.71611)
    else:
        base.update(N_ur=p.N_eff, N_ncdm=0)
    M = Class(); M.set(base); M.compute()
    cl = M.raw_cl(LMAX); l = cl["ell"]
    fac = l * (l + 1) / (2 * np.pi) * (p.T_CMB * 1e6) ** 2
    return np.interp(ells, l, fac * cl["tt"])


def main():
    ells = default_ells()
    core = {S: cmbcore_TT(S, ells) for S in MASSES}
    print("cmbcore done")
    cls_ = {}
    try:
        cls_ = {S: class_TT(S, ells) for S in MASSES}
    except Exception as e:  # noqa: BLE001
        print("classy unavailable:", e)

    np.savez(OUT / "massive_nu.npz", ells=ells,
             **{f"core_{S}": core[S] for S in MASSES},
             **{f"class_{S}": cls_[S] for S in cls_})

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cmbcore.plotstyle import use_style
    use_style()
    Sm = MASSES[-1]
    fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)
    ax[0].plot(ells, core[0.0], label="cmbcore massless")
    ax[0].plot(ells, core[Sm], label=f"cmbcore Sm={Sm} eV")
    if cls_:
        ax[0].plot(ells, cls_[Sm], "--", label=f"CLASS Sm={Sm} eV")
    ax[0].set_ylabel(r"$D_\ell^{TT}\ [\mu K^2]$"); ax[0].legend()
    ax[0].set_title(r"Massive-neutrino effect on TT ($h$ fixed)")
    ax[1].plot(ells, core[Sm] / core[0.0], label="cmbcore ratio")
    if cls_:
        ax[1].plot(ells, cls_[Sm] / cls_[0.0], "--", label="CLASS ratio")
    ax[1].axhline(1.0, c="k", lw=0.6)
    ax[1].set_xlabel(r"$\ell$"); ax[1].set_ylabel(rf"$C_\ell(Sm={Sm})/C_\ell(0)$")
    ax[1].legend()
    fig.savefig(OUT / "massive_nu.png", dpi=150, bbox_inches="tight")

    summary = {}
    if cls_:
        m = (ells >= 2) & (ells <= 1500)
        r_core = core[Sm][m] / core[0.0][m]
        r_class = cls_[Sm][m] / cls_[0.0][m]
        summary["Sigma_mnu"] = Sm
        summary["ratio_maxdiff_pct"] = round(float(np.max(np.abs(r_core - r_class)) * 100), 2)
        summary["TT_abs_median_diff_pct"] = round(
            float(np.median(np.abs(core[Sm][m] - cls_[Sm][m]) / cls_[Sm][m]) * 100), 2)
        (OUT / "massive_nu_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("Wrote massive_nu.npz, massive_nu.png")


if __name__ == "__main__":
    main()
