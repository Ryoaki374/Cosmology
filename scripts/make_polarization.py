#!/usr/bin/env python3
"""CMB polarization spectra TT/EE/TE (cmbcore) with CLASS validation (NB-pol).

Computes the temperature and E-mode spectra from cmbcore's polarized Boltzmann
hierarchy and, if ``classy`` is available, overlays the matched CLASS setting A
(TT+pol, no lensing, reionization off, massless nu, YHe minimized). Caches the
spectra and writes a 3-panel figure (TT, EE, TE).
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
OUT.mkdir(exist_ok=True)
LMAX = 1500


def class_spectra(p, ells):
    from classy import Class
    base = dict(output="tCl,pCl", lensing="no", modes="s",
                h=p.h, omega_b=p.Omega_b * p.h ** 2,
                omega_cdm=p.Omega_c * p.h ** 2, T_cmb=p.T_CMB,
                n_s=p.n_s, A_s=p.A_s, k_pivot=p.k_pivot,
                N_ur=p.N_eff, N_ncdm=0,
                reio_parametrization="reio_none", l_max_scalars=LMAX + 200)
    for yhe in (0.01, 0.24):
        try:
            M = Class(); M.set({**base, "YHe": yhe}); M.compute(); break
        except Exception:  # noqa: BLE001
            continue
    cl = M.raw_cl(LMAX); l = cl["ell"]
    fac = l * (l + 1) / (2 * np.pi) * (p.T_CMB * 1e6) ** 2
    return {k: np.interp(ells, l, fac * cl[k]) for k in ("tt", "ee", "te")}


def main(nk=200):
    p = Params.fiducial()
    ps = PowerSpectrum(p)
    ells = default_ells()
    print(f"cmbcore TT/EE/TE (nk={nk}) ...")
    ells, D = ps.dls_all(ks_si=k_grid(nk=nk) / const.Mpc)
    np.savez(OUT / "cls_polarization.npz", ells=ells,
             TT=D["TT"], EE=D["EE"], TE=D["TE"])

    cls_c = None
    try:
        print("CLASS TT/EE/TE ...")
        cls_c = class_spectra(p, ells)
    except Exception as e:  # noqa: BLE001
        print("classy unavailable:", e)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cmbcore.plotstyle import use_style
    use_style()
    fig, ax = plt.subplots(3, 1, figsize=(7, 9))
    for a, key, ck, ttl in [(ax[0], "TT", "tt", r"$D_\ell^{TT}$"),
                            (ax[1], "EE", "ee", r"$D_\ell^{EE}$"),
                            (ax[2], "TE", "te", r"$D_\ell^{TE}$")]:
        a.plot(ells, D[key], label="cmbcore")
        if cls_c is not None:
            a.plot(ells, cls_c[ck], "--", label="CLASS (setting A)")
        a.set_ylabel(ttl + r"$\ [\mu K^2]$"); a.legend()
    ax[2].set_xlabel(r"$\ell$")
    ax[0].set_title("CMB temperature & E-mode polarization: cmbcore vs CLASS")
    fig.savefig(OUT / "polarization.png", dpi=150, bbox_inches="tight")

    summary = {}
    if cls_c is not None:
        m = (ells >= 30) & (ells <= 1000)
        for key, ck in [("EE", "ee"), ("TE", "te")]:
            denom = np.max(np.abs(cls_c[ck][m]))
            summary[f"{key}_maxdiff_l30-1000_pct_of_peak"] = round(
                float(np.max(np.abs(D[key][m] - cls_c[ck][m])) / denom * 100), 2)
        (OUT / "polarization_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("Wrote cls_polarization.npz, polarization.png")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 200)
