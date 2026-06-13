#!/usr/bin/env python3
"""CLASS reference comparison (03_numerics_spec.md §5.3, NB8 / ch.14).

Runs CLASS with settings matched as closely as possible to cmbcore's minimal
construction (setting A: temperature only, no lensing, reionization off, massless
neutrinos, helium minimized) and reports the relative difference of D_ell^TT.
Caches the CLASS spectrum to figures/class_tt.csv so the comparison can be
reproduced without classy.

Tolerances (§5.3): rel.diff < 3% for l<=1000, < 5% for l<=1500.
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


def run_class(p: Params):
    from classy import Class
    base = {
        "output": "tCl", "lensing": "no", "modes": "s",
        "h": p.h, "omega_b": p.Omega_b * p.h ** 2,
        "omega_cdm": p.Omega_c * p.h ** 2, "T_cmb": p.T_CMB,
        "n_s": p.n_s, "A_s": p.A_s, "k_pivot": p.k_pivot,
        "N_ur": p.N_eff, "N_ncdm": 0,
        "reio_parametrization": "reio_none",
        "l_max_scalars": LMAX + 200,
    }
    # cmbcore uses Y_p=0; CLASS allows 0.01<=YHe<=0.5, so use the minimum (0.01)
    # to match as closely as possible; fall back to the standard BBN value.
    for yhe in (0.01, 0.24):
        try:
            M = Class()
            M.set({**base, "YHe": yhe})
            M.compute()
            return M, yhe
        except Exception as e:  # noqa: BLE001
            last = e
    raise RuntimeError(f"CLASS failed for all YHe values: {last}")


def class_dl(M, ells, T_cmb):
    cl = M.raw_cl(LMAX)
    l = cl["ell"]
    fac = l * (l + 1) / (2 * np.pi) * (T_cmb * 1e6) ** 2
    Dl = fac * cl["tt"]
    return np.interp(ells, l, Dl)


def main():
    p = Params.fiducial()
    cache = OUT / "cls_fiducial.npz"
    if cache.exists():
        d = np.load(cache)
        ells, Dl_core = d["ells"], d["cls"]
    else:
        ps = PowerSpectrum(p)
        ells, Dl_core = ps.dls(ks_si=k_grid(nk=200) / const.Mpc)

    M, yhe = run_class(p)
    Dl_class = class_dl(M, ells, p.T_CMB)
    np.savetxt(OUT / "class_tt.csv",
               np.column_stack([ells, Dl_class]), delimiter=",",
               header=f"l,Dl_class_uK2 (setting A, YHe={yhe})")

    m = ells >= 2
    rel = (Dl_core[m] - Dl_class[m]) / Dl_class[m]
    e = ells[m]
    def maxrel(lmax):
        s = e <= lmax
        return float(np.max(np.abs(rel[s])) * 100)
    summary = {
        "YHe_class": yhe,
        "max_reldiff_l<=1000_pct": round(maxrel(1000), 2),
        "max_reldiff_l<=1500_pct": round(maxrel(1500), 2),
        "median_reldiff_pct": round(float(np.median(np.abs(rel)) * 100), 2),
        "l1_core": int(e[np.argmax(np.where(e >= 150, Dl_core[m], 0))]),
    }
    (OUT / "class_comparison_summary.json").write_text(json.dumps(summary, indent=2))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from cmbcore.plotstyle import use_style
    use_style()
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7, 6.5), sharex=True,
                                 gridspec_kw={"height_ratios": [2, 1]})
    a1.plot(ells, Dl_core, label="cmbcore")
    a1.plot(ells, Dl_class, "--", label=f"CLASS (setting A, YHe={yhe:g})")
    a1.set_ylabel(r"$D_\ell\ [\mu K^2]$"); a1.legend(); a1.set_title("TT: cmbcore vs CLASS")
    a2.plot(e, rel * 100)
    a2.axhspan(-3, 3, alpha=0.15, color="green")
    a2.axhspan(-5, 5, alpha=0.07, color="orange")
    a2.set_xlabel(r"$\ell$"); a2.set_ylabel("rel. diff [%]"); a2.set_ylim(-30, 30)
    fig.savefig(OUT / "class_comparison.png", dpi=150, bbox_inches="tight")
    print(json.dumps(summary, indent=2))
    print("Wrote class_tt.csv, class_comparison.png, class_comparison_summary.json")


if __name__ == "__main__":
    main()
