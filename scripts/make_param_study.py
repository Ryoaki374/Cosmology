#!/usr/bin/env python3
"""Parameter-dependence study (README §1.4 required topics; NB6 / ch.13).

Computes D_ell for the fiducial cosmology and for the reionization-optical-depth
(tau) and recombination-shift (z_*) variations, and caches them so NB6 and the
figures can be regenerated without recomputing. Massive-neutrino (Sigma m_nu) is
CLASS-deferred in this version (04 §7) and is documented, not computed here.
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


def dl_for(params, nk=150):
    ps = PowerSpectrum(params)
    ks = k_grid(nk=nk) / const.Mpc
    ells, Dl = ps.dls(ks_si=ks)
    return np.asarray(ells), np.asarray(Dl), ps


def main():
    ells = None
    cases = {}
    summary = {}

    # Reionization optical depth: off, and tanh reionization at z_re=8.
    print("fiducial (no reionization) ...")
    e, Dl0, ps0 = dl_for(Params.fiducial())
    ells = e
    cases["fid"] = Dl0
    summary["tau_reio_fid"] = float(ps0.rec.tau(ps0.rec.x[0]) * 0 + 0.0)

    print("reionization z_re=8 ...")
    p_re = Params(tau_reio=1.0, z_reio=8.0)   # tau_reio!=0 enables the hook
    e, Dl_re, ps_re = dl_for(p_re)
    cases["reio_z8"] = Dl_re
    # Reionization optical depth = tau accumulated from now to just above the
    # reionized epoch (z~30); evaluating deeper would pick up the pre-
    # recombination plasma (tau ~ 1e5) and is not the reionization tau.
    x30 = float(np.log(1.0 / 31.0))
    tau_re = float(ps_re.rec.tau(x30))
    summary["tau_reio_z8"] = tau_re

    print("recombination shift +5% ...")
    e, Dl_rs, ps_rs = dl_for(Params(recomb_shift=0.05))
    cases["recomb_p5"] = Dl_rs
    summary["z_star_fid"] = float(ps0.rec.z_star())
    summary["z_star_shift_p5"] = float(ps_rs.rec.z_star())

    np.savez(OUT / "param_study.npz", ells=ells, **cases)

    # Quantitative effects (the "which direction, how many %, why" requirement).
    def pct(a, b):  # percent change of a vs b at l=500
        i = np.argmin(np.abs(ells - 500))
        return 100.0 * (a[i] - b[i]) / b[i]
    summary["reio_suppression_l500_pct"] = pct(Dl_re, Dl0)
    summary["reio_predicted_exp_-2tau_pct"] = 100.0 * (np.exp(-2 * tau_re) - 1)
    i1 = lambda D: int(ells[np.argmax(np.where(ells >= 150, D, 0))])
    summary["l1_fid"] = i1(Dl0)
    summary["l1_recomb_p5"] = i1(Dl_rs)

    (OUT / "param_study_summary.json").write_text(
        json.dumps({k: round(v, 5) if isinstance(v, float) else v
                    for k, v in summary.items()}, indent=2))
    print("Wrote param_study.npz and param_study_summary.json")
    for k, v in summary.items():
        print(f"  {k} = {v}")


if __name__ == "__main__":
    main()
