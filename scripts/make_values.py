#!/usr/bin/env python3
"""Generate ``figures/values.json`` — the numeric registry for the textbook.

Every number the prose cites must come from this file (``04`` §5). Run after
any change to ``cmbcore``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cmbcore import Params, BackgroundCosmology, Recombination, PowerSpectrum
from cmbcore import constants as const, analytic
from cmbcore.io import dump_values, save_cls

OUT = Path(__file__).resolve().parent.parent / "figures"
OUT.mkdir(exist_ok=True)


def main(nk: int = 200) -> None:
    p = Params.fiducial()
    bg = BackgroundCosmology(p)
    rec = Recombination(bg, p)

    x_star = rec.x_star()
    values = {
        "h": p.h,
        "Omega_b_h2": p.Omega_b * p.h ** 2,
        "Omega_c_h2": p.Omega_c * p.h ** 2,
        "Omega_gamma": p.Omega_gamma,
        "Omega_nu": p.Omega_nu,
        "Omega_Lambda": p.Omega_Lambda,
        "T_CMB_K": p.T_CMB,
        "z_eq": bg.z_eq(),
        "age_Gyr": bg.t0 / const.Gyr,
        "eta0_Mpc": bg.eta0 / const.Mpc,
        "z_star": rec.z_star(),
        "r_s_star_Mpc": float(rec.r_s(x_star)) / const.Mpc,
        "chi_star_Gpc": float(bg.comoving_distance(x_star)) / const.Gpc,
        "theta_star_100": analytic.theta_star(bg, rec),
        "visibility_norm": rec.visibility_norm(),
        "l_acoustic": analytic.acoustic_scale(bg, rec),
        "k_D_star_invMpc": analytic.silk_scale(bg, rec) * const.Mpc,
    }

    # Full spectrum for peak diagnostics (quadratic k-grid, §4.3).
    print(f"Computing spectrum (nk={nk}) ...")
    from cmbcore.spectrum import k_grid
    ps = PowerSpectrum(p)
    ks = k_grid(nk=nk) / const.Mpc
    ells, Dl = ps.dls(ks_si=ks, verbose=True)
    save_cls(OUT / "cls_fiducial.npz", ells, Dl)

    # Acoustic-peak finding via linear interpolation (avoids cubic-spline
    # overshoot between sparsely sampled high-l points), restricted to l >= 120
    # with a prominence cut so the SW/ISW rise is not mistaken for a peak.
    from scipy.signal import find_peaks
    lfine = np.arange(120, int(ells[-1]))
    Dfine = np.interp(lfine, ells, Dl)
    idx, _ = find_peaks(Dfine, prominence=0.10 * Dfine.max(), distance=120)
    peaks = lfine[idx]
    sp = lambda L: float(np.interp(L, ells, Dl))
    if len(peaks) >= 1:
        values["l1_peak"] = int(peaks[0])
        values["Dl_peak1_uK2"] = float(sp(peaks[0]))
    if len(peaks) >= 2:
        values["l2_peak"] = int(peaks[1])
        values["l2_over_l1"] = float(peaks[1] / peaks[0])
        values["Dl_peak2_uK2"] = float(sp(peaks[1]))
    if len(peaks) >= 3:
        values["l3_peak"] = int(peaks[2])
        values["Dl_peak3_uK2"] = float(sp(peaks[2]))

    dump_values(OUT / "values.json", values)
    print(f"Wrote {OUT / 'values.json'}")
    for k, v in sorted(values.items()):
        print(f"  {k:24s} = {v}")


if __name__ == "__main__":
    nk = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    main(nk)
