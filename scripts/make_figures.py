#!/usr/bin/env python3
"""Generate the textbook figures from ``cmbcore`` output (WP4).

Each figure is written to ``figures/`` and registered in
``figures/manifest.yaml`` (figure id -> filename). This script produces the
core spectrum/background/recombination figures; the full ~60-figure set is the
WP4 deliverable and extends this module.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cmbcore import Params, BackgroundCosmology, Recombination, PowerSpectrum
from cmbcore import constants as const
from cmbcore.plotstyle import use_style
from cmbcore.io import load_cls

OUT = Path(__file__).resolve().parent.parent / "figures"
OUT.mkdir(exist_ok=True)
MANIFEST = OUT / "manifest.yaml"


def _save(fig, fid, fname, manifest):
    path = OUT / fname
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    manifest[fid] = fname
    print(f"  {fid}: {fname}")


def fig_background(p, bg, manifest):
    x = np.linspace(-15, 0, 400)
    fig, ax = plt.subplots()
    Om = bg.Omega(x)
    for key, lbl in [("gamma", r"$\gamma$"), ("nu", r"$\nu$"),
                     ("b", "b"), ("c", "c"), ("Lambda", r"$\Lambda$")]:
        ax.plot(x, Om[key], label=lbl)
    ax.axvline(np.log(1 / (1 + bg.z_eq())), ls="--", c="k", lw=1,
               label=r"$z_{\rm eq}$")
    ax.set_xlabel(r"$x=\ln a$")
    ax.set_ylabel(r"$\Omega_i(x)$")
    ax.set_title("Background density fractions")
    ax.legend(ncol=3)
    _save(fig, "F2.1", "bg_density_fractions.png", manifest)


def fig_recombination(p, bg, rec, manifest):
    x = np.linspace(-9, -5, 600)
    z = np.exp(-x) - 1
    fig, ax = plt.subplots()
    ax.semilogy(z, rec.Xe(x))
    ax.axvline(rec.z_star(), ls="--", c="r", lw=1,
               label=fr"$z_*={rec.z_star():.0f}$")
    ax.set_xlabel("redshift $z$")
    ax.set_ylabel(r"$X_e(z)$")
    ax.set_title("Free-electron fraction (Saha + Peebles)")
    ax.invert_xaxis()
    ax.legend()
    _save(fig, "F2.2", "recomb_Xe.png", manifest)

    fig, ax = plt.subplots()
    ax.plot(x, rec.g_tilde(x), label=r"$\tilde g(x)$")
    ax2 = ax.twinx()
    ax2.semilogy(x, rec.tau(x), c="C1", label=r"$\tau(x)$")
    ax.set_xlabel(r"$x=\ln a$")
    ax.set_ylabel(r"visibility $\tilde g$")
    ax2.set_ylabel(r"$\tau$")
    ax.set_title("Optical depth and visibility function")
    _save(fig, "F2.3", "recomb_visibility.png", manifest)


def fig_spectrum(manifest):
    cache = OUT / "cls_fiducial.npz"
    if not cache.exists():
        print("  (run scripts/make_values.py first for cls_fiducial.npz)")
        return
    d = load_cls(cache)
    ells, Dl = d["ells"], d["cls"]  # cache stores D_ell under the 'cls' key
    fig, ax = plt.subplots()
    ax.plot(ells, Dl)
    ax.set_xlabel(r"multipole $\ell$")
    ax.set_ylabel(r"$D_\ell=\ell(\ell+1)C_\ell\,T_{CMB}^2/2\pi\ \ [\mu K^2]$")
    ax.set_title(r"CMB temperature power spectrum $D_\ell^{TT}$")
    ax.set_xlim(2, ells.max())
    _save(fig, "F11.1", "cl_tt_fiducial.png", manifest)


def fig_perturbations(p, bg, rec, manifest):
    """F3.x: single-mode evolution of Theta_0, delta_c, Phi (NB3)."""
    from cmbcore.perturbations import PerturbationSolver
    solver = PerturbationSolver(bg, rec, p)
    H0_c = p.H0 / const.c
    x = np.linspace(-12, 0, 500)
    a = np.exp(x)
    res = {n: solver.solve(n * H0_c) for n in (10, 100, 1000)}

    fig, ax = plt.subplots()
    r = res[1000]
    for l in (0, 1, 2):
        ax.plot(a, r[f"Theta{l}"](x), label=fr"$\Theta_{l}$")
    ax.set_xscale("log")
    ax.set_xlabel("a"); ax.set_ylabel(r"$\Theta_\ell$")
    ax.set_title(r"Photon multipoles, $k=1000\,H_0/c$ (acoustic oscillations)")
    ax.legend()
    _save(fig, "F3.6", "pert_theta_multipoles.png", manifest)

    fig, ax = plt.subplots()
    for n in (10, 100, 1000):
        ax.loglog(a, np.abs(res[n]["delta_c"](x)), label=fr"$k={n}H_0/c$")
    ax.set_xlabel("a"); ax.set_ylabel(r"$|\delta_c|$")
    ax.set_title("CDM growth: scale-dependent onset")
    ax.legend()
    _save(fig, "F3.5", "pert_delta_c.png", manifest)


def fig_param_study(manifest):
    """F6.4/F6.6: tau (reionization) and z_* (recomb shift) effects (NB6)."""
    import json
    cache = OUT / "param_study.npz"
    if not cache.exists():
        print("  (run scripts/make_param_study.py first for param_study.npz)")
        return
    d = np.load(cache)
    ells = d["ells"]
    S = json.loads((OUT / "param_study_summary.json").read_text())

    fig, ax = plt.subplots()
    ax.plot(ells, d["fid"], label="no reionization")
    ax.plot(ells, d["reio_z8"], label=fr"reion $z_{{re}}=8$ ($\tau\approx{S['tau_reio_z8']:.3f}$)")
    ax.set_xlabel(r"$\ell$"); ax.set_ylabel(r"$D_\ell\ [\mu K^2]$")
    ax.set_title(r"Optical depth: $e^{-2\tau}$ suppression at $\ell\gtrsim30$")
    ax.legend()
    _save(fig, "F6.4", "param_tau.png", manifest)

    fig, ax = plt.subplots()
    ax.plot(ells, d["fid"], label=fr"fiducial ($z_*={S['z_star_fid']:.0f}$)")
    ax.plot(ells, d["recomb_p5"], label=fr"$T_b\times1.05$ ($z_*={S['z_star_shift_p5']:.0f}$)")
    ax.set_xlim(2, 800)
    ax.set_xlabel(r"$\ell$"); ax.set_ylabel(r"$D_\ell\ [\mu K^2]$")
    ax.set_title(r"Recombination shift: $\ell_1$ moves with $z_*$")
    ax.legend()
    _save(fig, "F6.6", "param_zstar.png", manifest)


def fig_polarization(manifest):
    """F-pol: TT/EE/TE from the cached polarization spectra (NB-pol)."""
    cache = OUT / "cls_polarization.npz"
    if not cache.exists():
        print("  (run scripts/make_polarization.py first for cls_polarization.npz)")
        return
    d = np.load(cache)
    ells = d["ells"]
    fig, ax = plt.subplots(3, 1, figsize=(7, 8), sharex=True)
    for a, key, lab in [(ax[0], "TT", r"$D_\ell^{TT}$"),
                        (ax[1], "EE", r"$D_\ell^{EE}$"),
                        (ax[2], "TE", r"$D_\ell^{TE}$")]:
        a.plot(ells, d[key])
        a.set_ylabel(lab + r"$\ [\mu K^2]$")
    ax[2].set_xlabel(r"$\ell$")
    ax[0].set_title("Temperature and E-mode polarization (cmbcore)")
    _save(fig, "F-pol", "polarization_TTEETE.png", manifest)


def fig_massive_nu(manifest):
    """F6.7: massive-neutrino TT response vs CLASS (ch.13.8)."""
    cache = OUT / "massive_nu.npz"
    if not cache.exists():
        print("  (run scripts/make_massive_nu.py first for massive_nu.npz)")
        return
    d = np.load(cache)
    ells = d["ells"]
    keys = [k for k in d.files if k.startswith("core_")]
    Sm = sorted(float(k.split("_")[1]) for k in keys)[-1]
    fig, ax = plt.subplots()
    ax.plot(ells, d[f"core_{Sm}"] / d["core_0.0"], label="cmbcore")
    if f"class_{Sm}" in d.files:
        ax.plot(ells, d[f"class_{Sm}"] / d["class_0.0"], "--", label="CLASS")
    ax.axhline(1.0, c="k", lw=0.6)
    ax.set_xlabel(r"$\ell$")
    ax.set_ylabel(rf"$C_\ell(\Sigma m_\nu={Sm})/C_\ell(0)$")
    ax.set_title(r"Massive-neutrino TT response ($h$ fixed)")
    ax.legend()
    _save(fig, "F6.7", "massive_nu_response.png", manifest)


def main():
    use_style()
    p = Params.fiducial()
    bg = BackgroundCosmology(p)
    rec = Recombination(bg, p)
    manifest = {}
    print("Generating figures:")
    fig_background(p, bg, manifest)
    fig_recombination(p, bg, rec, manifest)
    fig_perturbations(p, bg, rec, manifest)
    fig_spectrum(manifest)
    fig_param_study(manifest)
    fig_polarization(manifest)
    fig_massive_nu(manifest)

    lines = ["# figure id -> filename (auto-generated by make_figures.py)"]
    for fid in sorted(manifest):
        lines.append(f"{fid}: {manifest[fid]}")
    MANIFEST.write_text("\n".join(lines) + "\n")
    print(f"Wrote {MANIFEST}")


if __name__ == "__main__":
    main()
