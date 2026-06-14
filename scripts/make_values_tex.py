#!/usr/bin/env python3
"""Generate textbook/values.tex (\\newcommand macros) from figures/values.json.

Lets the prose cite numbers as macros (e.g. \\valzeq) instead of bare literals,
so every number traces to the computed values.json (WP5/WP6, 04 §5).
"""
from __future__ import annotations
import json, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
vals = json.loads((ROOT / "figures" / "values.json").read_text())

# json key -> (LaTeX macro name, format string)
MACROS = {
    "z_eq": ("valzeq", "{:.0f}"),
    "z_star": ("valzstar", "{:.0f}"),
    "r_s_star_Mpc": ("valrs", "{:.1f}"),
    "chi_star_Gpc": ("valchistar", "{:.2f}"),
    "theta_star_100": ("valthetastar", "{:.3f}"),
    "k_D_star_invMpc": ("valkD", "{:.3f}"),
    "age_Gyr": ("valage", "{:.2f}"),
    "eta0_Mpc": ("valetazero", "{:.0f}"),
    "l1_peak": ("vallone", "{:.0f}"),
    "l2_peak": ("valltwo", "{:.0f}"),
    "l_acoustic": ("vallA", "{:.0f}"),
}

lines = ["% AUTO-GENERATED from figures/values.json by make_values_tex.py.",
         "% Do not edit by hand; run `python scripts/make_values_tex.py`."]
for key, (macro, fmt) in MACROS.items():
    if key in vals:
        lines.append(f"\\newcommand{{\\{macro}}}{{{fmt.format(vals[key])}}}")
out = ROOT / "textbook" / "values.tex"
out.write_text("\n".join(lines) + "\n")
print("wrote", out)
for ln in lines[2:]:
    print(" ", ln)
