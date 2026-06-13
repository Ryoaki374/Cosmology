"""Shared matplotlib styling for figures (``03_numerics_spec.md`` §6)."""

from __future__ import annotations


def use_style() -> None:
    """Apply the project's default plotting style."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.figsize": (7.0, 4.5),
        "figure.dpi": 110,
        "savefig.dpi": 150,
        "font.size": 12,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "axes.labelsize": 13,
        "legend.fontsize": 10,
        "lines.linewidth": 1.8,
    })


# Consistent colours for the four source-function terms (NB4/NB5 decomposition).
SOURCE_COLORS = {
    "sw": "#1f77b4",       # Sachs-Wolfe  g(Theta0+Psi)
    "doppler": "#ff7f0e",  # Doppler      d/dx(Hp g v_b)
    "isw": "#2ca02c",      # ISW          e^{-tau}(Psi'-Phi')
    "pol": "#d62728",      # polarization Theta2 correction
}
