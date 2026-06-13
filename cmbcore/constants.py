"""Physical constants in SI units.

All of ``cmbcore`` works in SI with the speed of light ``c`` written
explicitly, as mandated by ``03_numerics_spec.md`` §0. Cosmological inputs are
given as dimensionless density parameters and ``h``; the conversions to SI live
in :mod:`cmbcore.params` and :mod:`cmbcore.background`.
"""

from __future__ import annotations

import math

# Fundamental constants (CODATA-ish, sufficient for sub-percent cosmology).
c = 2.99792458e8            # speed of light            [m/s]
G = 6.67430e-11             # gravitational constant     [m^3 kg^-1 s^-2]
hbar = 1.054571817e-34      # reduced Planck constant    [J s]
h_planck = 2 * math.pi * hbar
k_B = 1.380649e-23          # Boltzmann constant         [J/K]

m_e = 9.1093837015e-31      # electron mass              [kg]
m_H = 1.67262192369e-27     # hydrogen (proton) mass     [kg]
sigma_T = 6.6524587321e-29  # Thomson cross section       [m^2]

# Hydrogen binding energy and 2s->1s two-photon rate.
eps_0_eV = 13.605693        # hydrogen ground-state binding energy [eV]
eV = 1.602176634e-19        # electron volt              [J]
eps_0 = eps_0_eV * eV       # binding energy             [J]
Lambda_2s1s = 8.227         # 2s -> 1s two-photon decay rate [1/s]

# Unit helpers.
Mpc = 3.085677581491367e22  # megaparsec                 [m]
Gpc = 1.0e3 * Mpc
year = 365.25 * 24 * 3600.0
Gyr = 1.0e9 * year

# Radiation: Stefan-Boltzmann a_rad = pi^2 k_B^4 / (15 hbar^3 c^3).
a_rad = (math.pi ** 2) * (k_B ** 4) / (15.0 * (hbar ** 3) * (c ** 3))
