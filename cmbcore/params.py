"""Cosmological parameter container.

The :class:`Params` dataclass holds the *input* cosmology. Derived radiation
densities (:math:`\\Omega_{\\gamma 0}, \\Omega_{\\nu 0}`) and :math:`\\Omega_{\\Lambda 0}`
are computed here so that every downstream module sees a single, consistent set
of density parameters (``00_README.md`` §3.2).
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math

from . import constants as const


@dataclass
class Params:
    """Input cosmological parameters (SI-friendly, dimensionless densities).

    Density parameters are *today's* values (subscript 0). Radiation densities
    and :math:`\\Omega_{\\Lambda 0}` are filled in by :meth:`__post_init__` for a
    flat universe unless ``Omega_k`` is set.
    """

    h: float = 0.674
    Omega_b: float = 0.0224 / 0.674 ** 2     # Omega_b   = (Omega_b h^2)/h^2
    Omega_c: float = 0.1200 / 0.674 ** 2     # Omega_c   = (Omega_c h^2)/h^2
    Omega_k: float = 0.0
    T_CMB: float = 2.7255                     # [K]
    N_eff: float = 3.046
    Sigma_mnu: float = 0.0                    # [eV]; minimal run uses 0
    n_s: float = 0.965
    A_s: float = 2.1e-9
    k_pivot: float = 0.05                     # [1/Mpc]
    tau_reio: float = 0.0                     # 0 disables reionization
    Yp: float = 0.0                           # helium fraction (default off)

    # Optional toy hooks for chapter 13 / NB6.
    z_reio: float = 7.7
    delta_z_reio: float = 0.5
    recomb_shift: float = 0.0                 # fractional T_b shift

    # Derived (filled in post-init).
    Omega_gamma: float = field(init=False, default=0.0)
    Omega_nu: float = field(init=False, default=0.0)
    Omega_Lambda: float = field(init=False, default=0.0)
    H0: float = field(init=False, default=0.0)

    def __post_init__(self) -> None:
        # H0 in SI [1/s]: 100 h km/s/Mpc.
        self.H0 = 100.0 * self.h * 1.0e3 / const.Mpc

        # Photon density today: Omega_gamma0 = (a_rad T^4 / c^2) / rho_crit,
        # with rho_crit = 3 H0^2 / (8 pi G).
        rho_crit = 3.0 * self.H0 ** 2 / (8.0 * math.pi * const.G)
        rho_gamma = const.a_rad * self.T_CMB ** 4 / const.c ** 2
        self.Omega_gamma = rho_gamma / rho_crit

        # Massless neutrinos: 7/8 (4/11)^{4/3} per effective species.
        self.Omega_nu = (
            self.N_eff * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0)
            * self.Omega_gamma
        )

        # Flat closure (or honour Omega_k if supplied).
        self.Omega_Lambda = (
            1.0 - self.Omega_k
            - (self.Omega_b + self.Omega_c + self.Omega_gamma + self.Omega_nu)
        )

    # --- convenience accessors -------------------------------------------------
    @property
    def Omega_m(self) -> float:
        """Total non-relativistic matter density today."""
        return self.Omega_b + self.Omega_c

    @property
    def Omega_r(self) -> float:
        """Total radiation density today (photons + massless neutrinos)."""
        return self.Omega_gamma + self.Omega_nu

    @property
    def f_nu(self) -> float:
        """Neutrino fraction of radiation, :math:`\\Omega_\\nu/(\\Omega_\\gamma+\\Omega_\\nu)`."""
        return self.Omega_nu / (self.Omega_gamma + self.Omega_nu)

    # --- factory presets -------------------------------------------------------
    @classmethod
    def fiducial(cls) -> "Params":
        """Planck-2018-near fiducial cosmology (``00_README.md`` §3.2).

        Minimal-construction defaults: no helium, massless neutrinos, no
        reionization, no polarization.
        """
        return cls()

    @classmethod
    def callin2006(cls) -> "Params":
        """Old-WMAP preset used to reproduce the attached slides (``04`` §3.4).

        ``Omega_b=0.046, Omega_m=0.224, h=0.7, massless nu, n_s=1``.
        """
        h = 0.7
        Omega_b = 0.046
        Omega_m = 0.224
        return cls(
            h=h,
            Omega_b=Omega_b,
            Omega_c=Omega_m - Omega_b,
            n_s=1.0,
            A_s=2.1e-9,
        )
