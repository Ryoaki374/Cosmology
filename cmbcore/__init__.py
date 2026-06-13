"""``cmbcore`` — a minimal-construction CMB temperature power-spectrum engine.

Pipeline (``03_numerics_spec.md``): background -> recombination ->
perturbations -> spectrum. See ``00_README.md`` for the project charter.
"""

from .params import Params
from .background import BackgroundCosmology
from .recombination import Recombination
from .perturbations import PerturbationSolver
from .spectrum import PowerSpectrum, default_ells, k_grid

__all__ = [
    "Params",
    "BackgroundCosmology",
    "Recombination",
    "PerturbationSolver",
    "PowerSpectrum",
    "default_ells",
    "k_grid",
]

__version__ = "0.1.0"
