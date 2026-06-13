"""Result caching and ``values.json`` export (``03_numerics_spec.md`` §6).

``values.json`` is the single source of numbers that the textbook (WP5) splices
in — no bare numeric literals are allowed in the prose (``04`` §5).
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


def save_cls(path, ells, cls, dls=None) -> None:
    """Cache a computed spectrum to an ``.npz`` file."""
    data = {"ells": np.asarray(ells), "cls": np.asarray(cls)}
    if dls is not None:
        data["dls"] = np.asarray(dls)
    np.savez(path, **data)


def load_cls(path):
    """Load a cached spectrum; returns a dict of arrays."""
    d = np.load(path)
    return {k: d[k] for k in d.files}


def dump_values(path, values: dict) -> None:
    """Write the ``values.json`` registry (sorted keys, 6 sig figs)."""
    def _round(v):
        if isinstance(v, float):
            return float(f"{v:.6g}")
        return v
    clean = {k: _round(v) for k, v in values.items()}
    Path(path).write_text(json.dumps(clean, indent=2, sort_keys=True,
                                     ensure_ascii=False))


def load_values(path) -> dict:
    return json.loads(Path(path).read_text())
