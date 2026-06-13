#!/usr/bin/env python3
"""Generate the starter notebooks (WP2 skeleton) with nbformat.

Produces NB0 (environment check) and NB1 (background cosmology) as a concrete,
executable demonstration of the ``cmbcore`` API. The full NB0-NB8 set is the
WP2 deliverable and extends this generator (``02_notebook_spec.md``).
"""

from __future__ import annotations

from pathlib import Path

import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT = Path(__file__).resolve().parent.parent / "notebooks"
OUT.mkdir(exist_ok=True)


def nb0():
    nb = new_notebook()
    nb.cells = [
        new_markdown_cell(
            "# NB0 — 環境チェックと `cmbcore` 入門\n\n"
            "このノートブックは依存関係を確認し、`cmbcore` の基本オブジェクト"
            "（`Params`, `BackgroundCosmology`）を生成する。"),
        new_code_cell(
            "import numpy as np, matplotlib.pyplot as plt\n"
            "from cmbcore import Params, BackgroundCosmology\n"
            "p = Params.fiducial()\n"
            "print('H0 [1/s] =', p.H0)\n"
            "print('Omega_gamma =', p.Omega_gamma, ' f_nu =', p.f_nu)"),
        new_code_cell(
            "bg = BackgroundCosmology(p)\n"
            "print('z_eq = %.0f' % bg.z_eq())\n"
            "print('age  = %.2f Gyr' % (bg.t0/3.1557e16))"),
    ]
    return nb


def nb1():
    nb = new_notebook()
    nb.cells = [
        new_markdown_cell(
            "# NB1 — 背景宇宙\n\n"
            "本文第1章に対応。$\\mathcal H(x)$、密度分率 $\\Omega_i(x)$、"
            "共形時間 $\\eta(x)$ を計算し、等密度 $z_{\\rm eq}$ を読み取る。\n\n"
            "**説明 → 実行 → 読み取り問い → 課題** の構造（`00_README` §5-E）。"),
        new_code_cell(
            "import numpy as np, matplotlib.pyplot as plt\n"
            "from cmbcore import Params, BackgroundCosmology\n"
            "from cmbcore.plotstyle import use_style; use_style()\n"
            "p = Params.fiducial(); bg = BackgroundCosmology(p)"),
        new_markdown_cell("## 密度分率の進化"),
        new_code_cell(
            "x = np.linspace(-15, 0, 400)\n"
            "Om = bg.Omega(x)\n"
            "for key in ['gamma','nu','b','c','Lambda']:\n"
            "    plt.plot(x, Om[key], label=key)\n"
            "plt.axvline(np.log(1/(1+bg.z_eq())), ls='--', c='k')\n"
            "plt.xlabel('x = ln a'); plt.ylabel('Omega_i(x)')\n"
            "plt.legend(); plt.title('density fractions'); plt.show()"),
        new_markdown_cell(
            "## 読み取り問い\n"
            "1. 放射優勢から物質優勢への遷移はどの $x$ で起きるか？\n"
            "2. $\\Omega_\\Lambda$ が 1 に近づくのはいつか？\n\n"
            "## 課題\n"
            "`Params(h=...)` を変えて $z_{\\rm eq}$ がどう動くか調べよ。"),
    ]
    return nb


def main():
    for name, builder in [("NB0_intro.ipynb", nb0), ("NB1_background.ipynb", nb1)]:
        nb = builder()
        nbf.write(nb, OUT / name)
        print("wrote", OUT / name)


if __name__ == "__main__":
    main()
