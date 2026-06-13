#!/usr/bin/env python3
"""Generate the textbook chapter/appendix skeletons (WP5 scaffold).

ch01/ch02 are written in full; this emits structural skeletons for ch03-ch14 and
appendices A-F from 01_textbook_spec.md: each has the chapter goal, the section
headings, derivation-contract placeholders, and a summary/exercise stub. Prose
and full derivations are the remaining WP5 work; the structure here lets `main.tex`
build the complete table of contents.
"""
from __future__ import annotations
from pathlib import Path

TB = Path(__file__).resolve().parent.parent / "textbook"

# (filename, chapter title, goal, [(sec-label, sec-title)...], [contract ids])
CHAPTERS = [
    ("ch03_recombination", "再結合とCMBの誕生",
     "電離度 $X_e(a)$ を Saha→Peebles で計算し、$\\tau,\\tilde g$ を定義して"
     "最終散乱面を数式として確立する。",
     [("3.1","水素再結合の素過程"),("3.2","Saha方程式（契約 S3.1）"),
      ("3.3","Peebles方程式（契約 S3.2）"),("3.4","光学的厚みと可視度関数（契約 S3.3）"),
      ("3.5","数値結果（Saha vs Peebles, $z_*\\approx1090$, $\\Delta z\\approx80$）"),
      ("3.6","再電離の予告（$z_{\\rm re}\\sim8$, $\\tau\\approx0.05$）"),
      ("3.7","光子平均自由行程の発散と脱結合")],
     "S3.1–S3.3"),
    ("ch04_gauge", "計量摂動とゲージ",
     "SVT分解とゲージ自由度を理解し、conformal Newtonian ゲージを選択理由つきで固定する。",
     [("4.1","なぜ摂動論か"),("4.2","計量摂動とSVT分解"),
      ("4.3","ゲージ変換（契約 S4.1）"),("4.4","conformal Newtonian ゲージの固定"),
      ("4.5","ゲージ不変量 $\\mathcal R$"),("4.6","流体変数の定義（速度・符号規約）")],
     "S4.1"),
    ("ch05_einstein", "アインシュタイン方程式の線形化",
     "Newtonianゲージで $\\delta G^\\mu{}_\\nu=8\\pi G\\,\\delta T^\\mu{}_\\nu$ を成分計算し、"
     "実装に使う2本（$\\Phi$ 発展式・$\\Psi$ 代数式）に整理する。",
     [("5.1","摂動 Christoffel・Ricci"),("5.2","$\\delta T^\\mu{}_\\nu$ と異方性応力（契約 S5.2）"),
      ("5.3","4本の方程式と独立2本の選択"),("5.4","ニュートン極限の確認")],
     "S5.1–S5.2"),
    ("ch06_boltzmann", "ボルツマン方程式",
     "$df/d\\lambda=C[f]$ から光子・CDM・バリオン・ニュートリノの発展方程式を全て導出する"
     "（本書の心臓部）。",
     [("6.1","位相空間と Liouville 演算子（契約 S6.1）"),("6.2","光子分布の摂動と輝度関数"),
      ("6.3","無衝突部分"),("6.4","トムソン衝突項"),
      ("6.5","多重極展開と階層方程式（契約 S6.3）"),("6.6","CDM（契約 S6.4）"),
      ("6.7","バリオンと運動量交換 $R=4\\Omega_{\\gamma0}/(3\\Omega_{b0}a)$（契約 S6.5）"),
      ("6.8","ニュートリノ（質量ゼロ／質量あり）")],
     "S6.1–S6.5"),
    ("ch07_initial", "初期条件と原始ゆらぎ",
     "超地平線極限で断熱初期条件一式を導き、$\\mathcal R$ 保存と規格化、原始スペクトルを確立する。",
     [("7.1","超地平線極限での縮約"),("7.2","断熱条件と全初期値（契約 S7.1）"),
      ("7.3","$\\mathcal R$ の保存と規格化（契約 S7.2）"),("7.4","インフレーションの最小限"),
      ("7.5","ガウシアン確率場とパワースペクトル")],
     "S7.1–S7.2"),
    ("ch08_acoustic", "強結合極限と音響振動",
     "強結合極限で階層を $\\Theta_0,\\Theta_1$ に縮約し、強制振動子方程式と WKB 解を導く。",
     [("8.1","強結合展開"),("8.2","光子バリオン流体の振動子方程式（契約 S8.1）"),
      ("8.3","WKB解と音地平線（契約 S8.2）"),("8.4","重力駆動とポテンシャル減衰（契約 S8.3）"),
      ("8.5","数値との突き合わせ")],
     "S8.1–S8.3"),
    ("ch09_damping", "光子拡散と減衰",
     "強結合1次の展開からせん断粘性・熱伝導による Silk 減衰を導き、$k_D(\\eta)$ を定量化する。",
     [("9.1","平均自由行程とランダムウォーク"),("9.2","強結合1次展開と $k_D$（契約 S9.1）"),
      ("9.3","可視度幅によるぼかし"),("9.4","数値確認（$k_D\\approx0.14$ Mpc$^{-1}$）")],
     "S9.1"),
    ("ch10_los", "自由伝播と視線積分",
     "line-of-sight 形式解を導き、ソース関数を SW・Doppler・ISW・四重極補正に分解する。",
     [("10.1","形式解の導出（契約 S10.1）"),("10.2","$\\mu$ 部分積分と $j_\\ell$ の出現（契約 S10.2）"),
      ("10.3","ソース関数の物理（4項）"),("10.4","瞬間再結合極限と Sachs-Wolfe（契約 S10.3）")],
     "S10.1–S10.3"),
    ("ch11_cl", "$C_\\ell$ の構成",
     "天球温度場の統計から $C_\\ell=4\\pi\\int\\frac{dk}{k}\\mathcal P_\\mathcal R|\\Theta_\\ell/\\mathcal R|^2$ を導出し、"
     "規格化監査と cosmic variance まで含める。",
     [("11.1","球面調和展開と $C_\\ell$（契約 S11.1）"),("11.2","平面波展開との接続（契約 S11.2）"),
      ("11.3","規格化監査"),("11.4","cosmic variance（契約 S11.3）"),
      ("11.5","結果（$\\mathcal D_\\ell$ vs Planck）"),("11.6","$\\theta\\leftrightarrow\\ell$ 対応")],
     "S11.1–S11.3"),
    ("ch12_peaks", "音響ピークの解剖学",
     "解析理解と射影を合成し、SWプラトー・ピーク位置・奇偶非対称・減衰尾・ISW を指差し確認する。",
     [("12.1","ピーク位置の公式（$\\ell_1\\approx220$）"),("12.2","奇数/偶数ピークの非対称"),
      ("12.3","SWプラトーと積分SW（契約 S12.1）"),("12.4","減衰尾"),("12.5","総合図解")],
     "S12.1"),
    ("ch13_params", "パラメータ依存性の定量分析",
     "各パラメータの $C_\\ell$ 応答を (i)機構 (ii)解析スケーリング (iii)数値（NB6）の3点セットで示す"
     "（要求仕様の中核）。",
     [("13.1","方法論: 応答の測り方"),("13.2","$A_s, n_s$"),
      ("13.3","$\\Omega_bh^2$（奇偶比・ピーク移動・減衰尾）"),("13.4","$\\Omega_mh^2$（$a_{\\rm eq}$・driving・ISW）"),
      ("13.5","$\\Omega_\\Lambda$・曲率（幾何縮退）"),
      ("13.6","光学的厚み $\\tau$ と $z_{\\rm re}$（$e^{-2\\tau}$ 抑制; 契約 S13.1）"),
      ("13.7","再結合期 $z_*$ のずれ（$\\partial\\ln\\ell_1/\\partial\\ln z_*$）"),
      ("13.8","ニュートリノ質量 $\\sum m_\\nu$（契約 S13.2; CLASS）"),
      ("13.9","まとめ: 応答行列と縮退")],
     "S13.1–S13.3"),
    ("ch14_numerics", "数値実装の設計",
     "`cmbcore` の設計を本文として解説する（`03_numerics_spec.md` の人間可読版）。"
     "ソースを粗 k で解き k方向スプライン補間→細 k で視線積分する手法を含む。",
     [("14.1","モジュール構成と単位系"),("14.2","ODEソルバと tight coupling 切替"),
      ("14.3","視線積分の収束（粗解＋細k）"),("14.4","検証プロトコルと CLASS 比較")],
     "—"),
]

APPENDICES = [
    ("appA_gr", "A", "一般相対論の計算詳細と規約対応表",
     "FLRW・摂動 Christoffel/Ricci 全成分表、Ma\\&Bertschinger・Dodelson・本書の変数対応表（符号込み）。"),
    ("appB_special", "B", "特殊関数公式集",
     "Legendre展開・直交性、Rayleigh展開、$j_\\ell$ 漸化式・漸近形、"
     "$\\int_0^\\infty j_\\ell^2(u)\\,du/u=1/[2\\ell(\\ell+1)]$ の証明。"),
    ("appC_gaussian", "C", "ガウシアン確率場の統計",
     "実空間/調和空間、Wickの定理、推定量 $\\hat C_\\ell$ の分散。"),
    ("appD_peebles", "D", "Peebles方程式の導出詳細",
     "3準位原子、case-B 再結合係数のフィット式、Lyα 脱出（Sobolev）の導出。"),
    ("appE_polarization", "E", "偏光の最小限",
     "トムソン散乱の偏光依存性、四重極→直線偏光、E/B の定義、$\\Pi=\\Theta_2+\\Theta_{P0}+\\Theta_{P2}$ の由来、"
     "再電離 EE バンプ、TT 近似誤差（減衰尾係数 $16/15\\to8/9$）。"),
    ("appF_notebooks", "F", "ノートブックガイド",
     "NB0–NB8 と本文章の対応表、実行環境。"),
]

CH_TMPL = """% 第{n}章 {title}（WP5 skeleton; 本文・導出は執筆予定）
\\section{{{title}}}
\\label{{ch:{label}}}

\\paragraph{{この章で導くこと（入力→出力）}}
{goal}

% --- 節構成（01_textbook_spec.md より）-------------------------------------
{sections}

\\paragraph{{導出契約}}
本章の導出契約: {contracts}。各契約は「入力式→出力式」を中間ステップ省略なし
（1ステップ＝学部4年生が5分で検算可能）で履行する。\\emph{{本文プローズは WP5 で執筆。}}

\\paragraph{{まとめ・チェックリスト・演習}}
% TODO(WP5): 章末まとめ、チェックリスト、演習 3–6 問（うち1問は対応 NB の課題）。
"""

SEC_TMPL = "\\subsection{{{title}}}\n% TODO(WP5): {label} の本文。\n"

APP_TMPL = """% 付録{letter} {title}（WP5 skeleton）
\\section*{{付録{letter}\\quad {title}}}
\\addcontentsline{{toc}}{{section}}{{付録{letter}\\quad {title}}}

{desc}

% TODO(WP5): 本文。
"""


def main():
    inputs = []
    for i, (label, title, goal, secs, contracts) in enumerate(CHAPTERS, start=3):
        body = "".join(SEC_TMPL.format(title=t, label=s) for s, t in secs)
        text = CH_TMPL.format(n=i, title=title, label=label, goal=goal,
                              sections=body, contracts=contracts)
        (TB / f"{label}.tex").write_text(text)
        inputs.append(label)
        print("wrote", label)
    for fname, letter, title, desc in APPENDICES:
        (TB / f"{fname}.tex").write_text(
            APP_TMPL.format(letter=letter, title=title, desc=desc))
        inputs.append(fname)
        print("wrote", fname)
    print("\nadd to main.tex:")
    for label, *_ in CHAPTERS:
        print(f"\\input{{{label}}}")


if __name__ == "__main__":
    main()
