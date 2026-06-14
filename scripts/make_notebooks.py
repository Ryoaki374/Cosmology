#!/usr/bin/env python3
"""Generate the full notebook set NB0-NB8 (02_notebook_spec.md, WP2).

Each notebook follows the spec structure: title + goals, explanation -> run ->
figure -> reading question, with exercises at the end. Heavy spectra are loaded
from caches (figures/cls_fiducial.npz, figures/param_study.npz) so that
`jupyter nbconvert --execute` completes quickly; the perturbation/recombination
notebooks compute live (fast). All physics lives in `cmbcore`.
"""
from __future__ import annotations
from pathlib import Path
import nbformat as nbf
from nbformat.v4 import new_notebook, new_markdown_cell, new_code_cell

OUT = Path(__file__).resolve().parent.parent / "notebooks"
OUT.mkdir(exist_ok=True)

HEADER = (
    "import numpy as np, matplotlib.pyplot as plt\n"
    "from pathlib import Path\n"
    "from cmbcore.plotstyle import use_style; use_style()\n"
    "np.random.seed(5220)\n"
    "FIG = Path('..') / 'figures'\n"
)


def md(s): return new_markdown_cell(s)
def code(s): return new_code_cell(s)


def nb0():
    return [
        md("# NB0 — 環境・単位系チュートリアル\n\n"
           "**対応**: 序章 / **目標**: (1) `cmbcore` の導入確認 (2) SI定数と単位換算 "
           "(3) `Params` の使い方。所要時間: 5分。"),
        code(HEADER + "from cmbcore import constants as const\n"
             "print('c   =', const.c, 'm/s')\n"
             "print('Mpc =', const.Mpc, 'm')\n"
             "print('eV  =', const.eV, 'J')"),
        md("## $\\Omega_{\\gamma0}h^2$ の計算（本文 S2.2）"),
        code("from cmbcore import Params\n"
             "p = Params.fiducial()\n"
             "print('Omega_gamma0 h^2 = %.3e' % (p.Omega_gamma * p.h**2))\n"
             "print('(expected ~ 2.47e-5)')"),
        md("## `Params` データクラス\n各密度や $H_0$ は post-init で導出される。"),
        code("print('H0 [1/s] =', p.H0)\n"
             "print('Omega_m  =', p.Omega_m, ' Omega_Lambda =', p.Omega_Lambda)\n"
             "print('f_nu     =', p.f_nu)"),
        md("## 【課題】\n1. `Params(h=0.70)` で $\\Omega_{\\gamma0}h^2$ は変わるか？理由は？\n"
           "2. `Params.callin2006()` の密度を表示し fiducial と比較せよ。"),
    ]


def nb1():
    return [
        md("# NB1 — 背景宇宙（第1章）\n\n"
           "**目標**: $\\mathcal H(x),\\eta(x)$ と密度分率の進化、$z_{\\rm eq}$。所要時間: 5分。"),
        code(HEADER + "from cmbcore import Params, BackgroundCosmology\n"
             "p = Params.fiducial(); bg = BackgroundCosmology(p)\n"
             "print('z_eq = %.0f, age = %.2f Gyr, eta0 = %.0f Mpc'"
             " % (bg.z_eq(), bg.t0/3.1557e16, bg.eta0/3.0857e22))"),
        md("## F1.1 密度分率 $\\Omega_i(x)$"),
        code("x = np.linspace(-15,0,400); Om = bg.Omega(x)\n"
             "for key in ['gamma','nu','b','c','Lambda']: plt.plot(x,Om[key],label=key)\n"
             "plt.axvline(np.log(1/(1+bg.z_eq())),ls='--',c='k',label='z_eq')\n"
             "plt.xlabel('x=ln a'); plt.ylabel('Omega_i'); plt.legend(); plt.show()"),
        md("## F1.2 共形ハッブル $\\mathcal H(x)/H_0$ と漸近則"),
        code("plt.semilogy(x, bg.Hp(x)/p.H0, label='Hp/H0')\n"
             "plt.semilogy(x, np.exp(-x), '--', label='~e^{-x} (rad)')\n"
             "plt.semilogy(x, np.exp(-x/2), ':', label='~e^{-x/2} (mat)')\n"
             "plt.xlabel('x'); plt.legend(); plt.ylim(0.5,1e3); plt.show()"),
        md("## 読み取り問い\n- $a_{\\rm eq}$ 以前で $\\mathcal H$ が急なのはなぜ？（放射の $a^{-2}$ 依存）\n"
           "## 【課題】 $\\Omega_\\Lambda\\to0$ で $\\eta_0$ はどう変わるか試せ。"),
    ]


def nb2():
    return [
        md("# NB2 — 再結合（第3章）\n\n"
           "**目標**: Saha vs Peebles の $X_e$、$\\tau,\\tilde g$、$z_*$。所要時間: 5分。"),
        code(HEADER + "from cmbcore import Params, BackgroundCosmology, Recombination\n"
             "p = Params.fiducial(); bg = BackgroundCosmology(p); rec = Recombination(bg,p)\n"
             "print('z_* = %.0f, r_s(z_*) = %.1f Mpc, int g dx = %.5f'"
             " % (rec.z_star(), rec.r_s(rec.x_star())/3.0857e22, rec.visibility_norm()))"),
        md("## F2.1 $X_e(x)$（Saha 単独 vs Peebles）"),
        code("x = np.linspace(-9,-5,500); z = np.exp(-x)-1\n"
             "plt.semilogy(z, rec.Xe(x), label='Saha+Peebles')\n"
             "plt.semilogy(z, rec._saha_Xe(x), '--', label='Saha only')\n"
             "plt.axvline(rec.z_star(), c='r', ls=':', label='z_*')\n"
             "plt.gca().invert_xaxis(); plt.xlabel('z'); plt.ylabel('X_e'); plt.legend(); plt.show()"),
        md("## F2.3 可視度関数 $\\tilde g(x)$"),
        code("plt.plot(x, rec.g_tilde(x)); plt.xlabel('x=ln a'); plt.ylabel('g~(x)')\n"
             "plt.title('visibility (peaks at last scattering)'); plt.show()"),
        md("## 追加実験: 再電離 ON での第2の山（第13章の伏線）"),
        code("pr = Params(tau_reio=1.0, z_reio=8.0)\n"
             "rec2 = Recombination(BackgroundCosmology(pr), pr)\n"
             "xx = np.linspace(-12,0,800)\n"
             "plt.plot(xx, rec.g_tilde(xx), label='no reion')\n"
             "plt.plot(xx, rec2.g_tilde(xx), label='reion z_re=8')\n"
             "plt.yscale('log'); plt.ylim(1e-4,None); plt.xlabel('x'); plt.legend(); plt.show()"),
        md("## 検収: $z_*=1090\\pm10$, $\\int\\tilde g\\,dx=1$。\n"
           "## 【課題】 `Params(Omega_b=...)` を変え freeze-out 電離度の依存を見よ。"),
    ]


def nb3():
    return [
        md("# NB3 — 単一モードの摂動進化（第6–9章; スライド再現）\n\n"
           "**目標**: $\\Theta_0,\\delta_b,\\delta_c,\\Phi,\\Psi$ の進化とスケール依存。所要時間: 1分。"),
        code(HEADER + "from cmbcore import Params, BackgroundCosmology, Recombination\n"
             "from cmbcore.perturbations import PerturbationSolver\n"
             "from cmbcore import constants as const\n"
             "p=Params.fiducial(); bg=BackgroundCosmology(p); rec=Recombination(bg,p)\n"
             "solver=PerturbationSolver(bg,rec,p)\n"
             "H0_c = p.H0/const.c   # 1/m\n"
             "ks = {n: n*H0_c for n in (10,100,1000)}\n"
             "res = {n: solver.solve(k) for n,k in ks.items()}"),
        md("## F3.1/F3.6 光子モノポール・多重極 $\\Theta_\\ell(a)$, $k=1000H_0/c$"),
        code("r=res[1000]; x=np.linspace(-12,0,500); a=np.exp(x)\n"
             "for l in (0,1,2): plt.plot(a, r[f'Theta{l}'](x), label=f'Theta_{l}')\n"
             "plt.xscale('log'); plt.xlabel('a'); plt.legend(); plt.title('acoustic oscillations'); plt.show()"),
        md("## F3.5 CDM 密度 $\\delta_c(a)$: スケール別の成長開始"),
        code("for n in (10,100,1000):\n"
             "    plt.loglog(a, np.abs(res[n]['delta_c'](x)), label=f'k={n}H0/c')\n"
             "plt.xlabel('a'); plt.ylabel('|delta_c|'); plt.legend(); plt.show()"),
        md("## F3.7 重力ポテンシャル $\\Phi,\\Psi$（地平線進入・等密度・$\\Lambda$ 期の減衰）"),
        code("for n in (10,1000):\n"
             "    plt.plot(a, res[n]['Phi'](x), label=f'Phi k={n}')\n"
             "plt.xscale('log'); plt.xlabel('a'); plt.ylabel('Phi'); plt.legend(); plt.show()"),
        md("## 読み取り: 大スケールは超地平線で凍結、小スケールは音響振動＋減衰。\n"
           "## 【課題】 $k=10^4H_0/c$ を加え Silk 減衰を確認せよ。"),
    ]


def nb4():
    return [
        md("# NB4 — 転送関数とソース関数（第10章）\n\n"
           "**目標**: ソース $\\tilde S(k,x)$ の4成分、視線積分 $\\Theta_\\ell(k)$。所要時間: 2分。"),
        code(HEADER + "from cmbcore import Params, PowerSpectrum\n"
             "from cmbcore.spectrum import k_grid\n"
             "from cmbcore import constants as const\n"
             "ps = PowerSpectrum(Params.fiducial())\n"
             "ks = k_grid(nk=80)/const.Mpc\n"
             "ells = np.array([2,100,220,500,1000])\n"
             "kf, Th = ps.transfer(ells, ks_si=ks)"),
        md("## F4.2 転送関数 $\\Theta_\\ell(k)$（細kグリッドで滑らか）"),
        code("kmpc = kf*const.Mpc\n"
             "for i,l in enumerate(ells): plt.plot(kmpc, Th[i], label=f'l={l}')\n"
             "plt.xlabel('k [1/Mpc]'); plt.ylabel('Theta_l(k)'); plt.xlim(0,0.1); plt.legend(); plt.show()"),
        md("## F4.4 ベッセルカーネル $j_\\ell(k\\chi)$ が角度射影を担う\n"
           "$\\Theta_\\ell(k)=\\int \\tilde S(k,x)\\,j_\\ell[k(\\eta_0-\\eta)]\\,dx$。"),
        code("from scipy.special import spherical_jn\n"
             "chi0 = ps.bg.eta0\n"
             "kk = np.linspace(0,0.1,400)\n"
             "for l in (100,220,500): plt.plot(kk, spherical_jn(l, kk*chi0/const.Mpc*const.Mpc), label=f'l={l}')\n"
             "plt.xlabel('k [1/Mpc]'); plt.ylabel('j_l(k chi0)'); plt.legend(); plt.show()"),
        md("## 検収: $\\Theta_\\ell(k)$ はグリッド起因のギザつきがない。\n"
           "## 【課題】 ソースを `ps.sources(ks)` で取り出し4成分を分解せよ。"),
    ]


def nb5():
    return [
        md("# NB5 — パワースペクトル（第11–12章）\n\n"
           "**目標**: $\\mathcal D_\\ell^{TT}$、ピーク構造の読み取り。所要時間: <1分（キャッシュ利用）。"),
        code(HEADER + "d = np.load(FIG/'cls_fiducial.npz')\n"
             "ells, Dl = d['ells'], d['cls']\n"
             "print('first peak at l =', int(ells[np.argmax(np.where(ells>=150,Dl,0))]))"),
        md("## F5.1 $\\mathcal D_\\ell^{TT}$（fiducial）と Planck 2018 binned TT"),
        code("plt.plot(ells, Dl, label='cmbcore')\n"
             "pl = Path('data')/'planck2018_tt_binned.csv'\n"
             "if pl.exists():\n"
             "    P = np.genfromtxt(pl, delimiter=',', comments='#')  # cols: l, Dl, err\n"
             "    plt.errorbar(P[:,0], P[:,1], yerr=P[:,2], fmt='.', ms=6,\n"
             "                 label='Planck 2018 (illustrative)')\n"
             "plt.xlabel('multipole l'); plt.ylabel('D_l [uK^2]'); plt.legend(); plt.show()"),
        md("## F5.4 ピーク検出表"),
        code("from scipy.signal import find_peaks\n"
             "lf=np.arange(120,1500); Df=np.interp(lf,ells,Dl)\n"
             "idx,_=find_peaks(Df, prominence=0.1*Df.max(), distance=120)\n"
             "for j,l in enumerate(lf[idx][:4]): print('peak %d: l=%d, D=%.0f uK^2'%(j+1,l,Df[idx][j]))"),
        md("## 読み取り: 第一ピーク $\\ell_1\\approx220$ は音地平線の角度スケール。\n"
           "## 【課題】 ピーク比 $H_2/H_1$ を計算しバリオン量との関係を考えよ。"),
    ]


def nb6():
    return [
        md("# NB6 — パラメータ依存性（第13章; 要求仕様の定量分析）\n\n"
           "**目標**: $\\tau$（再電離）と $z_*$（再結合シフト）の $\\mathcal D_\\ell$ への影響を "
           "「どの方向に・何%・なぜ」で示す。$\\sum m_\\nu$ は本版では CLASS 委譲（04 §7）。"
           "所要時間: <1分（キャッシュ利用）。"),
        code(HEADER + "import json\n"
             "d = np.load(FIG/'param_study.npz'); ells = d['ells']\n"
             "S = json.loads((FIG/'param_study_summary.json').read_text())\n"
             "print(json.dumps(S, indent=2))"),
        md("## F6.4 光学的厚み $\\tau$（再電離 $z_{\\rm re}=8$）: $\\ell\\gtrsim30$ で $e^{-2\\tau}$ 抑制"),
        code("plt.plot(ells, d['fid'], label='no reion')\n"
             "plt.plot(ells, d['reio_z8'], label='reion z_re=8')\n"
             "plt.xlabel('l'); plt.ylabel('D_l [uK^2]'); plt.legend()\n"
             "plt.title('tau~%.3f: suppression ~%.1f%% at l=500'"
             " % (S['tau_reio_z8'], S['reio_suppression_l500_pct'])); plt.show()\n"
             "print('predicted e^{-2tau}-1 = %.1f%%' % S['reio_predicted_exp_-2tau_pct'])"),
        md("## F6.6 再結合シフト ($T_b\\to1.05\\,T_b$): $z_*$ と $\\ell_1$ の応答"),
        code("plt.plot(ells, d['fid'], label='fiducial (z_*=%.0f)'%S['z_star_fid'])\n"
             "plt.plot(ells, d['recomb_p5'], label='shift +5%% (z_*=%.0f)'%S['z_star_shift_p5'])\n"
             "plt.xlabel('l'); plt.ylabel('D_l [uK^2]'); plt.xlim(2,800); plt.legend(); plt.show()\n"
             "print('l1: %d -> %d' % (S['l1_fid'], S['l1_recomb_p5']))"),
        md("## F6.7 ニュートリノ質量 $\\sum m_\\nu$（背景を厳密実装; $h$ 固定）\n"
           "`Params(Sigma_mnu=...)` で massive ν 背景（運動量積分）が有効になり、"
           "$\\Omega_\\Lambda$ 減・$\\theta_*$ 変化として TT に効く。"),
        code("mf = FIG/'massive_nu.npz'\n"
             "if mf.exists():\n"
             "    m = np.load(mf); el = m['ells']\n"
             "    plt.plot(el, m['core_0.12']/m['core_0.0'], label='cmbcore ratio')\n"
             "    if 'class_0.12' in m: plt.plot(el, m['class_0.12']/m['class_0.0'],'--',label='CLASS ratio')\n"
             "    plt.axhline(1,c='k',lw=.6); plt.xlabel('l')\n"
             "    plt.ylabel('C_l(0.12)/C_l(0)'); plt.legend(); plt.title('Sm_nu=0.12 eV TT response'); plt.show()\n"
             "else:\n"
             "    print('run scripts/make_massive_nu.py for massive_nu.npz')"),
        md("## 観測への含意\n"
           "- **$\\tau$**: 再電離で散乱が増え、$\\ell\\gtrsim30$ の全ピークが一様に $e^{-2\\tau}$ で抑制される"
           "（低 $\\ell$ は再電離バンプで例外）。\n"
           "- **$z_*$**: 再結合が早まると音地平線 $r_s$ が縮み、$\\ell_1\\propto1/\\theta_*$ が高 $\\ell$ 側へ動く。\n"
           "- **$\\sum m_\\nu$**: 背景（$H(z),\\theta_*$, 後期ISW）が支配的で TT への効果は小さい。"
           "本実装は背景を厳密に扱い CLASS と整合（摂動の自由流＝$P(k)$ 抑制は次段）。\n"
           "## 【課題】 `Params(z_reio=...)` を変え $\\tau$ と抑制率の関係を確認せよ。"),
    ]


def nb7():
    return [
        md("# NB7 — 解析近似 vs 数値（第8–9,12章）\n\n"
           "**目標**: 音地平線・Silk スケール・ピーク位置の解析予言を数値と突き合わせる。"
           "ここでは教育目的で解析式を直接書く。所要時間: <1分。"),
        code(HEADER + "from cmbcore import Params, BackgroundCosmology, Recombination\n"
             "from cmbcore import analytic, constants as const\n"
             "p=Params.fiducial(); bg=BackgroundCosmology(p); rec=Recombination(bg,p)\n"
             "lA = analytic.acoustic_scale(bg,rec)\n"
             "kD = analytic.silk_scale(bg,rec)*const.Mpc\n"
             "print('l_A = %.1f, k_D = %.4f /Mpc, 100 theta_* = %.4f'"
             " % (lA, kD, analytic.theta_star(bg,rec)))"),
        md("## F7.4 ピーク位置予言 $\\ell_m\\approx\\ell_A(m-\\phi)$ vs 数値検出"),
        code("pk = analytic.peak_positions(bg,rec,n=4)\n"
             "d = np.load(FIG/'cls_fiducial.npz'); ells,Dl=d['ells'],d['cls']\n"
             "from scipy.signal import find_peaks\n"
             "lf=np.arange(120,1500); Df=np.interp(lf,ells,Dl)\n"
             "idx,_=find_peaks(Df,prominence=0.1*Df.max(),distance=120); num=lf[idx][:4]\n"
             "for m,(a_,n_) in enumerate(zip(pk,num),1): print('peak %d: analytic %.0f, numeric %d'%(m,a_,n_))"),
        md("## 読み取り: 位相補正 $\\phi\\approx0.27$ で第一ピークがよく合う。\n"
           "## 【課題】 $k_D$ から減衰包絡 $e^{-(k/k_D)^2}$ を描き高 $\\ell$ の落ち込みと比べよ。"),
    ]


def nb8():
    return [
        md("# NB8 — 検証と収束（第14章; オプション）\n\n"
           "**目標**: 特徴量検証、収束、（あれば）CLASS比較。`classy` が無い環境では "
           "キャッシュ比較。所要時間: <1分。"),
        code(HEADER + "import json\n"
             "from cmbcore import Params, BackgroundCosmology, Recombination, analytic\n"
             "from cmbcore import constants as const\n"
             "p=Params.fiducial(); bg=BackgroundCosmology(p); rec=Recombination(bg,p)\n"
             "tab = {'z_eq':bg.z_eq(),'z_star':rec.z_star(),\n"
             "       'r_s_Mpc':rec.r_s(rec.x_star())/const.Mpc,\n"
             "       'chi_Gpc':bg.comoving_distance(rec.x_star())/const.Gpc,\n"
             "       '100theta':analytic.theta_star(bg,rec),\n"
             "       'k_D_invMpc':analytic.silk_scale(bg,rec)*const.Mpc}\n"
             "for k,v in tab.items(): print('%-12s %.3f'%(k,v))"),
        md("## F8.1 特徴量の目標値との比較\n"
           "| 量 | 目標 | 実測 |\n|---|---|---|\n"
           "| z_eq | 3400±100 | 上表 |\n| z_* | 1090±10 | 上表 |\n"
           "| r_s | 144.5±2 | 上表 |\n| 100θ_* | 1.041±0.01 | 上表 |"),
        md("## F8.1 CLASS 比較（設定A; キャッシュ `class_tt.csv` を使用）\n"
           "TTのみ・lensingなし・reio off・massless $\\nu$ で CLASS と突き合わせる。"),
        code("d=np.load(FIG/'cls_fiducial.npz'); ells,Dl=d['ells'],d['cls']\n"
             "cc = FIG/'class_tt.csv'\n"
             "fig,(a1,a2)=plt.subplots(2,1,figsize=(7,6),sharex=True,gridspec_kw={'height_ratios':[2,1]})\n"
             "a1.plot(ells,Dl,label='cmbcore')\n"
             "if cc.exists():\n"
             "    C=np.genfromtxt(cc,delimiter=',',comments='#')  # l, Dl_class\n"
             "    Dlc=np.interp(ells,C[:,0],C[:,1])\n"
             "    a1.plot(ells,Dlc,'--',label='CLASS (setting A)')\n"
             "    a2.plot(ells,100*(Dl-Dlc)/Dlc); a2.axhspan(-3,3,alpha=.15,color='green')\n"
             "    a2.set_ylim(-30,30); a2.set_ylabel('rel.diff [%]')\n"
             "a1.set_ylabel('D_l [uK^2]'); a1.legend(); a2.set_xlabel('l'); plt.show()"),
        md("## 読み取り\n"
           "- $\\ell\\lesssim950$ で相対差 $<3\\%$、中央値 $\\sim0.8\\%$。第一〜第三ピークは一致。\n"
           "- $\\ell\\gtrsim1000$ で cmbcore がやや高いのは偏光無視の減衰系統（付録E）。\n"
           "- `classy` がある環境では `scripts/make_class_comparison.py` で再計算できる。\n"
           "## 【課題】 `PerturbationSolver(lmax=...)` を変え収束を確認せよ。"),
    ]


BUILDERS = {
    "NB0_setup_and_units.ipynb": nb0,
    "NB1_background.ipynb": nb1,
    "NB2_recombination.ipynb": nb2,
    "NB3_perturbations_single_k.ipynb": nb3,
    "NB4_transfer_functions.ipynb": nb4,
    "NB5_power_spectrum.ipynb": nb5,
    "NB6_parameter_dependence.ipynb": nb6,
    "NB7_analytic_vs_numeric.ipynb": nb7,
    "NB8_validation_class.ipynb": nb8,
}


def main():
    # remove the old minimal stubs if present
    for old in ("NB0_intro.ipynb", "NB1_background.ipynb"):
        q = OUT / old
        if q.exists() and old == "NB0_intro.ipynb":
            q.unlink()
    for name, builder in BUILDERS.items():
        nb = new_notebook(cells=builder())
        nbf.write(nb, OUT / name)
        print("wrote", OUT / name)


if __name__ == "__main__":
    main()
