# 03 数値計算カーネル `cmbcore` 仕様書

実装言語: Python 3.11+。依存: numpy, scipy（必須）, matplotlib, （オプション）classy。SI単位、$c$ 明示。独立変数 $x=\ln a$、$' \equiv d/dx$。$\mathcal{H}=aH$（コード名 `Hp`）。

**式の正準典拠**: 本仕様の式は Callin (astro-ph/0606683) および Winther AST5220 講義ノート（https://cmb.wintherscoming.no/ ; GitHub HAWinther/AST5220-Cosmology）の規約に合わせてある。**実装エージェントは実装前に両典拠の該当式と本仕様を突き合わせ、符号・係数の不一致を発見したら実装を止めて `04` §6 の手順で差し戻すこと**（本仕様の転記ミスの可能性を常に疑う）。単体テストには「式レベルのチェック」（§5.1）を必ず含める。

## 1. モジュール `background.py`

クラス `BackgroundCosmology(params)`:
- $\mathcal{H}(x) = H_0\sqrt{(\Omega_{b0}+\Omega_{c0})e^{-x} + (\Omega_{\gamma0}+\Omega_{\nu0})e^{-2x} + \Omega_{k0} + \Omega_{\Lambda0}e^{2x}}$
- 派生量: $\mathcal{H}'$, $\mathcal{H}''$（解析微分で実装）, $H(x)$, $\Omega_i(x)$。
- $\eta(x)$: ODE $\frac{d\eta}{dx} = \frac{c}{\mathcal{H}}$ を $x_{\rm start}=-20$ から積分（初期値は放射優勢解析解 $\eta(x_{\rm start})=c/\mathcal{H}(x_{\rm start})$）。$t(x)$ 同様（$dt/dx = 1/H$, 初期値 $1/(2H)$）。
- スプライン化（scipy CubicSpline）して $x\in[-20,5]$ で評価可能に。
- 距離: 共動距離 $\chi(x) = \eta_0-\eta(x)$、$\Omega_k\ne0$ の一般式（$\sinh$/$\sin$）も実装（第13章用）。
- $\Omega_{\gamma0} = \frac{2\pi^2}{30}\frac{(k_BT_{\rm CMB})^4}{\hbar^3c^5}\cdot\frac{8\pi G}{3H_0^2}$、$\Omega_{\nu0} = N_{\rm eff}\cdot\frac{7}{8}(4/11)^{4/3}\Omega_{\gamma0}$（質量ゼロ既定）。

## 2. モジュール `recombination.py`

クラス `Recombination(bg, params)`:

### 2.1 Saha（$X_e > 0.99$ の間）
$$\frac{X_e^2}{1-X_e} = \frac{1}{n_b}\left(\frac{m_ek_BT_b}{2\pi\hbar^2}\right)^{3/2}e^{-\epsilon_0/(k_BT_b)},\quad \epsilon_0=13.605693\,\mathrm{eV},\ T_b=T_{\rm CMB}/a$$
$n_b = n_H = \frac{3H_0^2\Omega_{b0}}{8\pi G m_H}a^{-3}$（$Y_p=0$ 既定。$Y_p\ne0$ 対応はパラメータで用意するが既定OFF）。二次方程式の数値的に安定な解（桁落ち回避: $X_e = 2/(1+\sqrt{1+4/y})$ 形式）を用いる。

### 2.2 Peebles（$X_e \le 0.99$ で切替、以後ODE）
$$\frac{dX_e}{dx} = \frac{C_r(T_b)}{H}\left[\beta(T_b)(1-X_e) - n_H\,\alpha^{(2)}(T_b)\,X_e^2\right]$$
$$C_r = \frac{\Lambda_{2s\to1s} + \Lambda_\alpha}{\Lambda_{2s\to1s} + \Lambda_\alpha + \beta^{(2)}},\quad \Lambda_{2s\to1s}=8.227\,\mathrm{s}^{-1}$$
$$\Lambda_\alpha = H\,\frac{(3\epsilon_0)^3}{(8\pi)^2\,(\hbar c)^3\, n_{1s}},\quad n_{1s}=(1-X_e)n_H$$
$$\beta^{(2)} = \beta\,e^{3\epsilon_0/(4k_BT_b)},\qquad \beta = \alpha^{(2)}\left(\frac{m_ek_BT_b}{2\pi\hbar^2}\right)^{3/2}e^{-\epsilon_0/(k_BT_b)}$$
$$\alpha^{(2)} = \frac{8}{\sqrt{3\pi}}\,\sigma_T c\,\sqrt{\frac{\epsilon_0}{k_BT_b}}\,\phi_2,\qquad \phi_2 = 0.448\,\ln\!\frac{\epsilon_0}{k_BT_b}$$
注意: $\beta^{(2)}$ は $\beta$ の式に $e^{3\epsilon_0/4k_BT_b}$ を掛けた解析合成形で実装（指数のオーバーフロー回避）。$T_b=T_\gamma$ 近似（差は本文【数値の現実】で言及）。

### 2.3 光学的厚みと可視度
$$\tau' \equiv \frac{d\tau}{dx} = -\frac{c\,n_e\sigma_T}{H},\qquad n_e = X_e n_H$$
$\tau(x)$ は $x_0=0$ で $\tau=0$ の境界条件で後ろ向き積分。$\tilde g(x) = -\tau'e^{-\tau}$。$\tau'', \tilde g', \tilde g''$ はスプライン微分。検証: $\int\tilde g\,dx = 1\pm10^{-4}$。
$z_*$: $\tau(x_*)=1$ で定義。音地平線 $r_s(x) = \int_{-\infty}^{x} \frac{c_s\,c\,dx'}{\mathcal{H}}$, $c_s = c\sqrt{\frac{1}{3(1+R_s)}}$, $R_s = \frac{3\Omega_{b0}}{4\Omega_{\gamma0}}a$。（**注**: 摂動モジュールの結合定数 `R = 4Ωγ0/(3Ωb0 a)` は $1/R_s$。命名を `Rs_baryon_photon` / `R_coupling` と分けて混同を防ぐこと。）

### 2.4 拡張フック（第13章/NB6用）
- `reionization`: tanhモデル $X_e^{\rm re}(z) = \frac{1+f_{\rm He}}{2}\left[1+\tanh\frac{(1+z_{\rm re})^{3/2}-(1+z)^{3/2}}{\Delta_y}\right]$（$f_{\rm He}=0$ 既定、$\Delta z=0.5$）。Peebles解に加算合成。
- `recomb_shift`: Peebles/Saha評価時の $T_b \to (1+\epsilon_{\rm shift})T_b$ 一括スケール（$z_*$ を擬似的に動かすトイ実装; ±5%程度を想定）。

## 3. モジュール `perturbations.py`

クラス `PerturbationSolver(bg, rec, params)`。各 $k$ について ODE系を $x_{\rm start}=-18$（要件: $ck/\mathcal{H} < 10^{-3}$ かつ $|\tau'| > 10^{3}\,ck/\mathcal H$ を満たすこと; 満たさなければ自動で前倒し）から $x=0$ まで積分。

### 3.1 完全系（tight coupling 終了後）
光子（偏光なし, $\Pi = \Theta_2$）:
- (P1) $\Theta_0' = -\frac{ck}{\mathcal{H}}\Theta_1 - \Phi'$
- (P2) $\Theta_1' = \frac{ck}{3\mathcal{H}}\Theta_0 - \frac{2ck}{3\mathcal{H}}\Theta_2 + \frac{ck}{3\mathcal{H}}\Psi + \tau'\left[\Theta_1 + \frac{1}{3}v_b\right]$
- (P3) $\Theta_\ell' = \frac{\ell\,ck}{(2\ell+1)\mathcal{H}}\Theta_{\ell-1} - \frac{(\ell+1)ck}{(2\ell+1)\mathcal{H}}\Theta_{\ell+1} + \tau'\left[\Theta_\ell - \frac{1}{10}\Theta_2\,\delta_{\ell,2}\right]$ （$2\le\ell<\ell_{\max}$）
- (P4) 打ち切り: $\Theta_{\ell_{\max}}' = \frac{ck}{\mathcal{H}}\Theta_{\ell_{\max}-1} - \frac{c(\ell_{\max}+1)}{\mathcal{H}\,\eta(x)}\Theta_{\ell_{\max}} + \tau'\Theta_{\ell_{\max}}$

ニュートリノ（質量ゼロ）: (N1)–(N4) は (P1)–(P4) で $\tau'\to0$、$\Theta\to\mathcal{N}$、$\ell_{\max}^\nu$ 使用。

CDM・バリオン:
- (M1) $\delta_c' = \frac{ck}{\mathcal{H}}v_c - 3\Phi'$
- (M2) $v_c' = -v_c - \frac{ck}{\mathcal{H}}\Psi$
- (M3) $\delta_b' = \frac{ck}{\mathcal{H}}v_b - 3\Phi'$
- (M4) $v_b' = -v_b - \frac{ck}{\mathcal{H}}\Psi + \tau' R\,(3\Theta_1 + v_b)$、$R = \frac{4\Omega_{\gamma0}}{3\Omega_{b0}}e^{-x}$

計量:
- (G1) $\Phi' = \Psi - \frac{c^2k^2}{3\mathcal{H}^2}\Phi + \frac{H_0^2}{2\mathcal{H}^2}\left[\Omega_{c0}e^{-x}\delta_c + \Omega_{b0}e^{-x}\delta_b + 4\Omega_{\gamma0}e^{-2x}\Theta_0 + 4\Omega_{\nu0}e^{-2x}\mathcal{N}_0\right]$
- (G2) $\Psi = -\Phi - \frac{12H_0^2}{c^2k^2e^{2x}}\left[\Omega_{\gamma0}\Theta_2 + \Omega_{\nu0}\mathcal{N}_2\right]$（代数式; 右辺評価に使う）

### 3.2 初期条件（断熱; $f_\nu = \Omega_{\nu0}/(\Omega_{\gamma0}+\Omega_{\nu0})$）
$$\Psi = -\frac{1}{\frac{3}{2}+\frac{2f_\nu}{5}},\quad \Phi = -\left(1+\frac{2f_\nu}{5}\right)\Psi,\quad \delta_c=\delta_b=-\frac{3}{2}\Psi,\quad v_c=v_b=-\frac{ck}{2\mathcal{H}}\Psi$$
$$\Theta_0=-\frac{1}{2}\Psi,\quad \Theta_1=\frac{ck}{6\mathcal{H}}\Psi,\quad \Theta_2=-\frac{20ck}{45\mathcal{H}\tau'}\Theta_1,\quad \Theta_\ell=-\frac{\ell}{2\ell+1}\frac{ck}{\mathcal{H}\tau'}\Theta_{\ell-1}$$
$$\mathcal{N}_0=-\frac{1}{2}\Psi,\quad \mathcal{N}_1=\frac{ck}{6\mathcal{H}}\Psi,\quad \mathcal{N}_2=-\frac{c^2k^2e^{2x}(\Phi+\Psi)}{12H_0^2\,\Omega_{\nu0}},\quad \mathcal{N}_\ell=\frac{ck}{(2\ell+1)\mathcal{H}}\mathcal{N}_{\ell-1}$$
規格化は $\mathcal{R}$ 換算係数を `spectrum.py` で適用（§4.3）。

### 3.3 tight coupling（TC）レジーム
- 適用条件（全て満たす間）: $|\tau'| > 10$、$|\tau'| > 10\,\frac{ck}{\mathcal{H}}$、$x < x(X_e=0.99)$（再結合開始前）。
- TC中の積分変数: $\Theta_0, \Theta_1, \delta_c, v_c, \delta_b, v_b, \Phi$, $\mathcal{N}_\ell$。$\Theta_2 = -\frac{20ck}{45\mathcal{H}\tau'}\Theta_1$、$\Theta_{\ell\ge3}$ は §3.2 最終式の漸化で代数評価。
- TC中の修正方程式（Callin Eqs.(70)–(72)相当; **転写検証必須**）:
$$q = \frac{-\left[(1-R)\tau' + (1+R)\tau''\right](3\Theta_1+v_b) - \frac{ck}{\mathcal{H}}\Psi + \left(1-\frac{\mathcal{H}'}{\mathcal{H}}\right)\frac{ck}{\mathcal{H}}(-\Theta_0+2\Theta_2) - \frac{ck}{\mathcal{H}}\Theta_0'}{(1+R)\tau' + \frac{\mathcal{H}'}{\mathcal{H}} - 1}$$
$$v_b' = \frac{1}{1+R}\left[-v_b - \frac{ck}{\mathcal{H}}\Psi + R\left(q + \frac{ck}{\mathcal{H}}(-\Theta_0+2\Theta_2) - \frac{ck}{\mathcal{H}}\Psi\right)\right],\qquad \Theta_1' = \frac{1}{3}\left(q - v_b'\right)$$
- 切替時: TC終了時点の状態をフル系の初期値にコピー（$\Theta_{\ell\ge2}$ は代数値で初期化）。切替前後で $\Theta_0,\Theta_1$ の連続性を assert。

### 3.4 ソルバ
- `scipy.integrate.solve_ivp`, method=`LSODA`（fallback `BDF`）, `rtol=1e-8, atol=1e-10`（収束テストで根拠付け; NB8）。
- $\ell_{\max}=8$（光子）, $\ell_{\max}^\nu=10$ 既定（収束テストで±を確認）。
- 出力: 密な $x$ グリッド（$n_x=2000$, 再結合期 $x\in[-7.5,-6.5]$ を局所増点）上の全変数＋ $\Phi', \Psi', v_b'$ 等の導関数（RHS再評価で取得）。

## 4. モジュール `spectrum.py`

### 4.1 ソース関数（偏光なし, $\Pi=\Theta_2$）
$$\tilde S(k,x) = \tilde g\left[\Theta_0+\Psi+\tfrac{1}{4}\Theta_2\right] + e^{-\tau}\left[\Psi'-\Phi'\right] - \frac{1}{ck}\frac{d}{dx}\!\left(\mathcal{H}\tilde g\,v_b\right) + \frac{3}{4c^2k^2}\frac{d}{dx}\!\left[\mathcal{H}\frac{d}{dx}\!\left(\mathcal{H}\tilde g\,\Theta_2\right)\right]$$
微分はスプライン経由（端点処理に注意）。4項を個別に保持（NB4/NB5の分解図用）。

### 4.2 視線積分
$$\Theta_\ell(k) = \int_{x_{\rm start}}^{0}\tilde S(k,x)\,j_\ell\!\left[k\left(\eta_0-\eta(x)\right)\right]dx$$
- $j_\ell$: `scipy.special.spherical_jn` を $z$ グリッド（$\Delta z \le 2\pi/16$）で前計算しスプライン補間。$z<\ell/2$ 領域はゼロ近似可（指数的小ささ; 閾値は収束テスト）。
- 積分: 台形則、$x$ 刻みは再結合期で $\Delta x\le10^{-3}$、以降 $\Delta x\le0.05$（ISW用）。

### 4.3 $C_\ell$ と規格化
$$C_\ell = 4\pi\int_0^\infty \mathcal{P}_\mathcal{R}(k)\,\left|\Theta_\ell^{(\mathcal{R}=1)}(k)\right|^2\frac{dk}{k},\qquad \mathcal{P}_\mathcal{R}(k)=A_s\left(\frac{k}{k_{\rm pivot}}\right)^{n_s-1}$$
- **規格化監査**（本文11.3対応）: §3.2 の初期条件は特定の $\Psi$ 規格化なので、超地平線初期での $\mathcal{R} = \Phi - \frac{\mathcal H}{k}v$ …（本書符号規約で評価）を初期スナップショットから数値計算し、$\Theta_\ell^{(\mathcal{R}=1)} = \Theta_\ell/\mathcal{R}_{\rm ini}$ と換算する。換算係数は**初期スナップショットからの数値評価を正**とし、第7章S7.1/S7.2から導かれる超地平線解析値（$\mathcal{R}_{\rm ini}$ を $\Psi_{\rm ini}, f_\nu$ で表した閉形式。執筆エージェントが導出し、その式を単体テストの参照値として実装エージェントに引き渡す）との一致を単体テストで確認する。係数を取り違えると $C_\ell$ 全体が定数倍ずれるため、最終防衛線は CLASS との振幅比較（§5.3）。
- $k$ グリッド: $k_{\min}=5\times10^{-5}\,{\rm Mpc}^{-1}$, $k_{\max}=0.35\,{\rm Mpc}^{-1}$, $n_k=250$、二次スペーシング $k_i = k_{\min} + (k_{\max}-k_{\min})(i/n_k)^2$（低 $k$ 密）。積分はスプライン後 $\ln k$ で合成シンプソン。
- $\ell$ サンプリング: $\ell \in \{2,3,...,10,12,...,30(2刻み),35,...,100(5刻み),110,...,300(10刻み),325,...,1500(25刻み)\}$ で計算しスプライン。
- 単位変換: $\mathcal{D}_\ell = \frac{\ell(\ell+1)}{2\pi}C_\ell\,T_{\rm CMB}^2$ [μK²]。

### 4.4 拡張（第13章）
- `tau_reionization`: rec の reionization フックを通じ自動反映（ソース関数の $e^{-\tau}, \tilde g$ が変わるだけで構造は不変であることをテストで確認）。
- 過去の観測者モード（NB5 F5.5）: `eta_obs = η(a_obs)` で $j_\ell[k(\eta_{\rm obs}-\eta)]$、積分上限 $x_{\rm obs}$。

## 5. 検証仕様（単体・結合・参照比較）

### 5.1 式レベル単体テスト（pytest）
- 次元チェック: 全RHSをSI数値で評価し無次元性を assert（$c, H_0$ を10倍して不変）。
- 極限チェック: (a) $f_\nu\to0$ で $\Phi=-\Psi$ 初期、(b) 物質優勢解析解 $\Phi=$const, $\delta_c\propto a$ を $k=10H_0/c$ で1%再現、(c) 超地平線 $\mathcal{R}$ 保存（ドリフト<0.5%）、(d) TC接続連続性、(e) $\int\tilde g\,dx=1$、(f) Saha–Peebles接続点の滑らかさ。
- 光子バリオン運動量保存: $\frac{d}{dx}[\,$全運動量$\,]$ が衝突項相殺後に重力項のみとなることの数値確認。

### 5.2 物理特徴量テスト

| 量 | 目標値 | 許容 |
|---|---|---|
| $z_{\rm eq}$ | 3400 | ±100 |
| $z_*$ ($\tau=1$) | 1090 | ±10 |
| $r_s(z_*)$ | 144.5 Mpc | ±2 |
| $\chi(z_*)$ | 13.9 Gpc | ±0.2 |
| $100\,\theta_* = 100 r_s/\chi$ | 1.041 | ±0.01 |
| $\ell_1$ | 220 | ±5 |
| $\ell_2/\ell_1$ | ≈ 2.45 | ±0.05 |
| 可視度幅 $\Delta z$ | 80 | ±15 |
| $k_D(z_*)$ | ≈0.14 Mpc⁻¹ | ±20% |

（表の目標値は実装後に CLASS同等設定で再較正し、`04` §6 のレビューで確定値に更新すること。）

### 5.3 CLASS参照比較（CI; classy が無い環境ではキャッシュCSV比較）
- 設定A（フェア比較）: CLASSで `temperature contributions` をTTのみ・lensingなし、$Y_p\to$ 最小、$N_{\rm eff}=3.046$ 質量ゼロ、再電離off（`reio_parametrization=reio_none`）。許容: $\ell\le1000$ で3%、$\ell\le1500$ で5%。
- 設定B（フル）: CLASS既定（偏光込み計算のTT, lensed）との差を**図示のみ**（本書近似の誤差の可視化; 第14章 Fig 14.3）。

## 6. パッケージ構成

```
cmbcore/
  __init__.py  constants.py  params.py        # dataclass Params(fiducialファクトリ付き)
  background.py  recombination.py  perturbations.py  spectrum.py
  analytic.py    # WKB, k_D, r_s, ピーク位置公式（NB7用）
  io.py          # 結果キャッシュ(npz), values.json書き出し
  plotstyle.py
tests/           # §5の全テスト
scripts/make_figures.py  scripts/make_values.py
```
- docstring必須（式番号は本文の章.式番号で引用）。型ヒント必須。1関数50行以内目安。
- 性能目標: fiducial $C_\ell$ 一式（250 k × 階層）をシングルコア10分以内、`multiprocessing` で $k$ 並列化（既定: コア数）。
