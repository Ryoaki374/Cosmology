# 02 Jupyterノートブック仕様書

## 0. 共通要件

- 場所: リポジトリ `notebooks/`。命名 `NB1_background.ipynb` 等（下表）。
- 全NBの先頭セル: タイトル、対応する本文章、学習目標（3点以内）、所要時間目安。
- 計算の重い処理は全て `cmbcore` パッケージの関数呼び出しとし、**ノートブック内に物理計算ロジックを重複実装しない**（例外: NB7の解析近似式と、教育目的で「あえて素朴に書く」と明示したセル）。
- 各NBは「説明セル（Markdown, 数式は本文と同記法）→ 実行セル → 図 → 読み取り問い（Q&A形式、答えは折りたたみ）」の繰り返し構造。
- 各NBの末尾に【課題】2–3問（パラメータを変えて図を作り直す系）と、本文演習との対応。
- 実行要件: `jupyter nbconvert --execute` で完走 ≤ 10分/冊（NB4, NB5は ≤ 30分、結果キャッシュ機構を `cmbcore.io` に用意）。乱数は全て `seed=5220` 固定。
- 図は `figures/NBx_yy_<name>.pdf` と `.png` の両方を保存。本文が参照する図IDは下表の **太字**。
- スタイル: `cmbcore.plotstyle` を全NBで適用（serif、カラーユニバーサル配色、図サイズ既定 (7,4.5)）。x軸は原則 $x=\ln a$ または $a$（log）、スライド再現図は $a$（log）でスライドの軸範囲に合わせる。

## NB0_setup_and_units.ipynb（環境・単位系チュートリアル）
- 内容: インストール確認、`cmbcore.constants`（SI; CODATA値）と単位換算（Mpc, eV, K）、$\Omega_{\gamma0}h^2=2.47\times10^{-5}$ の計算演習、`Params` データクラスの使い方。
- 検収: 全定数の出典コメント、$\Omega_{\gamma0}$ 計算が本文S2.2と一致。

## NB1_background.ipynb（第1章対応）
- 計算: `cmbcore.background.BackgroundCosmology` で $\mathcal{H}(x), \eta(x), t(x)$, 密度パラメータ進化、距離。
- 図: **F1.1** $\Omega_i(a)$ と $\rho_i(a)$ / **F1.2** $\mathcal{H}(x)$ と漸近則 / **F1.3** $\eta(a)$ と $c/\mathcal{H}$、$k=10,100,1000\,H_0/c$ の水平線交差 / F1.4 $\eta_0, t_0, z_{\rm eq}, \chi(z_*)$ の数値表。
- 読み取り問い例: 「$a_{\rm eq}$ 以前に $\mathcal H$ が急になる理由」「$\Omega_\Lambda$ を0にすると $\eta_0$ は?」
- 検収値: $t_0 = 13.8\pm0.1$ Gyr、$z_{\rm eq}=3400\pm100$、$\chi(z_*)=13.9\pm0.2$ Gpc（その他は `03` §5.2 / `04` §3 の検証表に従う）。

## NB2_recombination.ipynb（第3章対応）
- 計算: `cmbcore.recombination.Recombination`。Saha単独、Peebles、両者比較。$\tau(x), \tau'(x), \tilde g(x)$。
- 図: **F2.1** $X_e(x)$（Saha vs Peebles、log）/ **F2.2** $\tau, |\tau'|$ / **F2.3** $\tilde g(x)$（全域＋再結合期拡大、$\int\tilde g\,dx=1$ の数値検証を出力）/ F2.4 $\lambda_{\rm mfp}/(c/\mathcal H)$。
- 追加実験セル: (a) freeze-out 残留電離度の $\Omega_b$ 依存、(b) 再電離ON（tanhモデル）時の $\tilde g$ の第2の山 — 第13章の伏線。
- 検収値: $z_*(\tau=1) = 1090\pm10$、$X_e$ freeze-out $\sim 2\times10^{-4}$ オーダー、$\tilde g$ ピーク位置と半値幅 $\Delta z = 80\pm15$。

## NB3_perturbations_single_k.ipynb（第6–9章対応; スライド再現の中核）
- 計算: `cmbcore.perturbations.PerturbationSolver.solve(k)` 単一モード。
- 図（スライド対応を明記）:
  - **F3.1** $\Theta_0(a)$, $k=1000H_0/c$（スライドp19再現）
  - **F3.2** $\delta_b(a)$ 同（p20再現; 振動後の単調成長まで）
  - **F3.3** $\Theta_0$ と $\delta_b$ の二軸重ね描き（p21）＋ $X_e$ 重ね（p22）
  - **F3.4** $k=10H_0/c$ と $k=10^4H_0/c$ の対比（p23: 超地平線凍結とSilk減衰）
  - **F3.5** $\delta_c(a)$, $k=10,100,1000H_0/c$（p25: スケール別成長開始）
  - **F3.6** $\Theta_{0,1,2}(a)$（p26: 多重極）
  - **F3.7** $\Phi, \Psi$ の進化（地平線進入・等密度期・Λ期の減衰; 第8章8.4と第12章ISWの数値的根拠）
- 各図に「何を読み取るか」のMarkdown解説（本文該当節を引用）。
- 検収: F3.1–F3.6 がスライドの定性的特徴（振動回数、位相同期、成長開始順序）を再現。tight coupling 切替点での解の連続性（跳び < 1%）。

## NB4_transfer_functions.ipynb（第10章対応）
- 計算: $k$ グリッド全体の解、ソース関数 $\tilde S(k,x)$、視線積分 → $\Theta_\ell(k)$。
- 図: **F4.1** $\tilde S$ の4成分比較（$k$ 2例）/ **F4.2** $\Theta_\ell(k)$（$\ell=2,100,220,500,1000$）/ F4.3 被積分関数 $\tilde S\,j_\ell$ の振動相殺の可視化（数値積分の難所の体感）/ F4.4 $j_\ell(k\chi)$ カーネル。
- 検収: $\Theta_\ell(k)$ が滑らか（グリッド起因のギザつきなし）、SW極限チェック（低 $k$ で $\Theta_\ell \approx \frac{1}{3}\Psi j_\ell$ に漸近、相対差<5%）。

## NB5_power_spectrum.ipynb（第11–12章対応）
- 計算: $C_\ell$（$2\le\ell\le1500$）、規格化、ソース項別分解、Planck 2018 binned TT データ（リポジトリに同梱するCSV; 出典明記）比較。
- 図: **F5.1** $\mathcal{D}_\ell$ vs Planck / **F5.2** ソース項別分解（SW/Doppler/eISW/lISW）/ **F5.3** cosmic variance帯 / **F5.4** ピーク位置・高さの自動検出表（$\ell_1,\ell_2,\ell_3$, $H_2/H_1$ 等）/ F5.5 スライドp33–34の「過去の宇宙でのスペクトル」再現（$a_{\rm obs}=0.1, 0.01$ で観測した場合; `eta0` を差し替えて計算）。
- 検収: `04` §3 の表（$\ell_1=220\pm5$ ほか）。F5.5はピークが低 $\ell$（大角度）へ移動することを示す（p33: ~3.6°, p34: ~15°）。

## NB6_parameter_dependence.ipynb（第13章対応; 要求仕様の定量分析）
- 構成: 13.2–13.8 の各節に対応するセクション。各セクションは (i) `cmbcore` で計算（可能なもの）、(ii) `classy` で計算（レンズ・偏光・質量ニュートリノが要るもの）、の2系統。`HAS_CLASSY` フラグで(ii)をスキップ可能にし、スキップ時はキャッシュ済みCSV（リポジトリ同梱）から図を再生。
- 必須実験と図:
  - **F6.1** $\Omega_bh^2 \in \{0.020, 0.0224, 0.025\}$: $\mathcal{D}_\ell$ と比、ピーク特徴量表
  - **F6.2** $\Omega_ch^2 \pm 10\%$: 同上＋ early ISW分解の変化
  - **F6.3** 幾何縮退: $(\Omega_\Lambda, \Omega_k)$ を $100\theta_*=$const 線上で動かして比プロット（$\ell>30$ でほぼ平坦、低 $\ell$ のlate ISWで割れる）
  - **F6.4** $\tau \in \{0.02, 0.054, 0.10, 0.15\}$（$A_se^{-2\tau}$ 固定版と $A_s$ 固定版の両方）: $\ell\gtrsim30$ の $e^{-2\tau}$ 抑制の実測、低 $\ell$ の例外。＋ `classy` でEE低 $\ell$ バンプ（$z_{\rm re}$ 依存、振幅 $\propto\tau^2$ のフィット確認）
  - **F6.5** $z_{\rm re} \in \{6,8,10,15\}$ ↔ $\tau(z_{\rm re})$ 解析式（S13.3）との照合表
  - **F6.6** 再結合期シフト: `Recombination(recomb_shift=±5%)`（Peebles係数の温度引数を一括スケールして $z_*$ を動かす実装; `03` §2.4）→ $r_s^*, \chi_*, \ell_1$, 減衰尾の応答。解析予測 $\partial\ln\ell_1/\partial\ln z_*$ と数値の比較表
  - **F6.7** $\sum m_\nu \in \{0, 0.06, 0.12, 0.3, 0.6\}$ eV（classy; $\theta_*$固定と$h$固定の2通り、lensed/unlensed両方）: TT比プロット、レンズ平滑化の変化、$P(k)$ 抑制 $-8f_\nu$ 線との比較
  - **F6.8** 応答行列ヒートマップ $\partial\ln\mathcal{D}_\ell/\partial\ln p$（5点ステンシル数値微分; $p \in \{A_s, n_s, \Omega_bh^2, \Omega_ch^2, \tau, h\}$）
- 各実験の末尾に「観測への含意」Markdown（本文13章の該当文と同期）。
- 検収: F6.4で $\ell=500$ の抑制率が $e^{-2\Delta\tau}$ と1%以内で一致 / F6.6で解析予測と数値の $\ell_1$ シフトが10%以内で整合 / F6.7で $\theta_*$ 固定時の unlensed TT 変化が $\sum m_\nu=0.3$ eV で数%以内に留まることを確認（「CMB単独では質量に鈍い→レンズ・LSSとの組合せ」の論点を数値で立証）。

## NB7_analytic_vs_numeric.ipynb（第8–9, 12章対応）
- 内容: 解析近似と `cmbcore` フル数値の突き合わせ。
  - WKB解（S8.2）vs $\Theta_0(x)$ 数値: **F7.1**
  - 放射優勢期 $\Phi$ 解析解 $3[\sin y - y\cos y]/y^3$ 型 vs 数値: F7.2
  - 減衰包絡線 $e^{-k^2/k_D^2}$ vs 数値振幅: **F7.3**
  - $r_s(x), k_D(x)$ の計算と $\ell_m$ 予言 vs NB5のピーク検出: F7.4（表）
  - SWプラトーの解析値 vs 数値低 $\ell$: F7.5
- ここだけはノートブック内に解析式を直接実装してよい（教育目的）。
- 検収: ピーク位置予言（位相補正込み）が数値と $\Delta\ell<10$、減衰包絡が $k<0.2\,{\rm Mpc}^{-1}$ で20%以内。

## NB8_validation_class.ipynb（第14章対応; オプション実行）
- 内容: 同一物理仮定（無偏光近似は不可能なので「CLASSフル vs cmbcore」と「CLASS設定をできる限り寄せた版」の2比較）での $C_\ell$ 相対差、収束テスト（$\ell_{\max}$, $n_k$, rtol を振る）、計算時間プロファイル。
- 図: **F8.1** 相対差 vs $\ell$（許容帯付き）/ F8.2 収束テスト3面 / F8.3 処理時間内訳。
- 検収: `00_README.md` §5 のDoD数値。

## 成果物一覧（実装エージェントの出力）
```
notebooks/NB0...NB8.ipynb
notebooks/data/planck2018_tt_binned.csv   # 出典URLをヘッダに
notebooks/data/class_cache/*.csv          # NB6/NB8のキャッシュ
figures/                                   # 全図（本文ビルドが参照）
```
