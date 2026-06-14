# プロジェクト進捗（成果物ステータス）

`00_README.md`–`04_task_plan_and_review.md` の仕様に基づく実装の現況。
仕様が要求する成果物は (a) 数値カーネル `cmbcore`、(b) ノートブック NB0–NB8、
(c) 本文 LaTeX（約200p）、(d) 検証 CI／図版／`values.json` パイプライン。
**実装順序（`04` §1: cmbcore → 図版 → 本文）に従い、基盤である `cmbcore` を
最優先で完成させた。**

## 完了（WP0–WP4 中核 ＋ WP5 構成）

### WP0 リポジトリ雛形
- `cmbcore/ tests/ notebooks/ scripts/ textbook/ figures/ specs/` 構成。
- `pyproject.toml`（依存固定, `pip install -e .` 可）、`Makefile`（`test/values/figures/notebooks/book/all`）。
- LaTeX 雛形（`textbook/main.tex` ＋ 3種ボックス環境）。

### WP1 `cmbcore`（`03_numerics_spec.md` 全節を実装）
| モジュール | 状態 | 備考 |
|---|---|---|
| `constants.py` `params.py` | ✅ | SI単位、fiducial / callin2006 プリセット、放射密度の自動導出 |
| `background.py` | ✅ | $\mathcal H,\mathcal H',\mathcal H''$（解析）、$\eta,t$ ODE、距離、$z_{\rm eq}$ |
| `recombination.py` | ✅ | Saha→Peebles、$\tau,\tilde g$、$r_s$、再電離/recomb_shift フック |
| `perturbations.py` | ✅ | TC＋フル階層、断熱IC、$\Theta_{\le8},\mathcal N_{\le10}$、ベクトル化記録 |
| `spectrum.py` | ✅ | 4項ソース関数、視線積分、$C_\ell$、k並列（multiprocessing） |
| `analytic.py` | ✅ | $\ell_A$、ピーク位置、$k_D$、$\theta_*$ |
| `io.py` `plotstyle.py` | ✅ | npzキャッシュ、`values.json` 書き出し、共通スタイル |

### 検証（`03` §5.2 特徴量表）— 実測値
| 量 | 目標 | 実測 | 判定 |
|---|---|---|---|
| $z_{\rm eq}$ | 3400±100 | 3403 | ✅ |
| $z_*(\tau=1)$ | 1090±10 | 1082 | ✅ |
| $r_s(z_*)$ | 144.5±2 Mpc | 143.4 | ✅ |
| $\chi(z_*)$ | 13.9±0.2 Gpc | 13.89 | ✅ |
| $100\theta_*$ | 1.041±0.01 | 1.033 | ⚠ ほぼ境界 |
| 可視度 $\int\tilde g\,dx$ | 1±1e-4 | 1.00000 | ✅ |
| $k_D(z_*)$ | 0.14±20% Mpc⁻¹ | 0.149 | ✅ |
| 音響ピーク $\ell_{1,2,3}$ | 220/540/810 | 220/525/800 | ✅ |
| 宇宙年齢 | ~13.8 Gyr | 13.81 | ✅ |

$C_\ell^{TT}$ は滑らかで実観測同様の音響ピーク系列（$\ell\approx220/525/800/1130/1430$、
SWプラトー $\mathcal D_\ell\sim10^3\,\mu K^2$、Silk減衰尾）を再現する
（`figures/cl_tt_fiducial.png`）。ピーク振幅は標準値より約 +12〜25% 高い
（偏光無視・$Y_p=0$・CLASS未較正の残差; §5.3 で詰める）。
**視線積分はソースを粗 k で解き k方向スプライン補間→細 k で積分する方式**
（CLASS/CAMB と同様）で、$\Theta_\ell(k)$ の急速 k 振動を解像しジッタを解消した。
振幅規格化は IC が $\mathcal R=1$ 正規化
（$\Psi_{\rm ini}=-\tfrac23\mathcal R/(1+\tfrac4{15}f_\nu)$）である事実から
$\mathcal R_{\rm ini}=1$ として確定（`spectrum._R_ini`）。
pytest は §5.1–5.2 の単体テスト全21件が緑（偏光 TT/EE/TE 検査を含む）。

### WP2 ノートブック NB0–NB8（`02_notebook_spec.md`）
- 全9冊を `scripts/make_notebooks.py` で生成。重いスペクトルはキャッシュ
  （`cls_fiducial.npz`, `param_study.npz`）を読み込み、`nbconvert --execute` が高速に完走。
- NB0 単位 / NB1 背景 / NB2 再結合 / NB3 単一モード摂動 / NB4 転送関数 /
  NB5 パワースペクトル（Planck重ね） / NB6 パラメータ依存性 / NB7 解析vs数値 /
  NB8 検証・収束。各冊「説明→実行→図→問い→課題」構造。

### WP3/WP4 検証・図版・values
- `scripts/make_values.py` → `figures/values.json`（本文差し込み用の全数値）。
- `scripts/make_figures.py` → 背景・再結合・摂動・スペクトル・パラメータ依存の図、
  `figures/manifest.yaml`。
- `scripts/make_param_study.py` → 必須トピック（$\tau$, $z_*$, $z_{\rm re}$）の
  定量計算と `param_study.npz` / `param_study_summary.json`。
- `.github/workflows/ci.yml`（pytest＋values 生成）。

### WP5 本文（全14章＋付録A–F 本文記述済み）✅
- 全14章を本文執筆: 第1（FLRW）・2（熱史）・3（再結合）・4（ゲージ）・5（Einstein線形化）・
  6（Boltzmann）・7（初期条件）・8（音響振動）・9（減衰）・10（視線積分）・11（$C_\ell$）・
  12（ピーク）・13（パラメータ依存）・14（数値実装）。各章に「入力→出力」宣言、導出契約、
  3種ボックス、章末まとめ・演習を備える。
- 付録 A（GR・規約対応表）・B（特殊関数）・C（ガウシアン統計）・D（Peebles導出）・
  E（偏光最小限）・F（NBガイド）も本文記述済み。
- 摂動方程式（G1/G2, P/N/M, 初期条件）は `cmbcore` 実装と文字単位で一致。
- `main.tex` が部・全14章・付録A–F を構成（LaTeX 環境のバランス検査済み；
  ビルドは toolchain のある環境で `make book`）。残るは図の最終差し込みと相互参照の
  通しビルド（WP6）。

### 偏光（Eモード）実装 — 完了 ✅
- 偏光階層 $\Theta^P_\ell$ を `perturbations.py` に追加し、温度の四重極ソースを
  $\Pi=\Theta_2+\Theta^P_0+\Theta^P_2$ に置換。tight coupling の四重極閉包も偏光込み
  （$8/15$）に更新、IC も Callin の偏光初期条件に。
- `spectrum.py` に E モードソース $\tilde S_E=\frac34\tilde g\Pi/(k\chi)^2$ と幾何因子
  $\sqrt{(\ell+2)!/(\ell-2)!}$ を実装、`cls_all/dls_all` で **TT/EE/TE** を出力。

### CLASS 比較（`03` §5.3）— 実施済み・**DoD 達成**
- `classy` 3.3.4、`scripts/make_class_comparison.py`。偏光実装により
  **全 $\ell\le1500$ で <1.5%（中央値 0.37%、$\ell\le1000$ 最大 1.43%）**一致し、
  DoD#1（全$\ell$5%・$\ell\le1000$ 3%）を達成。減衰尾の旧系統（+13%）は解消。
  EE/TE は `scripts/make_polarization.py` で CLASS と比較（`figures/polarization.png`）。

### 視線積分の性能 — 改善済み
- ベッセル表構築＋視線積分を ell 並列化（`spectrum.los_integral(nproc=...)`）。
  直列と bit 一致を確認、実測 2.4倍以上の高速化（高 $\ell$ ほど効く）。

## 既知の制約・残課題（次フェーズへ引き継ぎ）
- **WP6 通しビルド**: LaTeX toolchain のある環境での PDF 通しビルド、相互参照・索引解決、
  本文への図の最終差し込みと裸数値監査（`values.json` 由来の確認）。本文・図・数値は揃済み。
- **$\sum m_\nu$ 摂動**: 質量ニュートリノは**背景を厳密実装**済み（`massive_nu.py`、CLASS と
  unlensed TT で中央値0.47%一致）。摂動の自由流（$q$ グリッド階層, $P(k)$ 抑制）は次段
  （TT 寄与は小）。
- **B モード偏光**: スカラーは E のみ生成（B はテンソル/レンズ起源）。本版は EE/TE まで。

## 再現方法
```bash
pip install -e .            # 依存: numpy scipy matplotlib pytest nbformat nbconvert
make test                  # 単体テスト（§5.1-5.2）全21件
make values                # figures/values.json と cls_fiducial.npz
python scripts/make_param_study.py   # param_study.npz（τ, z_* 依存）
python scripts/make_polarization.py  # TT/EE/TE（偏光; classy 任意）
python scripts/make_class_comparison.py  # CLASS 比較（classy 任意）
make figures               # 全図 + figures/manifest.yaml
python scripts/make_notebooks.py     # NB0–NB8 生成
make notebooks             # nbconvert --execute で全NB完走確認
```
