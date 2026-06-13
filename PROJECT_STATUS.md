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
pytest は §5.1–5.2 の単体テスト全20件が緑。

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

### WP5 本文（10/14章 本文記述済み）
- 全14章＋付録A–F の `.tex` を `textbook/` に配置し `main.tex` で部・章・付録を構成。
- **本文記述済み（10章）**: 第1（FLRW）・2（熱史）・3（再結合）・8（音響振動）・
  9（減衰）・10（視線積分）・11（$C_\ell$）・12（ピーク）・13（パラメータ依存）・
  14（数値実装）。背景→ソース→$C_\ell$→数値→検証→パラメータ依存の本筋が通っている。
- **構造スケルトン（残4章）**: 第4–7章（ゲージ・Einstein線形化・Boltzmann・初期条件、
  第II部の GR/Boltzmann 導出）は節構成・導出契約ID・章末枠を持つ
  （`scripts/make_textbook_skeleton.py`）。これら formal な導出本文が残作業。

### CLASS 比較（`03` §5.3）— 実施済み
- `classy` 3.3.4 を導入し `scripts/make_class_comparison.py` で設定A比較を実行。
  **中央値 0.79%、$\ell\lesssim950$ で <3%**。第一〜第三ピークは CLASS とほぼ一致し
  振幅規格化の正しさを確認（旧「+12〜25%」は Planck フル物理との差で、較正不要だった）。
  $\ell\gtrsim1000$ の +最大13% は偏光無視の減衰系統（付録E）。`figures/class_tt.csv` に
  キャッシュ。

### 視線積分の性能 — 改善済み
- ベッセル表構築＋視線積分を ell 並列化（`spectrum.los_integral(nproc=...)`）。
  直列と bit 一致を確認、実測 2.4倍以上の高速化（高 $\ell$ ほど効く）。

## 既知の制約・残課題（次フェーズへ引き継ぎ）
- **偏光（減衰尾系統）**: $\ell\gtrsim1000$ の +最大13% を解消するには E/B 偏光の実装が必要
  （第2版; 付録E を昇格）。本版では系統として定量済み。
- **$\sum m_\nu$**: 質量ニュートリノは本版で CLASS 委譲（`04` §7）。$\tau$・$z_*$ は自前計算で
  定量提示済み（NB6/第13章）。
- **WP5 本文プローズ**: 第3–14章・付録の導出本文（構成・図・数値は完成、文章が残作業）。

## 再現方法
```bash
pip install -e .            # 依存: numpy scipy matplotlib pytest nbformat nbconvert
make test                  # 単体テスト（§5.1-5.2）全20件
make values                # figures/values.json と cls_fiducial.npz
python scripts/make_param_study.py   # param_study.npz（τ, z_* 依存）
make figures               # 全図 + figures/manifest.yaml
python scripts/make_notebooks.py     # NB0–NB8 生成
make notebooks             # nbconvert --execute で全NB完走確認
```
