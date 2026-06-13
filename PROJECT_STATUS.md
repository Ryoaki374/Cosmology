# プロジェクト進捗（成果物ステータス）

`00_README.md`–`04_task_plan_and_review.md` の仕様に基づく実装の現況。
仕様が要求する成果物は (a) 数値カーネル `cmbcore`、(b) ノートブック NB0–NB8、
(c) 本文 LaTeX（約200p）、(d) 検証 CI／図版／`values.json` パイプライン。
**実装順序（`04` §1: cmbcore → 図版 → 本文）に従い、基盤である `cmbcore` を
最優先で完成させた。**

## 完了（WP0, WP1 中核）

### WP0 リポジトリ雛形
- `cmbcore/ tests/ notebooks/ scripts/ textbook/ figures/ specs/` 構成。
- `pyproject.toml`（依存固定）、`Makefile`（`test/values/figures/notebooks/book`）。
- LaTeX 雛形（`textbook/main.tex` ＋ 3種ボックス環境）と第1章サンプル。

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

## 既知の制約・残課題（次フェーズへ引き継ぎ）
- **CLASS比較（`03` §5.3）**: 当環境に `classy` 未導入のため未実施。
  キャッシュCSV比較モードの参照データ生成が必要。ピーク振幅の +12〜25% は
  偏光無視・$Y_p=0$ の系統差とみられ、CLASS較正で確定する。
- **視線積分の性能**: ベッセル表の構築（高 $\ell$）が単一プロセスで律速。
  ell並列化が次の最適化候補（結果には影響しない）。
- **WP2 ノートブック**: NB0/NB1 のみ実装（`scripts/make_notebooks.py`）。
  NB2–NB8 は同ジェネレータを拡張。
- **WP5 本文**: 第1章のみ。第2–14章＋付録は `01_textbook_spec.md` の
  導出契約に沿って執筆。
- **WP3 CI**: GitHub Actions ワークフロー未作成（pytest＋NB＋CLASS比較）。

## 再現方法
```bash
pip install -e .            # または pip install numpy scipy matplotlib pytest
make test                  # 単体テスト（§5.1-5.2）
make values                # figures/values.json と cls_fiducial.npz を生成
make figures               # 背景・再結合・スペクトル図を生成
python scripts/make_notebooks.py
```
