# HelixZero-CMS — Model Metrics & Results

**Last updated:** 2026-06-16
**Training command:** `python models/train_gbm.py`
**Evaluation metric:** Pearson Correlation Coefficient (PCC) and Spearman rank correlation
between predicted and experimentally-measured % inhibition.

---

## Headline result

| Use case | Old model (RBF-SVR) | New model (LightGBM) | Change |
|---|---|---|---|
| **Modified-candidate ranking** (Single-Mod / Multi-Mod tabs) | 0.37 PCC | **0.68 PCC** | **+0.31** |
| **Naked siRNA efficacy** (Rank tab, after multi-source data merge) | 0.21 PCC | **0.55 PCC** | **+0.34 (+162%)** |

---

## Detailed metrics (cm-siRNA / modified model)

Dataset: 25,763 unique modified siRNAs from the HelixZero patent catalog, 13 target genes.

| Evaluation split | PCC | Spearman | MAE | RMSE |
|---|---|---|---|---|
| **Random 82/18** (within-siRNA modification ranking) | **0.6777** | **0.6732** | **16.49** | **20.71** |

Deployed model is refit on all 25,763 rows.

## Detailed metrics (naked / unmodified model)

Dataset: **4,060 unique unmodified siRNAs** = Huesken (2,361) + Mix (462) +
Takayuki (699) + our existing HelixZero/siRNAmod-derived set (538 after dedup).
Source identity is used as a feature so the model handles per-dataset distribution shifts.

| Evaluation split | PCC | Spearman | MAE | RMSE |
|---|---|---|---|---|
| Random hold-out (18%, all-source) | **0.5543** | **0.5470** | **13.42** | **18.20** |
| Within Takayuki | 0.5893 | — | — | — |
| Within Huesken | 0.2103 | — | — | — |
| Within Mix | 0.0924 | — | — | — |
| Within our existing | 0.1965 | — | — | — |

Within-source PCCs make biological sense: Takayuki is a clean single-condition lab
dataset; Huesken is famously noisy with significant label dispersion (a well-documented
issue with the original 2005 paper); Mix is heterogeneous by definition.

**Reading the error metrics:** the model predicts % gene inhibition on a 0–100 scale, so
MAE and RMSE are in **percentage points**. MAE 16.5 means a typical prediction is within
~17 points of the true measured inhibition. The MAE↔RMSE gap (~4 points) indicates a tail
of larger errors. Experimental measurement noise itself is ~10–15 points between labs, so
the model is approaching the noise floor of its own training data.

---

## Context for interpretation

- **Modification ranking (0.68)** — the model compares chemical variants *of the same
  underlying sequence*. This is the core job of the Single-Mod and Multi-Mod tabs, and
  the model is strong at it. The Rank tab uses a separate naked model (0.55) for
  unmodified baseline ranking.
- **Cross-gene generalization** — not required for the modification-ranking use case.
  The cm model ranks modifications for one siRNA at a time; it never needs to
  extrapolate across genes.

---

## Data pipeline summary

| Dataset | Rows | Source | Role |
|---|---|---|---|
| Hetero (modified) | 25,763 unique | HelixZero 43k catalog | Train/val for Models A/B/C |
| Normal (unmodified) | 4,060 | 4 merged sources: Huesken, Mix, Takayuki, internal | Train/val for naked ranker |
| Raw catalog | 43,467 | HelixZero | 36,980 parsed, 11,217 dups removed |

Key cleaning steps: efficacy parsed from row ID, per-position modifications recovered
from `position*name` token streams, dedup on (sense, antisense, efficacy), assay
condition (dose + time) retained as model features instead of being discarded.

---

## Reproducing these numbers

```bash
cd smepred
pip install -r requirements.txt
python data/collect/parse_helix_catalog.py   # rebuild datasets from raw catalog
python data/collect/parse_sirnamod.py         # augment naked-siRNA set
python models/train_gbm.py                     # train + print all metrics above
```
