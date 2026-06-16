# HelixZero-CMS — Technology Stack & Model Rationale

A plain-English guide to every technology used, the model training workflow, and **why
each choice was made**.

---

## 1. Overview

HelixZero-CMS predicts how effectively an siRNA molecule silences a gene, including the effect
of chemical modifications. It has three layers:

```
  Data layer        →  Model layer            →  Serving layer
  (parse + clean)      (LightGBM regressors)     (FastAPI + HTML UI)
```

---

## 2. Full technology stack

| Layer | Technology | Why this one |
|---|---|---|
| Language | **Python 3.11** | De-facto standard for ML + bioinformatics; rich ecosystem. |
| ML model | **LightGBM (gradient-boosted trees)** | Best accuracy/speed trade-off for tabular composition features (see §4). |
| Numerics | **NumPy** | Fast vectorised feature math. |
| Data wrangling | **pandas** | CSV parsing, dedup, train/val splits. |
| Eval / splits | **scikit-learn**, **SciPy** | `GroupShuffleSplit`, `pearsonr`, `spearmanr`. |
| Model persistence | **joblib** | Efficient pickling of fitted models (`.pkl`). |
| API | **FastAPI + Uvicorn** | Async, auto-generated Swagger docs, type-validated requests. |
| Frontend | **Single-file HTML/CSS/JS** | Zero build step; opens anywhere; talks to the API via `fetch`. |
| CLI | **Click** | Scriptable batch predictions. |

---

## 3. Models in the system

| File | Predicts | Trained on | Features |
|---|---|---|---|---|
| `model_normal.pkl` | Naked (unmodified) siRNA efficacy | 4,060 unmodified siRNAs (4 sources) | 152-d + 4 source one-hot = 156-d |
| `model_a.pkl` / `b` / `c` | Modified (cm-siRNA) efficacy | 25,763 modified siRNAs | 152-d composition + position + condition |

> Historically Models A/B/C were three different feature recipes (MNC, MNC+seed,
> MNC+tail) for an SVR. The LightGBM rebuild unifies them into one richer feature set,
> so A/B/C now share the same strong model while the UI selector still works.

---

## 4. Why LightGBM (and why we moved off SVR)

The original implementation followed the SMEpred paper: **Support Vector Regression
(RBF kernel) + mononucleotide composition (MNC)**. On the real 43k patent catalog it
plateaued at **0.37 PCC**. Three problems:

1. **Scalability** — RBF-SVR is ~O(n²–n³) in samples, so we could only train on a 6,000-row
   subsample of 25k+. LightGBM trains on the full set in seconds.
2. **Feature expressiveness** — MNC only counts *how many* of each base/modification appear,
   not *where*. Position matters enormously in siRNA biology (the seed region, the 3′ end).
3. **The condition confound** — patent data records the same sequence under different doses
   and timepoints, giving different inhibition values. Composition features can't see this.

**LightGBM (gradient-boosted decision trees)** solves all three:
- Handles the full dataset fast.
- Automatically learns non-linear interactions between composition, modification position,
  and assay condition — exactly the structure of this problem.
- Robust to mixed feature scales (fractions, counts, log-dose) without heavy preprocessing.

Result: **0.37 → 0.68 PCC** on the modification-ranking task.

Alternatives considered: XGBoost (comparable accuracy, heavier install), neural nets
(need far more data than 25k rows; harder to interpret), staying with SVR (rejected — see above).

---

## 5. Feature engineering

Each siRNA pair → a single feature vector (`src/features.py → extract_batch_gbm`):

| Group | Dims | What it captures |
|---|---|---|
| Base + modified MNC | 140 | Nucleotide & modification composition of both strands |
| Modification density | 8 | Where mods sit: overall, seed (pos 1–8), 3′ tail, count |
| GC content | 2 | Duplex stability proxy |
| Assay condition | 2 | log10(dose nM), time/24h — fixed to a 10 nM / 24 h reference at inference |
| **Total** | **152** | |

The condition features let us **keep all data** (instead of discarding everything but one
dose) while removing the confound — the model learns the dose-response and we query it at
a standard reference condition.

---

## 6. Training workflow

```
raw HelixZero 43k catalog
        │  parse_helix_catalog.py   (recover sequences, mods, efficacy, dose, time, gene)
        ▼
data/hetero_train.csv + hetero_val.csv      data/normal_siRNA.csv
        │                                            │ (+ parse_sirnamod.py augment)
        ▼                                            ▼
        └──────────────  train_gbm.py  ──────────────┘
                 │
                 └─ Refit on all rows → save model_a/b/c.pkl, model_normal.pkl
```

**Why two evaluations?** Different product features need different checks:
- The **Rank tab** uses the naked model to rank unmodified siRNAs for baseline selection.
- The **Single/Multi-Mod tabs** compare modifications of a known siRNA → random split is the
  faithful measure. Cross-gene generalization is not required because modification ranking
  happens within a single siRNA at a time.

---

## 7. Serving workflow (inference)

```
HTTP request (FastAPI)
   │
   ├─ /rank        → generate 21-mers → features → model_normal → ranked naked siRNAs
   ├─ /single-mod  → 1260 variants → features → model_a → ranked modifications (+Δ vs parent)
   └─ /multi-mod   → one custom design → features → model_a → single efficacy score
```

The model predicts inhibition directly on the 0–100 scale (clipped), so scores are real
percentages — not batch-relative rescalings.

---

## 8. Reproducibility

All randomness is seeded (`random_state=42`). Re-running `train_gbm.py` reproduces the
metrics in [METRICS.md](METRICS.md) exactly.
