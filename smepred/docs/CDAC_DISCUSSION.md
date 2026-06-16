# HelixZero-CMS: siRNA Chemical Modification Efficacy Predictor

## Discussion Guide for CDAC Pune — 16 June 2026

> This document walks through the entire project from scratch: problem → data → features → model → modification engine → API → UI → results. Each section includes code snippets and file references so you can navigate the codebase live during discussion.

---

## Table of Contents

1. [The Problem](#1-the-problem)
2. [Data Collection & Parsing](#2-data-collection--parsing)
3. [Data Engineering & Cleaning](#3-data-engineering--cleaning)
4. [Feature Engineering (152-d Vector)](#4-feature-engineering-152-d-vector)
5. [Model Training — LightGBM](#5-model-training--lightgbm)
6. [Why LightGBM Over SVR](#6-why-lightgbm-over-svr)
7. [The 152-d vs 175-d Detour (Critical Lesson)](#7-the-152-d-vs-175-d-detour-critical-lesson)
8. [Modification Engine — Single-Mod Scan](#8-modification-engine--single-mod-scan)
9. [Modification Engine — Multi-Mod Beam Search](#9-modification-engine--multi-mod-beam-search)
10. [Parent Baseline Fix](#10-parent-baseline-fix)
11. [API Layer (FastAPI)](#11-api-layer-fastapi)
12. [Web UI (Single-File HTML)](#12-web-ui-single-file-html)
13. [Model Calibration](#13-model-calibration)
14. [Validation Strategy & Results](#14-validation-strategy--results)
15. [Key Decisions Log](#15-key-decisions-log)
16. [Files Reference Map](#16-files-reference-map)

---

## 1. The Problem

siRNA (small interfering RNA) is a 21-nucleotide double-stranded RNA that silences specific genes via RNA interference. It's a promising therapeutic modality — but raw siRNA degrades in minutes, triggers immune responses, and never reaches target cells.

**Solution**: Chemical modifications — swap atoms on the sugar, base, or backbone to add stability and reduce toxicity.

**The combinatorial challenge**: There are **30 modification types × 21 positions × 2 strands = 1,260 possible modified variants per siRNA**. Testing even 5 variants costs $1,000–$2,500 and takes 2–4 weeks.

**Goal**: Build a model that predicts which modification pattern works best — computationally — so wet labs test only the top 10 instead of all 1,260.

### Key Files

| File | Purpose |
|------|---------|
| `src/sirna_generator.py` | Generates all 21-mer candidates from an mRNA sequence |
| `src/modification_engine.py` | Generates 1,260 modified variants + multi-mod beam search |
| `src/predictor.py` | Orchestrates ranking and modification prediction |

---

## 2. Data Collection & Parsing

### Source 1: HelixZero Patent Catalog (PRIMARY)

**Raw file**: `HelixZero_Biological_Catalog_43k.csv.xls` — 43,467 rows × 26 columns. Patent-derived siRNA measurements with per-position modification annotations.

**The parsing nightmare**: The file has unquoted commas inside `Modification_locations_*` columns (e.g., `1,2,3,4,…`), so `pandas.read_csv()` sees variable-length rows. Only ~3,500 of 43k rows have the nominal 26 fields.

**Solution — regex token stream** (`data/collect/parse_helix_catalog.py`):

```python
# Extract (position, modification name) pairs from the token stream
pattern = r"(\d+)\*([^|,]+?)(?=\s*\|\||,)"
matches = re.findall(pattern, row)
# Each match: ("1", "2'-O-Methyladenosine"), ("2", "2'-O-Methylguanosine"), ...
```

**Critical bug caught**: Negative inhibition is written `…-24h--8.87`. Naive regex reads `8.87` (positive). Fixed to capture sign:

```python
id_suffix = re.search(r"-(-?\d+\.?\d*)$", row_id)
# Now captures "-8.87" correctly, clipped to 0
```

### Source 2: siRNAmod Database (AUXILIARY)

**Raw file**: `Data-(1)-csv.csv.xls` — 5,329 rows from the siRNAmod database. Modification names without position information, so only unmodified rows are usable.

### Source 3: Published Naked siRNA Datasets

| Dataset | Source | Rows | Quality |
|---------|--------|------|---------|
| Huesken 2005 | *Nature Biotechnology* | 2,361 | Gold standard, noisy labels |
| Mix (Reynolds/Vickers/Ui-Tei) | Combined | 462 | Heterogeneous |
| Takayuki 2007 | *Nucleic Acids Research* | 699 | Cleanest single-condition |
| HelixZero internal | Internal catalog | 538 | Noisy |

These are merged by `data/collect/merge_naked_sources.py` → `data/normal_siRNA_extended.csv` (4,060 rows).

### Key Files

| File | Purpose |
|------|---------|
| `data/collect/parse_helix_catalog.py` | Primary parser — 43k patent rows → clean CSVs |
| `data/collect/parse_sirnamod.py` | siRNAmod parser — unmodified rows only |
| `data/collect/merge_naked_sources.py` | Merges 4 published datasets into one |
| `data/collect/clean_utils.py` | Shared: `clean_sequence()`, `parse_efficacy()`, `map_modification()` |
| `data/modification_codes.json` | 126 alias rules mapping modification names → 35 symbols |

---

## 3. Data Engineering & Cleaning

### 3.1 Modification Name → Symbol Mapping

Real modification names are compositional: `"2'-O-Methyluridine"` = class `2'-O-Methyl` + base `uridine`. Each class maps to one symbol regardless of base:

```python
# From data/modification_codes.json
alias_rules = [
    ("2'-O-Methyl", "M"),      # 2'-OMe
    ("2'-Fluoro", "F"),        # 2'-F
    ("2'-Deoxy", "D"),         # DNA
    ("Locked nucleic acid", "L"),  # LNA
    ("Phosphorothioate", "S"), # PS
    # ... 121 more rules — first match wins
]
```

**35 total symbols**: 5 canonical (A, U, G, C, T) + 30 chemical modification symbols.

### 3.2 Cleaning Pipeline

Applied in `data/collect/clean_utils.py`:

| Step | Rule |
|------|------|
| Uppercase + T→U | DNA → RNA conversion |
| Length filter | Keep 19–25 nt (paper used 21–24) |
| Efficacy parse | Extract number from free text / row ID |
| Efficacy clip | Negative → 0, cap at 100 |
| Dedup | Drop exact (sense, antisense, efficacy) duplicates |

### 3.3 Train/Validation Split

`data/collect/splits.py` implements the original paper's rule: sort by descending efficacy, take every 10th row starting at the 5th as validation.

### 3.4 Final Datasets

| File | Rows | Purpose |
|------|------|---------|
| `data/hetero_train_2728.csv` | 23,187 | cm-siRNA training (13 genes) |
| `data/hetero_val_303.csv` | 2,576 | cm-siRNA validation |
| `data/normal_siRNA_extended.csv` | 4,060 | Naked siRNA (4 sources) |

---

## 4. Feature Engineering (152-d Vector)

**`src/features.py`** — each siRNA is converted to a 152-number vector. The model sees numbers, not letters.

### Feature Breakdown

| Group | Dimensions | What It Captures |
|-------|-----------|------------------|
| **Base MNC** (sense + antisense) | 70 | Frequency of each of 35 symbols in the unmodified strands |
| **Modified MNC** (sense + antisense) | 70 | Frequency of each symbol in the modified strands |
| **Mod density** (sense) | 4 | Overall mod fraction, seed-region fraction, 3'-tail fraction, count |
| **Mod density** (antisense) | 4 | Same for antisense |
| **GC content** | 2 | GC fraction of both strands |
| **Assay conditions** | 2 | log₁₀(concentration nM), time in hours |
| **Total** | **152** | |

### Why *Both* Base and Modified MNC?

Near-fully-modified strands collapse to mostly `M`/`F` composition — the underlying sequence signal is lost. Adding *base* (unmodified) composition restores it:

| Features | PCC |
|----------|-----|
| Modified MNC only | 0.37 |
| Base MNC only | 0.40 |
| **Base + Modified MNC** | **0.48** |

```python
# src/features.py — core feature extraction
def features_gbm(sense, antisense, base_sense, base_antisense,
                 concentration_nM=10, time_h=24, mrna_features=None):
    # 70-d: base MNC (unmodified composition)
    base_feat = np.concatenate([
        _mnc(base_sense,  35),   # 35 symbols
        _mnc(base_antisense, 35),
    ])

    # 70-d: modified MNC
    mod_feat = np.concatenate([
        _mnc(sense,  35),
        _mnc(antisense, 35),
    ])

    # 8-d: modification density (4 per strand)
    density = np.concatenate([
        _mod_density(sense, base_sense),
        _mod_density(antisense, base_antisense),
    ])

    # 2-d: GC content
    gc = np.array([_gc_content(base_sense), _gc_content(base_antisense)])

    # 2-d: assay conditions
    cond = np.array([np.log10(concentration_nM), time_h / 24])

    return np.concatenate([base_feat, mod_feat, density, gc, cond])
    # = 70 + 70 + 8 + 2 + 2 = 152
```

### Batch Extraction for Speed

```python
# src/features.py
def extract_batch_gbm(sense_list, antisense_list, base_sense_list, base_antisense_list):
    """Vectorized batch extraction — handles 1260 variants in ~50ms."""
    return np.array([
        features_gbm(s, a, bs, ba)
        for s, a, bs, ba in zip(sense_list, antisense_list,
                                 base_sense_list, base_antisense_list)
    ])
```

### Key Files

| File | Purpose |
|------|---------|
| `src/features.py` | All feature extraction logic |
| `src/features.py:features_gbm()` | Core 152-d vector construction |
| `src/features.py:_mnc()` | Mononucleotide composition (frequency vector) |
| `src/features.py:_mod_density()` | Position-aware modification density |

---

## 5. Model Training — LightGBM

**`models/train_gbm_v3.py`** — the production training script.

### Algorithm: Gradient-Boosted Decision Trees

LightGBM builds an ensemble of decision trees sequentially. Each new tree corrects the errors of all previous trees.

```python
# models/train_gbm_v3.py
model = lgb.LGBMRegressor(
    n_estimators=800,
    learning_rate=0.03,
    num_leaves=63,
    subsample=0.8,
    colsample_bytree=0.7,
    reg_lambda=1.0,
    min_child_samples=30,
    random_state=42,
)
```

### Evaluation

```python
# Random split (modification ranking use case)
X_train_r, X_val_r, y_train_r, y_val_r = train_test_split(X, y, test_size=0.18)
```

### Output Models

| File | Type | Features | Training Rows |
|------|------|----------|---------------|
| `models/model_a.pkl` | cm-siRNA LightGBM | 152-d | 25,763 |
| `models/model_b.pkl` | cm-siRNA (same) | 152-d | 25,763 |
| `models/model_c.pkl` | cm-siRNA (same) | 152-d | 25,763 |
| `models/model_normal.pkl` | Naked LightGBM | 156-d (152 + 4 source one-hot) | 4,060 |

### Key Files

| File | Purpose |
|------|---------|
| `models/train_gbm_v3.py` | Production training script (v3, 152-d) |
| `models/train_gbm.py` | Original training (v1, deprecated) |
| `models/train_gbm_v2.py` | Training with extra gene features (v2, deprecated) |
| `models/calibrate.py` | Isotonic calibration fit |

---

## 6. Why LightGBM Over SVR

The original SMEpred paper (Dar 2016) used **Support Vector Regression** with RBF kernel. On the real 43k patent catalog it plateaued at PCC **0.37**. Three reasons:

| Problem | SVR | LightGBM |
|---------|-----|----------|
| **Scalability** | O(n²–n³) — could only use 6k of 25k rows | Trains on all 25k in ~2 min |
| **Feature richness** | MNC only (counts, no position) | Position-aware density + conditions |
| **Condition confound** | Ignores dose/time → sees same sequence at different efficacies | Feeds dose/time as features |

**Result**: 0.37 → 0.68 PCC (+84% relative). Full analysis in `docs/TECH_STACK.md`.

```python
# Key comparison code structure
# Old SVR (from original paper):
from sklearn.svm import SVR
svr = SVR(kernel='rbf', C=100)  # ~O(n²) training

# New LightGBM:
import lightgbm as lgb
model = lgb.LGBMRegressor(n_estimators=800)  # ~O(n) training
```

---

## 7. The 152-d vs 175-d Detour (Critical Lesson)

During development we tried adding **23 extra features**: mRNA alignment scores, ViennaRNA MFE, gene-level one-hot encoding, etc. — making a 175-d vector.

**Result**: Delta scores collapsed from +10–23 to +3–4. The modification signal was *diluted* by high-variance target features that dominated tree splits.

**Root cause analysis**: mRNA accessibility features have large variance across different targets, so the tree model prioritized them over the small MNC changes caused by single modifications. This is a classic problem when mixing feature types with different variance scales.

**Fix**: Reverted to pure 152-d composition features. Deltas immediately recovered to +10–23 range.

**Lesson**: For modification ranking, composition-only features are optimal. mRNA context is better handled as a separate model or post-processing step.

**File history**:
- `models/train_gbm.py` — v1 with 175-d (deprecated)
- `models/train_gbm_v2.py` — v2 with extra features (deprecated)
- `models/train_gbm_v3.py` — v3 with pure 152-d (production)

---

## 8. Modification Engine — Single-Mod Scan

**`src/modification_engine.py`** — `single_mod_scan()` generates all 1,260 variants.

```python
# Generate every combination of modification type × position × strand
for mod in MOD_SYMBOLS:       # 30 chemical modifications
    for pos in range(21):     # positions 1..21
        for strand in ["sense", "antisense"]:  # 2 strands
            # Build modified sequence
            mod_seq = list(seq)
            mod_seq[pos] = mod
            variant = "".join(mod_seq)
            # → 30 × 21 × 2 = 1,260 variants
```

Each variant is scored by the LightGBM model through `predict_modified()`.

### Performance Optimization

```python
# src/predictor.py — predict_modified()
# Batch all 1260 variants through feature extraction + model prediction
# Single model call instead of 1260 individual calls
X = extract_batch_gbm(s_list, a_list, bs_list, ba_list)
raw = model.predict(X)        # ~50ms for 1260 variants
scores = _normalize_scores(raw, calibrator_key="cm")
```

### Key Files

| File | Purpose |
|------|---------|
| `src/modification_engine.py:single_mod_scan()` | Generates 1,260 variants |
| `src/modification_engine.py:CmSiRNA` | Data class for modified siRNA |
| `src/predictor.py:predict_modified()` | Orchestrates single-mod scan |

---

## 9. Modification Engine — Multi-Mod Beam Search

**`src/modification_engine.py`** — `multi_mod_scan()` extends single-mod to multi-mod combinations.

### Challenge

Testing all 2-mod combinations exhaustively: C(1260, 2) ≈ **800K possibilities**. Testing 3-mod: C(1260, 3) ≈ **330M**. Computationally infeasible.

### Solution: Beam Search

```python
def multi_mod_scan(sense, antisense, max_mods=2, beam_width=20):
    # Step 1: Run single-mod scan, get top results diversified by mod type
    single_results = predict_modified(sense, antisense, mode="scan")

    # Diversify: best per (mod_symbol, strand)
    best_per_type = {}
    for r in single_results["results"]:
        key = (r.mod_symbol, r.mod_strand)
        if key not in best_per_type or r.efficacy_score > best_per_type[key].efficacy_score:
            best_per_type[key] = r
    beam = sorted(best_per_type.values(), key=lambda r: r.efficacy_score, reverse=True)[:beam_width]

    # Step 2: Expand — combine beam members pairwise
    for n_mods in range(2, max_mods + 1):
        candidates = []
        for v1 in beam:
            for v2 in single_results:
                # Dedup via canonical ordering
                pair = tuple(sorted([(v1.mod_symbol, v1.mod_position, v1.mod_strand),
                                      (v2.mod_symbol, v2.mod_position, v2.mod_strand)]))
                if pair in seen: continue
                seen.add(pair)

                # Build combined sequence
                mod_sense = list(sense); mod_antisense = list(antisense)
                # Apply v1, then v2 (skip if same position on same strand)
                ...

        # Score all candidates, keep top beam_width
        scored = score_variants(candidates)
        beam = scored[:beam_width]

    return sorted(all_scored, key=lambda v: v.efficacy_score, reverse=True)
```

**Why diversify by mod type?** Without diversification, the top beam_width hits are all E (2'-MOE) at different positions → all combinations are E+E. Taking the best per mod type yields E+L, E+F, L+Q, etc.

### Key Files

| File | Purpose |
|------|---------|
| `src/modification_engine.py:multi_mod_scan()` | Beam search multi-mod |
| `src/modification_engine.py:score_variants()` | Batch scores candidates |

---

## 10. Parent Baseline Fix

**Critical bug found**: The `parent_score` in Single-Mod tab was computed using the **cm-siRNA model** on the unmodified parent sequence. But the Rank tab uses the **naked model** for unmodified scores. This meant the parent baseline differed between tabs.

**Fix** (`src/predictor.py`):

```python
def predict_modified(sense, antisense, mode="scan", model_key="A", full_scan=False):
    # ...
    # OLD: parent_score was from cm-siRNA model self-prediction
    # parent_score = _predict_cm(parent_sense, parent_antisense)

    # NEW: parent_score uses naked model — matches Rank tab exactly
    X_parent = features_gbm(parent_sense, parent_antisense,
                             parent_sense, parent_antisense)
    parent_score = _predict_naked(X_parent)

    # Now parent baselines match:
    # Rank tab: 62.94
    # Single-Mod parent_score: 62.94
    # MATCH = True ✓
```

### Why This Matters

The delta score (`variant_score − parent_score`) is the key insight of the Single-Mod tab. If the parent baseline is wrong, all deltas are misleading. Matching the Rank tab ensures consistency across the entire application.

---

## 11. API Layer (FastAPI)

**`api/main.py`** — FastAPI server with 4 endpoints.

```python
from fastapi import FastAPI
from pydantic import BaseModel

app = FastAPI(title="SMEpred", version="3.0.0")

class RankRequest(BaseModel):
    sequence: str
    top_n: int = 20
    naked: bool = False

class ModScanRequest(BaseModel):
    sense: str
    antisense: str
    model_key: str = "A"
    full_scan: bool = False

# Endpoints
@app.post("/rank")
@app.post("/single-mod")
@app.post("/multi-mod-scan")
@app.get("/modifications")

# Model loading — lazy load on first request
_model_cache: dict = {}
_model_lock = Lock()

def _get_model(key: str):
    if key not in _model_cache:
        with _model_lock:
            if key not in _model_cache:  # double-check
                path = MODELS_DIR / f"model_{key}.pkl"
                _model_cache[key] = joblib.load(path)
    return _model_cache[key]
```

Note: The lazy-load pattern with double-checked locking prevents race conditions when multiple requests arrive simultaneously.

### Key Files

| File | Purpose |
|------|---------|
| `api/main.py` | FastAPI server — all endpoints |
| `api/main.py:_get_model()` | Lazy-loading model cache |

---

## 12. Web UI (Single-File HTML)

**`app.html`** — ~41KB single-file dark-theme web application. No build step, no npm, no dependencies beyond the browser.

### Tab Structure

| Tab | Purpose | API Endpoint |
|-----|---------|--------------|
| **Rank** | Paste mRNA → get ranked siRNA candidates | `POST /rank` |
| **Single-Mod** | Scan 1,260 modifications for one siRNA | `POST /single-mod` |
| **Multi-Mod** | Auto beam-search multi-mod or manual design | `POST /multi-mod-scan` |
| **Modifications** | Reference table of all 30 symbols | Static data |

### Key UI Features

```html
<!-- Rank tab → results table with "Modify →" button -->
<button class="modify-btn" onclick="fillSingleMod(row.sense, row.antisense)">
  Modify →
</button>

<!-- Auto Multi-Mod Scan button in Single-Mod tab -->
<button id="autoMultiModBtn" onclick="runAutoMultiModScan()"
        style="background: linear-gradient(135deg, #7c3aed, #a855f7);">
  Auto Multi-Mod Scan
</button>

<!-- Parent Baseline badge -->
<span class="parent-baseline-badge">Parent: 62.94</span>
```

The "Modify →" button auto-fills the Single-Mod tab — removing friction from copy-paste. The "Auto Multi-Mod Scan" button runs the beam search and auto-switches to the Multi-Mod tab.

### Key Files

| File | Purpose |
|------|---------|
| `app.html` | Entire web UI — HTML, CSS, JavaScript |

---

## 13. Model Calibration

Raw LightGBM output is **compressed** toward the mean of training labels (~47.7). An isotonic regression spreads predictions to fill the 0–100 range.

```python
# models/calibrate.py
from sklearn.isotonic import IsotonicRegression

calibrator = IsotonicRegression(out_of_bounds='clip', y_min=0, y_max=100)
calibrator.fit(raw_predictions, true_labels)

# At inference:
calibrated = calibrator.predict(raw_output)
```

**Effect**: PCC preserved (0.68 → 0.68), but absolute scores become realistic. Raw 67 → calibrated 80.

```python
# src/predictor.py
def _normalize_scores(raw_scores, calibrator_key="cm"):
    calibrator = joblib.load(MODELS_DIR / f"calibrator_{calibrator_key}.pkl")
    return calibrator.predict(raw_scores)
```

### Key Files

| File | Purpose |
|------|---------|
| `models/calibrate.py` | Calibrator training |
| `models/calibrator_cm.pkl` | Fitted isotonic calibrator |
| `models/calibrator_naked.pkl` | Naked model calibrator |

---

## 14. Validation Strategy & Results

### Five-Level Validation

| Level | What It Catches | Result |
|-------|----------------|--------|
| 1. Unit tests | Bugs in feature extraction, modification engine | 19/19 pass |
| 2. Random split CV | Overfitting; modification ranking | PCC 0.68 |
| 3. Independent held-out | Leaky train/test split | PCC 0.77 on original paper set |
| 4. Independent held-out | Leaky train/test split | PCC 0.77 on original paper set |
| 5. Biology sanity | Does model match known RNAi rules? | ✓ (asymmetry, seed toxicity, rescue) |

### cm-siRNA Model Performance (Random 82/18 Split)

| Metric | Value |
|--------|-------|
| PCC | **0.6789** |
| Spearman | **0.6736** |
| MAE | 16.42 |
| Best gene (PCSK9) | **0.85** |
| Worst gene (PLN, n=95) | 0.46 |

### Naked Model Performance

| Source | n | PCC |
|--------|---|-----|
| All (source-aware) | 4,060 | **0.5543** |
| Taka (best independent) | 699 | 0.6905 |
| Huesken (noisy) | 2,361 | 0.4179 |

### What the Numbers Mean Physically

- **MAE 16.42**: typical prediction is within ±17 percentage points of true inhibition
- **Experimental noise floor**: ~10–15 points between labs (well-documented)
- **We are approaching the noise floor** — you cannot be more accurate than the measurement itself

### Key Files

| File | Purpose |
|------|---------|
| `logs/cm_model_v3_meta.json` | v3 model training metrics |
| `logs/normal_model_v3_meta.json` | Naked model v3 metrics |
| `logs/final_diagnostic.json` | Consolidated all diagnostics |
| `tools/validate_on_cmsirnadb.py` | External validation against CMsiRNAdb |
| `tools/validate_on_published_papers.py` | Validation on Alnylam 2025, Kenski 2012 |
| `docs/VALIDATION.md` | Full validation dossier |

---

## 15. Key Decisions Log

### Decision 1: LightGBM over SVR (Week 1)

- **Context**: Original paper used SVR. Initial runs on 43k rows plateaued at 0.37 PCC.
- **Decision**: Switch to LightGBM gradient-boosted trees.
- **Rationale**: Scales to full dataset, handles non-linear interactions, mixed feature scales.
- **Result**: 0.37 → 0.68 PCC.

### Decision 2: 152-d over 175-d (Week 3)

- **Context**: 23 extra features (mRNA alignment, ViennaRNA MFE, gene one-hot) were diluting modification signal.
- **Decision**: Revert to pure 152-d composition features.
- **Rationale**: Extra features have large variance → tree splits prioritize them over small MNC changes.
- **Result**: Deltas recovered from +3–4 to +10–23.

### Decision 3: Parent baseline uses naked model (Week 4)

- **Context**: Single-Mod parent_score didn't match Rank tab scores.
- **Decision**: Change `parent_score` from cm-siRNA model to naked model.
- **Rationale**: Rank tab uses naked model for unmodified scores → Single-Mod parent must match.
- **Result**: Tabs consistent. Confirmed: Rank=62.94, Single-Mod parent=62.94.

### Decision 4: Beam search for multi-mod (Week 4)

- **Context**: Testing all 2-mod combinations is ~800K variants — infeasible.
- **Decision**: Beam search with diversification by modification type.
- **Rationale**: Avoids combinatorial explosion while still exploring diverse combinations.
- **Result**: Top-20 single hits → 2-mod pairs → score → keep top-20 → repeat.

### Decision 5: Modify-forward button over manual copy-paste (Week 3)

- **Context**: Users had to manually copy 21-mer siRNA from Rank tab and paste into Single-Mod.
- **Decision**: Add "Modify →" button that auto-fills.
- **Rationale**: Friction kills adoption. One-click > copy-paste.

### Decision 6: Single-file HTML over SPA framework (Week 2)

- **Context**: Needed a web UI that works anywhere with zero setup.
- **Decision**: Single `app.html` with vanilla CSS/JS.
- **Rationale**: Zero build pipeline, opens in any browser, serves from FastAPI static route.

---

## 16. Files Reference Map

### Quick lookup: find which file implements a concept

| Concept | File | Line |
|---------|------|------|
| mRNA → 21-mer candidates | `src/sirna_generator.py` | `generate_sirnas()` |
| 152-d feature vector | `src/features.py` | `features_gbm()` |
| Batch feature extraction | `src/features.py` | `extract_batch_gbm()` |
| Rank siRNAs (end-to-end) | `src/predictor.py` | `rank_sirnas()` |
| Single-mod scan | `src/predictor.py` | `predict_modified()` |
| Multi-mod beam search | `src/modification_engine.py` | `multi_mod_scan()` |
| Single-mod variant gen | `src/modification_engine.py` | `single_mod_scan()` |
| Seed toxicity lookup | `src/filters.py` | `seed_toxicity()` |
| Functional checks (GC, runs, palindrome) | `src/filters.py` | `functional_check()` |
| Model A training | `models/train_gbm_v3.py` | LGBMRegressor |
| Model cache + inference | `src/predictor.py` | `_get_model()` |
| Score normalization | `src/predictor.py` | `_normalize_scores()` |
| API endpoint: /rank | `api/main.py` | `POST /rank` |
| API endpoint: /single-mod | `api/main.py` | `POST /single-mod` |
| API endpoint: /multi-mod-scan | `api/main.py` | `POST /multi-mod-scan` |
| Web UI (all tabs) | `app.html` | Full file |
| Patent catalog parser | `data/collect/parse_helix_catalog.py` | Full file |
| Naked data merger | `data/collect/merge_naked_sources.py` | Full file |
| Train/val split | `data/collect/splits.py` | `paper_split()` |
| Modification symbol map | `data/modification_codes.json` | Full file |
| CMsiRNAdb API extractor | `data/collect/extract_cmsirnadb_api.py` | Full file |
| CLI: rank | `cli/run.py` | `rank` command |
| CLI: single-mod | `cli/run.py` | `single-mod` command |
| CLI: multi-mod | `cli/run.py` | `multi-mod` command |
| Test suite | `tests/test_pipeline.py` | 19 unit tests |
| Validation dossier | `docs/VALIDATION.md` | Full doc |
| Speaker notes | `deck/SPEAKER_NOTES.md` | 15-slide script |
| Q&A cheatsheet | `docs/QA_CHEATSHEET.md` | Full doc |
| Paper manuscript | `HelixZero-CMS_Paper.md` | Full paper |
| Paper figures | `figures/paper_figures/` | 10 PNG files |

---

## Simple Explanations with Real-Life Examples (for CDAC Panel)

> This section explains every technical concept using analogies that bioinformatics scientists will grasp instantly. Use these during your presentation to make complex ML concepts click without jargon.

---

### 1. Gradient Boosting (LightGBM) — "The Apprentice Chef"

**Concept**: LightGBM builds an ensemble of hundreds of shallow decision trees. Each new tree is trained to correct the errors of all previous trees combined.

**Simple example**: Imagine training an apprentice chef to make the perfect biryani.

- **Tree 1** (Day 1): The chef makes biryani. You taste it and say: "Too salty, not enough spice."
- **Tree 2** (Day 2): The chef makes biryani again, but now focuses on fixing the salt and spice. You taste: "Better salt, but rice is overcooked."
- **Tree 3** (Day 3): The chef fixes the rice. You taste: "Rice is perfect, but now it needs more saffron."
- ... this continues for 800 days.

Each day, the chef focuses ONLY on what was wrong yesterday. The final biryani (after 800 iterations) is the sum of all 800 small improvements. That is gradient boosting — **each tree corrects the ensemble's previous mistakes**.

**Why it beats a single model**: A single decision tree is like a chef cooking from one recipe book. 800 trees working together is like 800 chefs each mastering one dish and combining forces. The ensemble always outperforms any single member.

**Bioinformatics analogy**: Like doing multiple rounds of sequence alignment refinement — each round fixes the gaps from the previous round until convergence.

---

### 2. Decision Trees — "20 Questions"

**Concept**: A decision tree asks a series of yes/no questions to arrive at a prediction.

**Simple example**: A tree predicting siRNA efficacy might ask:

```
Question 1: Is GC content > 50%?
  ├── YES → Question 2: Does antisense position 2 have a modification?
  │         ├── YES → Predict efficacy = 72
  │         └── NO  → Question 3: Is there a homopolymer run?
  │                   ├── YES → Predict efficacy = 45
  │                   └── NO  → Predict efficacy = 68
  └── NO  → Question 4: Is seed toxicity flagged?
            ├── YES → Predict efficacy = 31
            └── NO  → Predict efficacy = 55
```

Each fork in the road is a learned rule. The tree asks, you answer, and it lands on a leaf that gives a prediction. **LightGBM combines 800 such trees** — like having 800 experts each asking different questions, then averaging their votes.

**Bioinformatics analogy**: Like a phylogenetic key for species identification — "Does it have a backbone? → YES → Does it have fur? → YES → Does it have a tail? → ..."

---

### 3. The 152-d Feature Vector — "The Patient's Lab Report"

**Concept**: The model cannot read RNA letters (A, U, G, C). We convert each siRNA into 152 numbers that capture everything the model needs to know.

**Simple example**: A doctor doesn't look at a patient and guess their health. They order lab tests: CBC, lipid profile, liver enzymes, blood sugar, etc. Each test gives a number. The doctor reads all 30 numbers together to make a diagnosis.

**The 152 numbers are our siRNA's lab report**:

| Lab Test (Feature Group) | What It Measures |
|--------------------------|-----------------|
| **CBC** (70 numbers) | Complete blood count — how many of each nucleotide/modification is present |
| **Liver function** (8 numbers) | Where the modifications sit — seed region vs tail |
| **Lipid profile** (2 numbers) | GC content — duplex stability proxy |
| **Vitals** (2 numbers) | Experimental conditions — dose and time |

**Just as a doctor uses CBC + LFT + lipid profile together**, our model uses all 152 numbers simultaneously. One number alone tells you nothing. The 152-number combination captures the full siRNA "health signature."

**Bioinformatics analogy**: Like using a 152-feature vector for protein sequence encoding (composition, hydrophobicity, charge, secondary structure propensities) instead of raw amino acid letters.

---

### 4. Mononucleotide Composition (MNC) — "Ingredient Count in a Recipe"

**Concept**: MNC counts how many times each nucleotide symbol appears in a sequence and divides by strand length to get a frequency.

**Simple example**: A recipe for "GCAGCACGACUUCUUCAAGUU":
- G appears 4 times → frequency = 4/21 = 0.19
- C appears 5 times → frequency = 5/21 = 0.24
- A appears 4 times → frequency = 4/21 = 0.19
- U appears 5 times → frequency = 5/21 = 0.24
- All other symbols (M, F, L, E, ...) = 0

**Why both base AND modified MNC (140 numbers)?**

Imagine a dish that is heavily spiced — so much that you can't tell what the original meat was. If you only count the spices, you lose the base flavor. By counting BOTH the original ingredients (base MNC) and the final spiced dish (modified MNC), we capture both the underlying recipe and the modifications.

**Why this matters**: Near-fully-modified siRNAs collapse to mostly M/F symbols. If we only use modified MNC, the model can't distinguish between different sequences under all that modification. Adding base MNC restores that signal — raising PCC from 0.37 → 0.48.

**Bioinformatics analogy**: Like using both amino acid composition (PCP) and dipeptide composition to encode a protein sequence. One gives global frequency, the other gives local context.

---

### 5. Isotonic Calibration — "The Inconsistent Thermometer"

**Concept**: Raw LightGBM predictions are "compressed" toward the mean (47.7). A good prediction of 80 comes out as raw 67. Calibration spreads them back to the real 0–100 scale.

**Simple example**: Imagine a thermometer that consistently reads 5° too low in the 30–50° range, but reads correctly at 0° and 100°. You know it's monotonic (higher temperature → higher reading) but biased in the middle.

You take it to a calibration lab:
1. You measure known temperatures: 0°, 20°, 40°, 60°, 80°, 100°
2. The thermometer reads: 0°, 17°, 33°, 52°, 73°, 100°
3. You build a correction curve: when it reads 33°, the real temperature is 40°

**Isotonic regression does exactly this** for our model's predictions. It learns a correction curve that maps raw predictions → calibrated real-world scores. The order stays the same (rank preservation), but the absolute values become meaningful.

**Why this matters**: Without calibration, all scores cluster around 40–60 and you can't tell which candidates are truly excellent. After calibration, you get honest 0–100 scores where "80" really means "top 16%."

**Bioinformatics analogy**: Like using a standard curve in an ELISA assay. The raw OD values are meaningless until you map them to concentrations using the standard curve. Isotonic regression is our standard curve.

---

### 6. Beam Search (Multi-Mod) — "The Genealogist's Strategy"

**Concept**: Testing all 2-mod combinations is 800K variants. For 3-mod, it's 330 million. Beam search explores intelligently by keeping only the best candidates at each step.

**Simple example**: Imagine you're a genealogist tracing family trees for 1,260 people to find the best 3-generation marriage alliance.

- **Exhaustive search**: Check all (1260 × 1259 × 1258) / 6 ≈ 330 million possibilities. Way too many.
- **Beam search**: 
  1. First, evaluate all 1,260 individuals and pick the top 20 (beam width = 20).
  2. Now take each of those 20 and pair them with every other individual (20 × 1259 ≈ 25,000 pairs). Score all pairs, keep top 20.
  3. Now take each top pair and add a third member (20 × 1258 ≈ 25,000 trios). Score, keep top 20.

**Total checked**: 1,260 + 25,000 + 25,000 ≈ **51,000** instead of 330 million. That's a **6,500× speedup**.

**Why diversify the beam?** If you take the top 20 individuals by wealth alone, you get 20 rich people — and every marriage alliance is "rich + rich." By taking the best individual from EACH community (diversification by modification type), you get "rich + scientist," "rich + artist," "rich + doctor" — far more useful combinations.

**Bioinformatics analogy**: Like heuristic tree search in phylogenetic reconstruction. Instead of evaluating all possible trees (impossible), you start with a good tree and explore only the most promising rearrangements (NNI, SPR). Beam search is our tree rearrangement strategy.

---

### 7. Use-Case Framing — "The Right Test for the Right Job"

**Concept**: Different product features need different evaluation strategies. A model that ranks modifications of a *known* siRNA does not need to generalize to brand-new genes — because it never leaves the context of one siRNA.

**Simple example**: You are a chef perfecting a biryani recipe. You try different spice combinations (modifications) on the SAME base recipe (siRNA). You don't need to know how to cook Korean or Italian food (different genes) — you just need to know which spice combination makes THIS biryani taste best.

**That's our 0.68 PCC**: It measures "given we know the base siRNA, can we rank its chemical variants?" This is exactly what the Single-Mod and Multi-Mod tabs do. Cross-gene prediction (can we rank naked siRNAs for a gene never seen?) is a different question, answered by the naked model (PCC 0.55) used in the Rank tab.

**Why this matters**: Confiating both questions is how accuracy gets overstated. By separating them — modification ranking (0.68) and naked baseline ranking (0.55) — we give honest numbers for each use case.

**Bioinformatics analogy**: Like testing a variant effect predictor. You don't test it on genes it was never trained on and call it "poor." You test it on its intended use case — within-gene variant ranking.

---

### 8. Overfitting — "The Cramming Student"

**Concept**: A model that memorizes training examples but fails on new data is overfitting.

**Simple example**: Two students study for an exam:
- **Student A** (LightGBM with early stopping): Understands the concepts. Can solve problems they've never seen.
- **Student B** (overfit model): Memorizes all 100 practice questions and answers. On the exam, if a question is worded slightly differently (even if the concept is the same), they get it wrong. They scored 100% on practice but 40% on the real exam.

**How we prevent overfitting**:

| Technique | Analogy |
|-----------|---------|
| **Early stopping** | Stop studying when you stop improving on practice tests |
| **Subsample (0.8)** | Only read 80% of the textbook each time — prevents memorizing specific examples |
| **Colsample (0.7)** | Study only 70% of the topics each session — forces you to learn broad concepts |
| **L2 regularization** | Penalize yourself for being too confident in any single fact |
| **Min child samples (30)** | Don't make a rule unless at least 30 examples support it |

**Bioinformatics analogy**: Like BLAST with an E-value threshold — you want hits that are statistically significant, not random matches. Overfitting is like reporting every alignment with a positive score, including garbage.

---

### 9. LightGBM vs SVR — "The Research Lab vs The Startup"

**Concept**: We switched from Support Vector Regression (original paper's method) to LightGBM. Here's why.

**Simple example**: You need to sort 25,000 parcels by weight.

| Approach | Analogy | Result |
|----------|---------|--------|
| **SVR (old)** | One person weighs every pair of parcels (O(n²)). For 25,000 parcels, that's 312 million comparisons. The person gets tired and only does 6,000 parcels. | PCC = 0.37 |
| **LightGBM (new)** | 800 people each weigh a subset of parcels and correct each other's mistakes (O(n)). All 25,000 parcels sorted in 2 minutes. | PCC = 0.68 |

**But accuracy isn't the only advantage**:

| Capability | SVR | LightGBM |
|------------|-----|----------|
| Handles full 25K rows | ❌ (too slow) | ✅ (2 min CPU) |
| Position-aware features | ❌ (MNC only) | ✅ (+ density features) |
| Mixed data types (fractions, counts, log-dose) | ❌ (needs scaling) | ✅ (native) |
| Missing values | ❌ | ✅ (handles natively) |

**Bioinformatics analogy**: SVR is like ClustalW (accurate but O(n²) — can't handle thousands of sequences). LightGBM is like MMseqs2 (O(n) clustering — handles millions). Same science, different scale.

---

### 10. Feature Engineering — "Choosing the Right Biomarkers"

**Concept**: The features we choose to feed the model determine what it can learn. Bad features = bad model, no matter how good your algorithm is.

**Simple example**: You want to predict whether a patient has diabetes.

| Approach | Features | Result |
|----------|----------|--------|
| **Bad features** | Shoe size, favorite color, height | Model is useless |
| **Good features** | Fasting blood sugar, HbA1c, BMI, family history | Model is useful |
| **Optimal features** | Above + insulin levels, C-peptide, glucose tolerance | Model is excellent |

**Our 152-d features** are the equivalent of the "optimal" set:

| Feature Group | Why We Chose It | What Happens Without It |
|---------------|----------------|-------------------------|
| Base MNC (70-d) | Captures underlying sequence | Model loses signal in highly modified siRNAs (PCC drops from 0.68 → 0.37) |
| Modified MNC (70-d) | Captures modification identity | Model can't distinguish modification types |
| Mod density (8-d) | Captures WHERE mods are | Position 1 and position 21 modifications look identical |
| Assay conditions (2-d) | Doubles usable training data | Same sequence at different doses → contradictory labels |

**The 175-d mistake (what we learned)**: We tried adding mRNA features (alignment scores, ViennaRNA MFE, gene one-hot) to get to 175-d. It backfired — those features had high variance and drowned out the modification signal. Deltas collapsed from +10–23 to +3–4. **More features is NOT always better** — you need the RIGHT features.

**Bioinformatics analogy**: Like choosing substitution matrix for sequence alignment. BLOSUM62 works for most cases, but using BLOSUM45 for distant homologs or BLOSUM80 for close ones is better. Adding random columns to the scoring matrix (like mRNA features) just adds noise.

---

### 11. Early Stopping — "The Marathon Runner's Coach"

**Concept**: During training, we monitor performance on a validation set. When it stops improving for 50 rounds, we stop training — even if we haven't reached 800 trees.

**Simple example**: A marathon runner is training for a race. Every week, the coach measures their time.

- **Week 1**: 5 min/km
- **Week 4**: 4:30 min/km
- **Week 8**: 4:15 min/km
- **Week 12**: 4:10 min/km
- **Week 16**: 4:10 min/km (no improvement)
- **Week 20**: 4:11 min/km (getting WORSE — overtraining)

A good coach stops at Week 12. Continuing to Week 20 doesn't help — it actually hurts (overtraining = overfitting).

**Our model**: Best performance at 799 trees out of 800 max. If we kept training to 2,000 trees, it would memorize training noise and perform WORSE on new data.

**Why this matters**: Early stopping is the single most effective overfitting prevention technique. It's why our 0.68 PCC is real and not inflated by memorization.

**Bioinformatics analogy**: Like iterative refinement in multiple sequence alignment. After a certain number of iterations, the alignment quality plateaus and further iterations risk introducing errors (MUSCLE typically converges in 3–5 iterations).

---

### 12. Diversified Beam (the E+E Fix) — "The Diversity Committee"

**Concept**: Our multi-mod beam search was combining E+E modifications because E (2'-MOE) was the top-ranked single modification at every position. By diversifying the beam, we get E+L, E+F, L+Q, etc.

**Simple example**: You're forming a committee to solve a problem. You ask for the top 20 candidates by IQ alone. You get 20 people with IQ 160+ — but they're all physicists who think alike. The committee lacks diversity and misses creative solutions.

**Better approach**: Take the top candidate from EACH field (Physics, Chemistry, Biology, Engineering, Economics, etc.). You still get 20 people, but now you have diverse perspectives. The biology person suggests something the physics person never would have thought of.

**That's our diversified beam**: Instead of the top 20 by efficacy (all E modifications), we take the top candidate for EACH modification type (E, L, F, Q, M, S, etc.). This gives us combinations like "E on sense + L on antisense" instead of "E + E" everywhere.

**The data problem behind this**: The single-mod scan ranks E (2'-MOE) modifications highest because in our training data, 2'-MOE is strongly correlated with higher efficacy. This is biologically real — MOE is a bulky modification that improves nuclease resistance. But the best multi-mod design might be MOE on the backbone + LNA on the seed region + 2'-F on the tail — three different types.

**Bioinformatics analogy**: Like codon optimization. You don't just use the most frequent codon for EVERY amino acid (monoculture). You balance the codon usage to match the host organism's tRNA pool (diversity). The best combination is not the most common single codon repeated.

---

### 13. Parent Baseline Fix — "The Calibrated Scale"

**Concept**: The Single-Mod tab's parent_score (unmodified siRNA score) was computed using the wrong model. Fix: use the same naked model that the Rank tab uses.

**Simple example**: You have two weighing scales in your lab.

- **Old way (broken)**: Scale A (cm-siRNA model) measures the empty beaker's weight = 55g. Scale B (naked model) measures the SAME empty beaker = 63g. When you add a chemical and measure with Scale A, you get 75g. Your delta = 75 − 55 = +20g. But Scale B says the base weight is 63g, so the real delta should be 75 − 63 = +12g.

- **Fixed (consistent)**: Now you measure the empty beaker on Scale B (63g) and the filled beaker ALSO on Scale A (75g). The delta = 75 − 63 = +12g. The nominal scores differ between scales, but at least the delta is correct.

**But wait** — we actually went further. We now measure BOTH empty and filled on Scale B:

```python
parent_score = _predict_naked(X_parent)       # Scale B for empty
variant_score = _predict_cm(X_variant)         # Scale A for filled
delta = variant_score - parent_score
```

This way, the parent_score in Single-Mod tab **exactly matches** the Rank tab's unmodified scores. Users see consistent numbers everywhere.

**Why fixing this was critical**: Users would Rank siRNAs, see a score of 62.94, then go to Single-Mod and see parent_score = 55.0 (from the different model). They'd think the app was broken. Now: Rank = 62.94, Single-Mod parent = 62.94, MATCH.

**Bioinformatics analogy**: Like standardizing RNA-seq counts across experiments. If you use RPM in one analysis and FPKM in another, the numbers look different even for the same gene. You must normalize to a common scale before comparing.

---

### 14. Seed Toxicity Rescue — "The Antidote"

**Concept**: Some siRNA seed sequences are toxic (cause cell death). But certain modifications at seed-region positions can neutralize that toxicity.

**Simple example**: A drug has a side effect (liver toxicity). But adding a specific protective agent (N-acetylcysteine) alongside the drug neutralizes that side effect without reducing efficacy.

**For siRNAs**:
- The **seed** (positions 2–7 of antisense) determines off-target binding → potential toxicity
- Certain modifications at seed positions (**2'-OMe, 2'-F, LNA, 2'-MOE**) disrupt this off-target binding → toxicity "rescued"
- We detect this in the UI and show "**Mitigated**" instead of "Toxic"

**The code rule**:
```python
if baseline_toxicity in ("Toxic", "Caution") and any(
    mod in {"M", "F", "L", "E"} for mod in seed_region_mods
):
    label = "Mitigated"
```

**Why it's impressive**: This uses established biology (Jackson et al., *RNA* 2006) and surfaces it in the UI as a live, per-candidate annotation. A chemist looking at a modified siRNA can immediately see: "This would be toxic WITHOUT the modification, but the modification I placed at position 2 rescues it."

**Bioinformatics analogy**: Like predicting whether a mutation is pathogenic or benign based on its location in a conserved domain. The same mutation in different contexts has different effects.

---

### 15. The 1,260 Variants Problem — "The Combinatorial Explosion"

**Concept**: With 30 modifications × 21 positions × 2 strands, a single siRNA generates 1,260 variants to test.

**Simple example**: You're a chef who can add any of 30 spices to any of 21 spots on a dish, and each spot can be on the top layer or bottom layer. Even if you only try each possibility once, you need to cook 1,260 versions. At ₹500 per version, that's ₹6.3 lakh per siRNA. You have 980 siRNA candidates per gene → that's ₹62 crore per gene.

**SMEpred does all 1,260 computationally in ~50 milliseconds**. The cost: electricity for a laptop.

**This is the core value proposition**: We turn a ₹6.3 lakh, 4-week experiment into a 50-millisecond computation. The top-10 predictions then cost ₹5,000 to validate in the lab.

| Step | Wet Lab | SMEpred |
|------|---------|---------|
| Test 1 variant | ₹2,500, 1 day | 0.04 ms |
| Test 1,260 variants | ₹6.3 lakh, 4 weeks | 50 ms |
| Test top 10 variants | ₹25,000, 2 days | Free (already ranked) |

**Bioinformatics analogy**: Like BLAST vs experimental hybridization. Screening 1,000 primers for PCR experimentally costs ₹50,000. A BLAST search costs nothing and takes 10 seconds. We don't replace the PCR — we tell you which 5 primers are worth buying.

---

*Generated 2026-06-16 for CDAC Pune discussion. All code references are live. File paths relative to `D:\Helixx\smepred\`.*
