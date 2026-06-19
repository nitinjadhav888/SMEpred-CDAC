# HelixZero-CMS: Complete Workflow, Architecture & Tech Stack

## Table of Contents

1. [System Architecture Overview](#1-system-architecture-overview)
2. [End-to-End Data Flow](#2-end-to-end-data-flow)
3. [Module-by-Module Walkthrough](#3-module-by-module-walkthrough)
4. [Model Selection Rationale](#4-model-selection-rationale)
5. [Feature Engineering](#5-feature-engineering)
6. [Training Pipeline](#6-training-pipeline)
7. [Inference Pipeline](#7-inference-pipeline)
8. [Modification Engine Algorithms](#8-modification-engine-algorithms)
9. [API Layer](#9-api-layer)
10. [Frontend Architecture](#10-frontend-architecture)
11. [Pseudocode for Key Algorithms](#11-pseudocode-for-key-algorithms)
12. [Tech Stack Summary](#12-tech-stack-summary)
13. [File Structure Reference](#13-file-structure-reference)

---

## 1. System Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────┐
│                        USER INTERFACES                                   │
│                                                                          │
│  ┌──────────────┐    ┌──────────────┐    ┌───────────────────────────┐  │
│  │  Web Browser  │    │    CLI/Term   │    │  REST Client (curl, etc) │  │
│  │  (app.html)   │    │  (cli/run.py) │    │                          │  │
│  └──────┬───────┘    └──────┬───────┘    └─────────────┬─────────────┘  │
│         │                  │                          │                  │
└─────────┼──────────────────┼──────────────────────────┼──────────────────┘
          │          HTTP    │  CLI invocation          │  HTTP
          ▼                  ▼                          ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                           API LAYER                                      │
│                                                                          │
│              ┌──────────────────────────────────────┐                   │
│              │           FastAPI Server              │                   │
│              │  ┌──────────────────────────────┐    │                   │
│              │  │  /rank           → Predictor  │    │                   │
│              │  │  /single-mod     → Predictor  │    │                   │
│              │  │  /multi-mod      → Predictor  │    │                   │
│              │  │  /multi-mod-scan → ModEngine  │    │                   │
│              │  │  /modifications  → JSON file  │    │                   │
│              │  └──────────────────────────────┘    │                   │
│              └──────────────────────────────────────┘                   │
└─────────────────────────────┬───────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                        PREDICTION ENGINE                                 │
│                                                                          │
│  ┌────────────┐  ┌────────────┐  ┌────────────────┐  ┌──────────────┐  │
│  │  Predictor  │  │  Features  │  │  Modification  │  │   Filters    │  │
│  │  (predict  │──│ (extract   │──│  Engine (scan, │──│ (toxicity,   │  │
│  │  .py)      │  │  .py)      │  │  beam, etc)    │  │  func check) │  │
│  └──────┬─────┘  └────────────┘  └────────────────┘  └──────────────┘  │
│         │                                                               │
│         ▼                                                               │
│  ┌──────────────┐  ┌────────────────┐                                  │
│  │  Model B v4  │  │  Naked Model   │                                  │
│  │  (LightGBM   │  │  (LightGBM     │                                  │
│  │   .pkl)      │  │   .pkl)        │                                  │
│  └──────────────┘  └────────────────┘                                  │
│                                                                          │
└─────────────────────────────┬───────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────────────┐
│                       DATA & CONFIGURATION                               │
│                                                                          │
│  ┌─────────────────┐  ┌────────────────────┐  ┌─────────────────────┐  │
│  │  Trained Models  │  │  Calibrators       │  │  Toxicity Table     │  │
│  │  (model_*.pkl)   │  │  (calibrator_*)    │  │  (cell_viability)   │  │
│  ├─────────────────┤  ├────────────────────┤  ├─────────────────────┤  │
│  │  model_b.pkl    │  │  calibrator_naked  │  │  4,097 seed → %     │  │
│  │  model_normal   │  │  .pkl              │  │  viability mappings │  │
│  │  .pkl           │  │                    │  │                     │  │
│  └─────────────────┘  └────────────────────┘  └─────────────────────┘  │
│                                                                          │
└─────────────────────────────────────────────────────────────────────────┘
```

---

## 2. End-to-End Data Flow

### Workflow 1: Rank siRNA Candidates (unmodified)

```
User Input (mRNA sequence or FASTA file)
    │
    ▼
┌────────────────────────┐
│     parser.py          │  ← Normalize T→U, validate A/U/G/C only
│  load_sequence()       │     Accepts inline, FASTA, or file path
└───────────┬────────────┘
            │ 21+ nt RNA sequence
            ▼
┌────────────────────────┐
│  sirna_generator.py    │  ← Sliding window: position → (sense, antisense)
│  generate_candidates() │      Window size = 21, step = 1
└───────────┬────────────┘
            │ List[SiRNACandidate] (N − 20 candidates)
            ▼
┌────────────────────────┐
│  features.py           │  ← Extract V4 features (214-d per candidate)
│  extract_batch_v4()    │     Position one-hot + TNC + GC content
└───────────┬────────────┘
            │ np.array shape: (n_candidates, 214)
            ▼
┌────────────────────────┐
│  predictor.py          │  ← Load Naked LightGBM model
│  _predict_naked()      │     Pad source one-hot, predict, normalize
└───────────┬────────────┘
            │ Raw scores → [0, 100] scaled
            ▼
┌────────────────────────┐
│  filters.py            │  ← Seed toxicity lookup (4,097 seeds)
│  annotate_candidates() │     Functional checks (GC, runs, palindrome)
└───────────┬────────────┘
            │ RankedSiRNA list, sorted by efficacy (high → low)
            ▼
    Return to API/CLI/UI
```

### Workflow 2: Single-Mod Scan (modified)

```
User Input (sense, antisense, options)
    │
    ▼
┌──────────────────────────────┐
│  modification_engine.py      │  ← Generate 1,260 single-mod variants
│  single_mod_scan()           │     30 mod symbols × 21 positions × 2 strands
└──────────────┬───────────────┘
               │ List[CmSiRNA]
               ▼
┌──────────────────────────────┐
│  features.py                 │  ← Extract positional features (1,467-d)
│  extract_positional_features │     Per-variant: 33 flags/pos × 42 pos
│  _batch()                    │     + 70 global counts + 16 summary + 1 log_conc
└──────────────┬───────────────┘
               │ np.array shape: (1260, 1467)
               ▼
┌──────────────────────────────┐
│  predictor.py                │  ← Load Model B LightGBM (1115 trees)
│  _get_model("B") → predict() │     Identity normalization → [0, 100]
└──────────────┬───────────────┘
               │ Raw scores → scores[0, 100]
               ▼
┌──────────────────────────────┐
│  filters.py                  │  ← Modification-aware seed toxicity
│  toxicity_for_modified()     │     "Mitigated" if rescue mod present
└──────────────┬───────────────┘
               │ RankedCmSiRNA list, sorted by efficacy
               ▼
     Return to API: top-N variants + parent score + toxicity
```

### Workflow 3: Multi-Mod Beam Search

```
User Input (sense, antisense, max_mods=auto, beam_width)
    │
    ▼
┌─────────────────────────────────────┐
│  Step 1: Run single-mod scan        │
│  (same as Workflow 2, full 1,260)   │
└──────────────┬──────────────────────┘
               │ List[RankedCmSiRNA] + parent_score
               ▼
┌─────────────────────────────────────┐
│  Step 2: Diversify initial beam     │
│  Round-robin through all 30 mod     │
│  types, pick top from each          │
└──────────────┬──────────────────────┘
               │ beam = [variant_1, variant_2, ..., variant_K]
               ▼
┌─────────────────────────────────────┐
│  Step 3: Score all beam candidates  │
│  (extract features + predict)       │
└──────────────┬──────────────────────┘
               │ Scored + sorted by efficacy
               ▼
┌─────────────────────────────────────┐
│  Step 4: Expand (iterative)         │
│  for n = 2 to max_mods:             │
│    for each beam_candidate:         │
│      for each single_result:        │
│        combine if no position       │
│        conflict                     │
│    score new candidates             │
│    keep top-K as new beam           │
│    if score plateau → EARLY STOP   │
└──────────────┬──────────────────────┘
               │ All scored variants (1-mod through N-mod)
               ▼
┌─────────────────────────────────────┐
│  Step 5: Compute composite score    │
│  for each variant:                  │
│    composite = 0.5×efficacy_norm    │
│               + 0.12×nuclease       │
│               + 0.12×immunogenicity  │
│               + 0.12×RISC_loading    │
│               + 0.14×thermo_stability│
└──────────────┬──────────────────────┘
               │ Variants sorted by composite score
               ▼
┌─────────────────────────────────────┐
│  Step 6: Return top-N variants      │
│  with full biophysics breakdown     │
└─────────────────────────────────────┘
```

---

## 3. Module-by-Module Walkthrough

### 3.1 `src/parser.py` — Input Parsing

**Responsibilities:**
- Accept mRNA sequence in plain text or FASTA format
- Normalize DNA→RNA (T→U)
- Validate only A/U/G/C characters
- Raise on empty or too-short sequences (<21 nt)

**Key function:** `load_sequence(source)` returns a single RNA string.

**Design decision:** Single-file input keeps the API simple. Multi-gene inputs are handled by the user running the app on one gene at a time.

### 3.2 `src/sirna_generator.py` — Candidate Generation

**Responsibilities:**
- Slide a 21-nt window across the mRNA
- For each window, generate:
  - **Sense strand:** identical to mRNA segment (T→U)
  - **Antisense strand:** reverse complement of sense
- Return `List[SiRNACandidate]`

**Key function:** `generate_candidates(seq)` returns candidates with 0-based position, sense, antisense.

**Edge cases:**
- Short sequences (<21 nt): zero candidates
- Non-A/U/G/C characters: flagged by parser upstream

### 3.3 `src/features.py` — Feature Extraction

#### Model B v4 Features (1,467-d)

Used for modified siRNA prediction.

```
Per position (21 pos × 2 strands × 33 flags):
  ├── 31 × is_mod_type_X     (binary: is this modification type at this position?)
  ├── 1 × is_canonical       (binary: no modification here)
  └── 1 × is_modified        (binary: any modification here)

Per strand (×2):
  ├── 31 × global_count_X    (how many of each mod type on this strand)
  ├── 5 × summary_stats     (frac_mod, seed_2F_frac, seed_2OMe_frac,
  │                          cleave_2F_count, cleave_2OMe_count, cleave_LNA_count)
  ├── 1 × gc_content         (of base sequence)
  ├── 1 × term_5_ps          (is position 1 PS?)
  └── 1 × term_3_ps          (is position 21 PS?)

Global:
  └── 1 × log_conc           (log concentration + 1, default 10 nM)

Total: 33×42 + 31×2 + 8×2 + 1 = 1,467 features
```

#### Naked Model V4 Features (214-d)

Used for unmodified siRNA ranking.

```
Sense strand:
  ├── 84 = 21 pos × 4 one-hot (A/U/G/C)
  ├── 64 = tri-nucleotide composition (4×4×4 / 19, normalized)
  └──  1 = GC content

Antisense strand:
  ├── 64 = tri-nucleotide composition (4×4×4 / 19, normalized)
  └──  1 = GC content

Total: 84 + 64 + 1 + 64 + 1 = 214 features
```

**Why two feature extractors?**
- The naked model does not use modification information (there are none)
- The modified model needs position-aware modification flags to distinguish "2'-OMe at position 5" from "2'-OMe at position 15"

### 3.4 `src/predictor.py` — Prediction Orchestration

**Responsibilities:**
- Load and cache LightGBM models (lazy loading)
- Route workflow: rank naked siRNAs or predict modified variants
- Normalize scores to 0–100
- Coordinate feature extraction, model inference, and annotation

**Key functions:**
- `rank_by_naked_score(source)` — rank unmodified candidates
- `predict_modified(sense, antisense, mode, ...)` — predict modified variants
- `_get_model(key)` — lazy-load models from disk
- `_normalize_scores(raw, mode)` — identity/clip/rescale normalization

### 3.5 `src/modification_engine.py` — Chemical Modification Generator

**Responsibilities:**
- Generate single-modification scan variants (1,260)
- Generate user-specified multi-modification designs
- Run beam-search multi-mod scan (combinatorial optimization)
- Track conflicts (same position modified twice)

**Three modes:**

| Mode | Function | Output |
|------|----------|--------|
| Single-mod scan | `single_mod_scan()` | 1,260 variants |
| MultiModGen | `multimod_gen()` | 1 custom variant |
| Beam search | `multi_mod_scan()` | K candidates (K grows with expansion) |

### 3.6 `src/filters.py` — Safety & Quality Filters

**Responsibilities:**
- Seed toxicity lookup (Janas et al. 2018 table)
- Modification-aware toxicity (rescue mod override)
- Functional checks (GC content, homopolymers, palindromes)

**Data dependency:** `data/oligoformer/cell_viability.tsv` — 4,097 seed→viability mappings.

---

## 4. Model Selection Rationale

### Why LightGBM Over Alternatives?

| Model | Pros | Cons | Decision |
|-------|------|------|----------|
| **SVR (RBF kernel)** | — Good on small data | — O(n²) scaling → intractable at 83k rows | ❌ Rejected |
| **Random Forest** | — Good baseline | — Worse accuracy than GBDT on tabular data | ❌ Rejected |
| **XGBoost** | — Industry standard | — Level-wise growth slower for high-dim features | ❌ Rejected |
| **LightGBM** | — Leaf-wise, histogram-based, fast training on 83k × 1,467; native categorical support; built-in CV | — Can overfit on very small data | ✅ **Selected** |
| **Neural Network** | — High capacity | — Requires more data; less interpretable; harder to deploy | ❌ Rejected |

### LightGBM Advantages for This Problem

1. **Scalability**: Leaf-wise tree growth with histogram-based splitting trains in minutes on 83,535 × 1,467 data
2. **Feature importance**: Built-in permutation importance helps identify which position/mod combinations matter most
3. **Regularization**: `num_leaves`, `min_data_in_leaf`, `lambda_l1`, `lambda_l2` control overfitting
4. **Handling heterogeneous data**: Robust to the label variance from mixed experimental conditions (different doses, cell types)
5. **Calibration**: Gradient-boosted trees produce reasonably calibrated probabilities/scores without post-processing

### Model Specifications

| Parameter | Naked Model | Model B v4 |
|-----------|-------------|------------|
| **Algorithm** | LightGBM | LightGBM |
| **Training rows** | 4,060 | 83,535 |
| **Features** | 214-d (source-aware) | 1,467-d |
| **Trees** | 1,000 | 1,115 |
| **Objective** | regression_l2 | regression_l2 |
| **Learning rate** | 0.05 | 0.05 |
| **Test PCC** | 0.48 (Takayuki) | 0.650 (hetero_val_303) |
| **Calibration** | Isotonic (calibrator_naked.pkl) | Identity (no calibrator) |

### Why No Neural Network?

- 83,535 rows is moderate for deep learning
- LightGBM matches or exceeds NN performance on tabular data with fewer hyperparameters
- Deployment is simpler: pickle file vs PyTorch/TF serving
- Interpretability is critical for scientific acceptance

---

## 5. Feature Engineering

### Evolution of Feature Space

| Version | Features | Description | Status |
|---------|----------|-------------|--------|
| SMEpred (Dar 2016) | 70-d | MNC only | Historical reference |
| Model A (v1–v2) | 140-d | MNC base + MNC modified | ❌ Deleted |
| Model B v1–v3 | 595-d | 140 base + 13 feature groups | ❓ Legacy |
| Model B v4 (current) | **1,467-d** | Position flags (31 mods × 21 pos × 2 strands) | ✅ **Active** |
| Naked V4 (current) | **214-d** | One-hot + TNC + GC | ✅ **Active** |

### Why Position-Aware Features?

The key insight: modification position matters biologically. A 2'-F at position 5 of the antisense strand affects seed-region stability differently than a 2'-F at position 20. The 31 per-position binary flags capture this spcificity.

```
Old approach (MNC/DNC):
  "2'-OMe count = 5"    → loses which positions have 2'-OMe

New approach (positional):
  "is_2OMe at pos 5 AS" → captures exact modification pattern
  "is_2OMe at pos 20 AS" → separately flagged
```

### Feature Extraction Time

| Extractor | Time per batch (1,260 variants) |
|-----------|--------------------------------|
| `extract_positional_features_batch` | ~50 ms |
| `extract_batch_v4` | ~10 ms |

Both are optimized with vectorized NumPy operations (no per-row Python loops for the heavy computation).

---

## 6. Training Pipeline

### Data Sources (Model B v4)

| Source | Rows | Description |
|--------|------|-------------|
| Position-aware dataset | 55,730 | Sequences with per-position mod annotations |
| Hetero_train (SMEpred) | 23,187 | Modified siRNA efficacy benchmark |
| CMsiRNAdb | 4,618 | External patent-derived dataset |
| **Total** | **83,535** | |

### Training Script

Located at `models/train_model_b_v4.py`.

**Algorithm:**
```
1. Load all 3 data sources
2. Extract 1,467-d features for all rows
3. Merge into single training array
4. Split: every 10th row starting at index 5 → validation (≈8,353)
5. Train LightGBM:
   - objective: regression_l2
   - num_leaves: 31
   - learning_rate: 0.05
   - feature_fraction: 0.8
   - early_stopping_rounds: 50
   - max_rounds: 10,000
6. Save: model_b.pkl + model_b_meta.json
```

### Training Data for Naked Model

- **Source**: 4,060 unmodified siRNA sequences from 4 independent sources
  - Huesken (Huesken et al., 2005): 2,182 sequences
  - Takayuki (Ui-Tei et al.): 594 sequences
  - Mix (multiple labs): 539 sequences
- **Source encoding**: Each training example is one-hot encoded with its source ID
- **Calibration**: Isotonic regression calibrator maps raw → 0–100

---

## 7. Inference Pipeline

### Score Normalization

```
Raw output from LightGBM
    │
    ▼
Normalize mode:
  ┌─────────────────────────────────────────────────────────────┐
  │  "identity"  → clip(raw, 0, 100)     [Model B — default]   │
  │  "clip"      → clip(raw, 0, 100)     [safety net]          │
  │  "rescale"   → clip(raw/113.8*100)   [legacy]              │
  │  "calibrate" → isotonic_transform    [naked model only]     │
  └─────────────────────────────────────────────────────────────┘
```

### Label Thresholds

| Score Range | Label |
|-------------|-------|
| ≥ 80 | Very High |
| 70–79 | High |
| 55–69 | Moderate |
| < 55 | Low |

Thresholds are based on training-data percentiles (P50 = 48, P75 = 72, P84 = 80, P94 = 90).

---

## 8. Modification Engine Algorithms

### 8.1 Single-Mod Scan

```
Algorithm SingleModScan(sense, antisense):
    variants = []
    for symbol in MODIFICATION_SYMBOLS (30):
        for position in 1..21:
            variant = apply_mod(sense, position, symbol)
            variants.add(variant)
        for position in 1..21:
            variant = apply_mod(antisense, position, symbol)
            variants.add(variant)
    return variants  // 1,260 total
```

### 8.2 MultiModGen

```
Algorithm MultiModGen(sense, antisense, sense_mods, sense_positions,
                       antisense_mods, antisense_positions):
    mod_sense = list(sense)
    mod_antisense = list(antisense)

    // Parse "F,,M" and "2,5,,10,12" → [(F, [2,5]), (M, [10,12])]
    sense_groups = parse_multimod_input(sense_mods, sense_positions)
    for (symbol, positions) in sense_groups:
        for pos in positions:
            mod_sense[pos - 1] = symbol

    antisense_groups = parse_multimod_input(antisense_mods, antisense_positions)
    for (symbol, positions) in antisense_groups:
        for pos in positions:
            mod_antisense[pos - 1] = symbol

    return CmSiRNA(modified sense, modified antisense, ...)
```

### 8.3 Beam Search Multi-Mod

```
Algorithm MultiModBeamSearch(sense, antisense):
    // Step 1: Single-mod scan
    single_results = predict_modified(sense, antisense, scan=true)
    parent_score = single_results.parent_score

    // Step 2: Diversify initial beam
    per_mod = group_by(single_results, mod_symbol)
    beam = round_robin_select(per_mod, beam_width)

    // Step 3: Initial scoring
    beam = predict_and_score(beam)

    all_results = copy(beam)

    // Step 4: Iterative expansion
    for n_mods = 2 to MAX_POSSIBLE (21):
        candidates = []

        for v1 in beam:
            for v2 in single_results:
                conflict = (v2.position already modified in v1)
                if conflict: continue

                combined = apply(v1.modifications + v2.modification)
                candidates.add(combined)

        scored = predict_and_score(candidates)
        beam = top_K(scored, beam_width)
        all_results.extend(scored)

        // Early stop: if top-3 composite scores plateau
        if score_improvement < threshold:
            break

    // Step 5: Compute composite scores
    for variant in all_results:
        variant.composite = compute_composite(variant)

    // Step 6: Sort by composite and deduplicate
    return sort(all_results, by=composite, descending)
```

### Beam Search Complexity

```
Let B = beam_width (~20), S = single_results (1,260)

Per expansion round:
    candidates = B × S = 25,200
    predictions = 25,200

Total for 5 rounds:         1260 + 5 × 25,200 = 127,260 predictions
Total for early-stop (3):    1260 + 3 × 25,200 = 76,860 predictions
Time at ~10ms/1260:          ≈ 0.6–1.0 seconds per round
Total time (5 rounds):       ≈ 3–5 seconds
```

---

## 9. API Layer

### Endpoint Reference

| Method | Path | Purpose |
|--------|------|---------|
| GET | `/` | Serve app.html frontend |
| GET | `/app.html` | Serve frontend (alternative path) |
| POST | `/rank` | Rank unmodified siRNA candidates |
| POST | `/rank/upload` | Rank from FASTA file upload |
| POST | `/single-mod` | Single-mod scan (1,260 variants) |
| POST | `/multi-mod` | Custom multi-modification design |
| POST | `/multi-mod-scan` | Auto beam search (legacy, max_mods=2–3) |
| POST | `/multi-mod-from-single` | Auto beam search (full, all mod counts) |
| GET | `/modifications` | List 30 supported modification symbols |

### Request/Response Examples

**POST /single-mod**
```json
{
    "sense": "GCAGCACGACUUCUUCAAGUU",
    "antisense": "CUUGAAGAAGUCGUGCUGCUU",
    "model": "B",
    "top_n": 50,
    "full_scan": true
}
```
```json
{
    "parent_sense": "GCAGCACGACUUCUUCAAGUU",
    "parent_antisense": "...",
    "parent_score": 62.05,
    "model": "B",
    "total_variants": 1260,
    "results": [
        {
            "rank": 1, "sense": "...", "antisense": "...",
            "mod_symbol": "8", "mod_position": 1, "mod_strand": "sense",
            "efficacy_score": 88.5, "delta_score": 26.45,
            "efficacy_label": "Very High",
            "toxicity_score": null, "toxicity_label": "Unknown",
            "toxicity_note": ""
        }
    ]
}
```

### Error Handling

All endpoints return structured HTTP errors:
- `422` — Validation error (bad input)
- `503` — Model file not found (run training first)
- `500` — Unexpected internal error

---

## 10. Frontend Architecture

### Design

- **Single-file HTML** (app.html) — no build step, no npm, no framework
- Dark theme with CSS variables for consistency
- Vanilla JavaScript (ES6 async/await)
- Communicates with FastAPI backend via `fetch()`

### Tabs

| Tab | Function | API Endpoint |
|-----|----------|-------------|
| Rank siRNAs | Paste mRNA, rank unmodified candidates | POST /rank |
| Single-Mod Scan | Choose siRNA, scan 1,260 mod variants | POST /single-mod |
| Multi-Mod Design | Custom modifications at specific positions | POST /multi-mod |
| Modifications | Reference guide for 30 mod symbols | GET /modifications |

### Key UI Features

- **Multi-Mod from Rank tab**: One-click beam search from any ranked candidate → sends to `/multi-mod-from-single`
- **Single-Mod → Multi-Mod pipeline**: Load single-mod results into beam search with seed variant
- **Modification-aware toxicity**: "Mitigated" label when rescue mod present
- **Score bars**: Visual efficacy representation with color coding
- **Unified model label**: "HelixZero (unified)" displayed throughout

### Upcoming: Biophysics Display

The biophysics parameter balance will be displayed in the Multi-Mod results as a compact breakdown:

```
  Composite: ████████░░ 82.3
  ─────────────────────────
  Efficacy:  ██████████ 100.0  (raw)
  Nuclease:  ████████░░  75.0
  Immuno:    ████████░░  82.0
  RISC:      █████████░  90.0
  Thermo:    ███████░░░  71.0
```

---

## 11. Pseudocode for Key Algorithms

### 11.1 Positional Feature Extraction

```
FUNCTION extract_positional_features_batch(sense_list, antisense_list,
                                            base_sense_list, base_antisense_list,
                                            conc_list):

    FEATURE_DIM = 33 flags × 21 pos × 2 strands    // 1,386
                 + 31 counts × 2 strands            // 62
                 + 8 summary × 2 strands            // 16
                 + 1 log_conc                       // 1
                 = 1,467

    result = zeros(len(sense_list), FEATURE_DIM)

    FOR i = 0 TO len(sense_list):
        feats = []

        // ── Per-position flags for each strand ──
        FOR strand in [sense, antisense]:
            FOR pos = 1 TO 21:
                nt = strand[pos]
                base = base_strand[pos]
                is_mod = (nt != base)

                // 31 mod-type flags
                FOR each mod_type in MOD_CHAR_MAP:
                    feats.append(nt == mod_type)
                // 1 canonical flag
                feats.append(not is_mod)
                // 1 modified flag
                feats.append(is_mod)

        // ── Global per-strand summary ──
        FOR strand in [sense, antisense]:
            counts = count_mods(strand)
            feats.extend(counts_array)

            feats.extend([
                frac_mod / 21,
                seed_2F_frac,
                seed_2OMe_frac,
                cleave_2F, cleave_2OMe, cleave_LNA,
                gc_content,
                term_5_ps, term_3_ps
            ])

        // ── Log concentration ──
        feats.append(log(conc + 1))

        result[i] = array(feats)

    RETURN result
```

### 11.2 Beam Search with Early Stop

```
FUNCTION multi_mod_scan(sense, antisense, beam_width=20):

    // Phase 1: single-mod results
    single_out = predict_modified(sense, antisense, full_scan=True)
    parent_score = single_out.parent_score
    single_results = single_out.results

    // Phase 2: diversify initial beam
    per_mod_type = {}
    FOR r IN single_results:
        per_mod_type[r.mod_symbol].append(r)

    beam = []
    FOR rank = 0 TO max_len(per_mod_type):
        FOR sym IN sorted(per_mod_type.keys()):
            IF len(beam) >= beam_width: BREAK
            IF rank < len(per_mod_type[sym]):
                beam.append(per_mod_type[sym][rank])

    // Phase 3: iterative expansion with early stop
    all_variants = beam.copy()
    prev_best_score = -INF

    FOR n_mods = 2 TO 21:  // upper bound
        candidates = []

        FOR b IN beam:
            FOR s IN single_results:
                IF position_conflict(b.modifications, s.modification):
                    CONTINUE
                combined = merge_modifications(b, s)
                candidates.add(combined)

        scored = predict_batch(candidates)
        beam = top_K(scored, beam_width)
        all_variants.extend(scored)

        // Early stop condition
        current_best = mean(top_3_composites(scored))
        IF current_best - prev_best_score < 1.0:
            BREAK  // plateau detected
        prev_best_score = current_best

    // Phase 4: composite scoring
    FOR v IN all_variants:
        v.composite = WEIGHTS × [
            v.efficacy_normalized,
            nuclease_score(v),
            immuno_score(v),
            risc_score(v),
            thermo_score(v)
        ]

    RETURN sort(all_variants, by=composite, descending)
```

### 11.3 Nuclease Resistance Score

```
FUNCTION nuclease_score(variant):
    score = 0

    // PS protection at termini
    IF PS_at(variant.sense, [1, 20, 21])
        OR PS_at(variant.antisense, [1, 20, 21]):
        score += 25

    // Multiple PS linkages
    IF total_PS_count(variant) >= 3:
        score += 20

    // 2'-modification density
    pct_2mod = percent_modified(variant, modifiers=["F","M","L","E"])
    IF pct_2mod >= 0.3: score += 15
    IF pct_2mod >= 0.5: score += 10

    // LNA at 3' terminus
    IF LNA_at(variant.sense, 21) OR LNA_at(variant.antisense, 21):
        score += 10

    // No exposed termini
    IF all_termini_modified(variant):
        score += 10

    RETURN min(score, 100)
```

### 11.4 Immunogenicity Score

```
FUNCTION immuno_score(variant):
    score = 50  // start at neutral

    antisense_U_pos = get_positions(variant.antisense, base="U")
    sense_U_pos = get_positions(variant.sense, base="U")

    // Every U modified with 2'-OMe → strong TLR suppression
    IF all_have_mod(antisense_U_pos, ["M", "W"]):
        score += 30
    IF all_have_mod(sense_U_pos, ["M", "W"]):
        score += 20

    // Any pseudouridine → immune silent
    IF any_have_mod(variant.antisense + variant.sense, ["W"]):
        score += 15

    // Partial U modification
    pct_U_mod = percent_U_modified(variant, ["M","W","E"])
    IF pct_U_mod >= 0.5:
        score += 10

    // PS backbone reduces immune recognition
    IF total_PS_count(variant) > 0:
        score += 10

    // Penalties for unmodified Uridines
    FOR strand in [sense, antisense]:
        FOR each unmodified U:
            score -= 5

    // Penalties for GU-rich motifs
    FOR each GU-rich stretch (≥4 nt, unmodified):
        score -= 10

    RETURN clamp(score, 0, 100)
```

---

## 12. Tech Stack Summary

| Layer | Technology | Version | Purpose |
|-------|-----------|---------|---------|
| **Language** | Python | ≥3.10 | Primary development language |
| **ML framework** | LightGBM | ≥4.0 | Gradient-boosted trees (Model B, Naked) |
| **Web server** | FastAPI | ≥0.104 | REST API backend |
| **ASGI server** | Uvicorn | ≥0.24 | Production server |
| **CLI framework** | Click | ≥8.1 | Command-line interface |
| **Data processing** | NumPy | ≥1.24 | Feature extraction (vectorized) |
| **Data processing** | Pandas | ≥2.0 | CSV/TSV data handling |
| **Model I/O** | Joblib | ≥1.3 | Model serialization |
| **SciKit-Learn** | sklearn | ≥1.3 | Isotonic calibration |
| **Frontend** | Vanilla HTML/CSS/JS | — | Single-file web UI |
| **API docs** | Swagger (auto) | — | Built-in at /docs |
| **Deployment** | Railway | — | Cloud deployment (railway.json) |

### Dependencies (requirements.txt)

```
biopython>=1.81
numpy>=1.24
pandas>=2.0
scikit-learn>=1.3
lightgbm>=4.0
scipy>=1.11
fastapi>=0.104
uvicorn[standard]>=0.24
click>=8.1
joblib>=1.3
requests>=2.31
python-multipart>=0.0.6
```

---

## 13. File Structure Reference

```
HelixZero-CMS/
├── smepred/
│   ├── api/
│   │   └── main.py              ★ FastAPI server (8 endpoints)
│   ├── cli/
│   │   └── run.py               ★ Click CLI (rank, single-mod, multi-mod)
│   ├── data/
│   │   ├── modification_codes.json   ★ 30 mod symbols + alias rules
│   │   ├── hetero_train_2728.csv     SMEpred training set (23,187 rows)
│   │   ├── hetero_val_303.csv        SMEpred validation set (2,576 rows)
│   │   ├── cmsirnadb_full.csv        CMsiRNAdb external data (4,618 rows)
│   │   ├── normal_siRNA.csv          Naked siRNA dataset
│   │   ├── mrna_sequences.json       mRNA reference sequences
│   │   └── oligoformer/
│   │       ├── cell_viability.tsv    ★ Seed toxicity table (4,097 seeds)
│   │       ├── Hu.csv, Mix.csv, Taka.csv   Naked model training sources
│   ├── models/
│   │   ├── model_b.pkl           ★ Model B v4 (LightGBM, 1,115 trees)
│   │   ├── model_b_meta.json     Training metadata (date, rows, PCC)
│   │   ├── model_normal.pkl      ★ Naked model (LightGBM)
│   │   ├── calibrator_naked.pkl  Isotonic calibrator
│   │   ├── train_model_b_v4.py   Model B training script
│   │   └── backup/               Legacy model backups
│   ├── src/
│   │   ├── features.py           ★ Feature extraction (1,467-d + 214-d)
│   │   ├── predictor.py          ★ Prediction orchestration
│   │   ├── modification_engine.py ★ Modification generation + beam search
│   │   ├── filters.py            Seed toxicity + functional filters
│   │   ├── parser.py             FASTA/sequence parsing
│   │   ├── sirna_generator.py    21-mer sliding window
│   │   └── mrna_features.py      mRNA target-site features
│   ├── tests/
│   │   └── test_pipeline.py      ★ 18 unit tests
│   ├── app.html                  ★ Single-file web UI
│   ├── requirements.txt          Python dependencies
│   ├── pyproject.toml             Package config
│   ├── EXPLANATION.md            Parameter explanations
│   └── README.md                 Quick-start guide
│
├── HelixZero-CMS_Paper.md        Full paper manuscript
├── DATASET_README.md             Dataset documentation
└── test_record/                  Functional test results
```

**Legend**: ★ = key files for understanding the system

---

## Appendix: Performance Benchmarks

| Operation | Time | Notes |
|-----------|------|-------|
| Feature extraction (1,260 variants) | ~50 ms | Vectorized NumPy |
| Model B predict (1,260 variants) | ~30 ms | LightGBM predict |
| Single-mod scan + predict | ~200 ms | 1,260 variants |
| Beam search (max_mods=5, beam=20) | ~3–8 s | Full pipeline |
| Rank 100 candidates | ~50 ms | Naked model |
| Server startup (first request) | ~2 s | Lazy-loads models |
| Subsequent requests | <100 ms | Cached models |
