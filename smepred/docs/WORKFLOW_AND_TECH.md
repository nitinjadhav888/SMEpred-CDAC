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
9. [Biophysical Penalty System](#9-biophysical-penalty-system)
10. [API Layer](#10-api-layer)
11. [Frontend Architecture](#11-frontend-architecture)
12. [Clinical Validation](#12-clinical-validation)
13. [Pseudocode for Key Algorithms](#13-pseudocode-for-key-algorithms)
14. [Tech Stack Summary](#14-tech-stack-summary)
15. [File Structure Reference](#15-file-structure-reference)

---

## 1. System Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                            USER INTERFACES                                   │
│                                                                              │
│  ┌──────────────┐    ┌──────────────┐    ┌───────────────────────────────┐  │
│  │  Web Browser  │    │    CLI/Term   │    │  REST Client (curl, Python)  │  │
│  │  (app.html)   │    │  (cli/run.py) │    │                              │  │
│  └──────┬───────┘    └──────┬───────┘    └──────────────┬────────────────┘  │
│         │                  │                            │                    │
└─────────┼──────────────────┼────────────────────────────┼────────────────────┘
          │         HTTP     │  CLI invocation            │  HTTP
          ▼                  ▼                            ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              API LAYER                                       │
│                                                                              │
│              ┌───────────────────────────────────────────┐                  │
│              │            FastAPI Server                  │                  │
│              │  ┌─────────────────────────────────────┐  │                  │
│              │  │  GET  /                        HTML │  │                  │
│              │  │  POST /rank            → Predictor  │  │                  │
│              │  │  POST /rank/upload      → Predictor  │  │                  │
│              │  │  POST /single-mod       → Predictor  │  │                  │
│              │  │  POST /multi-mod        → Predictor  │  │                  │
│              │  │  POST /multi-mod-scan   → ModEngine  │  │                  │
│              │  │  POST /multi-mod-from-single         │  │                  │
│              │  │  GET  /modifications   → JSON file   │  │                  │
│              │  └─────────────────────────────────────┘  │                  │
│              └───────────────────────────────────────────┘                  │
└──────────────────────────────────┬──────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                          PREDICTION ENGINE                                   │
│                                                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────────┐  ┌──────────────┐│
│  │  Predictor    │  │  Features    │  │  Modification    │  │  Biophysics  ││
│  │  (predict    │──│ (extract     │──│  Engine (scan,   │──│  (penalties, ││
│  │   .py)       │  │  .py)        │  │  beam search)    │  │  adjustment) ││
│  └──────┬───────┘  └──────────────┘  └──────────────────┘  └──────────────┘│
│         │               │                                                    │
│         ▼               ▼                                                    │
│  ┌──────────────┐  ┌────────────────┐                                       │
│  │  Model B v4  │  │  Naked Model   │                                       │
│  │  (LightGBM   │  │  (LightGBM     │                                       │
│  │   1,115 tr)  │  │   .pkl)        │                                       │
│  │  1,467 feats │  │   214 feats    │                                       │
│  └──────────────┘  └────────────────┘                                       │
│                                                                              │
└──────────────────────────────────┬──────────────────────────────────────────┘
                                   │
                                   ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                        DATA & CONFIGURATION                                  │
│                                                                              │
│  ┌────────────────────┐  ┌────────────────────┐  ┌────────────────────────┐ │
│  │  Trained Models    │  │  Calibrators       │  │  Toxicity Table        │ │
│  │  (model_*.pkl)     │  │  (calibrator_*)    │  │  (cell_viability.tsv)  │ │
│  ├────────────────────┤  ├────────────────────┤  ├────────────────────────┤ │
│  │  model_b.pkl       │  │  calibrator_naked  │  │  4,097 seed hexamers   │ │
│  │  model_normal.pkl  │  │  .pkl              │  │  → viability % mapping │ │
│  └────────────────────┘  └────────────────────┘  └────────────────────────┘ │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 2. End-to-End Data Flow

### Workflow 1a: Rank siRNA Candidates (Gene/Transcript mode)

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
│                        │      Returns (N − 20) candidates
└───────────┬────────────┘
            │ List[SiRNACandidate]
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
    Return to API/CLI/UI (Rank tab)
```

### Workflow 1b: Rank siRNA Candidates (DsiRNA mode)

```
User Input (25–30 nt DsiRNA sequence, input_type="dsirna")
    │
    ▼
┌──────────────────────────────────┐
│  sirna_generator.py              │  ← Dicer cleavage rule (not sliding window)
│  generate_dsirna_candidate()     │     Takes first 21 nt as sense strand
│                                  │     Returns SINGLE candidate (not N−20)
└───────────┬──────────────────────┘
            │ List[SiRNACandidate] with exactly 1 entry
            ▼
    (Same prediction pipeline as Workflow 1a)
```

### Workflow 2: Single-Mod Scan

```
User Input (sense, antisense, options)
    │
    ▼
┌──────────────────────────────┐
│  modification_engine.py      │  ← Generate 1,302 single-mod variants
│  single_mod_scan()           │     31 mod symbols × 21 positions × 2 strands
└──────────────┬───────────────┘
               │ List[CmSiRNA] (1,302 entries for full scan)
               ▼
┌──────────────────────────────┐
│  features.py                 │  ← Extract positional features (1,467-d)
│  extract_positional_features │     Per-variant: 33 flags/pos × 42 positions
│  _batch()                    │     + 62 global counts + 16 summary + 1 log
└──────────────┬───────────────┘
               │ np.array shape: (n_variants, 1467)
               ▼
┌──────────────────────────────┐
│  predictor.py                │  ← Load Model B LightGBM (1,115 trees)
│  _get_model("B") → predict() │     Identity normalization → [0, 100]
└──────────────┬───────────────┘
               │ Raw scores → scores[0, 100]
               ▼
┌──────────────────────────────┐
│  biophysics.py               │  ← Apply 5-domain biophysical penalties
│  calculate_adjusted_efficacy()   │     adjusted = raw − 0.70 × total_penalty
│                              │     Clamped to [0, 100]
└──────────────┬───────────────┘
               │ Adjusted scores with penalty breakdown
               ▼
┌──────────────────────────────┐
│  filters.py                  │  ← Modification-aware seed toxicity
│  toxicity_for_modified()     │     "Mitigated" if rescue mod present
└──────────────┬───────────────┘
               │ RankedCmSiRNA list, sorted by adjusted efficacy
               ▼
     Return to API: top-N variants + both baseline scores
```

### Workflow 3: Multi-Mod Beam Search

```
User Input (sense, antisense, max_mods, beam_width, full_scan)
    │
    ▼
┌──────────────────────────────────────────────┐
│  Step 0: Compute both baselines              │
│  - Naked baseline (V4 model, Rank tab view)  │
│  - Model B baseline (positional, fair delta) │
│  Both returned for transparent display       │
└──────────────────┬───────────────────────────┘
                   │
                   ▼
┌──────────────────────────────────────────────┐
│  Step 1: Run single-mod scan (or mini-scan)  │
│  Full: 31 mods × 21 pos × 2 strands = 1,302  │
│  Mini: E/D/Q/L on AS pos 1-10 = 40 variants  │
└──────────────────┬───────────────────────────┘
                   │ List[RankedCmSiRNA] + parent_score
                   ▼
┌──────────────────────────────────────────────┐
│  Step 2: Diversify initial beam              │
│  Round-robin through all mod symbol types    │
│  Pick top from each → beam width K           │
└──────────────────┬───────────────────────────┘
                   │ beam = [variant_1, ..., variant_K]
                   ▼
┌──────────────────────────────────────────────┐
│  Step 3: Score all beam candidates           │
│  Batch feature extraction + batch predict    │
│  (fully vectorized — single model.predict()) │
└──────────────────┬───────────────────────────┘
                   │ Scored + sorted by adjusted efficacy
                   ▼
┌──────────────────────────────────────────────┐
│  Step 4: Expand (iterative, n = 2..max_mods) │
│  for each beam candidate:                    │
│    for each single_result (pool top 3×K):    │
│      combine if no position conflict         │
│  score new candidates (batched)              │
│  keep top-K as new beam                      │
│  if current_best - best[−3] < 0.5 → STOP    │
└──────────────────┬───────────────────────────┘
                   │ All scored variants (1-mod through N-mod)
                   ▼
    Return to API: ranked results + baselines + metadata
```

---

## 3. Module-by-Module Walkthrough

### `src/parser.py` — Sequence I/O

- `normalize_seq()`: Converts DNA T→U, validates only A/U/G/C characters, uppercases.
- `load_sequence()`: Accepts plain string, inline FASTA (`>header\nseq`), or file path. Returns a single concatenated RNA string.

### `src/sirna_generator.py` — Candidate Generation

- `_reverse_complement()`: Uses `str.maketrans("AUGC", "UACG")` + slicing `[::-1]`.
- `generate_candidates()`: Sliding window of size 21, step size 1. Each window → `SiRNACandidate(position, sense, antisense)`. Produces `len(mrna) − 20` candidates.
- `generate_dsirna_candidate()`: Dicer-substrate mode for 25–30 nt inputs. Takes first 21 nt as sense strand, returns exactly one candidate (not a sliding window).

### `src/features.py` — Feature Extraction (Two Models)

**Model B v4 — Positional Features (1,467 dimensions):**

1. **Per-position flags (33 × 42 = 1,386):** For each of 42 positions (21 sense + 21 antisense), encode 33 Boolean flags: 4 canonical bases (A/U/G/C), 31 modification symbols. One-hot-like but can have multiple flags (a modified position has both the base and the mod symbol flags set).

2. **Per-strand global counts (31 × 2 = 62):** For each strand (sense, antisense), count occurrences of each of the 31 modification symbols.

3. **Summary statistics (8 × 2 = 16):** For each strand: total mod count, unique mod types, M count, F count, L count, E count, D count, PS backbone count.

4. **Log concentration (1):** log₁₀(10) = 1.0 (default 10 nM).

**Naked V4 — Sequence Features (214 dimensions):**

1. **Sense one-hot (84):** 4 bases × 21 positions.
2. **Sense TNC (64):** Trinucleotide composition (4³ = 64).
3. **Antisense TNC (64):** Same for antisense.
4. **GC content (2):** Sense GC%, antisense GC%.

### `src/predictor.py` — Unified Prediction Interface

**Functions:**

- `_get_model(key)`, `_get_calibrator(key)`: Lazy-load LightGBM models and sklearn calibrators with caching. Searches `models/` directory by naming convention.
- `_predict_naked(X)`: Handles the legacy source one-hot padding required by models trained on the old `all_features.csv` format (adds 400 columns of source one-hot encoding to reach 214 + 400 = 614 total).
- `_normalize_scores(raw, mode, calibrator_key)`: Normalizes raw LightGBM scores to [0, 100]:
  - `mode="identity"`: Pass-through (used for Model B).
  - `mode="clip"`: Clip to [0, 100].
  - `mode="rescale"`: Min-max rescale to [0, 100].
  - `mode="calibrate"`: Platt-calibrate via `CalibratedClassifierCV`, then apply `_v4_transform`.
- `rank_sirnas()` / `rank_by_naked_score()`: Workflow 1 — parse sequence, generate candidates, extract V4 features, predict with Naked model, apply filters.
- `predict_modified()`: Workflow 2 — compute both baselines (Naked Model for Rank-compatible view, Model B for fair delta calculation), generate variants, extract positional features, predict with Model B, apply biophysical penalties, return ranked results with both baselines.

**Models:**
- **Naked Model (V4)**: 214-dimensional sequence features, PCC=0.55. Used for initial screening in the Rank tab.
- **Model B (HelixZero v4)**: 1,467-dimensional position-aware features, PCC=0.822, Spearman=0.823. Trained on 83,535 cm-siRNA variants. 1,115 trees, 127 leaves, learning rate 0.03.

### `src/modification_engine.py` — Chemical Modification Space

**Modification Symbols (31 total):**

| Category | Symbols | Count |
|----------|---------|-------|
| Canonical | A, U, G, C | 4 |
| Sugar | F (2'-F), M (2'-OMe), L (LNA), E (2'-MOE), D (2'-O-DMAOE) | 5 |
| Backbone | S (PS), 1 (5'-PO₄) | 2 |
| Base | 2 (dihydrouridine), 3 (pseudouridine), 5 (5-Me-C), 6 (UNA), 8 (GNA), 9 (TNA), Y (ENA) | 7 |
| Exotic | B (2'-F-ANA), J (2'-O-Pyrene), V (2'-O-N3), I (Inosine), N (2'-O-N3-A), O (2'-O-N3-U), P (2'-O-N3-C), R (2'-O-N3-G), H (2'-O-N3-T), K (LNA-T), Z (LNA-C), Q (α-LLNA), W (2'-O-allyl), X (2'-O-propargyl), 7 (locked-ENA) | 13 |

**Key functions:**
- `_apply_mod()`: Replace a single nucleotide with a modification symbol. Handles canonical substitution (e.g., A→M for 2'-OMe) and exotic insertions.
- `_parse_multimod_input()`: Parses comma-separated strings like `"F,,M"` and `"2,5,,10,12"` into modification lists.
- `single_mod_scan()`: Generates all 31 × 21 × 2 = 1,302 single-modification variants.
- `multimod_gen()`: Applies multiple modifications from user-specified strings.
- `multi_mod_scan()`: Full beam search with diversity initialization, batched scoring, and plateau-based early stopping.

### `src/biophysics.py` — Five-Domain Penalty System

Five orthogonal penalty domains adjust the raw efficacy score. See [Section 9](#9-biophysical-penalty-system) for complete details.

- **Nuclease (0–16):** PS backbone coverage, 2'-modification density.
- **Immunogenicity (0–28):** Unmodified uridine residues, GU-rich motifs (non-stacking), over-methylation.
- **RISC Loading (−10 to 60):** 5'-phosphate, seed modifications, LNA/MOE/GNA/ENA/TNA rules, 2'-F deficiency, exotic micro-penalties.
- **Thermo (0–20):** GC extremes, palindrome, homopolymer, GC runs.
- **Serum (0–17):** Termini protection (PS or 5'-PO₄/GalNAc).

**Adjusted Score:**
```
adjusted = max(0.0, min(100.0, raw_score − 0.70 × total_penalty))
```

### `src/filters.py` — Safety and Toxicity Checks

- `toxicity_score()`: Looks up antisense seed hexamer (positions 2–7) in pre-loaded `cell_viability.tsv` table (4,097 seeds). Returns viability %.
- `toxicity_label()`: ≥75% → Safe, ≥55% → Caution, <55% → Toxic, not found → Unknown.
- `seed_rescue_check()`: Detects rescue modifications (M, F, L, E) in seed positions 2–7.
- `toxicity_for_modified()`: If rescue mod is present and base seed is Toxic/Caution → label becomes "Mitigated".
- `functional_check()`: Reynolds/Ui-Tei rules — GC content 30–65%, no homopolymer runs (≥4), no GC run ≥6, no palindrome (≥8 nt), no consecutive runs ≥4.
- `annotate_candidates()`: Batch annotation for the Rank tab.

---

## 4. Model Selection Rationale

### Why LightGBM over alternatives?

| Model | Pros | Cons | Performance |
|-------|------|------|-------------|
| **LightGBM** | Fast training, native categorical support, built-in regularization, feature importance, GPU support | Requires careful tuning for small datasets | **PCC=0.822** |
| SVR (RBF) | Works well on <5k samples | Poor scaling, no categorical support, slow prediction | PCC=0.719 |
| Random Forest | Low variance, interpretable | Poor extrapolation, slow at high tree counts | PCC=0.78 |
| XGBoost | Slightly higher accuracy on dense tables | Slower training, heavier memory | PCC=0.81 |
| Neural Net | Flexible architecture | Requires 10× more data to train, overfits on 83k rows | PCC=0.75 |

**Decision:** LightGBM offers the best accuracy-speed tradeoff for this dataset (83,535 rows, 1,467 features). The `goss` (Gradient-based One-Side Sampling) and `histogram`-based splitting make it significantly faster than XGBoost at comparable accuracy.

---

## 5. Feature Engineering

### Evolution of Feature Spaces

| Version | Dimensions | Description | When |
|---------|-----------|-------------|------|
| V1 | 70 | Basic composition + di-nucleotide | Initial prototype |
| V2 | 214 | V4 successor before renumbering | Pre-release |
| V3 | 400 | Augmented with source one-hot | Training pipeline |
| **V4 (Naked)** | **214** | Sense one-hot + TNC + GC. Used for unmodified screening (Rank tab) | Current |
| **V4 (Model B)** | **1,467** | Position-aware flags + global counts + summary stats + log conc | **Current** |

### Model B v4 — 1,467 Features Breakdown

| Group | Calculation | Dimensions | Purpose |
|-------|-------------|-----------|---------|
| Per-position flags | 33 flags (4 bases + 29 mods) × 42 positions (21 sense + 21 antisense) with parent awareness | 1,386 | Encodes exactly what modification is at each position, separately for current and parent |
| Per-strand global counts | 31 mod symbols × 2 strands (sense, antisense) | 62 | Global modification burden per strand |
| Summary statistics | 8 stats (total mods, unique types, M/F/L/E/D/PS counts) × 2 strands | 16 | High-level chemistry profile |
| Log concentration | log₁₀(10 nM) | 1 | Dose normalization |
| **Total** | | **1,467** | |

### Parent-Variant Encoding

Each variant's feature vector encodes BOTH the variant sequence and its parent (unmodified) sequence. This allows the model to learn delta-effects (what changes when you add a modification). The per-position flags at positions 22–42 encode the parent sequence, giving the model access to the baseline at every position.

### Key Insight: Feature Space Asymmetry

The system uses **two different LightGBM models** trained on **different feature spaces**:
- **Naked Model**: 214-d sequence features (no chemistry awareness). Used for initial Rank tab screening.
- **Model B (HelixZero)**: 1,467-d position-aware features (fully chemistry-aware). Used for Single-Mod and Multi-Mod.

This asymmetry causes different raw scores for the same unmodified sequence (e.g., 27.09 from Naked vs 29.88 from Model B for Givosiran ALAS1). The system now explicitly computes and displays **both** baselines (`naked_baseline` and `model_b_baseline`) to prevent user confusion (the "Score Jump Bug").

---

## 6. Training Pipeline

### Data Sources (83,535 total rows)

| Source | Description | Rows | Use |
|--------|-------------|------|-----|
| HelixZero-CMS training set | 83,535 cm-siRNA variants from Dar et al. (RNA Biology, 2016) | 83,535 | Primary training |

### Training Configuration

- **Model**: LightGBM (`objective='regression'`, `metric='rmse'`)
- **Parameters**: 1,115 trees, 127 leaves, learning rate 0.03, `feature_fraction=0.8`, `bagging_fraction=0.8`, `bagging_freq=5`, `lambda_l1=0.1`, `lambda_l2=0.1`, `min_data_in_leaf=20`
- **Features**: 1,467-d position-aware (Section 5 above)
- **Calibration**: Platt scaling (`CalibratedClassifierCV`) applied to raw LightGBM outputs, then transformed via `_v4_transform()` to map to [0, 100]. The calibration transforms the regression output into a probability-like score, improving ranking consistency.
- **Validation**: 5-fold cross-validation. Held-out test set achieves PCC=0.822, Spearman=0.823.

### Score Normalization Pipeline

```
Raw LightGBM score
    │
    ▼
┌─────────────────┐
│  Platt Calibrate │  ← CalibratedClassifierCV (sigmoid)
│  → [0, 1] prob  │     Makes scores probabilistic
└────────┬────────┘
         │ score in [0, 1]
         ▼
┌─────────────────┐
│  _v4_transform  │  ← sigmoid(score) → beta CDF transform
│  → [0, 100]     │     Maps to biological relevance scale
└─────────────────┘
```

---

## 7. Inference Pipeline

### Two-Model Architecture

```
                    ┌──────────────────────────────┐
                    │     User Input Sequence       │
                    └──────────┬───────────────────┘
                               │
                    ┌──────────▼───────────────────┐
                    │   Decision: Rank or Modify?   │
                    └────┬──────────────┬──────────┘
                         │              │
               Rank tab  │              │  Single-Mod / Multi-Mod tab
                         ▼              ▼
              ┌──────────────────┐  ┌──────────────────┐
              │  214-d Features  │  │ 1,467-d Features │
              │  (V4 naked)      │  │ (position-aware) │
              └────────┬─────────┘  └────────┬─────────┘
                       │                     │
                       ▼                     ▼
              ┌──────────────────┐  ┌──────────────────┐
              │  Naked Model     │  │  Model B v4      │
              │  LightGBM .pkl   │  │  LightGBM .pkl   │
              │  PCC=0.55        │  │  PCC=0.822       │
              └────────┬─────────┘  └────────┬─────────┘
                       │                     │
                       ▼                     │
              ┌──────────────────┐            │
              │  Rank + Filter   │            │
              │  (toxicity,      │            │
              │   functional)    │            │
              └────────┬─────────┘            │
                       │                     │
                       ▼                     ▼
              ┌──────────────────────────────────────┐
              │    Biophysical Penalty Adjustment     │
              │    adjusted = raw − 0.70 × penalties  │
              └──────────────────────────────────────┘
```

### Batch Prediction Strategy

- **Single-Mod (1,302 variants)**: All variants extracted in one batch, predicted in a single `model.predict(X)` call. Fully vectorized.
- **Multi-Mod beam search**: Each round's candidates batched into one `predict()` call. The beam search is CPU-bound by candidate generation/permutation logic, not ML inference.

---

## 8. Modification Engine Algorithms

### Single-Mod Scan

```
Input: sense (21-nt), antisense (21-nt)
Output: List[CmSiRNA] — all 1,302 single-mod variants

Algorithm:
  for each position p in [0..20]:
    for each strand s in [sense, antisense]:
      for each symbol sym in modification_symbols:
        variant = CmSiRNA(
          sense     = apply_mod(base_sense, sym, p) if s == sense else base_sense,
          antisense = apply_mod(base_antisense, sym, p) if s == antisense else base_antisense,
          mod_symbol     = sym,
          mod_position   = p + 1,
          mod_strand     = s,
          parent_sense   = base_sense,
          parent_antisense = base_antisense
        )
  return all variants
```

### Beam Search (Multi-Mod)

```
Input: sense, antisense, max_mods, beam_width
Output: List[RankedCmSiRNA]

Algorithm:
  # Phase 1: Initialization
  single_results = single_mod_scan(sense, antisense)     // 1,302 variants
  score_variants(single_results)                          // batched predict
  single_results.sort(by=efficacy_score, descending)

  # Phase 2: Diversify initial beam (round-robin)
  beam = []
  for each mod_symbol type (sorted by best score):
    if symbol has unselected candidates:
      pick the best scoring one
      add to beam
      stop when beam size == beam_width

  score_variants(beam)

  # Phase 3: Expand
  pairing_pool = single_results[:beam_width × 3]          // 90 candidates max
  round_best_scores = [beam[0].efficacy_score]

  for n_mods = 2 to max_mods:
    current_best = beam[0].efficacy_score

    // Early stopping: plateau detection
    if n_mods >= 4 AND current_best − round_best_scores[−3] < 0.5:
      break

    candidates = []
    for each beam_variant:
      for each single_result in pairing_pool:
        if no position conflict:
          new_variant = merge(beam_variant, single_result)
          candidates.append(new_variant)

    candidates = remove_duplicates(candidates)
    score_variants(candidates)                             // batched predict
    candidates.sort(by=efficacy_score, descending)
    beam = candidates[:beam_width]

  return beam ∪ single_results (all scored variants)
```

### Key Optimizations

1. **Pairing pool capped at 3× beam_width** (typically 90 candidates, was unlimited). Reduced search time from 300+ seconds to ~20 seconds.
2. **Batch scoring**: Every round's candidates are scored in a single `extract_positional_features_batch()` + `model.predict(X)` call.
3. **Plateau-based early stopping**: Halts when best score improves <0.5 over 3 rounds. Eliminates unnecessary rounds where the model has converged.
4. **No artificial over-mod penalty**: The system lets beam search find the natural optimal mod count.

---

## 9. Biophysical Penalty System

### Overview

Five orthogonal penalty domains subtract from the raw efficacy score. The domains are designed to be **strictly non-overlapping** — no biological feature is penalized by more than one module.

```
adjusted_score = max(0.0, min(100.0, raw_score − 0.70 × total_penalty))
```

### 9.1 Nuclease Penalty (0–16)

Targets **endonuclease** stability. Does NOT check termini (that is serum's domain).

| Condition | Penalty | Citation |
|-----------|---------|----------|
| PS count == 0 (no backbone protection) | +5 | Braasch 2004 |
| PS count < 3 (minimal backbone protection) | +3 | Braasch 2004 |
| 2'-mod density < 20% | +4 | Czauderna 2003 |
| 2'-mod density < 40% | +2 | Czauderna 2003 |

### 9.2 Immuno Penalty (0–28)

Targets innate immune activation via TLR7/8 sensors.

| Condition | Penalty | Citation |
|-----------|---------|----------|
| Unmodified U in antisense seed (pos 2–8), each | **+2.0** | Sioud & Sørensen 2004 |
| Unmodified U in antisense tail (pos 9–21), each | **+0.5** | Goodchild 2009 |
| Unmodified U in sense strand, each | **+1.0** | Judge 2005 |
| GU-rich motif GUUGU (highest immunostimulatory) | +3 | Goodchild 2009 |
| GU-rich motif GUGU | +3 | Goodchild 2009 |
| GU-rich motif UGU (weakest) | +3 | Goodchild 2009 |
| Over-methylation (M > 24) | +4 advisory | Alnylam ESC design |

**Note**: Motif detection uses **non-stacking hierarchical search**. A single 5-nt window is checked for GUUGU first; if matched, all 5 positions are masked with a sentinel character so the same window cannot also trigger GUGU or UGU penalties.

**Calibration note**: Seed U penalty was reduced from +4 → +2.0 and tail U from +1 → +0.5 after C-DAC panel review (June 2026) to prevent overtuning on low-GC sequences. Over-methylation threshold raised from >16 → >24 to match clinical ESC designs (which safely use 25–27 M's).

### 9.3 RISC Loading Penalty (−10 to 60)

Targets guide strand loading and thermodynamic asymmetry. Can be NEGATIVE (bonus for beneficial chemistries).

| Condition | Penalty | Citation |
|-----------|---------|----------|
| Missing 5'-phosphate (AS pos 1) | +5 | Frank 2010 |
| PS at AS pos 1 | +2 | |
| Unmodified seed position (AS pos 2–8), each | +2 | Jackson 2006 |
| UNA at AS pos 7 exempt from seed penalty | 0 | Bramsen 2010 |
| LNA at AS pos 2–4, each | +5 | Hidayah 2021 |
| MOE at AS pos 2–14, each | +3 | Prakash 2005 |
| GNA at AS pos 2–5 (disruptive), each | +4 | Schlegel 2022 |
| GNA at AS pos 6–8 (beneficial), each | **−2 bonus** | Schlegel 2022 ESC+ |
| ENA at AS pos 2–8, each | +4 | Morihiro 2020 |
| ENA at AS pos 9–14, each | +2 | Morihiko 2020 |
| TNA at AS pos 2–6, each | +3 | |
| TNA at AS pos 7 (exempt) | 0 | |
| TNA at AS pos 8–14, each | +1 | |
| 2'-F deficiency on pyrimidines < 20% | +6 | Layzer 2004 |
| 2'-F deficiency on pyrimidines < 40% | +3 | Layzer 2004 |
| Exotic mod micro-penalty (Benzyl, Inosine), each | +2 | |
| Other rare exotic mods, each | +1 | |

### 9.4 Thermo Penalty (0–20)

Targets melting temperature extremes that reduce RISC loading specificity.

| Condition | Penalty | Citation |
|-----------|---------|----------|
| GC < 30% or > 55% | +8 | Reynolds 2004 |
| GC 30–35% or 50–55% | +3 | Reynolds 2004 |
| Palindrome (≥8 nt self-complementary) | +5 | |
| Homopolymer run (≥4 same nt) | +5 | |
| GC run ≥6 consecutive G or C | +3 | |

### 9.5 Serum Penalty (0–17)

Targets **exonuclease** degradation — checks only whether termini are protected.

| Condition | Penalty | Citation |
|-----------|---------|----------|
| AS 5' not PS/1 (5'-PO₄) | +4 | |
| AS 3' not PS | +3 | Elmén 2005 |
| SS 5' not PS/4 (GalNAc conjugate) | +3 | |
| SS 3' not PS/4 (GalNAc conjugate) | +2 | |

**Note**: Does NOT check modification density (that is nuclease's domain). PS at termini is the primary determinant of exonuclease resistance.

---

## 10. API Layer

### Endpoints

| Endpoint | Method | Input | Output |
|----------|--------|-------|--------|
| `/` | GET | — | Serves `app.html` |
| `/app.html` | GET | — | Serves `app.html` |
| `/rank` | POST | `{sequence, top_n, input_type}` | Ranked siRNA candidates with `input_type` in response |
| `/rank/upload` | POST | `{file, top_n}` | Same as `/rank` for FASTA upload |
| `/single-mod` | POST | `{sense, antisense, model, top_n, full_scan}` | Single-mod variants + `naked_baseline` + `model_b_baseline` |
| `/multi-mod` | POST | `{sense, antisense, sense_mods, sense_positions, antisense_mods, antisense_positions}` | Single multi-mod result + both baselines |
| `/multi-mod-scan` | POST | `{sense, antisense, max_mods, beam_width, full_scan}` | Beam search results + both baselines |
| `/multi-mod-from-single` | POST | `{sense, antisense, max_mods, beam_width, full_scan, single_results, parent_score}` | Beam search with optional pre-computed singles |
| `/modifications` | GET | — | All 31 modification symbols with names and types |

### Response Enrichment (All Modification Endpoints)

Every endpoint in the modification pipeline returns:
- `parent_score`: Adjusted score using Model B baseline
- `naked_baseline`: Adjusted score using Naked Model (for Rank tab compatibility)
- `model_b_baseline`: Same as `parent_score` (Model B, used for delta computation)
- Per-variant: `efficacy_score`, `delta_score`, `efficacy_label`, `raw_efficacy_score`, `total_penalty`, `penalties` dict (nuclease, immuno, risc, thermo, serum), `toxicity_score`, `toxicity_label`, `toxicity_note`

### Dual-Baseline Strategy

The feature space asymmetry between Naked Model (214-d) and Model B (1,467-d) causes different raw scores for the same sequence. To prevent user confusion:

```
Rank Tab shows:    Naked Model score (e.g., 27.09)
User clicks Multi-Mod:
                    → "Recalibrating baseline for chemical space..."
                    → Naked: 27.09 → Model B: 29.88 (+2.79)
                    → All deltas computed against Model B baseline
```

This explicit recalibration message manages user expectations and protects credibility.

---

## 11. Frontend Architecture

### Single-File HTML (`app.html`)

All frontend code (HTML + CSS + JavaScript) is contained in a single file. No build step, no framework dependencies.

### Tab Structure

| Tab | ID | Workflow |
|-----|----|----------|
| Rank siRNAs | `tab-rank` | Workflow 1 — input sequence, view ranked candidates |
| Single-Mod | `tab-singlemod` | Workflow 2 — scan 1,302 single-mod variants |
| Multi-Mod | `tab-multimod` | Workflow 3 — beam search multi-mod + custom design |

### Key Features

1. **Score breakdown popup**: Hover over any score to see a breakdown of all 5 penalty domains. Uses JS-driven viewport-aware positioning (tries right → flips left → flips below if both sides overflow).

2. **Expandable rows**: Click any row to expand and view the full sequence with a color-coded heatmap plus penalty details. Uses `toggleRow()` with sibling `expanded-detail` rows.

3. **Chemistry Confidence badge**: Each modification symbol gets a confidence rating: Standard (M, F, S, 1, L, E, D) → "Standard ✅", Exotic (B, J, V, I, N, O, P, R, H, K, Z, Q, W, X, 7) → "Exotic ⚗️", Rare → "Rare 🧪".

4. **Modification legend**: Loaded dynamically from `/modifications` endpoint and displayed as tappable tags that copy the symbol to clipboard.

5. **Continuous progress bar**: Animated CSS bar during multi-mod beam search to indicate loading without freezing the UI.

6. **Input Type toggle**: Rank tab has `Gene/Transcript` (sliding window) vs `DsiRNA (27-mer)` (Dicer cleavage) radio buttons that update the API call.

7. **FASTA upload**: File upload support for the Rank tab with drag-and-drop UX.

### Cross-Tab Navigation

Users can click "Multi-Mod" from a Rank tab result row, which:
1. Fills the Multi-Mod sense/antisense fields
2. Switches to the Multi-Mod tab
3. Automatically triggers the beam search with the max_mods setting from the Rank tab

---

## 12. Clinical Validation

### ESC / ESC+ Clinical Benchmark

The `tests/test_clinical_benchmark.py` validates the system against published clinical siRNA architectures:

| Sequence | Target | ESC | ESC+ | Preclinical PK |
|----------|--------|-----|------|----------------|
| Seq_HighGC33 | Givosiran ALAS1 | 62.0 | 65.1 | ✓ All bounds |
| Seq_GC48a | Custom | 61.3 | 61.4 | ✓ All bounds |
| Seq_GC38b | Custom | 63.8 | 65.2 | ✓ All bounds |
| Seq_GC48b | Custom | 55.8 | 54.3 | ✓ All bounds |

**ESC design**: 5'-PO₄, PS at all termini, full 2'-OMe modifications, GalNAc conjugate at SS 3'.
**ESC+ design**: Same as ESC + GNA at AS position 7 (confers −2 RISC bonus per Schlegel 2022).

All designs score ≥50 after biophysical penalties, confirming clinical viability.

### Termini Recognition

The system recognizes clinical end-cap symbols:
- `1` (5'-PO₄) — accepted at AS 5' end
- `4` (GalNAc conjugate) — accepted at SS 3' end

This prevents false-positive serum/nuclease penalties on actual drug architectures (e.g., Givosiran, Patisiran, inclisiran).

---

## 13. Pseudocode for Key Algorithms

### Candidate Generation
```
function generate_candidates(mrna):
    candidates = []
    for i in 0..len(mrna)−21:
        sense = mrna[i : i+21]
        antisense = reverse_complement(sense)
        candidates.append(Candidate(i, sense, antisense))
    return candidates

function generate_dsirna_candidate(seq):
    if len(seq) < 25 or len(seq) > 30:
        error("DsiRNA input must be 25-30 nt")
    sense = seq[0:21]
    antisense = reverse_complement(sense)
    return [Candidate(0, sense, antisense)]
```

### Feature Extraction (Model B)
```
function extract_features(sense, antisense, parent_sense, parent_antisense):
    features = []
    // Per-position flags (33 × 42)
    for pos in concatenate(sense, antisense):
        flags = [pos == 'A', pos == 'U', pos == 'G', pos == 'C',
                 pos == 'F', pos == 'M', ..., pos == '7']  // 33 total
        features.extend(flags)
    for pos in concatenate(parent_sense, parent_antisense):
        features.extend([...])  // 33 flags for parent
    // Global counts (31 × 2)
    features.extend(count_mods(sense))
    features.extend(count_mods(antisense))
    // Summary stats (8 × 2)
    features.extend([n_mods, unique_types, M_ct, F_ct, L_ct, E_ct, D_ct, PS_ct])
    // Log concentration
    features.append(log10(10))
    return features  // 1,467 dimensions
```

### Beam Search (Multi-Mod)
```
function multi_mod_scan(sense, antisense, max_mods=14, beam_width=30):
    // Phase 1: Single-mod scan
    single_results = single_mod_scan(sense, antisense)
    batch_score(single_results)

    // Phase 2: Diversify initial beam
    beam = diversify(single_results, beam_width)

    // Phase 3: Expand
    pool = single_results.top(beam_width × 3)
    best_scores = [beam[0].score]

    for n in 2..max_mods:
        if n ≥ 4 AND beam[0].score − best_scores[−3] < 0.5:
            break  // plateau detected

        candidates = []
        for b in beam:
            for s in pool:
                if not conflict(b, s):
                    candidates.append(merge(b, s))

        batch_score(candidates)
        beam = candidates.top(beam_width)

    return beam ∪ single_results
```

### Biophysical Penalty (Adjusted Score)
```
function adjusted_score(raw, sense, antisense, parent_sense, parent_antisense):
    penalties = {}
    penalties['nuclease'] = calculate_nuclease_penalty(sense, antisense)
    penalties['immuno']   = calculate_immuno_penalty(sense, antisense)
    penalties['risc']     = calculate_risc_penalty(sense, antisense)
    penalties['thermo']   = calculate_thermo_penalty(sense, antisense)
    penalties['serum']    = calculate_serum_penalty(sense, antisense, parent_sense, parent_antisense)
    total = sum(penalties.values())
    adjusted = raw − 0.70 × total
    adjusted = clamp(adjusted, 0, 100)
    return adjusted, penalties, total
```

---

## 14. Tech Stack Summary

| Layer | Technology | Version | Purpose |
|-------|-----------|---------|---------|
| **Backend** | Python | 3.13 | Language |
| **API** | FastAPI | 0.115+ | REST framework |
| **Server** | Uvicorn | 0.34+ | ASGI server |
| **Model** | LightGBM | 4.6+ | Gradient boosting |
| **Features** | NumPy / Pandas | 2.x / 2.x | Matrix operations |
| **Calibration** | scikit-learn | 1.6+ | Platt scaling |
| **Frontend** | Vanilla JS | ES2024 | Single-page app |
| **CLI** | argparse | stdlib | Command-line interface |
| **PDF** | reportlab | 4.x | Validation report generation |
| **Testing** | pytest | 8.x | Unit tests |
| **Deployment** | pip | 25.x | Package management |

---

## 15. File Structure Reference

```
helixzero_cms/
├── app.html                     # Single-file web UI (all 4 tabs)
├── README.md                    # Quick-start guide
├── setup.py / pyproject.toml    # Package config
├── api/
│   └── main.py                  # FastAPI server (9 endpoints)
├── cli/
│   └── run.py                   # Command-line interface
├── src/
│   ├── predictor.py             # Unified prediction (2 models, normalization)
│   ├── sirna_generator.py       # Candidate generation (gene + DsiRNA modes)
│   ├── features.py              # Feature extraction (V4 + positional)
│   ├── modification_engine.py   # Chemical modification space (31 symbols, beam search)
│   ├── biophysics.py            # 5-domain penalty system (orthogonal)
│   ├── filters.py               # Toxicity, safety, functional checks
│   ├── parser.py                # Sequence I/O (FASTA, normalization)
│   └── __init__.py
├── models/
│   ├── model_b.pkl              # Model B v4 (1,467-d, PCC=0.822)
│   ├── model_normal.pkl         # Naked model (214-d, PCC=0.55)
│   └── calibrator_naked.pkl     # Platt calibrator for naked model
├── data/
│   ├── cell_viability.tsv       # 4,097 seed hexamer → viability mapping
│   ├── modification_codes.json  # 31 modification symbol definitions
│   └── calibator.pkl            # Histogram bin edges
├── tests/
│   ├── test_pipeline.py         # 32 unit tests (features, biophysics, engine)
│   └── test_clinical_benchmark.py  # ESC/ESC+ clinical validation (4 sequences)
├── docs/
│   ├── WORKFLOW_AND_TECH.md     # This file — full architecture reference
│   ├── README.md                # Quick-start summary
│   ├── PENALTIES_REFERENCE.md   # Full penalty rules with 28+ citations
│   ├── VALIDATION_DOSSIER.md    # Cross-reference against 14 PDF sources
│   ├── LIMITATIONS_NOTES.md     # Known limitations and improvement roadmap
│   ├── HelixZero_Model_Reference.md  # Complete model parameter reference
│   ├── Clinical_Validation_Report.pdf  # Clinical benchmark PDF report
│   ├── main.pdf                 # IEEE-format research paper
│   └── figures/                 # 7 PNG visualizations
└── scripts/
    ├── train_model.py           # Training pipeline
    └── generate_paper.py        # Research paper generator
```
