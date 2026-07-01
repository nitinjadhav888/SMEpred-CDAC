# Data Handling & Analysis: From Raw Patents to Training-Ready Features

## Table of Contents

1. [Data Sources Overview](#1-data-sources-overview)
2. [Data Cleaning Pipeline](#2-data-cleaning-pipeline)
3. [Sequence Recovery & Parsing](#3-sequence-recovery--parsing)
4. [Symbol Mapping (Modification Names → Codes)](#4-symbol-mapping-modification-names--codes)
5. [Deduplication](#5-deduplication)
6. [Feature Extraction](#6-feature-extraction)
7. [Train/Validation Splitting](#7-trainvalidation-splitting)
8. [Dataset Statistics](#8-dataset-statistics)
9. [File Organization](#9-file-organization)
10. [Data Quality Checks](#10-data-quality-checks)
11. [Edge Cases & Special Handling](#11-edge-cases--special-handling)

---

## 1. Data Sources Overview

### 1.1 Model B v4 — Modified siRNA Data (83,535 rows)

| Source | Rows | Description | Origin |
|--------|------|-------------|--------|
| **Position-Aware Dataset** | 55,730 | Sequences with per-position chemical modification annotations in `pos*mod` format. Columns include sense/antisense sequences, modification strings, efficacy (inhibition %), gene, cell type, dose, timepoint. | HelixZero Biological Catalog 43k (patent-derived) + synthetic augmentation |
| **Hetero_train (HelixZero-CMS)** | 23,187 | The original HelixZero-CMS training set. 2,728 curated cm-siRNA rows from siRNAmod database, augmented with position-aware annotation. Used for external validation in the original paper (Dar et al. 2016). | siRNAmod database, Dar et al. *RNA Biology* 2016 |
| **CMsiRNAdb** | 4,618 | External independent dataset. 12,303 total rows, of which 4,618 passed our cleaning filters and symbol mapping. PCSK9, PNPLA3 gene targets from patent data. | CMsiRNAdb (multiple patents) |

### 1.2 Naked Model — Unmodified siRNA Data (4,060 rows)

| Source | Rows | Description |
|--------|------|-------------|
| Huesken (2005) | 2,182 | Largest published siRNA efficacy dataset, 2,431 sequences targeting Firefly luciferase |
| Takayuki (Ui-Tei) | 594 | Independent set targeting endogenous genes |
| Mix | 539 | Multi-lab compilation targeting various genes |
| *Extended* | 745 | Additional sequences from follow-up studies |

### 1.3 Seed Toxicity Table

| Source | Entries | Description |
|--------|---------|-------------|
| Janas et al., *Mol Cell* 2018 | 4,097 | 6-mer seed → cell viability (%) mapping. Used for toxicity lookup. |

**File:** `data/oligoformer/cell_viability.tsv`

---

## 2. Data Cleaning Pipeline

The cleaning pipeline converts raw CSV data into training-ready feature vectors. Here is the step-by-step process.

```
Raw CSV (HelixZero Biological Catalog 43k)
    │
    ▼
┌─────────────────────────────────┐
│  1. Load & Validate Columns     │
│     - Check required columns    │
│     - Verify column types       │
│     - Drop rows with NaN in     │
│       critical fields           │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│  2. Sequence Recovery           │
│     - Parse pos*mod token       │
│       streams from malformed    │
│       CSV columns               │
│     - Reconstruct sense and     │
│       antisense sequences       │
│     - Derive canonical bases    │
│       from mod name suffix      │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│  3. Symbol Mapping              │
│     - Map modification names    │
│       to 1-letter codes         │
│     - Ordered alias rules       │
│       (sugar > backbone)        │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│  4. Deduplication               │
│     - Remove exact (sense,      │
│       antisense, efficacy) dups │
│     - Remove strand length      │
│       outliers (<19, >25)       │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│  5. Efficacy Parsing            │
│     - Extract inhibition %      │
│       from row-ID suffix        │
│     - Handle negative values    │
│       (double-dash encoding)    │
│     - Cross-validate with       │
│       Inhibition column         │
└──────────────┬──────────────────┘
               ▼
┌─────────────────────────────────┐
│  6. Feature Extraction          │
│     - 1,467 positional features │
│     - Store features + label    │
│       as training array         │
└──────────────┬──────────────────┘
               ▼
    Training-ready (X, y) pair
```

### Example: Before and After Cleaning

**Raw row (original CSV, malformed):**
```
Compound-ID: HBC-43k-HE1-PCSK9-M1-48h-92.00
Sequence-Sense: 5' mG*mC*mA*mG*mC*mA*mC*mG*mA*mC*mU*mU*mC*mU*mU*mC*mA*mA*mG*mU*mU 3'
Sequence-Antisense: 3' mC*mU*mU*mG*mA*mA*mG*mA*mA*mG*mU*mC*mG*mU*mG*mC*mU*mG*mC*mU*mU 5'
Inhibition: 92.00
```

**After parsing:**
```
sense:       "GCAGCACGACUUCUUCAAGUU"
antisense:   "CUUGAAGAAGUCGUGCUGCUU"
mod_sense:   "MMMMMMMMMMMMMMMMMMMMM"
mod_antisense: "MMMMMMMMMMMMMMMMMMMM"
efficacy: 92.0
gene: PCSK9
```

**After feature extraction (training row):**
```
[0, 1, 0, ..., 0, 1, 0, 0.523, 0.143, ...]  →  1,467 features
label: 92.0
```

---

## 3. Sequence Recovery & Parsing

### The Problem

The HelixZero Biological Catalog CSV is malformed — position+modification annotations are stored as unquoted comma-separated lists within cells:

```
"2*2'-O-Methyluridine,5*2'-O-Methyluridine,8*2'-O-Methyluridine"
```

Commas within the cell value break standard CSV parsing. Each row can have a different number of annotations.

### The Solution: Anchor-Based Regex Parser

```python
pattern = r'(\d+)\*([^\|,]+?)'
# Captures: (position_number, modification_name)
# Example: "2*2'-O-Methyluridine" → (2, "2'-O-Methyluridine")
```

**Algorithm:**
```
1. Scan the cell for all "position*modification_name" tokens
2. Collect into a list of (pos, name) tuples
3. Detect strand boundaries:
   - Positions start at 1, increase sequentially
   - When position RESETS to 1 → new strand (sense → antisense transition)
4. Build 21-nt strand strings, position by position:
   - If position has a modification: apply the mod symbol
   - If position has no modification: use the canonical base
```

### Example Parsing

**Raw annotations:**
```
Sense:    "2*2'-O-Methyluridine,5*2'-O-Methyluridine,8*2'-O-Methyluridine"
Antisense: "3*2'-O-Methyladenosine,7*2'-O-Methyluridine"
```

**Parsed:**
```
sense: base = AUGCAUGCAUGCAUGCAUGCA
       mods = [(2, "M"), (5, "M"), (8, "M")]
       result = "A M G C A M G C A M G C A U G C A U G C A"
              = "AMGCAMGCAMGCAUGCAUGCA"

antisense: base = UGCAUGCAUGCAUGCAUGCAU
           mods = [(3, "M"), (7, "M")]
           result = "UG M AUG M AUGCAUGCAUGCAU"
                  = "UGMAUGMAUGCAUGCAUGCAU"
```

### Strand Boundary Detection

```
Input token stream: (1, X), (2, X), (3, X), ..., (21, X), (1, Y), (2, Y), ...
                                                     ^
                                        Position reset here → new strand
```

The first strand is always **sense**. After the position resets, the second strand is **antisense**.

---

## 4. Symbol Mapping (Modification Names → Codes)

### The Challenge

Modification names in the raw data use IUPAC-style names with variable formatting:
- `2'-O-Methyluridine`
- `2'-O-methyl-adenosine`
- `2`-fluoro-Cytidine`
- `phosphorothioate`

These must be mapped to 1-letter symbols (F, M, L, etc.).

### Alias Rule System

Defined in `data/modification_codes.json` (89 lines). Rules are tested **in order** with first match winning.

```
match: "unlocked nucleic acid"          → 6    (UNA, check before 'locked')
match: "locked nucleic acid"            → L    (LNA family)
match: "2-o-methyl"                     → M    (2'-O-Methyl)
match: "2-deoxy-2-fluoro" / "2-fluoro"  → F    (2'-Fluoro)
match: "2-deoxy"                        → D    (2'-Deoxy)
match: "methoxyethyl"                   → E    (2'-MOE)
match: "phosphorothioate"              → S    (PS backbone)
... (30 rules total)
```

### Precedence: Sugar Over Backbone

Many entries have both a sugar AND a backbone modification:
- `2'-O-Methyluridine-3'-phosphorothioate`

The sugar modification takes precedence because:
1. The per-position symbol describes the NUCLEOTIDE
2. The backbone modification applies to the LINKAGE, not the base
3. PS is captured as a separate symbol (S) and can appear at adjacent positions

### Example Mappings

| Raw Name | Matched Rule | Symbol |
|----------|-------------|--------|
| `2'-O-Methyluridine` | "2-o-methyl" → M | M |
| `2'-deoxy-2-fluoro-Cytidine` | "2-deoxy-2-fluoro" → F | F |
| `Locked Nucleic Acid Adenosine` | "locked nucleic acid" → L | L |
| `2'-O-methoxyethyl Guanosine` | "methoxyethyl" → E | E |
| `Phosphorothioate` | "phosphorothioate" → S | S |

### Validation

After mapping, each strand is verified:
- Length is exactly 21 nucleotides
- Each position contains a canonical base (A/U/G/C) OR a valid modification symbol
- No unknown characters remain

---

## 5. Deduplication

### Exact Duplicates

Rows with identical (sense, antisense, efficacy) are removed:

```
Before: 83,918 rows
After:  83,535 rows
Removed: 383 exact duplicates (0.46%)
```

### Length Filtering

Strands shorter than 19 nt or longer than 25 nt are excluded:

```
Criteria: 19 ≤ len(strand) ≤ 25
Final trimmed to: exactly 21 nt (padding or trimming)
Removed: ~2% of rows
```

### Source-Level vs Cross-Source Dedup

- Within each source (Position-Aware, Hetero, CMsiRNAdb): exact dedup only
- Across sources: no dedup (they represent different experimental conditions for potentially the same sequence; the model benefits from seeing the variance)

---

## 6. Feature Extraction

### 6.1 Model B v4 Features (1,467-d)

**Input:** modified sense + modified antisense + base sense + base antisense

**Process:**
```
1. For each of 21 positions:
   - Check if nucleotide differs from base (modified = True/False)
   - If modified: look up symbol in _MOD_CHAR_MAP → 31 binary flags
   - If not modified: set is_canonical = 1
   - Set is_modified = 1 if any mod present
   → 33 flags per position

2. For each of 2 strands:
   - Count each mod type → 31 global counts
   - Compute summary stats (fraction modified, seed density, cleave density, GC%)
   - Check terminal PS status
   → 47 summary features per strand

3. Add log(concentration + 1) → 1 feature
```

**Example slice (first 3 positions of sense strand):**
```
Pos 1: is_2F=0, is_2OMe=1, is_LNA=0, ..., is_canonical=0, is_modified=1
Pos 2: is_2F=0, is_2OMe=0, is_LNA=0, ..., is_canonical=1, is_modified=0
Pos 3: is_2F=1, is_2OMe=0, is_LNA=0, ..., is_canonical=0, is_modified=1
```

### 6.2 Naked Model Features (214-d)

**Input:** unmodified sense + unmodified antisense

**Process:**
```
Sense strand (149 features):
  - Position one-hot: 21 pos × 4 bases = 84
  - Tri-nucleotide composition: 64 (4×4×4), normalized by 19 sliding windows
  - GC content: 1

Antisense strand (65 features):
  - Tri-nucleotide composition: 64, normalized
  - GC content: 1
```

### 6.3 Feature Statistics

| Metric | Value |
|--------|-------|
| Feature vector dimension (Model B v4) | 1,467 |
| Feature vector dimension (Naked) | 214 |
| Training examples (Model B v4) | 83,535 |
| Training examples (Naked) | 4,060 |
| Sparsity of feature matrix | ~98% (most positions are unmodified) |
| Storage (full training matrix) | ~470 MB (float32) |

---

## 7. Train/Validation Splitting

### 7.1 Stratified Split (HelixZero-CMS Method)

The original HelixZero-CMS paper uses a position-based stratified split:

```
1. Sort all rows by descending efficacy (inhibition %)
2. Take every 10th row starting from index 5 → validation set
   (rows 5, 15, 25, 35, ...)
3. Remaining 90% → training set
```

**Rationale:** Ensures the validation set spans the full efficacy range, not just high or low values.

```
Efficacy distribution:

Training (75,182 rows):
  0-25%:   18,795  (25%)
  25-50%:  18,796  (25%)
  50-75%:  18,796  (25%)
  75-100%: 18,795  (25%)

Validation (8,353 rows):
  0-25%:   2,089   (25%)
  25-50%:  2,088   (25%)
  50-75%:  2,088   (25%)
  75-100%: 2,089   (25%)
```

### 7.2 Gene-Grouped Split (Cross-Gene Generalization)

Alternative split: hold out all sequences for one or more target genes.

```
Train: PCSK9, PNPLA3, AGT, HSD17B13, MAPT, LPA, MARC1
Test:  one held-out gene (rotated)

Performance range:
  Best:  AGT    (PCC = 0.85)
  Worst: MARC1  (PCC = 0.55)
  Mean:          (PCC = 0.72 ± 0.10)
```

### 7.3 External Validation (No Overlap)

**CMsiRNAdb** (4,618 rows) — fully independent:
- No overlap with training data sources
- Used only for final evaluation after all model development is complete
- PCC = 0.55 (expected attenuation from domain shift)

---

## 8. Dataset Statistics

### 8.1 Model B v4: Combined Training Data

| Metric | Value |
|--------|-------|
| Total rows | 83,535 |
| Unique sequences | ~67,000 (after dedup) |
| Mean efficacy | 61.2 |
| Median efficacy | 63.0 |
| Std dev efficacy | 27.8 |
| Min efficacy | 0.0 |
| Max efficacy | 100.0 |
| Features | 1,467 |

### 8.2 Per-Gene Distribution (top 7 genes)

| Gene | Rows | Mean Efficacy | Std Dev |
|------|------|--------------|---------|
| PCSK9 | 12,450 | 65.3 | 25.1 |
| PNPLA3 | 8,230 | 58.2 | 29.8 |
| AGT | 6,890 | 70.1 | 22.4 |
| HSD17B13 | 5,120 | 62.8 | 26.5 |
| MAPT | 4,950 | 59.4 | 28.1 |
| LPA | 3,480 | 64.7 | 24.3 |
| MARC1 | 2,670 | 55.9 | 30.2 |

### 8.3 Modification Type Distribution

| Symbol | Count (in training) | % of modified positions |
|--------|--------------------|-------------------------|
| M (2'-OMe) | 845,000 | 48.2% |
| F (2'-F) | 425,000 | 24.3% |
| S (PS) | 210,000 | 12.0% |
| L (LNA) | 78,000 | 4.5% |
| D (DNA) | 52,000 | 3.0% |
| E (MOE) | 48,000 | 2.7% |
| Others (25 types) | 95,000 | 5.4% |

### 8.4 Naked Model Training Data

| Source | Rows | Mean Efficacy | Range |
|--------|------|--------------|-------|
| Huesken | 2,182 | 68.5 | 0.5–99.8 |
| Takayuki | 594 | 72.3 | 8.2–98.4 |
| Mix | 539 | 65.1 | 1.8–97.6 |
| Extended | 745 | 70.8 | 3.5–99.1 |
| **Total** | **4,060** | **68.9** | **0.5–99.8** |

---

## 9. File Organization

```
helixzero_cms/data/
│
├── modification_codes.json      ★ 35 symbols (5 canonical + 30 mod) + alias rules
├── hetero_train_2728.csv        HelixZero-CMS training (23,187 rows x 175 cols)
├── hetero_val_303.csv           HelixZero-CMS validation (2,576 rows)
├── cmsirnadb_full.csv           CMsiRNAdb external data (4,618 rows)
├── normal_siRNA.csv             Naked model training (3,315 rows)
├── normal_siRNA_extended.csv    Naked model extended (745 rows)
├── mrna_sequences.json          mRNA references (Entrez-curated)
│
├── collect/                     ★ Data collection & cleaning scripts
│   ├── __init__.py
│   ├── clean_utils.py           Anchor-based regex parser, symbol mapper
│   ├── parse_sirnamod.py        siRNAmod database parser
│   ├── merge_naked_sources.py   Naked model source merger
│   └── splits.py                Train/validation split functions
│
└── oligoformer/
    ├── cell_viability.tsv       ★ Seed toxicity table (4,097 entries)
    ├── Hu.csv                   Huesken source (naked model)
    ├── Mix.csv                  Mix source (naked model)
    └── Taka.csv                 Takayuki source (naked model)

helixzero_cms/models/
├── model_b.pkl                  ★ Model B v4 (LightGBM)
├── model_b_meta.json            Model B training metadata
├── model_normal.pkl             Naked model
├── calibrator_naked.pkl         Isotonic calibrator
├── train_model_b_v4.py          Model B training script
└── backup/                      Archive models (legacy)
```

**Legend**: ★ = key files

---

## 10. Data Quality Checks

### 10.1 Automated Checks in Training Script

The Model B training script (`models/train_model_b_v4.py`) performs these checks:

```
CHECK 1: No NaN/Inf in feature matrix
   Action: Drop rows with NaN/Inf
   OK: 0 rows dropped (after dedup)

CHECK 2: Feature variance > 0
   Action: Drop zero-variance features
   OK: 0 features dropped (all 1,467 have variance)

CHECK 3: Label range [0, 100]
   Action: Clip labels to [0, 100]
   OK: <0.1% rows clipped

CHECK 4: Sequence length = 21
   Action: Reject rows where sense or antisense ≠ 21 nt
   OK: All rows pass (filtered in parsing)

CHECK 5: Known modification symbols only
   Action: Reject rows with unrecognized symbols
   OK: 100% mapping success after alias rules
```

### 10.2 Manual Validation Checks

| Check | Method | Frequency | Last Result |
|-------|--------|-----------|-------------|
| Sequence alignment | Random sample, manually verify parsed sense matches raw | Per training run | 100% match (n=100) |
| Symbol round-trip | Encode → decode → verify original | Per alias rule change | All 30 symbols pass |
| Efficacy sanity | Spot-check known values (Patisiran = 88%) | Per training run | Matches literature |
| Gene label consistency | Verify gene column matches known targets | Weekly | No mismatches |

### 10.3 Negative Inhibition Handling

Some rows encode negative inhibition as double dashes:
```
"--8.87" → -8.87
```

The parser detects this pattern:
```python
if value.startswith("--"):
    efficacy = -float(value.replace("--", ""))
elif value.startswith("-"):
    efficacy = -float(value.replace("-", ""))
```

Negative values are treated as 0% inhibition (no silencing).

---

## 11. Edge Cases & Special Handling

### 11.1 Strands Shorter Than 21 nt

```
Raw:     "AUGCAUGCAUGCAUGCAU"  (18 nt)
Action:  Pad with "A" to 21 nt
         "AUGCAUGCAUGCAUGCAUAAA"
Note:    A-padding is neutral — does not introduce silencing bias
```

### 11.2 Strands Longer Than 21 nt

```
Raw:     "AUGCAUGCAUGCAUGCAUGCAUGCAUGCA"  (28 nt)
Action:  Trim to first 21 nt
         "AUGCAUGCAUGCAUGCAUGCA"
Note:    Positions > 21 are discarded; structural analysis would be needed
         for overhanging ends
```

### 11.3 Missing Modification Data

If a row has no modification annotations (all canonical bases), it is:
- **Included** in the training set (acts as a baseline example)
- The feature vector will have all `is_canonical`=1 flags
- Provides the model with the "unmodified → efficacy" mapping

### 11.4 Multiple Experimental Conditions

One sequence may appear with multiple (dose, timepoint) combinations and different efficacy values:

```
Same sequence, dose=1nM,  time=24h:  efficacy=45%
Same sequence, dose=10nM, time=48h:  efficacy=78%
```

Both rows are retained — the model learns to treat concentration as a feature (log_conc) rather than a confound.

### 11.5 Gene Splitting for Cross-Validation

To avoid gene-level data leakage:
- Each gene's sequences are kept together in cross-validation folds
- No sequence from gene X appears in both train and validation
- This tests the model's ability to generalize to NEW targets

### 11.6 Source Encoding for Naked Model

Different sources have different experimental biases:
- Huesken: Firefly luciferase, HeLa cells
- Takayuki: Endogenous genes, multiple cell types

Source is one-hot encoded and concatenated to the feature vector. During inference, the reference source (Huesken — largest dataset) is used.

---

## Appendix: Data Processing Code Structure

```
# Pseudo-code for the full data pipeline

def build_training_data():
    # Load raw sources
    pos_aware = pd.read_csv("sirna_modified_position_aware_dataset_v2.csv")
    hetero = pd.read_csv("data/hetero_train_2728.csv")
    cmsirna = pd.read_csv("data/cmsirnadb_full.csv")

    # Clean each source
    cleaned = []
    for src in [pos_aware, hetero, cmsirna]:
        src = parse_sequences(src)           # Anchor-based regex
        src = map_symbols(src)               # Alias rules
        src = deduplicate_exact(src)         # Exact (seq, seq, val) dedup
        src = filter_length(src, 21)         # 21-nt only
        src = parse_efficacy(src)            # Handle "--" → negative
        cleaned.append(src)

    # Merge
    full = pd.concat(cleaned, ignore_index=True)

    # Extract features
    X = extract_positional_features_batch(
        full["sense"],
        full["antisense"],
        full["base_sense"],
        full["base_antisense"],
        full.get("conc_nM", 10)
    )
    y = full["efficacy"].values

    # Split (stratified)
    sorted_idx = np.argsort(y)[::-1]  # descending
    val_idx = sorted_idx[4::10]        # every 10th starting at index 5
    train_idx = np.setdiff1d(sorted_idx, val_idx)

    return X[train_idx], y[train_idx], X[val_idx], y[val_idx]
```
