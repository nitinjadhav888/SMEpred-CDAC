# HelixZero-CMS — Implementation Guide

**Plain English guide to every algorithm, concept, and design decision.**

> **Production models use LightGBM** with position-aware + assay-condition features.
> Train with `python models/train_gbm.py`.

---

## What This System Does (One Paragraph)

You give it an mRNA or gene sequence. It chops that sequence into every possible
21-nucleotide fragment and predicts how well each one would work as an siRNA
(a molecule that silences the gene). It ranks those candidates. Then, for any
candidate you pick, it systematically applies 30 types of chemical modifications
at every position (1260 combinations), predicts which modified version works best,
and returns a ranked list. You can also manually specify custom modification
patterns and get a score for that exact design.

---

## Background — What is siRNA and Why Modify It?

### siRNA (Small Interfering RNA)
A siRNA is a short double-stranded RNA molecule, 21 nucleotides long. One strand
(the "antisense" or "guide" strand) is loaded into a protein complex called RISC.
RISC then searches for matching mRNA in the cell, binds to it, and destroys it,
preventing the gene from being expressed. This is called RNA interference (RNAi).

### The Problem with Unmodified siRNAs
Natural siRNAs get rapidly degraded in the bloodstream by enzymes called nucleases.
They also trigger the immune system, get filtered out before reaching target cells,
and sometimes hit unintended mRNAs (off-targets). These problems block their use
as medicines.

### Chemical Modifications
Chemists can swap out or alter individual atoms in the siRNA backbone, sugar, or
bases to make the molecule more stable, less immunogenic, and more cell-specific.
Examples: replacing the 2'-OH group with Fluorine (2'-F), adding a methyl group
(2'-OMe), or linking ribose rings together (LNA). The challenge is that not every
modification helps — some reduce efficacy. There are 30 commonly used modifications,
and each can be placed at 21 positions on 2 strands = 1260 options per siRNA. This
system predicts which ones will work before you spend money testing them in a lab.

---

## Algorithm 1 — Sequence Parsing (src/parser.py)

### What it does
Takes raw input (file path, FASTA text, or inline sequence string) and returns a
clean RNA sequence.

### FASTA Format
FASTA is the standard file format for biological sequences. It looks like:
```
>gene_name  ← header line starting with ">"
AUGCAUGCAUGCAUGCAUGCA  ← the sequence (can span multiple lines)
```

### DNA to RNA conversion
DNA uses the base Thymine (T). RNA uses Uracil (U). Biologically they are the same
position — when you read a gene file, it's often DNA format. We replace all T with U.

### Validation
Any character other than A, U, G, C (or T for DNA input) is rejected. Sequences
shorter than 21 nt are rejected (can't generate even one siRNA).

---

## Algorithm 2 — 21-mer siRNA Generation (src/sirna_generator.py)

### Sliding Window
We move a window of width 21 across the sequence one position at a time:
```
Position 0:  AUGCAUGCAUGCAUGCAUGCA  (nucleotides 0–20)
Position 1:  UGCAUGCAUGCAUGCAUGCAU  (nucleotides 1–21)
Position 2:  GCAUGCAUGCAUGCAUGCAUG  (nucleotides 2–22)
...
```
For a 1000-nt mRNA, this produces 980 candidates.

### Antisense Derivation (Reverse Complement)
The antisense strand is derived from the sense strand in two steps:
1. Complement each base: A↔U and G↔C
2. Reverse the result

Example:
```
Sense:     5' - GCAGCACGACUUCUUCAAGUU - 3'
Step 1 complement:   CGUCGUGCUGAAGAAGUUCAA
Step 2 reverse:      AACUUGAAGAAGUCGUGCUGC  ← this is the antisense (guide strand)
```
The antisense strand binds to the mRNA target and guides cleavage.

---

## Algorithm 3 — Feature Extraction (src/features.py)

The model needs numbers, not letters. We convert each siRNA into a numerical vector
using three strategies from the paper.

### Feature Type 1: MNC (Mononucleotide Composition)

**What it is:** The fraction of each nucleotide type in the sequence.

**Why it works:** The composition of a siRNA strongly affects its efficacy. For example:
- High AU content at the 5'-end of the antisense strand helps RISC load the right strand.
- Extreme GC content (>70% or <30%) tends to reduce efficacy.
- The type of chemical modification present tells the model about stability changes.

**How we compute it:**
1. Define an alphabet of 35 symbols: 5 canonical bases (A, U, G, C, T) + 30 modification symbols.
2. Count how many times each symbol appears in the sense strand (21 positions).
3. Divide by strand length to get frequency (0.0 to 1.0).
4. Do the same for the antisense strand.
5. Concatenate both frequency vectors → 70-dimensional vector.

**Example:**
```
Sense strand: GCAGCACGACUUCUUCAAGUU (21 nt)
G count = 4 → G frequency = 4/21 = 0.190
C count = 5 → C frequency = 5/21 = 0.238
A count = 4 → A frequency = 4/21 = 0.190
U count = 5 → U frequency = 5/21 = 0.238
(all other symbols = 0 for unmodified siRNA)
```

**Paper result:** MNC alone achieved PCC = 0.80 (the best single feature). This is Model-A.

---

### Feature Type 2: DNC (Dinucleotide Composition)

**What it is:** The fraction of each consecutive two-nucleotide pair in the sequence.

**Why it's useful:** Captures nearest-neighbor stacking energy effects — the stability
of any RNA duplex depends on how pairs of adjacent bases stack on top of each other.

**How we compute it:**
1. Slide a window of width 2 across the sequence.
2. Count each of the 35×35 = 1225 possible pairs.
3. Divide by (length − 1) for frequency.
4. Concatenate sense + antisense → 2450-dimensional vector.

**Paper result:** DNC alone gave PCC = 0.53 — worse than MNC, because absolute position
matters more than pair frequencies for modified siRNAs.

---

### Feature Type 3: BIN (Binary Position Pattern)

**What it is:** Encodes *where* each nucleotide type appears in the sequence,
not just how often.

**Why it matters:** Position is critical in siRNA biology:
- The "seed region" (positions 2–8 from the 5'-end of the antisense strand) is where
  RISC first contacts the mRNA — modifications here have the biggest impact.
- The 3'-end affects strand loading specificity.
- Chemical modifications at the first and last positions affect stability differently
  than the same modification in the middle.

**How we compute it:**
For each of the 35 nucleotide types, write a 1 at every position in the sequence
where that type appears, and 0 everywhere else.
Result: a binary vector of length 35 × sequence_length.

**Example for a 4-nt sequence "AUGC" with only canonical bases:**
```
A: [1, 0, 0, 0]  ← A is only at position 1
U: [0, 1, 0, 0]
G: [0, 0, 1, 0]
C: [0, 0, 0, 1]
... (all 31 modification symbols get [0,0,0,0])
Flat vector: [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1, 0,...0]  (35×4 = 140 values)
```

**Paper BIN variants used in models:**
- **Model-B**: MNC (70-d) + BIN of antisense seed region (first 13 nt from 5'-end) = 70 + 455 = 525-d
  → Best independent validation PCC = 0.86
- **Model-C**: MNC (70-d) + BIN of antisense last 8 nt from 3'-end = 70 + 280 = 350-d
  → Independent validation PCC = 0.78

---

## Algorithm 4 — LightGBM Training (models/train_gbm.py)

### What is LightGBM (Gradient Boosting Machine)?

LightGBM is a gradient-boosted decision tree algorithm. It builds an ensemble of
shallow decision trees sequentially — each new tree corrects the errors of all previous
trees. It handles non-linear relationships, mixed feature types, and scales to tens
of thousands of rows.

### Key hyperparameters
- **n_estimators=800**: number of trees in the ensemble
- **learning_rate=0.03**: how much each tree contributes
- **num_leaves=63**: maximum leaves per tree (controls model capacity)
- **subsample=0.8**: row sampling per tree (prevents overfitting)
- **colsample_bytree=0.7**: feature sampling per tree
- **reg_lambda=1.0**: L2 regularization

### Feature set (152-d)
The LightGBM uses richer features than the original MNC:
- **140-d**: base + modified MNC (70 per strand)
- **8-d**: position-aware modification density (4 per strand)
- **2-d**: GC content per strand
- **2-d**: experimental condition (dose log10, time/24)

### Gene-Grouped Split
The model is evaluated by holding out entire target genes — no gene appears in both
train and test. This gives an honest estimate of performance on **never-before-seen
targets**.

### Pearson Correlation Coefficient (PCC)
PCC measures the linear relationship between two lists of numbers:
- PCC = 1.0 → perfect positive correlation
- PCC = 0.0 → no correlation
- PCC = −1.0 → perfect inverse correlation

### Three Models Saved
- `model_a.pkl` — LightGBM (cm-siRNA efficacy, default)
- `model_b.pkl` — LightGBM (same model, for compatibility)
- `model_c.pkl` — LightGBM (same model, for compatibility)
- `model_normal.pkl` — LightGBM for ranking unmodified siRNA candidates

---

## Algorithm 5 — Modification Engine (src/modification_engine.py)

### Mode 1: Single-Modification Scan

For each siRNA, we generate every combination of:
- 30 modification types (F, M, L, D, E, B, ...)
- 21 positions (1 to 21)
- 2 strands (sense and antisense)

Total: 30 × 21 × 2 = **1260 cm-siRNAs**.

Each cm-siRNA is the parent sequence with exactly one character swapped at one position.

**Example:**
```
Parent sense: GCAGCACGACUUCUUCAAGUU
Apply F (2'-Fluoro) at position 3:
Modified:     GCFGCACGACUUCUUCAAGUU
               ↑ position 3 replaced with F
```

### Mode 2: MultiModGen

User specifies modifications directly. The input format uses `,,` to separate
different modification types:
```
Modification symbols:  "F,,M"        → type 1 = F,  type 2 = M
Positions:             "2,5,,10,12"  → F at 2,5  and  M at 10,12
```

Multiple modifications at different positions are applied simultaneously to produce
one custom cm-siRNA design.

---

## Algorithm 6 — Score Normalization (src/predictor.py)

The LightGBM models predict inhibition directly on the 0–100 scale, so we simply
clip values into range rather than min-max rescaling per batch.

**Efficacy labels:**
- ≥ 90 → "Very High"
- 80–89 → "High"
- 70–79 → "Moderate"
- < 70 → "Low"

**Delta score:** The difference between a cm-siRNA's score and its parent unmodified
siRNA's score. Positive delta = the modification improves efficacy. Negative delta =
the modification hurts efficacy.

---

## System Architecture (End-to-End Flow)

```
INPUT: mRNA FASTA file / sequence string
         │
         ▼
  [parser.py]  →  clean RNA sequence (uppercase, T→U)
         │
         ▼
  [sirna_generator.py]  →  list of (N-20) SiRNACandidate objects
    each has: position, sense (21 nt), antisense (21 nt)
         │
         ▼
  [features.py]  →  numerical feature vectors
    Model A: 70-d MNC
    Model B: 525-d MNC + BIN seed
    Model C: 350-d MNC + BIN tail
         │
         ▼
  [model_normal.pkl]  →  efficacy score for each unmodified siRNA
         │
         ▼
  [predictor.py rank_sirnas()]  →  RankedSiRNA list (best → worst)
         │
  USER SELECTS ONE siRNA
         │
         ▼
  [modification_engine.py]  →  1260 cm-siRNA variants
         │
         ▼
  [features.py]  →  feature vectors for each cm-siRNA
         │
         ▼
  [model_a/b/c.pkl]  →  predicted efficacy for each cm-siRNA
         │
         ▼
  [predictor.py predict_modified()]  →  RankedCmSiRNA list with delta scores

OUTPUT: ranked tables (CSV or JSON) or REST API response
```

---

## Data Requirements

The models need experimentally validated training data. The paper used:

| Dataset      | Sequences | Description                          |
|--------------|-----------|--------------------------------------|
| Hetero-T2728 | 2728      | cm-siRNAs for training (10-fold CV)  |
| Hetero-V303  | 303       | Held-out independent validation set  |
| Homo-T1900   | 1900      | Same-condition siRNAs (supplementary)|
| Normal-2182  | 2182      | Unmodified siRNAs for normal ranking |

**Source:** siRNAmod database — Dar et al., Scientific Reports 2016.
Download at: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4741272/

**CSV format required** (one row per siRNA):
```
sense,antisense,efficacy
GCAGCACGACUUCUUCAAGUU,CUUGAAGAAGUCGUGCUGCUU,85.3
...
```

If no real data files are present, the training script generates synthetic data as a
placeholder. Synthetic models will NOT have meaningful predictive power.

---

## File Structure Reference

```
smepred/
├── data/
│   ├── modification_codes.json    ← 35 nucleotide symbols with binary codes
│   ├── hetero_train_2728.csv      ← (you provide) cm-siRNA training data
│   ├── hetero_val_303.csv         ← (you provide) independent validation data
│   └── normal_siRNA_2182.csv      ← (you provide) unmodified siRNA data
│
├── src/
│   ├── parser.py                  ← Sequence input (FASTA, file, inline)
│   ├── sirna_generator.py         ← 21-mer sliding window + antisense derivation
│   ├── features.py                ← MNC, DNC, BIN feature extraction
│   ├── modification_engine.py     ← 1260-variant scan + MultiModGen
│   └── predictor.py               ← Orchestrates full pipeline
│
├── models/
│   ├── train_gbm.py               ← Training script
│   ├── model_a.pkl                ← Trained LightGBM Model-A (after running train)
│   ├── model_b.pkl                ← Trained LightGBM Model-B
│   ├── model_c.pkl                ← Trained LightGBM Model-C
│   ├── model_normal.pkl           ← Normal siRNA ranker
│   └── model_homo.pkl             ← Dose-controlled model
│
├── api/
│   └── main.py                    ← FastAPI REST server
│
├── cli/
│   └── run.py                     ← Command-line interface
│
├── tests/
│   └── test_pipeline.py           ← Unit tests
│
├── docs/
│   └── IMPLEMENTATION.md          ← This document
│
└── requirements.txt
```

---

## Quick Start Commands

```bash
# 1. Install dependencies
pip install -r requirements.txt

# 2. Train models (uses synthetic data if real CSVs not present)
python models/train_gbm.py

# 3. Run unit tests
python tests/test_pipeline.py

# 4. Rank siRNAs from a sequence
python cli/run.py rank --sequence AUGCAUGCAUGCAUGCAUGCAUGCAUGCAUG --top 10

# 5. Rank from FASTA file
python cli/run.py rank --input gene.fasta --top 20 --output results.csv

# 6. Scan 1260 modifications for a chosen siRNA
python cli/run.py single-mod \
  --sense GCAGCACGACUUCUUCAAGUU \
  --antisense CUUGAAGAAGUCGUGCUGCUU \
  --model A --top 20

# 7. Custom multi-modification prediction
python cli/run.py multi-mod \
  --sense GCAGCACGACUUCUUCAAGUU \
  --antisense CUUGAAGAAGUCGUGCUGCUU \
  --sense-mods F,,M --sense-positions 2,5,,10,12 \
  --model A

# 8. Start REST API server
uvicorn api.main:app --reload --port 8000
# Visit http://localhost:8000/docs for Swagger UI
```

---

## Model Performance Reference (from Paper)

| Model | Features                         | Training PCC | Independent PCC |
|-------|----------------------------------|-------------|----------------|
| A     | MNC (70-d)                       | 0.80        | 0.80           |
| B     | MNC + BIN antisense seed 13 nt   | 0.77        | 0.86           |
| C     | MNC + BIN antisense last 8 nt    | 0.76        | 0.78           |
| Normal| MNC hybrid (unmodified siRNAs)   | 0.72        | —              |

PCC = Pearson Correlation Coefficient. Values closer to 1.0 = better predictions.

---

## How to Add Real Training Data

1. Download cm-siRNA data from the siRNAmod database.
2. Format it as CSV with three columns: `sense,antisense,efficacy`.
   - `sense`    : 21-nt sense strand (RNA, uppercase)
   - `antisense`: 21-nt antisense strand
   - `efficacy` : percentage inhibition (0–100)
3. Save as:
   - `data/hetero_train_2728.csv` (training set, ~2728 rows)
   - `data/hetero_val_303.csv`    (validation set, ~303 rows)
   - `data/normal_siRNA_2182.csv` (unmodified siRNAs, ~2182 rows)
4. Run `python models/train_gbm.py` again.
5. The `.pkl` model files will be overwritten with properly trained models.
