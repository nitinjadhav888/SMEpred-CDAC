# SMEpred — Data Collection & Engineering (Plain English)

This document explains, step by step, how the two raw siRNA datasets were collected,
cleaned, transformed, normalized, and turned into the CSVs the models train on. It also
reports the honest model results and *why* they land where they do.

---

## 1. The Two Source Files

Both files arrived with a misleading `.csv.xls` extension but are **plain CSV text**
(verified by the file's magic bytes).

### File 1 — `Data-(1)-csv.csv.xls`  (siRNAmod export, AUXILIARY)
- 5,329 rows × 7 columns: `siRNAmodDB, PMID, sense seq, sense modification name,
  antisense seq, antisense modification name, inhibition %`.
- Efficacy is free text, e.g. `"69.5 percent target mRNA inhibition"`.
- Multiple modifications are separated by `*`. **No position information.**
- **Role:** because it lacks positions, it cannot build the positioned strands the models
  need. We use only its genuinely *unmodified* rows to enlarge the normal-siRNA dataset.

### File 2 — `HelixZero_Biological_Catalog_43k.csv.xls`  (PRIMARY)
- 43,467 data rows × 26 columns. Patent-derived siRNAs.
- Crucially contains **per-position modifications** as `position*name || position*name`,
  e.g. `1*2'-O-Methylcytidine || 2*2'-O-Methyladenosine || …`.
- Inhibition is encoded in the row ID suffix (e.g. `…-48h-88.00` → 88.0) and verified
  against the `Inhibition` column.
- **Role:** primary source for all chemically-modified siRNA models (A, B, C).

---

## 2. The Parsing Problem (and the fix)

File 2 is a **malformed CSV**: the `Modification_locations_*` / `position_*` columns hold
unquoted comma-lists (`1,2,3,4,…`), so the number of comma-fields varies per row (only
~3,500 of 43k rows have the nominal 26 fields). `pandas.read_csv` cannot parse it.

**Solution — an anchor/regex line parser** (`data/collect/parse_helix_catalog.py`):
- We verified across 857k tokens that **no modification name contains a comma**, so the
  token regex `(\d+)\*([^|,]+?)(?=\s*\|\||,)` cleanly extracts every `(position, name)`.
- Strands are separated by detecting the **position reset** (positions run 1..N for the
  antisense strand, then restart at 1 for the sense strand; antisense appears first).
- From each token we derive two things:
  - the **canonical base** (from the name's base suffix: `…uridine`→U, `…adenosine`→A, …)
  - the **35-alphabet symbol** (via the alias map; canonical positions keep their base)

This reconstructs both the unmodified base strand and the modified-symbol strand without
ever relying on the broken sequence columns.

### A subtle but critical bug we caught
Negative inhibition is written with a **double dash**: `…-24h--8.87`. A naive regex reads
`8.87` (positive). We fixed the regex to capture the sign (`-(-?\d+\.?\d*)$`); negatives
are then clipped to 0 per the cleaning rule. Without this, ~19k poor/non-functional
siRNAs were mislabeled as functional — which alone wrecks the model.

---

## 3. Modification Name → Symbol Mapping

Real names are *compositional* (`<modification-class><base>`, e.g. "2'-O-Methyluridine").
File 1 has 307 distinct names; File 2 has 27,924 (because each base gives a variant).

Per the paper, **each modification class maps to one symbol regardless of base**
("2'-O-Methyluridine" and "2'-O-Methyladenosine" both → `M`). The mapping lives in
`data/modification_codes.json` as an ordered `alias_rules` list (first substring match
wins; sugar modifications take precedence over backbone). Examples:
`2'-O-Methyl→M, 2'-Fluoro→F, 2'-Deoxy→D, Locked nucleic acid→L, GalNAc→4,
Phosphorothioate→S (fallback)`.

Names matching no rule are logged to `data/collect/unmapped_report.txt` (only **11**
distinct remain, all rare/junk) and that position is treated as canonical.

---

## 4. Cleaning, Transformation, Normalization, Standardization

Shared helpers live in `data/collect/clean_utils.py`.

| Step | Rule |
|------|------|
| **Sequence cleaning** | Uppercase; convert DNA `T`→RNA `U`; keep only `A/U/G/C`; drop junk. |
| **Length filter** | Keep strand lengths 19–25 nt (paper used 21–24). |
| **Efficacy parse** | Extract the number from free text / ID suffix. |
| **Efficacy clip** | Negative → 0 (no silencing); cap at 100. Bounds target to [0,100]. |
| **Drop invalid** | Rows with unparseable efficacy or invalid sequences are removed. |
| **Deduplication** | Drop exact duplicates on `(sense, antisense, efficacy)`. |
| **Standardization** | `StandardScaler` (mean 0 / std 1 per feature) inside every pipeline. |
| **Score normalization** | At inference, LightGBM output is clipped to [0, 100]. |

**Train / validation split** (`data/collect/splits.py`, `paper_split`): sort by
**descending efficacy**, then take **every 10th row starting at the 5th** as the
independent validation set. This is the paper's exact rule and guarantees the validation
set spans the full efficacy range.

---

## 5. Feature Representation (and why it changed)

The paper's MNC placed modification symbols at known positions to form a 35-symbol
sequence. That works when modifications are *sparse* (canonical bases still dominate).

**This patent data is near-fully modified**, so the 35-symbol composition collapses to
mostly `M`/`F` and **loses the underlying base sequence** — the strongest efficacy signal.
We measured this directly on the dose-controlled set:

| Features | PCC |
|----------|-----|
| modified-symbol MNC only (original) | 0.37 |
| base-sequence MNC only | 0.40 |
| **base + modified MNC (adopted)** | **0.48** |

So every model now uses **140-d features = MNC(base strands) ⊕ MNC(modified strands)**;
All LightGBM models (A/B/C) use the same 152-d feature set internally — the A/B/C
distinction is preserved in the UI only for backward compatibility. For unmodified
siRNAs the two halves coincide.

---

## 6. Output Datasets

Produced by running the two parsers (File 2 first, then File 1):

| File | Rows | Purpose |
|------|------|---------|
| `data/hetero_train_2728.csv` | 23,187 | cm-siRNA training (all conditions mixed) |
| `data/hetero_val_303.csv` | 2,576 | independent validation |
| `data/homo_train.csv` | 4,244 | **dose-controlled** (10 nM / 24 h) training |
| `data/homo_val.csv` | 472 | dose-controlled validation |
| `data/normal_siRNA_extended.csv` | 4,060 | unmodified siRNAs (4 merged sources: Huesken, Mix, Takayuki, internal) for the ranker |

---

## 7. Model Results (LightGBM v3 — current)

Metrics from the deployed models (trained by `python models/train_gbm.py`):

| Model | Split | PCC | Spearman | MAE |
|-------|-------|-----|----------|-----|
| **cm-siRNA** (25,763 rows) | Random 82/18 | **0.6789** | 0.6736 | 16.42 |
| **Naked** (4,060 rows, 4 sources) | Random 82/18 | **0.5543** | 0.5470 | 13.42 |
| **Naked per-source**: Taka (699) | Within-source | 0.6905 | — | — |
| **Naked per-source**: Huesken (2,361) | Within-source | 0.4179 | — | — |
| **Homo** (dose-controlled, 4,716 rows) | Random split | **0.7370** | — | — |

**Why below the original SMEpred paper's 0.80:**
1. **Dataset size** — the paper trained on 2,728 hand-curated rows (low noise). Our dataset
   is 25,763 rows from real pharma patents (heterogeneous doses, cell lines, timepoints).
2. **Dose confound** — the same siRNA has very different knockdown at 0.1 nM vs 100 nM.
   We mitigate by feeding dose + time as features, but residual noise remains.
3. **Full modification** — near-fully-modified strands collapse composition to mostly M/F,
   losing sequence signal. We solve this by feeding *both* base and modified MNC (140-d),
   which raised PCC from 0.37 → 0.48 on early tests, with LightGBM taking it further to 0.68.

**Recommendation:** within-siRNA modification ranking (PCC 0.68) is strong and production-ready. The naked model (PCC 0.55) is suitable for baseline ranking.

---

## 8. How to Reproduce

```bash
# 1. Parse the primary catalog (writes hetero/homo/normal CSVs + unmapped report)
python data/collect/parse_helix_catalog.py

# 2. Merge File 1's unmodified rows into the normal dataset
python data/collect/parse_sirnamod.py

# 3. Train all models on the real data
python models/train_gbm.py

# 4. Verify
python tests/test_pipeline.py
python cli/run.py rank --sequence <mRNA-or-FASTA>
```

File paths are currently hard-coded to `D:\Helixx\...`; edit the `CATALOG` / `SIRNAMOD`
constants at the top of the two parser scripts if the files move.
