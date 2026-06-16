# HelixZero-CMS — Full Explanation of All Parameters

## 1. Efficacy Score (LightGBM Column in Rank Tab)

### What It Is

The **Efficacy Score** (0–100) is a predicted **% target-gene silencing** (knockdown) for an siRNA candidate. A score of 80 means the model predicts that siRNA will reduce the target mRNA level by ~80% relative to a negative control.

### How It Is Calculated — The Full Pipeline

```
mRNA sequence
    ↓
Generate 21-mer siRNA candidates (sliding window)
    ↓
Extract 152-dimensional feature vector per candidate
    ↓
LightGBM model (gradient-boosted decision trees) predicts raw score
    ↓
Isotonic calibrator maps raw score → calibrated 0–100 score
    ↓
Efficacy label assigned based on thresholds
```

### The 152-D Feature Vector (What Goes Into the Model)

The model does **not** see sequence letters directly — it sees numerical features derived from the sequence. Every siRNA is converted to 152 numbers:

| Feature Group | Dimensions | What It Captures |
|---------------|-----------|------------------|
| **Base MNC (sense + antisense)** | 70 | Mononucleotide composition of the *unmodified* strands — the fraction of each of the 35 nucleotide symbols (A, U, G, C, T + 30 chemical modification symbols) |
| **Modified MNC (sense + antisense)** | 70 | Mononucleotide composition of the *modified* strands — same 35-bin frequency vector, computed on the actual modified sequence |
| **Mod density (sense)** | 4 | How heavily modified the sense strand is: overall mod fraction, seed-region (pos 1–8) mod fraction, 3'-tail mod fraction, count of modified positions |
| **Mod density (antisense)** | 4 | Same as above for the antisense strand |
| **GC content** | 2 | GC fraction of the unmodified sense and antisense strands |
| **Assay conditions** | 2 | log₁₀(concentration in nM) and time in hours (default: 10 nM, 24 h at inference) |

**Why both base and modified MNC?** Near-fully-modified strands collapse to mostly M or F composition, losing the underlying sequence. Adding the base (unmodified) composition restores that sequence signal. Empirically, this raises PCC from 0.37 (modified-only) to 0.48 (base + modified).

### The Model: LightGBM

**LightGBM** is a **gradient-boosted decision tree** algorithm — an ensemble of hundreds of decision trees where each tree corrects the errors of the previous ones.

- **Algorithm**: Gradient-boosted decision trees (GBDT) with leaf-wise tree growth
- **Training data**: 25,765 modified siRNA sequences with measured % inhibition values
- **Data source**: Patent data (HelixZero catalog), covering 13 target genes
- **Trees**: 799 (set by early stopping on a 5% holdout)
- **Performance**: PCC = **0.68** within genes, MAE = **16.4** points

**Why not a neural network / deep learning?** Gradient-boosted trees are the gold standard for tabular data with noisy, mixed-type features. They handle missing values natively, are resistant to outliers, and train quickly without GPU. The 152-d feature vector is engineered RNA chemistry knowledge — the tree model extracts the interactions.

### The Calibrator (Isotonic Regression)

The raw LightGBM output tends to be **compressed** toward the mean of the training labels (~47.7). This means most raw predictions fall in the 30–70 range even though true labels span 0–100.

An **isotonic regression** model learns a monotonic mapping from raw predictions to true labels:

```python
IsotonicRegression(out_of_bounds='clip', y_min=0, y_max=100)
```

It is fitted using **5-fold cross-validation** predictions (to prevent overfitting), then applied at inference time. This spreads the scores to fill the full 0–100 range while preserving rank order.

**Effect**: Raw PCC=0.68 → Calibrated PCC=0.68 (rank order preserved), but absolute scores become more realistic (e.g. raw 67 → calibrated 80).

---

## 2. Seed Toxicity

### What It Means

**Seed toxicity** predicts the risk of **off-target cell toxicity** caused by the siRNA's seed region (antisense nucleotides 2–7). This is a miRNA-like off-target effect: the seed can bind to unintended mRNAs and silence them, potentially causing cell death.

The score is **predicted cell viability (%)** — higher is safer.

| Label | Viability Range | Meaning |
|-------|----------------|---------|
| **Safe** | ≥70% | Low off-target toxicity risk |
| **Caution** | 50–69% | Moderate risk — consider seed-rescue modifications |
| **Toxic** | <50% | High risk — avoid this seed region |
| **Unknown** | N/A | Seed sequence not in the reference database |

### How It Is Calculated

```
Antisense strand (21 nt)
    ↓
Extract seed = positions 2–7 (6-mer)
    ↓
Look up in seed-toxicity table (4,097 entries)
    ↓
Return viability % or None if not found
```

**Step 1 — Seed extraction:**

```python
antisense = AUAUUCACUAAACGACUGCTT
               ↑↑↑↑↑↑
seed = UAUUCA    (positions 2–7, 0-indexed: [1:7])
```

The seed region (positions 2–7 of the guide strand) is the primary determinant of miRNA-like off-target binding (Janas et al., Mol Cell 2018; Burchard et al., 2009).

**Step 2 — Lookup table:**

The table (`data/oligoformer/cell_viability.tsv`) contains experimentally measured cell viability for 4,096 possible 6-mer RNA seeds (Janas et al., Mol Cell 2018):

```
Seed    cell_viability
GGUGGG  4.9     ← highly toxic
GCUAAC  5.2     ← highly toxic
...
CUGGGC  99.2    ← safe
```

The data comes from published siRNA off-target profiling experiments where each seed was transfected into cells and viability was measured by ATP-based assays (e.g. CellTiter-Glo).

### Modification-Aware Toxicity (for Modified siRNAs)

When evaluating **modified** siRNAs (in the Single-Mod tab), the baseline seed toxicity is computed from the *unmodified* (canonical) seed — since modifications change the chemistry but not the underlying base sequence that drives off-target binding.

However, certain modifications at **seed-region positions (2–7)** can **rescue** toxicity:

```python
# These modifications mitigate seed-driven off-target effects:
SEED_RESCUING_MODS = {"M", "F", "L", "E"}
# M = 2'-OMe, F = 2'-Fluoro, L = LNA, E = 2'-MOE
```

If a rescue modification is present in the seed region **and** the baseline was Toxic or Caution, the label is overridden to **"Mitigated"**.

**Why they work**: 2'-OMe (M) and 2'-Fluoro (F) at seed positions reduce the thermodynamic stability of seed-target duplexes, decreasing off-target silencing while preserving on-target activity (Janas et al., Nucleic Acids Res 2018).

---

## 3. Functional Checks

These are **structural filters** that flag siRNAs likely to have poor specificity or synthesis issues. They are applied to the **sense strand** only (the guide/antisense is the active strand):

| Check | Rule | Why It Matters |
|-------|------|---------------|
| **GC content** | Must be 30–65% | Low GC = weak target binding; high GC = too stable, risk of off-targets |
| **Homopolymer run** | No 5 consecutive identical bases (AAAAA, UUUUU, etc.) | Causes synthesis errors and secondary structure |
| **GC-rich run** | No 6-base run of G/C in any combination | Forms stable G-quadruplex structures, reduces siRNA activity |
| **Internal palindrome** | No 4-base reverse complement within the strand | Causes self-hairpinning, reducing effective concentration |

### GC Content

```python
gc_pct = (sense.count("G") + sense.count("C")) / 21 * 100
# Pass: 30% ≤ gc_pct ≤ 65%
```

**Science**: siRNA activity is optimal in this range. Below 30%, the duplex is too weak for RISC loading. Above 65%, it becomes too stable, increasing off-target seed effects and reducing strand bias.

### Homopolymer Run

```python
re.compile(r"A{5}|U{5}|G{5}|C{5}")
```

A run of 5+ identical nucleotides causes:
- **Synthesis errors**: phosphoramidite coupling efficiency drops for homopolymers
- **Secondary structure**: the strand can self-fold, reducing RISC loading

### GC-Rich Run

```python
# All 64 combinations of G/C of length 6
re.compile("GGGGGG"), re.compile("GGGGGC"), ..., re.compile("CCCCCC")
```

Six consecutive G/C nucleotides (in any combination) can form **G-quadruplex** or **GC-clamp** secondary structures that prevent RISC loading.

### Internal Palindrome

For each 4-base window, checks if the reverse complement appears anywhere downstream:

```python
# Example: AAUGCA → reverse complement of AAUG is CAUU
# If CAUU appears later in the strand → palindrome
```

A palindrome allows the strand to **self-hybridize** (hairpin), reducing the effective concentration available for RISC loading.

---

## 4. Efficacy Label Thresholds

| Label | Score Range | Meaning in Training Data |
|-------|------------|--------------------------|
| **Very High** | ≥80 | Top 16% of training examples |
| **High** | 70–79 | Top ~25% |
| **Moderate** | 55–69 | Above average (~median+) |
| **Low** | <55 | Below median |

**Why these thresholds?** The training data has a wide spread (mean=47.7, std=28.3, range 0–100). A calibrated score of 80 corresponds to the **84th percentile** of the training labels — meaning only 16% of measured siRNAs in the patent data achieved ≥80% silencing. The thresholds are set to make the labels meaningful: "Very High" genuinely means the model expects this candidate to outperform 84% of known sequences.

**Note**: The absolute score values are predictions, not measurements. The **ranking order** (candidate #1 vs #20) is more reliable than the specific number (80.0 vs 79.5). Rank ordering has Spearman ρ ≈ 0.67 against held-out data.

---

## 5. Delta Score (Single-Mod Tab)

```
Delta Score = Modified Variant Score − Parent (Unmodified) Score
```

A positive delta means the chemical modification **improves** efficacy relative to the unmodified parent. A negative delta means the modification **reduces** efficacy.

**What to look for**: The most useful delta scores are positive (+8, +12, etc.) — they identify modification placements that boost silencing. In typical scans, E (2'-MOE) modifications on the antisense strand at positions 1–10 give the largest deltas (+10 to +20 points).

---

## 6. Modification Symbols (All 35)

### Sugar Modifications (14 symbols)

These modify the **ribose sugar** at the 2' position — the most common class for therapeutic siRNA:

| Symbol | Name | Effect on siRNA |
|--------|------|----------------|
| **F** | 2'-Fluoro | Increases nuclease resistance; slight duplex stabilization; well-tolerated at most positions |
| **M** | 2'-OMethyl | Strong nuclease resistance; mild duplex stabilization; gold standard for therapeutic siRNAs |
| **E** | 2'-MOE | Bulky modification — strong nuclease resistance; mild duplex destabilization; good for seed-region toxicity rescue |
| **L** | LNA (Locked Nucleic Acid) | Locks sugar in C3'-endo conformation — strong duplex stabilization; use sparingly |
| **D** | DNA | Deoxyribose — natural DNA nucleotide; weakens RNA duplex; used in some overhang designs |
| **I** | FANA | 2'-Fluoroarabino NA — strong nuclease resistance; unique hybridization properties |
| **N** | 4'-thio | Sulfur replaces 4'-oxygen; very strong nuclease resistance |
| **Y** | ENA | Ethylene-bridge NA; similar to LNA but with 6-membered ring |
| **B** | 2'-O-Benzyl | Large aromatic group — strong nuclease resistance; steric bulk |
| **Z** | 2'-OMe-4'-thio | Double modification — maximum nuclease resistance |
| **X** | 2'-O-allyl | Unsaturated side chain; moderate nuclease resistance |
| **Q** | Abasic | Missing base entirely — disrupts base pairing; used as a spacer |
| **6** | UNA | Unlocked nucleic acid — flexible acyclic sugar; strong duplex destabilization |
| **7** | ANA | Arabino nucleic acid — 2' epimer of RNA |

### Backbone Modifications (4 symbols)

These modify the **phosphodiester linkage**:

| Symbol | Name | Effect |
|--------|------|--------|
| **S** | Phosphorothioate (PS) | Sulfur replaces non-bridging oxygen → nuclease resistance, but can cause toxicity at high density |
| **P** | Boranophosphate | Borane group → nuclease resistance, less toxic than PS |
| **R** | Methylphosphonate | Neutral backbone — no charge; improves cell penetration |
| **H** | Phosphoramidate | N replaces O in backbone; very stable |

### Base Modifications (5 symbols)

These modify the **nucleobase** (the A/U/G/C part):

| Symbol | Name | Effect |
|--------|------|--------|
| **V** | 5-Methyl Cytidine (m5C) | Increases base-pairing stability; reduces immune recognition |
| **W** | Pseudouridine | Natural RNA modification; reduces immune activation |
| **J** | Inosine | Pairs with A, U, C — reduces specificity; used in seed-region for off-target reduction |
| **K** | 2-thio Uridine | Increases duplex stability (stronger U-A pairing) |
| **O** | Dihydrouridine | Flexible base — destabilizes duplex |

### Terminus Modifications (2 symbols)

| Symbol | Name | Effect |
|--------|------|--------|
| **1** | 5'-Phosphate | Required for RISC loading; often added chemically |
| **2** | 3'-Phosphate | Blocks exonuclease degradation |
| **3** | 5'-OMe cap | Blocks 5'-phosphorylation → prevents sense-strand RISC loading (strand bias) |

### Conjugates (2 symbols)

| Symbol | Name | Effect |
|--------|------|--------|
| **4** | Cholesterol | Enables *in vivo* delivery without lipid nanoparticles |
| **5** | PEG | Increases circulation half-life |

---

## 7. Why Rank and Single-Mod Use Different Scores

A common question: *"Why does the Rank tab show score 55 for a candidate, but after running Single-Mod on it the best variant scores 86?"*

**They use different models:**

| Tab | Model | Training Data | Predicts |
|-----|-------|-------------|----------|
| **Rank** | Naked model (model_normal) | 4,060 unmodified siRNAs (4 sources) | **Unmodified (naked) siRNA efficacy** |
| **Single-Mod** | cm-siRNA model (model_a/b/c) | 25,765 modified siRNAs | **Individual modified variant scores** |

The **naked model** (PCC=0.55) was trained only on unmodified siRNAs. It ranks candidates by their predicted efficacy as bare, unmodified sequences — this is the starting point before chemical optimization.

The **cm-siRNA model** (PCC=0.68) was trained on modified siRNAs and knows how modifications affect silencing. When you run Single-Mod on a selected siRNA, it computes the delta: `modified_score - parent_score`. The parent_score uses the **same naked model** the Rank tab uses, so baselines are consistent (both show, e.g., 53.5).

**Rank**: Scores each 21-mer candidate with the naked model and ranks by baseline (unmodified) efficacy.

**Single-Mod**: Runs the **full 1260-variant scan** (all 30 modifications × 21 positions × 2 strands) on ONE selected candidate, showing you which modification/position combination gives the best delta over baseline.

---

## 8. Model Training Summary

### cm-siRNA Model (model_a/b/c)

| Aspect | Detail |
|--------|--------|
| **Algorithm** | LightGBM (gradient-boosted decision trees) |
| **Features** | 152-d (MNC + mod density + GC + conditions) |
| **Training data** | 25,765 modified siRNAs (13 genes, patent data) |
| **Validation** | 5% holdout with early stopping |
| **Trees** | 799 |
| **Performance** | Within-gene PCC = 0.68, MAE = 16.4 |
| **Calibrator** | Isotonic regression (fitted via 5-fold CV) |

### Naked Model (model_normal)

| Aspect | Detail |
|--------|--------|
| **Algorithm** | LightGBM |
| **Features** | 156-d (152 seq + 4 source one-hot) |
| **Training data** | 4,060 unmodified siRNAs from 4 published sources |
| **Sources** | Huesken (2005), Takayuki (2007), Mix (Reynolds/Ui-Tei/Vickers), internal |
| **Performance** | All-source PCC = 0.55, Best source (Taka) = 0.69 |

### Data Sources

- **Huesken et al. 2005**: "Design of a genome-wide siRNA library using an artificial neural network" — *Nature Biotechnology* 23, 995–1001. 2,361 siRNAs with measured % inhibition.
- **Reynolds et al. 2004**: "Rational siRNA design for RNA interference" — *Nature Biotechnology* 22, 326–330.
- **Ui-Tei et al. 2004**: "Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference" — *Nucleic Acids Research* 32, 936–948.
- **Takayuki (Naito) 2007**: "siDirect: highly effective, target-specific siRNA design software" — *Nucleic Acids Research* 35 (Web Server issue).
- **HelixZero patent data**: Internally measured modified siRNA efficacy values for therapeutic targets.
- **Seed-toxicity reference**: Janas et al., Mol Cell 2018 — 4,096-entry siRNA seed → cell viability table.

---

## Glossary

| Term | Definition |
|------|-----------|
| **MNC** | Mononucleotide composition — frequency of each nucleotide symbol in a sequence |
| **LightGBM** | Gradient-boosted decision tree framework from Microsoft |
| **PCC** | Pearson Correlation Coefficient — measures linear correlation between predicted and true values |
| **Spearman ρ** | Rank-order correlation coefficient — measures whether rankings match |
| **MAE** | Mean Absolute Error — average absolute difference between predicted and true values |
| **Seed region** | Antisense positions 2–7 (6-mer) — primary determinant of off-target binding |
| **RISC** | RNA-Induced Silencing Complex — the protein complex that uses the siRNA guide strand to find and cleave target mRNA |
| **2'-OMe** | 2'-O-Methyl modification — replaces the 2'-OH with a methyl group |
| **PS** | Phosphorothioate — replaces a non-bridging oxygen in the backbone with sulfur |
| **Isotonic regression** | Monotonic function fitting — maps predictions to a different scale while preserving order |
