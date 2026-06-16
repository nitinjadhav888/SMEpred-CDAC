# Scoring System — Detailed Explanation

## 1. LightGBM Score (Rank Tab)

The Rank tab shows the **naked (unmodified) siRNA efficacy score** — predicted percent mRNA
knockdown at standard assay conditions (10 nM, 24 h). This is computed by the LightGBM model
`model_normal.pkl` from a **152-dimensional feature vector**.

### The 152 Features — Four Groups

| Group | Dim | What it captures | Always varies? |
|---|---|---|---|
| **Mononucleotide Composition (MNC)** | 140 | Frequency of each of 35 nucleotide symbols in both strands | ✅ Yes — different siRNA sequences have different A/U/G/C composition |
| **Modification Density** | 8 | How many positions differ between modified and base strand | ❌ **No** — for naked siRNAs, modified = base, so all 8 features are **zero** |
| **GC Content** | 2 | GC fraction of the unmodified sense & antisense strands | ✅ Yes — varies with sequence |
| **Assay Condition** | 2 | log₁₀(dose) and normalized time (fixed at inference) | ❌ No — pinned to 10 nM / 24 h at inference |

### 1a. MNC Features (140-d)

Each nucleotide is one of 35 symbols — 4 canonical (A, U, G, C), and 31 chemical
modifications (F=2'-Fluoro, M=2'-OMe, L=LNA, E=2'-MOE, S=PS, etc.).

For **each strand** (sense + antisense), we count how many of the 21 positions are
occupied by each symbol and divide by 21 to get a **frequency vector** (35-d).

```
MNC(sense)       → 35 frequencies
MNC(antisense)   → 35 frequencies
                   ─────
Total per strand → 70-d
```

Then **both** the unmodified (base) and modified strands are encoded:

```
Base MNC      →  70-d
Modified MNC  →  70-d
                   ─────
Total MNC     → 140-d
```

For **naked** siRNAs, the base and modified strands are the same sequence, so
Modified MNC = Base MNC (identical copy). The MNC block captures the **sequence
composition** (A/U/G/C frequencies) — it is the primary source of signal for
the naked model.

For **modified** siRNAs, the Modified MNC captures which chemical symbols (F, M,
L, E, etc.) were introduced and at what density — this gives the cm-siRNA model
its predictive power.

### 1b. Modification Density (8-d)

Per strand, 4 position-aware features comparing modified vs base strand:

| Feature | What it measures |
|---|---|
| **Overall mod fraction** | Fraction of the 21 positions that are modified |
| **Seed-region mod fraction** | Fraction of positions 1–8 that are modified (seed = positions 2–7, expanded window) |
| **3′-tail mod fraction** | Fraction of the last 3 positions that are modified |
| **Normalized mod count** | Number of modified positions / 21 |

Sense → 4-d, Antisense → 4-d = **8-d total**.

> ** ⚠ Important for the Rank tab:** The Rank tab uses the **naked** model, which
> scores unmodified siRNA candidates. For naked siRNAs, modified = base, so all 8
> modification density features are **identically zero** for every candidate. These
> features only carry signal when scoring **modified** variants (Single-Mod / Multi-Mod
> tabs, where the cm-siRNA model is used).

### 1c. GC Content (2-d)

Simple GC fraction of each unmodified strand:

$$\text{GC\%} = \frac{G + C}{A + U + G + C}$$

- Low GC (<30%) → weak silencing
- High GC (>65%) → off-target risks
- Sense GC + Antisense GC = **2-d**

### 1d. Assay Condition (2-d)

Fixed at reference values during inference:

| Feature | Formula | Default |
|---|---|---|
| Dose | log₁₀(conc_nM + 0.01) | log₁₀(10.01) ≈ 1.00 |
| Time | time_h / 24.0 | 24/24 = 1.0 |

These two "nuisance" features allow the model to learn assay-variance during training;
at prediction time they are pinned to standard conditions so the score represents
"predicted activity at 10 nM, 24 h".

### Putting It Together

```
Input: sense strand (21-mer), antisense strand (21-mer)
         │
         ▼
   1. MNC (140-d) —— base MNC + modified MNC
   2. Mod density (8-d) —— 4 per strand  ⚠ 0 for naked candidates
   3. GC content (2-d) —— sense GC + antisense GC
   4. Condition (2-d) —— log₁₀(dose), normalized time  ⚠ fixed at inference
         │
         ▼
   152-d feature vector
         │
         ▼
   LightGBM regressor (84 trees, depth 9)
         │
         ▼
   Raw prediction ──► Isotonic calibrator ──► Score (0–100)

> **Key insight:** For naked siRNAs (Rank tab), features 1b (mod density = 0) and
> the modified half of 1a (modified MNC = base MNC) carry no signal. The model
> effectively runs on ~72 informative dimensions: 70-d base MNC + 2-d GC content.
> The remaining 80 dimensions are structural zeros/duplicates — harmless but
> uninformative for naked prediction.
```

### Efficacy Labels

| Score Range | Label | Meaning |
|---|---|---|
| ≥80 | Very High | Top 6% of training distribution |
| 70–79 | High | 84th–94th percentile |
| 55–69 | Moderate | 75th–84th percentile |
| <55 | Low | Below 75th percentile |

---

## 2. Seed Toxicity Score (Rank Tab Seed Tox Column)

The seed toxicity score is a **lookup-based** prediction of cell viability (%) based on
the 6-mer seed region of the antisense strand.

### How It Works

```
Antisense: 5'- N N N N N N N N N N N N N N N N N N N N N -3'
                 ↑1 2 3 4 5 6 7                         ↑21

Seed = positions 2–7 (1-based), i.e., antisense[1:7] (0-based)
```

**Step 1:** Extract the 6-mer seed from the antisense strand (T → U converted).

**Step 2:** Look up the seed in a precomputed table of 4,096 seeds × cell viability
(Janas et al., oligonucleotide seed-toxicity measurements).

**Step 3:** Classify into a label:

| Viability (%) | Label | Meaning |
|---|---|---|
| ≥70 | Safe | Acceptable off-target risk |
| 50–69 | Caution | Moderate risk — may need seed-rescuing modifications |
| <50 | Toxic | High risk — avoid or apply rescue modifications |
| Not in table | Unknown | No data for this 6-mer |

### Seed Toxicity for Modified siRNAs (Single-Mod / Multi-Mod tabs)

For chemically modified siRNAs, the **baseline** toxicity still comes from the
**unmodified antisense seed** (because the toxicity table is based on canonical
A/U/G/C sequences). Then a secondary check scans the modified antisense seed
region (positions 2–7) for known **off-target-rescuing modifications**:

| Modification | Name | Rescue effect |
|---|---|---|
| **M** | 2′-O-Methyl | Blocks seed-mediated off-target silencing |
| **F** | 2′-Fluoro | Partial rescue |
| **L** | LNA | Strong rescue |
| **E** | 2′-MOE | Strong rescue |

If a rescue modification is present in the seed region and the baseline label was
Toxic or Caution, the label is overridden to **Mitigated** — the underlying
liability exists but the modification strategy addresses it.

### Limitations

- Only 4,096 possible 6-mers are in the lookup table — seeds outside this set
  are marked **Unknown**
- The table is derived from oligonucleotide toxicity assays, not siRNA-specific
  RISC loading data
- Rescue modifications are heuristic (literature-based), not predicted by the model

---

## 3. Functional Checks (Rank Tab Func Column)

Three simple sequence filters applied to the sense strand:

| Check | Rule | Fail if |
|---|---|---|
| **GC content** | 30–65% | Sense GC outside range |
| **Homopolymer run** | No 5 identical bases | AAAAA, UUUUU, GGGGG, CCCCC |
| **GC-only run** | No 6 consecutive G/C | e.g., GGCGCC |
| **Palindrome** | Reverse-complement within strand | 4-mer self-complementarity |

These are standard siRNA design rules (Elbashir, Reynolds, Ui-Tei). siRNAs that
fail these checks are marked with func_ok=False and a reason string, but are still
scored and displayed.
