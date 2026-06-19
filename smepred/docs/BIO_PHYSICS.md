# Biological & Biophysical Parameters: Concepts, Literature & Implementation

## Table of Contents

1. [siRNA Mechanism of Action](#1-sirna-mechanism-of-action)
2. [Chemical Modifications Catalogue](#2-chemical-modifications-catalogue)
3. [HelixZero Efficacy Score](#3-helixzero-efficacy-score)
4. [Nuclease Degradation Resistance](#4-nuclease-degradation-resistance)
5. [Immunogenicity](#5-immunogenicity)
6. [RISC Loading & Ago2 Binding](#6-risc-loading--ago2-binding)
7. [Thermodynamic Stability](#7-thermodynamic-stability)
8. [Serum Stability (Half-Life Proxy)](#8-serum-stability-half-life-proxy)
9. [Seed Toxicity](#9-seed-toxicity)
10. [Functional Filters](#10-functional-filters)
11. [Composite Multi-Objective Score](#11-composite-multi-objective-score)
12. [Literature Reference Table](#12-literature-reference-table)
13. [Glossary](#13-glossary)

---

## 1. siRNA Mechanism of Action

### Concept

Small interfering RNA (siRNA) is a double-stranded RNA molecule, typically 21 nucleotides in length, that triggers RNA interference (RNAi) to silence target gene expression in a sequence-specific manner.

### Biological Pathway

```
Exogenous siRNA (21-nt dsRNA)
        │
        ▼
    Loading into RISC (RNA-Induced Silencing Complex)
        │
        ▼
    Strand Selection (guide/passenger discrimination)
        │  └── Argonaute 2 (Ago2) retains the guide strand
        ▼
    Guide-strand-directed mRNA recognition
        │  └── Seed region (positions 2–8 of guide) initiates pairing
        ▼
    Ago2-catalyzed mRNA cleavage
        │  └── Endonucleolytic cleavage between positions 10–11 relative to guide 5'
        ▼
    mRNA degradation → Gene silencing
```

**Key biological actors:**
- **Dicer**: Processes long dsRNA into 21-nt siRNA (exogenous siRNA bypasses Dicer)
- **RISC loading complex (RLC)**: Facilitates strand selection
- **Argonaute 2 (Ago2)**: The "slicer" endonuclease that cleaves target mRNA
- **Guide strand**: The antisense strand that binds target mRNA
- **Passenger strand**: The sense strand that is degraded or ejected

### Relevance to HelixZero-CMS

HelixZero-CMS predicts the **efficacy of chemically modified siRNAs** measured as percentage inhibition of target gene expression (0–100 scale). The system does not model the RNAi pathway explicitly; instead, it learns statistical correlations between modification patterns and measured silencing from 83,535 training examples. The biophysical parameters in this document serve as **post-hoc interpretability and ranking signals** when efficacy alone saturates.

### References

- Fire et al., *Nature* 1998 — RNAi discovery (Nobel Prize 2006)
- Elbashir et al., *Nature* 2001 — 21-nt siRNAs function in mammalian cells
- Rand et al., *Cell* 2005 — Ago2 is the slicer in mammalian RNAi
- Khvorova et al., *Cell* 2003 — Strand selection asymmetry

---

## 2. Chemical Modifications Catalogue

### Concept

Chemical modifications replace atoms or functional groups on the ribose sugar, phosphate backbone, or nucleobase of an siRNA to improve pharmacological properties: nuclease resistance, reduced immunogenicity, enhanced target binding, and enabled delivery.

### Modification Categories (31 total)

| Category | Count | Symbols | Modification Site |
|----------|-------|---------|-------------------|
| **Canonical** | 5 | A, U, G, C, T | Natural nucleotides (T = DNA) |
| **Sugar (2' & 4')** | 13 | F, M, L, E, D, B, N, I, Z, Y, X, Q, 6, 7, 8, 9 | Ribo or deoxyribose ring |
| **Backbone** | 4 | S, P, R, H | Phosphate linkage |
| **Base** | 5 | V, W, J, K, O | Nucleobase |
| **Terminus** | 2 | 1, 2, 3 | 5' or 3' end |
| **Conjugate** | 2 | 4, 5 | Covalent tethers |

### Key Modifications in Approved Drugs

| Drug (Company) | Target | Key Modifications |
|----------------|--------|-------------------|
| Patisiran (Alnylam) | TTR | 2'-OMe, 2'-F, PS termini |
| Givosiran (Alnylam) | ALAS1 | 2'-OMe, 2'-F, PS, 3'-dTdT |
| Inclisiran (Novartis) | PCSK9 | 2'-OMe, 2'-F, PS, GalNAc conjugate |
| Lumasiran (Alnylam) | HAO1 | 2'-OMe, 2'-F, PS termini |

### Relevance to HelixZero-CMS

Each position on the 21-nt sense and antisense strand can carry exactly one modification symbol. The 31 symbols (5 canonical + 31 mod types, excluding the 5 that overlap with canonical symbols actually used) are encoded into a 1,467-dimensional feature vector with per-position binary flags. The model learns which modification at which position on which strand maximally improves silencing efficacy.

### References

- Bramsen & Kjems, *Nucleic Acid Therapeutics* 2012 — Modification category review
- Khvorova & Watts, *Nature Biotechnology* 2017 — siRNA therapeutic chemistry
- Allerson et al., *J Med Chem* 2005 — 2'-OMe/2'-F structure-activity

---

## 3. HelixZero Efficacy Score

### Concept

The primary output of HelixZero-CMS. A 0–100 scaled score predicting the percentage inhibition of target gene expression achieved by a specific chemically modified siRNA (cm-siRNA).

### How It Is Computed

```
Input: (sense_strand, antisense_strand, base_sense, base_antisense)
    │
    ▼
Feature Extraction (extract_positional_features_batch)
    │  1467 features per variant:
    │  ├── 1386 = 33 flags/pos × 21 pos × 2 strands
    │  │    31 mod type flags + 1 is_canonical + 1 is_modified
    │  ├── 70 = global counts per mod type × 2 strands
    │  ├── 16 = aggregate summary stats × 2 strands
    │  └── 1 = log(concentration + 1)
    ▼
Model B (LightGBM gradient-boosted trees)
    │  1,115 trees, trained on 83,535 cm-siRNA rows
    ▼
Raw score → Identity normalization (clip to [0, 100])
    │
    ▼
Efficacy Score (0–100)
```

### Model Performance

| Metric | Value |
|--------|-------|
| PCC (held-out test, n=2,576) | 0.650 |
| R² (position-aware test set) | 0.72 |
| MAE (held-out test) | 16.90 |
| Hetero-V303 external validation | 0.650 |
| Training data | 83,535 rows |

### Known Limitation: Score Saturation

At approximately 5 modifications, the efficacy score frequently reaches 100. This is not a prediction failure but a genuine property of heavily modified siRNAs — clinical candidates routinely score at ceiling. The biophysical parameters (Sections 4–8) provide **differentiating signal** above the saturation point.

### References

- Dar et al., *RNA Biology* 2016 — Original SMEpred workbench
- HelixZero-CMS Paper (Section 3) — Model architecture and evaluation

---

## 4. Nuclease Degradation Resistance

### Concept

Unmodified siRNA is rapidly degraded by endonucleases and exonucleases present in serum and cellular compartments. Chemical modifications can sterically block nuclease access to the phosphodiester backbone and sugar moieties, extending the functional lifetime of the therapeutic.

### Key Findings from Literature

| Finding | Evidence Level | Source |
|---------|---------------|--------|
| PS (phosphorothioate) at 3′ and 5′ termini protects against exonucleases | **Clinical standard** — used in all 4 FDA-approved siRNAs | Alnylam clinical data; Bramsen & Kjems 2012 |
| 2′-OMe (M) provides broad endonuclease resistance proportional to density | **High** — comprehensive SAR study | Allerson et al., *J Med Chem* 2005 |
| 2′-F (F) provides moderate nuclease resistance, less than 2′-OMe | **High** | Layzer et al., *RNA* 2004 |
| LNA (L) at termini enhances exonuclease resistance | **Moderate** | Singh et al., *Chem Comm* 1998; Kurreck et al., *Nucleic Acids Res* 2002 |
| 4′-thio (N) improves nuclease stability | **Moderate** | Hoshika et al., *Nucleic Acids Res* 2004 |
| Full 2′ modification (all positions) maximizes nuclease resistance | **Clinical standard** | Patisiran (Alnylam) — fully modified |

### Implementation in HelixZero-CMS

Rule-based scoring function (0–100):

```
NucleaseScore = 
  + 25  if PS at positions 1, 20, or 21 on either strand
  + 20  if ≥3 PS linkages total
  + 15  if ≥30% of positions are 2'-modified (F/M/L/E)
  + 10  if ≥50% 2'-modified
  + 10  if LNA at 3' terminus
  + 10  if 4'-thio present
  + 10  if no unmodified 5' or 3' terminus
  (max 100, min 0)
```

### Caveats

- Nuclease resistance is **ordinal**, not cardinal — the scale ranks "better vs worse" but does not predict precise half-life in serum
- Different nucleases (endonucleases vs exonucleases) have different modification sensitivities
- Delivery system (LNP, GalNAc) contributes more to in vivo half-life than modification pattern alone

---

## 5. Immunogenicity

### Concept

Unmodified siRNA can activate the innate immune system through Toll-like receptors (TLR7, TLR8), triggering interferon responses and inflammatory cytokines that suppress RNAi and cause toxicity. Chemical modifications can mask immunostimulatory motifs.

### Key Findings from Literature

| Finding | Evidence Level | Source |
|---------|---------------|--------|
| 2′-OMe at **every Uridine** suppresses TLR7/8 activation | **Very High** — definitive study | Judge et al., *Nature Biotechnology* 2006 |
| GU-rich sequences are potent TLR7/8 agonists | **High** — structural basis | Heil et al., *Science* 2004 |
| Pseudouridine (W) substitution reduces immune recognition | **High** (mRNA vaccine data) | Karikó et al., *Immunity* 2005; *Mol Ther* 2008 |
| 5′ triphosphate triggers RIG-I activation | **High** | Hornung et al., *Science* 2006 |
| PS backbone reduces but does not eliminate TLR activation | **Moderate** | Robbins et al., *Oligonucleotides* 2007 |
| Unmodified U-rich 6-mers are TLR7 ligands | **High** — crystal structure | Tanji et al., *Nature* 2015 |
| 2'-F does NOT suppress TLR activation (unlike 2'-OMe) | **Moderate** | Judge et al., 2006 (differential effect) |

### Implementation in HelixZero-CMS

Rule-based scoring function (0–100, higher = less immunogenic):

```
ImmunoScore =
  + 30  if every U in antisense is modified with 2'-OMe (M) or pseudouridine (W)
  + 20  if every U in sense is modified with M or W
  + 15  if any U has pseudouridine (W)
  + 10  if ≥50% of U positions have any 2'-modification
  + 10  if PS backbone present
  + 10  if no terminal 5'-P unmodified (position 1 antisense has mod)
  +  5  if no 4+ nucleotide GU-rich stretch unmodified
  (max 100, min 0)

Penalties (subtracted):
  - 20  if any U is unmodified in antisense strand
  - 15  if any U is unmodified in sense strand
  - 10  per unmodified GU-rich motif (GUUGU, GUG, UGU motifs)
```

### Caveats

- TLR7/8 activation is sequence-specific, not just modification-specific — some sequences are inherently more immunogenic
- Delivery vehicle can mask or exacerbate immune recognition
- These rules are derived from *in vitro* PBMC assays; *in vivo* immunogenicity may differ
- Approved drugs (Patisiran, Inclisiran) have extensive 2'-OMe coverage and minimal immunogenicity

---

## 6. RISC Loading & Ago2 Binding

### Concept

For an siRNA to function, the correct strand (antisense) must be loaded into Argonaute 2 (Ago2) within the RISC complex. Chemical modifications at critical positions can enhance or disrupt this loading process.

### Key Findings from Literature

| Finding | Evidence Level | Source |
|---------|---------------|--------|
| 5′-phosphate (symbol 1) on antisense is essential for Ago2 binding | **Very High** — crystal structure | Ma et al., *Nature* 2005; Frank et al., *Nature* 2010 |
| Thermodynamic asymmetry guides strand selection — less stable 5′ end is preferred for loading | **Very High** — foundational discovery | Khvorova et al., *Cell* 2003; Schwarz et al., *Cell* 2003 |
| Modification at antisense position 1 can disrupt 5′-P recognition | **High** | Bramsen & Kjems, *Nucleic Acid Ther* 2012 |
| Heavy 2′-modification in seed (positions 2–8) can reduce target cleavage efficiency | **Moderate** | Jackson et al., *RNA* 2006 |
| LNA in seed region may enhance or disrupt depending on position | **Moderate** | Braasch et al., *Biochemistry* 2003 |
| PS at guide strand 5′ end is compatible with RISC loading | **Clinical** — used in all approved drugs | Alnylam development data |
| 2′-F is compatible with RISC and is well-tolerated at most positions | **High** | Allerson et al., *J Med Chem* 2005 |

### Implementation in HelixZero-CMS

Rule-based scoring function (0–100):

```
RISCScore =
  + 25  if 5'-P (symbol 1) present on antisense (position 1)
  + 20  if antisense positions 1–8 have ≤3 modifications (minimal seed disruption)
  + 15  if antisense 5' end (pos 1–3) has fewer modifications than sense 5' end
  + 10  if no LNA in antisense positions 2–4 (LNA in early seed can block loading)
  + 10  if sense 3' end has stability-enhancing mods (promotes correct strand bias)
  + 10  if PS at antisense position 1 is absent (PS at the 5'-P site can slightly reduce affinity)
  + 10  if overall modification density on antisense ≤ 80% (leaves RNAi-compatible surface)
  (max 100, min 0)
```

### Design Principle: Thermodynamic Asymmetry

```
                 5' ────────────────────── 3'  Sense (passenger)
                    |||||||||||||||||||||
                 3' ────────────────────── 5'  Antisense (guide)
                    
  AS 5' end (pos 1-4):         LESS stable = better loading
  Sense 5' end (pos 1-4):      MORE stable = better loading
  (Less base pairing at AS 5' → easier strand separation → guide loading)
```

Modifications that stabilize base pairing (LNA, 2'-F) can alter thermodynamic asymmetry and potentially reverse strand preference.

### Caveats

- The strand-biasing effect of modifications is cumulative and complex; the rules above are heuristics
- In fully-modified siRNAs (like clinical candidates), RISC loading efficiency is typically maintained through a balance of modifications
- Some of these rules derive from unmodified siRNA studies and may not fully generalize to heavily modified molecules

---

## 7. Thermodynamic Stability

### Concept

The thermodynamic properties of the siRNA duplex — its melting temperature (Tm), GC content, and internal structures — determine target binding affinity, strand separation energetics, and overall silencing efficiency.

### Key Findings from Literature

| Finding | Evidence Level | Source |
|---------|---------------|--------|
| GC content of 30–52% (optimal ≈ 35–45%) for maximal silencing | **Very High** — validated in 3 independent studies | Reynolds et al., *Nat Biotechnol* 2004; Ui-Tei et al., *Nucleic Acids Res* 2004 |
| Low GC (<30%) → weak target binding; High GC (>55%) → off-target, poor RISC loading | **Very High** | Reynolds eight-criteria |
| No runs of 5+ identical bases (homopolymers) | **High** | Reynolds et al. 2004 |
| No GC-rich runs (≥6 G/C in any combination) | **Moderate** | Reynolds criteria |
| Internal palindromes (4-base self-complementarity) reduce efficacy | **Moderate** | Ui-Tei criteria |
| Stable base pair at 3′ of antisense enhances RISC loading | **Moderate** | Khvorova et al., *Cell* 2003 |
| 2′-OMe and 2′-F increase Tm and duplex stability | **High** | Allerson et al. 2005 |
| LNA significantly increases Tm (+2–8°C per modification) | **High** | Kurreck et al. 2002 |
| PS slightly decreases Tm (−0.5°C per linkage) | **Moderate** | Bramsen & Kjems 2012 |

### Implementation in HelixZero-CMS

Rule-based scoring function (0–100):

```
ThermoScore =
  + 30  if GC% of base sequence is 30–52% (ideal range)
  + 15  if GC% = 35–45% (tight optimal)
  + 10  if terminal GC clamp (pos 20–21 have G or C)
  + 10  if no homopolymer runs (AAAAA, UUUUU, GGGGG, CCCCC)
  + 10  if no GC-only 6-mers
  + 10  if no internal 4-base palindrome
  + 15  if ≥50% modifications are Tm-stabilizing (F/M/L/E/D)
  - 15  if GC% < 30% or > 55%
  - 10  if PS content > 6 linkages (cumulative Tm reduction)
  (max 100, min 0)
```

### Why GC Content Matters

```
    GC% < 30%:        duplex too weak → poor target affinity
    GC% 35–45%:      optimal balance of affinity and specificity
    GC% 45–55%:      moderate — still functional, more off-target risk
    GC% > 55%:       duplex too stable → RISC loading impaired, off-target↑
```

### Caveats

- GC content rules were derived for unmodified siRNA; modifications can compensate for non-ideal GC
- The model (Model B) already accounts for GC indirectly through its feature set and training data
- These rules serve as a **secondary quality signal** when efficacy scores are indistinguishable

---

## 8. Serum Stability (Half-Life Proxy)

### Concept

Serum stability measures how long an siRNA resists degradation in biological fluids. Unlike nuclease resistance (Section 4), which focuses on *mechanisms* of protection, serum stability integrates modification pattern, terminal protection, and overall chemical robustness into a single proxy score.

### Why This Is a Proxy, Not a Direct Measure

True half-life prediction requires:
- Delivery system identity (LNP, GalNAc, lipid conjugate)
- Route of administration (IV, SC, local)
- Species-specific nuclease profiles
- Pharmacokinetic modeling (compartmental)

HelixZero-CMS cannot predict half-life because this information is absent from the training data. Instead, we compute a **serum stability proxy** based on modification-mediated protection features, which correlates directionally with in vitro serum stability.

### Implementation in HelixZero-CMS

Rule-based scoring function (0–100):

```
SerumStabilityScore =
  + 25  if PS linkages protect both 3' and 5' ends of antisense
  + 20  if PS linkages protect both 3' and 5' ends of sense
  + 15  if ≥60% of all positions are 2'-modified (F/M/L/E)
  + 10  if LNA at any terminus (pos 1, 21 on either strand)
  + 10  if any 4'-thio (N) or FANA (I) present
  + 10  if >50% of positions modified overall
  + 10  if at least one terminal conjugate (symbol 4 or 5)
  (max 100, min 0)
```

### Caveats

- **This score measures modification-mediated protection, not pharmacokinetic half-life**
- Clinical half-life is dominated by delivery system, clearance mechanisms, and tissue distribution
- This score should be interpreted as "how well the modification pattern protects the RNA backbone" — a useful but incomplete signal

---

## 9. Seed Toxicity

### Concept

The seed region (antisense positions 2–7) of an siRNA binds to complementary sequences in unintended mRNAs, causing off-target silencing (miRNA-like activity). This can lead to cellular toxicity. Certain seed sequences are more toxic than others, independent of the intended target.

### Data Source

- **Janas et al., *Molecular Cell* 2018**: Systematic measurement of cell viability for 4,097 unique 6-mer seed sequences
- Data is loaded from `data/oligoformer/cell_viability.tsv`

### Implementation: Exact Calculation

```python
# Step 1 — Extract seed
seed = antisense[1:7]     # positions 2–7 (0-indexed: 1–6)
                           # e.g. antisense = "CUUGAAGAAGUCGUGCUGCUU"
                           #       seed      = "UUGAAG"

# Step 2 — Look up cell viability
table = load_tsv("data/oligoformer/cell_viability.tsv")
viability = table.get(seed)   # returns float % or None
```

### Thresholds & Label Assignment

```python
def toxicity_label(viability):
    if viability is None:
        return "Unknown"         # seed not in table
    if viability >= 70.0:
        return "Safe"            # ≥ 70% cell viability
    if viability >= 50.0:
        return "Caution"         # 50–69% cell viability
    return "Toxic"               # < 50% cell viability
```

### UI Display Format

Each label is displayed with the numerical viability percentage (when available):

| Label | Example Display | Meaning |
|-------|----------------|---------|
| **Safe** | `Safe · 79.6%` | Cell viability ≥ 70% — low off-target risk |
| **Caution** | `Caution · 51.1%` | Cell viability 50–69% — moderate risk, use with care |
| **Toxic** | `Toxic · 37.1%` | Cell viability < 50% — high off-target toxicity risk |
| **Unknown** | `Unknown` | Seed not in lookup table — cannot assess |
| **Mitigated** | `Mitigated · 51.1%` | Baseline was Toxic/Caution, but a rescue mod (M/F/L/E) is present in seed positions 2–7 |

The viability percentage always refers to the **canonical seed** (unmodified bases). The percentage does not change with modifications — toxicity risk is assessed at the sequence level, and modifications either mitigate or do not change that risk.

### Modification-Aware Toxicity (Seed Rescue)

Certain modifications in the seed region mitigate off-target toxicity:

| Modification | Rescue Mechanism | Source |
|-------------|------------------|--------|
| **2′-OMe (M)** at pos 2 | Disrupts miRNA-like off-target silencing | Jackson et al., *RNA* 2006 |
| **2′-F (F)** in seed | Reduces seed-pairing stability | Bramsen & Kjems 2012 |
| **LNA (L)** in seed | Blocks seed-mediated repression | Obad et al., *Nat Genet* 2011 |
| **2′-MOE (E)** in seed | Same family as 2′-OMe | Alnylam data |

**Algorithm:**

```python
_SEED_RESCUING_MODS = {"M", "F", "L", "E"}

def toxicity_for_modified(modified_antisense, base_antisense):
    # 1. Get baseline from canonical (unmodified) seed
    base_viab = toxicity_score(base_antisense)
    base_label = toxicity_label(base_viab)

    # 2. Check for rescue modifications in seed positions 2–7
    rescues = []
    for i in range(1, 7):  # 0-indexed: 1..6
        if modified_antisense[i] in _SEED_RESCUING_MODS:
            rescues.append((i + 1, modified_antisense[i]))

    # 3. Override label if rescue present
    if rescues and base_label in {"Toxic", "Caution"}:
        return base_viab, "Mitigated", rescue_note
    if rescues and base_label == "Safe":
        return base_viab, "Safe", rescue_note  # still safe, note the rescue
    return base_viab, base_label, ""           # no change
```

### Relevant Files

- **Toxicity table**: `data/oligoformer/cell_viability.tsv` (4,097 seed→% mappings)
- **Implementation**: `src/filters.py` — `toxicity_score()`, `toxicity_label()`, `toxicity_for_modified()`, `seed_rescue_check()`

### References

- Janas et al., *Mol Cell* 2018 — Seed toxicity table
- Jackson et al., *RNA* 2006 — 2′-OMe seed rescue

---

## 10. Functional Filters

### Concept

Standard siRNA design filters that identify sequence liabilities independent of chemical modifications. These are the classic Reynolds/Ui-Tei rules applied to the canonical (unmodified) sense strand.

### Rules Implemented — Exact Calculations

The filters are applied **in order**. The first failure stops evaluation and returns an error message. All checks use the **sense** strand (unmodified canonical bases).

#### Rule 1: GC Content (30–65%)

```python
def _gc_pct(seq: str) -> float:
    """GC percentage = (G_count + C_count) / total_length * 100"""
    return (seq.count("G") + seq.count("C")) / len(seq) * 100.0
```

**Check:** `30.0 <= gc_pct <= 65.0`

| Condition | Result |
|-----------|--------|
| GC = 29% | `False, "GC 29% out of 30–65%"` |
| GC = 30% | `True` (passes) |
| GC = 52% | `True` (passes) |
| GC = 66% | `False, "GC 66% out of 30–65%"` |

#### Rule 2: No 5-Base Homopolymer Run

```python
import re
_FIVE_RUN = re.compile(r"A{5}|U{5}|G{5}|C{5}")
```

**Check:** No 5 consecutive identical bases anywhere in the sequence.

| Sequence | Result |
|----------|--------|
| `AUGCAUGCAUGCAUGCAUGCA` | Passes |
| `AUUUUUCAUGCAUGCAUGCA` | `False, "5-base homopolymer run"` (UUUUU at pos 3–7) |
| `AGGGGGCACGACUUCUUCA` | `False, "5-base homopolymer run"` (GGGGG at pos 2–6) |

#### Rule 3: No 6-Base GC-Only Run

```python
import itertools
_GC6 = [re.compile("".join(p)) for p in itertools.product("GC", repeat=6)]
# Generates 64 patterns: GGGGGG, GGGGGC, GGGGCG, ..., CCCCCC
```

**Check:** No stretch of 6 consecutive bases where every base is either G or C.

| Sequence | Result |
|----------|--------|
| `GCAGCACGACUUCUUCAAGUU` | Passes (longest GC run is 4: GCAGCA at pos 1–6 → has an A) |
| `GCGCGCGCGAUCUUCAAGUU` | `False, "6-base GC run"` (GCGCGCGCG at pos 1–9 → positions 1–6 are all G/C) |
| `AUGCAUGCGCGCGCGCGCGC` | `False, "6-base GC run"` (GCGCGC at pos 8–13) |

#### Rule 4: No Internal 4-Base Palindrome

```python
def _has_palindrome(seq: str, half: int = 4) -> bool:
    """Internal palindrome: a 4-mer whose reverse-complement appears downstream."""
    trans = str.maketrans("AUCG", "UAGC")
    for i in range(len(seq) - 2 * half + 1):
        rc = seq[i:i+half][::-1].translate(trans)
        if rc in seq[i+half:]:
            return True
    return False
```

**Check:** For each position, take a 4-mer, compute its reverse complement, and search for that exact reverse complement anywhere downstream (after the 4-mer).

| Sequence | Result |
|----------|--------|
| `AUGCAUGCAUGCAUGCAUGCA` | `False, "internal palindrome"` (AUGCAUGC... → AUGU at pos 1–4, its RC CAUA appears at pos 10–13) |
| `GCAGCACGACUUCUUCAAGUU` | Passes |

**What "reverse complement" means for RNA:**
```
A ↔ U
U ↔ A
G ↔ C
C ↔ G

Example: seq[0:4] = "AUGC"
  Reverse:      "CGUA"
  Complement:   "GCAU"
  Check: does "GCAU" appear anywhere in seq[4:]?
```

### Complete Algorithm

```python
def functional_check(siRNA_strand: str) -> tuple[bool, str]:
    """
    Returns (ok, reason). Checks are ORDERED:
      1. GC% in [30, 65]
      2. No 5-base homopolymer run
      3. No 6-base GC-only run
      4. No internal 4-base palindrome
    First failure → immediate return with reason.
    """

    s = siRNA_strand.upper().replace("T", "U")

    gc = _gc_pct(s)
    if not (30.0 <= gc <= 65.0):
        return False, f"GC {gc:.0f}% out of 30–65%"

    if _FIVE_RUN.search(s):
        return False, "5-base homopolymer run"

    for p in _GC6:
        if p.search(s):
            return False, "6-base GC run"

    if _has_palindrome(s):
        return False, "internal palindrome"

    return True, ""
```

### UI Display Format

In the Rank tab results table, each candidate shows a function (✓/✗) column:

| Display | Meaning |
|---------|---------|
| ✓ | All 4 checks pass |
| ✗ GC 29% out of 30–65% | GC content failure with exact value |
| ✗ 5-base homopolymer run | Homopolymer detected (AAAAA/UUUUU/GGGGG/CCCCC) |
| ✗ 6-base GC run | Six consecutive G/C bases |
| ✗ internal palindrome | 4-mer self-complementarity detected |

The user can still select and use the candidate — the flag is advisory, not blocking.

### Source

- Reynolds et al., *Nature Biotechnology* 2004 — Original eight-criteria paper
- Ui-Tei et al., *Nucleic Acids Research* 2004 — Independent validation

### Implementation

Defined in `src/filters.py` (`functional_check()`). Returns `(ok: bool, reason: str)`.

### Relevance in HelixZero-CMS

Applied during siRNA ranking (Rank tab). Candidates failing any filter are flagged with a reason. The user can still use the candidate, but the flag warns of potential experimental failure.

---

## 11. Composite Multi-Objective Score

### Purpose

When the HelixZero efficacy score saturates at 100 (common at ≥5 modifications), the composite score provides a secondary ranking signal that captures biophysical quality across five additional axes.

### Formula

```
Composite = 0.50 × Efficacy_norm + 0.12 × Nuclease + 0.12 × Immunogenicity
          + 0.12 × RISC_Loading + 0.14 × ThermoStability
```

Where:
- `Efficacy_norm` = raw efficacy score / 100 (0–1 range)
- All other scores are in 0–100 range, divided by 100 for normalization
- Serum stability is excluded from the composite to avoid over-weighting degradation protection (it correlates highly with nuclease resistance)

### Weight Rationale

| Parameter | Weight | Rationale |
|-----------|--------|-----------|
| Efficacy | 0.50 | Primary measured phenotype — what the model was trained to predict |
| Thermodynamic Stability | 0.14 | Best-validated rules (Reynolds/Ui-Tei, 2004) — independent replication |
| Nuclease Resistance | 0.12 | Clinically essential — every approved drug has it |
| Immunogenicity | 0.12 | Critical for in vivo use — failure here means drug candidate fails |
| RISC Loading | 0.12 | Mechanistically necessary — without it, no silencing occurs |

### Interpretation

| Composite Range | Interpretation |
|----------------|---------------|
| ≥ 85 | Excellent candidate — top-tier across all biophysical axes |
| 70–84 | Good candidate — minor trade-offs in one or two parameters |
| 50–69 | Adequate — may need optimizations for specific axes |
| < 50 | Poor — significant liabilities in multiple parameters |

### Caveats

- The composite score is **directionally correct but not quantitatively precise** — the weights are heuristic, not trained
- Users should examine individual parameter scores (displayed in the UI) to understand trade-offs
- The composite is best used as a **tiebreaker** when multiple candidates have similar efficacy scores

---

## 12. Literature Reference Table

| # | Parameter | Key Reference | Validation Level | Year |
|---|-----------|---------------|------------------|------|
| 1 | RNAi Mechanism | Fire et al., *Nature* | Nobel Prize | 1998 |
| 2 | 21-nt siRNA | Elbashir et al., *Nature* | Foundational | 2001 |
| 3 | Strand selection asymmetry | Khvorova et al., *Cell* | High (3000+ citations) | 2003 |
| 4 | Strand selection asymmetry | Schwarz et al., *Cell* | High (3000+ citations) | 2003 |
| 5 | Ago2 structure & 5′-P binding | Ma et al., *Nature* | Crystallographic | 2005 |
| 6 | Ago2 slicer mechanism | Rand et al., *Cell* | Biochemical | 2005 |
| 7 | GC content rules | Reynolds et al., *Nat Biotechnol* | High (2000+ citations) | 2004 |
| 8 | Sequence rules (independent) | Ui-Tei et al., *Nucleic Acids Res* | High | 2004 |
| 9 | 2′-OMe suppresses TLR | Judge et al., *Nat Biotechnol* | Very High | 2006 |
| 10 | TLR7 binds GU-rich RNA | Heil et al., *Science* | High | 2004 |
| 11 | 2′-OMe & 2′-F SAR | Allerson et al., *J Med Chem* | High | 2005 |
| 12 | 5′-triphosphate & RIG-I | Hornung et al., *Science* | High | 2006 |
| 13 | Pseudouridine immune evasion | Karikó et al., *Immunity* | High (mRNA vaccines) | 2005 |
| 14 | Chemical mod review | Bramsen & Kjems, *Nucleic Acid Ther* | Comprehensive review | 2012 |
| 15 | siRNA therapeutic chemistry | Khvorova & Watts, *Nat Biotechnol* | Review | 2017 |
| 16 | Seed toxicity | Janas et al., *Mol Cell* | High (4,097 seeds) | 2018 |
| 17 | Seed rescue (2′-OMe) | Jackson et al., *RNA* | High | 2006 |
| 18 | SMEpred workbench | Dar et al., *RNA Biology* | Original workbench | 2016 |
| 19 | LNA properties | Kurreck et al., *Nucleic Acids Res* | Biophysical | 2002 |
| 20 | PS modifications | Eckstein, *Antisense Nucleic Acid* | Review | 2014 |

---

## 13. Glossary

| Term | Definition |
|------|------------|
| **Ago2** | Argonaute 2, the endonuclease that cleaves target mRNA in RISC |
| **Antisense (AS)** | The guide strand that binds target mRNA (5′→3′, positions 1–21) |
| **Beam search** | Heuristic search algorithm that keeps the top-K candidates at each expansion step |
| **cm-siRNA** | Chemically modified small interfering RNA |
| **Composite score** | Weighted combination of efficacy + biophysical parameter scores |
| **Efficacy** | Percentage inhibition of target gene expression (0–100) |
| **GalNAc** | N-acetylgalactosamine conjugate for liver-targeted delivery |
| **Guide strand** | The antisense strand retained by Ago2 |
| **LightGBM** | Gradient-boosted decision tree framework used for Model B |
| **LNP** | Lipid nanoparticle delivery system |
| **Modification density** | Fraction of 21 positions carrying a non-canonical symbol |
| **Passenger strand** | The sense strand that is degraded during RISC loading |
| **PCC** | Pearson correlation coefficient |
| **PS** | Phosphorothioate backbone modification (symbol S) |
| **RISC** | RNA-Induced Silencing Complex |
| **Seed region** | Antisense positions 2–7 (sometimes 2–8), critical for target recognition |
| **Sense (SS)** | The passenger strand (5′→3′, positions 1–21) |
| **TLR** | Toll-like receptor (TLR7, TLR8 sense GU-rich RNA) |
| **Tm** | Melting temperature of the siRNA duplex |
| **2′-F** | 2′-Fluoro modification (symbol F) |
| **2′-OMe** | 2′-O-Methyl modification (symbol M) |
