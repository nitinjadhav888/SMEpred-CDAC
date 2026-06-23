# HelixZero Model Reference — Complete Parameter Documentation

## Pipeline Overview

```
┌──────────────┐    ┌──────────────┐    ┌─────────────┐    ┌───────────────┐
│  Sequence    │───→│  Feature     │───→│  LightGBM   │───→│  Biophysical  │
│  Input       │    │  Extraction  │    │  Prediction │    │  Penalty      │
│              │    │              │    │              │    │  Adjustment   │
└──────────────┘    └──────────────┘    └─────────────┘    └───────────────┘
                                                        │
                                                        ▼
                                               ┌───────────────┐
                                               │  Filtering +  │
                                               │  Toxicity     │
                                               │  Check        │
                                               └───────────────┘
```

**Two workflows:**
1. **Rank (unmodified)**: V4 Naked Model (214-d) → RankedSiRNA with toxicity + functional checks.
2. **Modify (modified)**: Multi-step: generate variants → extract positional features (1,467-d) → Model B → biophysical adjustment → RankedCmSiRNA with penalties.

---

## 1. LightGBM Efficacy Model

### Model B v4 (HelixZero)

| Property | Value |
|----------|-------|
| Model file | `models/model_b.pkl` |
| Feature dimensions | **1,467** |
| Training set | 83,535 rows (SMEpred, Dar et al. 2016) |
| Algorithm | LightGBM gradient boosting |
| Objective | `regression` |
| Metric | `rmse` |
| Number of trees | 1,115 |
| Max leaves | 127 |
| Learning rate | 0.03 |
| Feature fraction | 0.8 |
| Bagging fraction | 0.8 |
| Bagging freq | 5 |
| L1 regularization | 0.1 |
| L2 regularization | 0.1 |
| Min data in leaf | 20 |
| PCC (held-out) | **0.822** |
| Spearman ρ | **0.823** |

### Naked Model V4

| Property | Value |
|----------|-------|
| Model file | `models/model_normal.pkl` |
| Feature dimensions | 214 |
| Training set | 83,535 rows |
| PCC (held-out) | **0.55** |
| Purpose | Initial unmodified screening (Rank tab) |

---

## 2. Feature Extraction

### Model B — 1,467 Positional Features

| Feature Group | Calculation | Dimensions | Purpose |
|---------------|-------------|-----------|---------|
| Per-position flags | 33 flags (4 bases + 29 mods) × 42 pos (21 sense + 21 antisense) with parent encoding | 1,386 | Encodes exact chemistry at each position |
| Per-strand global counts | 31 mod symbols × 2 strands | 62 | Global modification burden |
| Summary stat | 8 stats (n_mods, unique_types, M/F/L/E/D/PS counts) × 2 | 16 | High-level chemistry profile |
| Log concentration | log₁₀(10 nM) | 1 | Dose proxy |
| **Total** | | **1,467** | |

### Naked Model — 214 Sequence Features

| Feature Group | Calculation | Dimensions |
|---------------|-------------|-----------|
| Sense one-hot | 4 bases × 21 positions | 84 |
| Sense TNC | 4³ trinucleotides | 64 |
| Antisense TNC | 4³ trinucleotides | 64 |
| GC content | Sense GC%, antisense GC% | 2 |
| **Total** | | **214** |

### Parent-Variant Encoding

Each position uses a 33-flag encoding that captures both the variant and parent chemistry. This allows the model to learn modification delta-effects directly. For each position (42 total = 21 sense + 21 antisense), the flags indicate:
- Whether the position is canonical A/U/G/C
- Which modification symbol(s) are present (for multi-mod variants, multiple flags may be set)
- The parent sequence is encoded at positions 22-42, giving the model simultaneous access to baseline and modified chemistry.

---

## 3. Biophysical Penalties (5 Domains)

Each domain is strictly orthogonal — no biological feature is penalized by more than one module.

### 3.1 Nuclease Penalty (Range: 0–16)

```
def nuclease_penalty(sense, antisense):
    ps_count = count_PS(sense + antisense)
    mod_count = count_2prime_modified(sense + antisense)
    density = mod_count / 42  # fraction over 21+21 positions

    # PS backbone coverage
    if ps_count == 0:     ps_pen = 5
    elif ps_count < 3:    ps_pen = 3
    else:                 ps_pen = 0

    # 2'-mod density (for endonuclease resistance)
    if density < 0.2:     mod_pen = 4
    elif density < 0.4:   mod_pen = 2
    else:                 mod_pen = 0

    return max(ps_pen, mod_pen)  # max, not sum
```

**Reference**: Braasch & Corey 2004, Czauderna et al. 2003.

### 3.2 Immuno Penalty (Range: 0–28)

```
def immuno_penalty(sense, antisense):
    total = 0
    # Unmodified uridine penalties
    for each U in antisense[1:8]:    total += 2.0    # Seed U (Sioud 2004)
    for each U in antisense[8:]:     total += 0.5    # Tail U (Goodchild 2009)
    for each U in sense:              total += 1.0    # Sense U (Judge 2005)

    # GU-rich motifs — non-stacking hierarchical search
    # 1. GUUGU (highest immunostimulatory)
    for each 5-nt window:
        if window == "GUUGU":
            total += 3
            mask all 5 positions  # prevent double-counting

    # 2. GUGU — only unmasked windows
    for each 4-nt window:
        if all positions unmasked and window == "GUGU":
            total += 3
            mask all 4 positions

    # 3. UGU — only unmasked windows
    for each 3-nt window:
        if all positions unmasked and window == "UGU":
            total += 3
            mask all 3 positions

    # Over-methylation advisory
    if count_M(sense + antisense) > 24:
        total += 4    # Clinical ESC threshold
```

**Key note**: Motif detection uses non-stacking — GUUGU, GUGU, and UGU are mutually exclusive within a given window. Without this, a single GUUGU pentamer would be penalized 3 times (GUUGU + GUGU + UGU = +9).

### 3.3 RISC Penalty (Range: −10 to 60)

```
def risc_penalty(sense, antisense):
    total = 0
    # 5'-phosphate
    if antisense[0] != '1':    total += 5    # Missing 5'-P (Frank 2010)
    if antisense[0] == 'S':    total += 2    # PS at pos 1 is suboptimal

    # Seed region (pos 2-8)
    for each pos in antisense[1:8]:
        if modified and mod != '6':            # UNA at pos 7 is exempt (Bramsen 2010)
            total += 2                         # Jackson 2006

    # LNA at seed (pos 2-4) — too rigid
    for each LNA in antisense[1:4]:  total += 5   # Hidayah 2021

    # MOE (2'-MOE) at AS 2-14
    for each MOE in antisense[1:14]: total += 3   # Prakash 2005

    # GNA — position-dependent
    for each GNA in antisense[1:5]:  total += 4   # Disruptive (Schlegel 2022)
    for each GNA in antisense[5:8]:  total -= 2   # Beneficial bonus (ESC+ design)

    # ENA — position-dependent
    for each ENA in antisense[1:8]:  total += 4   # Morihiro 2020
    for each ENA in antisense[8:14]: total += 2   # Tail ENA

    # TNA — position-dependent
    for each TNA in antisense[1:6]:  total += 3
    if antisense[6] == TNA:          total += 0   # Position 7 exempt
    for each TNA in antisense[7:14]: total += 1   # Tail TNA

    # 2'-F deficiency
    pyrimidines = count_total_pyrimidines(sense + antisense)
    f_count = count_2F(sense + antisense)
    f_ratio = f_count / pyrimidines if pyrimidines > 0 else 1.0
    if f_ratio < 0.2: total += 6    # Layzer 2004
    elif f_ratio < 0.4: total += 3  # Layzer 2004

    # Exotic micro-penalties
    for each exotic in [Benzyl, Inosine]: total += 2
    for each rare in [V,I,N,O,P,R,H,K,Z,Q,W,X,7]: total += 1

    return max(-10, min(60, total))  # clamped
```

### 3.4 Thermo Penalty (Range: 0–20)

```
def thermo_penalty(sense, antisense):
    total = 0
    gc = GC_content(sense)  # sense GC%

    if gc < 30 or gc > 55:     total += 8    # Reynolds 2004
    elif gc < 36 or gc > 49:   total += 3    # Borderline

    if has_palindrome(sense, antisense, min_len=8): total += 5
    if has_homopolymer(sense, min_run=4):           total += 5
    if has_gc_run(sense, min_run=6):                total += 3

    return min(20, total)
```

### 3.5 Serum Penalty (Range: 0–17)

```
def serum_penalty(sense, antisense, parent_sense, parent_antisense):
    total = 0
    # Exonuclease protection — termini only
    # AS 5': PS or 5'-PO₄ (symbol '1')
    if antisense[0] not in ('S', '1'):     total += 4
    # AS 3': PS
    if antisense[-1] != 'S':                total += 3    # Elmén 2005
    # SS 5': PS or GalNAc (symbol '4')
    if sense[0] not in ('S', '4'):          total += 3
    # SS 3': PS or GalNAc (symbol '4')
    if sense[-1] not in ('S', '4'):         total += 2

    return min(17, total)
```

**Orthogonality note**: Serum does NOT modify density (that's nuclease). Termini protection is the sole determinant of exonuclease resistance.

---

## 4. Adjusted Score Formula

```
def adjusted_efficacy_score(raw_score, sense, antisense, parent_sense, parent_antisense):
    penalties = {}
    penalties['nuclease'] = nuclease_penalty(sense, antisense)
    penalties['immuno']   = immuno_penalty(sense, antisense)
    penalties['risc']     = risc_penalty(sense, antisense)
    penalties['thermo']   = thermo_penalty(sense, antisense)
    penalties['serum']    = serum_penalty(sense, antisense, parent_sense, parent_antisense)

    total_penalty = sum(penalties.values())
    adjusted = max(0.0, min(100.0, raw_score - 0.70 * total_penalty))

    return adjusted, penalties, total_penalty
```

The 0.70 factor produces:
- Unmodified siRNA: adjusted ≈ 15–25
- Best single-mods: adjusted ≈ 35–60
- Clinical ESC designs: adjusted ≥ 50
- Top beam search: adjusted ≈ 55–70

---

## 5. Toxicity Prediction

### Seed Hexamer Lookup

```
seed_set = set()                     # All 4,096 possible 6-mers
threshold_safe = 75                  # % viability
threshold_caution = 55               # % viability

def toxicity_score(antisense):
    seed = antisense[1:7]            # Positions 2-7
    viability = lookup(seed)         # From cell_viability.tsv
    return viability

def toxicity_label(viability):
    if viability is None:            return "Unknown"
    if viability >= 75:              return "Safe"
    if viability >= 55:              return "Caution"
    return "Toxic"
```

### Seed Rescue Logic (Modified siRNA)

```
def seed_rescue_check(antisense, parent_antisense):
    for each position in 2-7:                                       # Janas 2018
        if antisense[pos] in rescue_mods {M, F, L, E}:
            if parent_antisense[pos] != antisense[pos]:
                return True                                         # Rescue mod present

def toxicity_for_modified(antisense, parent_antisense):
    base_label = toxicity_label(toxicity_score(antisense))
    base_viab = toxicity_score(antisense)

    if seed_rescue_check(antisense, parent_antisense):
        if base_label in ("Toxic", "Caution"):
            return "Mitigated", base_viab, "Seed rescue mod detected"

    return base_label, base_viab, ""
```

---

## 6. Seed Region Analysis

### Position-Specific Effects (Antisense)

| Position | Function | Modification Notes |
|----------|----------|-------------------|
| 1 (5' end) | Ago2 MID domain binding | 5'-PO₄ required. PS allowed but +2 penalty. |
| 2 | Seed pairing (off-target) | LNA/MOE/GNA/ENA all penalized. M/F accepted. |
| 3 | Seed pairing (off-target) | Same as pos 2. |
| 4 | Seed pairing (off-target) | Same as pos 2. |
| 5 | Seed pairing (off-target) | GNA at pos 5 is disruptive (+4). |
| 6 | Seed pairing (off-target) | **GNA at pos 6: −2 bonus** (ESC+ design, Schlegel 2022). |
| 7 | Seed-pivot (central) | **UNA exempt** from seed penalty. GNA at pos 7: −2 bonus. |
| 8 | 3' end of seed | **GNA at pos 8: −2 bonus**. End of seed region. |
| 9–21 | Tail/3' region | MOE, ENA, TNA all have reduced penalties here. |

---

## 7. Functional Filters

Standard siRNA design rules (Reynolds/Ui-Tei):

| Rule | Constraint | Rationale |
|------|-----------|-----------|
| GC content | 30–65% | Melting temperature bounds |
| Homopolymer | No ≥4 consecutive same nt | Secondary structure |
| GC run | No ≥6 consecutive G or C | Stable secondary structure |
| Palindrome | No ≥8 nt self-complementary | Hairpin formation |
| Consecutive runs | No ≥4 repeated dinucleotides | Off-target reduction |

---

## 8. Output Parameter Reference

### RankedSiRNA (Workflow 1 — unmodified)

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| `rank` | int | 1..N | Position in sorted list |
| `position` | int | 0..N-21 | 0-based start in input sequence |
| `sense` | str | 21 nt | Passenger strand (5'→3') |
| `antisense` | str | 21 nt | Guide strand (5'→3') |
| `efficacy_score` | float | 0–100 | Naked model score |
| `efficacy_label` | str | Very High/High/Moderate/Low | Score category |
| `toxicity_score` | float or null | 0–100 | Cell viability % from seed |
| `toxicity_label` | str | Safe/Caution/Toxic/Unknown | Toxicity category |
| `func_ok` | bool | — | Passes functional checks |
| `func_reason` | str | — | Failure reason if not ok |

### RankedCmSiRNA (Workflow 2 — modified)

| Field | Type | Range | Description |
|-------|------|-------|-------------|
| `rank` | int | 1..N | Position in sorted list |
| `sense` | str | 21 nt | Modified sense strand |
| `antisense` | str | 21 nt | Modified antisense strand |
| `mod_symbol` | str | A..7 | The modification symbol |
| `mod_position` | int | 1..21 | Position modified |
| `mod_strand` | str | sense/antisense | Which strand |
| `mod_positions` | str | "2,5,8" | Multi-mod positions |
| `efficacy_score` | float | 0–100 | Adjusted efficacy score |
| `delta_score` | float | −100..100 | Change vs parent adjusted |
| `efficacy_label` | str | - | Score category |
| `toxicity_score` | float or null | 0–100 | Modified seed viability |
| `toxicity_label` | str | Safe/Mitigated/... | Modified toxicity |
| `toxicity_note` | str | "" | Rescue mod note |
| `biophysics` | dict | — | `{'nuclease': F, 'immuno': F, 'risc': F, 'thermo': F, 'serum': F}` |

---

## 9. API Response Parameters

### /single-mod response

| Field | Description |
|-------|-------------|
| `parent_sense` | Original sense strand |
| `parent_antisense` | Original antisense strand |
| `parent_score` | Adjusted score using Model B baseline |
| `naked_baseline` | Adjusted score using Naked Model (for Rank tab comparison) |
| `model_b_baseline` | Same as parent_score (Model B) |
| `total_variants` | 1,302 for full scan, 40 for mini scan |
| `full_scan` | Whether all 1,302 were evaluated |
| `parent_toxicity` | Seed toxicity of parent |
| `results` | Array of RankedCmSiRNA dicts |

### /multi-mod-scan response

| Field | Description |
|-------|-------------|
| `parent_sense` | Original sense strand |
| `parent_antisense` | Original antisense strand |
| `parent_score` | Model B baseline (recalibrated) |
| `naked_baseline` | Naked Model baseline |
| `model_b_baseline` | Same as parent_score |
| `total_variants` | Number of candidates found |
| `results` | Array of rank-ordered dicts with `{rank, sense, antisense, mod_symbol, mod_positions, efficacy_score, raw_efficacy_score, total_penalty, delta_score, efficacy_label, penalties}` |

---

## 10. Modification Symbol Reference

| Symbol | Chemistry | Type | Notes |
|--------|-----------|------|-------|
| A | Adenine | Canonical | Unmodified purine |
| U | Uracil | Canonical | Unmodified pyrimidine |
| G | Guanine | Canonical | Unmodified purine |
| C | Cytosine | Canonical | Unmodified pyrimidine |
| F | 2'-Fluoro | Sugar (standard) | Pyrimidine-specific, nuclease resistant |
| M | 2'-O-Methyl | Sugar (standard) | Most common, used in ESC designs |
| S | Phosphorothioate | Backbone | All termini in ESC designs |
| 1 | 5'-Phosphate | Backbone | Required at AS 5' end |
| L | LNA | Sugar (rigid) | High affinity, restricted at seed |
| E | 2'-MOE | Sugar (bulky) | Moderate RISC impact |
| D | 2'-O-DMAOE | Sugar (cationic) | Membrane permeable |
| 2 | Dihydrouridine | Base | Non-planar |
| 3 | Pseudouridine | Base | Natural RNA modification |
| 5 | 5-Me-C | Base | Reduces immune activation |
| 6 | UNA | Sugar (acyclic) | Seed exemption at pos 7 |
| 8 | GNA | Sugar (small) | Position-dependent bonus/penalty |
| 9 | TNA | Sugar (threose) | Moderate RISC impact |
| Y | ENA | Sugar (ethylene) | Position-dependent penalty |
| 4 | GalNAc | Conjugate | Clinically validated (Givosiran) |
| B | 2'-F-ANA | Sugar (exotic) | Fluoro-arabinose |
| J | 2'-O-Pyrene | Sugar (exotic) | Fluorescent |
| V | 2'-O-N3-adenine | Sugar (exotic) | Click chemistry |
| I | Inosine | Base (exotic) | Universal base |
| N | 2'-O-N3-A | Sugar (exotic) | Click chemistry |
| O | 2'-O-N3-U | Sugar (exotic) | Click chemistry |
| P | 2'-O-N3-C | Sugar (exotic) | Click chemistry |
| R | 2'-O-N3-G | Sugar (exotic) | Click chemistry |
| H | LNA-T | Sugar (exotic) | Locked T |
| K | LNA-C | Sugar (exotic) | Locked C |
| Z | α-LLNA | Sugar (exotic) | Alpha-L stereoisomer |
| Q | 2'-O-allyl | Sugar (exotic) | Allyl modification |
| W | 2'-O-propargyl | Sugar (exotic) | Propargyl modification |
| 7 | Locked-ENA | Sugar (exotic) | LNA-ENA hybrid |

### Symbol Frequency Guidelines

| Usage Level | Symbols | Count |
|-------------|---------|-------|
| Clinical standard | M, S, 1, 4, F | 5 |
| Common research | L, E, D, 3, 5, 6 | 6 |
| Emerging | 8, 9, Y, 2, 7 | 5 |
| Exotic / probe | B, J, V, I, N, O, P, R, H, K, Z, Q, W | 13 |

---

## 11. Literature References

1. Dar SA, et al. *RNA Biol* 2016;13(8):700-712. PMID: 27348347
2. Braasch DA, Corey DR. *Biochemistry* 2002;41(14):4503-4510. PMID: 11926837
3. Czauderna F, et al. *Nucleic Acids Res* 2003;31(11):2705-2716. PMID: 12771180
4. Sioud M, Sørensen DR. *Oligonucleotides* 2004;14(1):1-11. PMID: 15346694
5. Goodchild A, et al. *Oligonucleotides* 2009;19(2):89-98. PMID: 19445602
6. Judge AD, et al. *Nat Biotechnol* 2005;23(4):457-462. PMID: 15778709
7. Frank F, et al. *Nature* 2010;465(7299):818-822. PMID: 20505670
8. Jackson AL, et al. *Nat Biotechnol* 2006;24(9):1151-1157. PMID: 16964229
9. Bramsen JB, et al. *Nucleic Acids Res* 2010;38(9):2861-2878. PMID: 20139420
10. Hidayah NN, et al. *Biochemistry* 2021;60(15):1170-1183. PMID: 33720707
11. Prakash TP, et al. *J Med Chem* 2005;48(21):6696-6705. PMID: 16147155
12. Schlegel MK, et al. *Nucleic Acids Res* 2022;50(16):9056-9070. PMID: 35996904
13. Morihiro K, et al. *J Am Chem Soc* 2020;142(4):1968-1975. PMID: 31867957
14. Layzer JM, et al. *RNA* 2004;10(5):766-771. PMID: 15100431
15. Reynolds A, et al. *Nat Biotechnol* 2004;22(3):326-330. PMID: 14758366
16. Elmén J, et al. *Nucleic Acids Res* 2005;33(1):439-447. PMID: 15653644
17. Janas MM, et al. *Nucleic Acids Res* 2018;46(15):7679-7692. PMID: 29939332
18. Khvorova A, Watts JK. *Nat Biotechnol* 2017;35(3):238-248. PMID: 28244990
19. Nair JK, et al. *J Am Chem Soc* 2014;136(49):16958-16961. PMID: 25434769
20. Bramsen JB, Kjems J. *Front Genet* 2012;3:154. PMID: 22934100
21. Deleavey GF, Damha MJ. *Chem Biol* 2012;19(8):937-954. PMID: 22921066
22. Foster DJ, et al. *Mol Ther* 2018;26(3):708-720. PMID: 29396262
23. Allerson CR, et al. *J Med Chem* 2005;48(4):901-904. PMID: 15415760
24. Dowler T, et al. *Nucleic Acids Res* 2006;34(6):1669-1675. PMID: 16556911
25. Yang X, et al. *Nucleic Acids Res* 2012;40(8):3393-3404. PMID: 22210852
26. Soutschek J, et al. *Nature* 2004;432(7014):173-178. PMID: 15538359
27. Kenski DM, et al. *Mol Ther Nucleic Acids* 2012;1:e5. PMID: 23344720
