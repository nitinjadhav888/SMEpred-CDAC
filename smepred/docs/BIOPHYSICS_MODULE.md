# Biophysics Scoring Module: Design & Implementation

## Status: Legacy Design Document (Superseded)

> **This document describes the *original planned* composite-score architecture.**
> The **current implementation** (`src/biophysics.py`) uses a **penalty-based** approach:
> `adjusted = raw - 0.70 × total_penalty` with 5 domains (nuclease 0-16, immuno 0-28,
> risc 0-50, thermo 0-20, serum 0-17). As of June 2026, RISC includes MOE, GNA, and
> 2'-F deficiency penalties. See [`PENALTIES_REFERENCE.md`](PENALTIES_REFERENCE.md).

---

## 1. Module Interface

```python
# src/biophysics.py

def compute_biophysics(sense: str, antisense: str,
                       base_sense: str, base_antisense: str
                       ) -> dict:
    """
    Compute all biophysical parameter scores for one cm-siRNA variant.

    Parameters
    ----------
    sense : str
        Modified sense strand (21-mer with modification symbols)
    antisense : str
        Modified antisense strand
    base_sense : str
        Unmodified (canonical) sense strand
    base_antisense : str
        Unmodified (canonical) antisense strand

    Returns
    -------
    dict with keys:
        nuclease_resistance : float (0-100)
        immunogenicity     : float (0-100, 100 = least immunogenic)
        risc_loading       : float (0-100)
        thermo_stability   : float (0-100)
        serum_stability    : float (0-100, proxy only)
    """

def compute_composite(biophysics: dict, efficacy: float) -> dict:
    """
    Compute composite score from efficacy + biophysics parameters.

    Returns dict with composite score + individual breakdown.

    Formula:
        composite = 0.50 * efficacy_norm
                  + 0.12 * nuclease_resistance
                  + 0.12 * immunogenicity
                  + 0.12 * risc_loading
                  + 0.14 * thermo_stability

    Where efficacy_norm = efficacy / 100.0, all others / 100.0.
    """
```

### Key Design Decisions

1. **Pure functions, no state** — each call is independent, no model loading needed
2. **Single variant per call** — caller batches if needed (typically called per-candidate during result formatting, not during bulk scoring)
3. **Returns dict, not dataclass** — easy JSON serialization for API response

---

## 2. Scoring Functions — Detailed Logic

### 2.1 `nuclease_resistance_score(sense, antisense, base_sense, base_antisense)`

| Factor | Points | Condition |
|--------|--------|-----------|
| PS at termini | +25 | Either strand has PS (S) at pos 1, 20, or 21 |
| Multiple PS | +20 | Total PS count across both strands ≥ 3 |
| 2'-mod density (moderate) | +15 | ≥30% of 42 positions are F/M/L/E/D |
| 2'-mod density (high) | +10 | ≥50% of 42 positions are F/M/L/E/D |
| LNA at 3' terminus | +10 | LNA at sense or antisense position 21 |
| 4'-thio present | +10 | Any N symbol on either strand |
| No exposed termini | +10 | Position 1 on both strands has a modification |

**Cap:** 100

**Pseudo-code:**
```python
def nuclease_resistance_score(sense, antisense, base_sense, base_antisense):
    score = 0.0
    combined = sense + antisense

    # PS at termini
    term_positions = [sense[0], sense[20], antisense[0], antisense[20]]
    if any(p == 'S' for p in term_positions if p != base_sense[0]):
        score += 25

    # Multiple PS
    ps_count = combined.count('S')
    if ps_count >= 3:
        score += 20
    elif ps_count >= 1:
        score += 10

    # 2'-mod density
    mod_count = sum(1 for c in combined if c in ('F','M','L','E','D'))
    density = mod_count / 42.0
    if density >= 0.5:
        score += 25  # cumulative: 15 + 10
    elif density >= 0.3:
        score += 15

    # LNA at 3'
    if sense[20] == 'L' or antisense[20] == 'L':
        score += 10

    # 4'-thio
    if 'N' in combined:
        score += 10

    # No exposed termini
    if (sense[0] != base_sense[0] and antisense[0] != base_antisense[0]):
        score += 10

    return min(score, 100.0)
```

### 2.2 `immunogenicity_score(sense, antisense, base_sense, base_antisense)`

| Factor | Points | Condition |
|--------|--------|-----------|
| All antisense U modified with M/W | +30 | Every U in antisense has M or W substitution |
| All sense U modified with M/W | +20 | Every U in sense has M or W substitution |
| Any pseudouridine | +15 | W present on either strand |
| Partial U modification | +10 | ≥50% of U positions have M, W, or E |
| PS backbone present | +10 | Any S on either strand |
| Antisense 5' end protected | +10 | Position 1 antisense has a mod (hides 5'-P from RIG-I) |
| **Penalty: unmodified U** | −10 per U | Any unmodified U in antisense seed (pos 2–8) |
| **Penalty: GU-rich motif** | −15 | Unmodified GUUGU, GUGU, or UGU motif |

**Start:** 50 (neutral baseline)

**Pseudo-code:**
```python
def immunogenicity_score(sense, antisense, base_sense, base_antisense):
    scores = 50.0  # neutral baseline

    # U modification analysis
    def u_status(strand, base_strand):
        u_positions = [i for i, b in enumerate(base_strand) if b == 'U']
        modified_u = [i for i in u_positions if strand[i] != 'U']
        u_mod_symbols = [strand[i] for i in modified_u]
        pct_mw = sum(1 for s in u_mod_symbols if s in ('M','W')) / len(u_positions) if u_positions else 1.0
        all_mw = all(s in ('M','W') for s in u_mod_symbols)
        return pct_mw, all_mw, u_positions

    as_pct, as_all, as_u = u_status(antisense, base_antisense)
    ss_pct, ss_all, ss_u = u_status(sense, base_sense)

    if as_all:
        scores += 30
    if ss_all:
        scores += 20

    if 'W' in sense + antisense:
        scores += 15

    if as_pct >= 0.5 or ss_pct >= 0.5:
        scores += 10

    if 'S' in sense + antisense:
        scores += 10

    # Penalties: unmodified U in antisense seed (positions 2-8, 0-indexed 1-7)
    seed_as = antisense[1:8]
    base_seed_as = base_antisense[1:8]
    for i in range(min(len(seed_as), len(base_seed_as))):
        if base_seed_as[i] == 'U' and seed_as[i] == 'U':
            scores -= 10

    # Penalty: GU-rich motifs
    combined_base = (base_sense + base_antisense).upper()
    for motif in ['GUUGU', 'GUGU', 'UGU']:
        if motif in combined_base:
            # Only penalize if the motif region is not modified
            idx = combined_base.find(motif)
            combined_mod = (sense + antisense)[idx:idx+len(motif)]
            if combined_mod == motif:  # no modifications
                scores -= 15

    return max(0.0, min(100.0, scores))
```

### 2.3 `risc_loading_score(sense, antisense, base_sense, base_antisense)`

| Factor | Points | Condition |
|--------|--------|-----------|
| 5'-P on antisense | +25 | Position 1 antisense = '1' (5'-phosphate symbol) |
| Minimal seed disruption | +20 | Antisense positions 2–8 have ≤3 modifications |
| Strand bias (AS 5' less modified) | +15 | AS positions 1–3 have fewer mods than SS positions 1–3 |
| No LNA in early seed | +10 | No LNA at antisense positions 2–4 |
| Sense 3' stabilized | +10 | Modifications at sense positions 19–21 |
| Overall mod density ≤ 80% | +10 | ≤80% of 42 positions modified on AS |
| PS at AS position 1 absent | +10 | Position 1 antisense is NOT PS (PS at 5'-P reduces Ago2 affinity) |

**Pseudo-code:**
```python
def risc_loading_score(sense, antisense, base_sense, base_antisense):
    score = 0.0

    # 5'-P on antisense
    if antisense[0] == '1':
        score += 25

    # Seed disruption (positions 2-8, 0-indexed 1-7)
    seed_mods = sum(1 for i in range(1, 8)
                    if antisense[i] != base_antisense[i])
    score += max(0, 20 - seed_mods * 5)  # 0 mods → +20, 4+ mods → +0

    # Strand bias
    as_5prime_mods = sum(1 for i in range(3) if antisense[i] != base_antisense[i])
    ss_5prime_mods = sum(1 for i in range(3) if sense[i] != base_sense[i])
    if as_5prime_mods <= ss_5prime_mods:
        score += 15

    # No LNA in early seed
    if not any(antisense[i] == 'L' for i in range(1, 4)):
        score += 10

    # Sense 3' stabilized
    if any(sense[i] != base_sense[i] for i in range(18, 21)):
        score += 10

    # Mod density on antisense ≤ 80%
    as_mod_count = sum(1 for i in range(21) if antisense[i] != base_antisense[i])
    if as_mod_count <= 16:  # 16/21 ≈ 76%
        score += 10

    # No PS at antisense position 1
    if antisense[0] != 'S':
        score += 10

    return min(score, 100.0)
```

### 2.4 `thermo_stability_score(sense, antisense, base_sense, base_antisense)`

| Factor | Points | Condition |
|--------|--------|-----------|
| GC% in ideal range | +30 | 30–52% GC in base sense |
| GC% in tight optimal | +15 | 35–45% GC (cumulative with above) |
| Terminal GC clamp | +10 | Positions 20–21 have G or C |
| No homopolymer run | +10 | No AAAA, UUUU, GGGG, CCCC (4+, not 5+) |
| No GC-only 6-mer | +10 | No stretch of 6 only G/C |
| No palindrome | +10 | No 4-base internal palindrome |
| Stabilizing mods dominate | +15 | ≥50% of modifications are Tm-stabilizing (F/M/L/E/D) |
| **Penalty: GC out of range** | −15 | GC < 30% or > 55% |
| **Penalty: heavy PS** | −10 | >6 PS linkages (PS lowers Tm) |

**Pseudo-code:**
```python
def thermo_stability_score(sense, antisense, base_sense, base_antisense):
    score = 0.0
    base = base_sense.upper().replace('T', 'U')

    # GC content
    gc_count = base.count('G') + base.count('C')
    gc_pct = gc_count / 21.0 * 100

    if 30 <= gc_pct <= 52:
        score += 30
    if 35 <= gc_pct <= 45:
        score += 15

    # Terminal GC clamp
    if base[19] in ('G','C') and base[20] in ('G','C'):
        score += 10

    # No homopolymer runs
    has_homopolymer = any(
        seq in base for seq in ['AAAA', 'UUUU', 'GGGG', 'CCCC']
    )
    if not has_homopolymer:
        score += 10

    # No GC-only 6-mer
    import re
    has_gc6 = bool(re.search(r'[GC]{6}', base))
    if not has_gc6:
        score += 10

    # No palindrome
    if not _has_palindrome(base):
        score += 10

    # Stabilizing modifications
    combined_mod = sense + antisense
    mod_positions = sum(1 for i, c in enumerate(combined_mod)
                        if c != (base_sense + base_antisense)[i])
    stabilizing = sum(1 for c in combined_mod if c in ('F','M','L','E','D'))
    if mod_positions > 0 and stabilizing / mod_positions >= 0.5:
        score += 15

    # Penalties
    if gc_pct < 30 or gc_pct > 55:
        score -= 15
    ps_count = combined_mod.count('S')
    if ps_count > 6:
        score -= 10

    return max(0.0, min(100.0, score))


def _has_palindrome(seq: str, half: int = 4) -> bool:
    """Check for internal 4-base palindrome."""
    trans = str.maketrans('AUCG', 'UAGC')
    for i in range(len(seq) - 2 * half + 1):
        rc = seq[i:i+half][::-1].translate(trans)
        if rc in seq[i+half:]:
            return True
    return False
```

### 2.5 `serum_stability_score(sense, antisense, base_sense, base_antisense)`

| Factor | Points | Condition |
|--------|--------|-----------|
| PS protects AS termini | +25 | PS at antisense positions 1, 20, or 21 |
| PS protects SS termini | +20 | PS at sense positions 1, 20, or 21 |
| High 2'-mod density | +15 | ≥60% of positions are 2'-modified |
| LNA at terminus | +10 | LNA at any terminal position |
| Exotic mods (N, I) | +10 | 4'-thio or FANA present |
| Majority modified | +10 | >50% of 42 positions modified |
| Conjugate present | +10 | Symbol 4 or 5 present |

**Cap:** 100

---

## 3. Composite Score

```python
WEIGHTS = {
    "efficacy": 0.50,
    "nuclease": 0.12,
    "immuno":   0.12,
    "risc":     0.12,
    "thermo":   0.14,
}


def compute_composite(bio: dict, efficacy_raw: float) -> dict:
    efficacy_norm = efficacy_raw / 100.0
    composite = (
        WEIGHTS["efficacy"] * efficacy_norm
        + WEIGHTS["nuclease"] * bio["nuclease_resistance"] / 100.0
        + WEIGHTS["immuno"]   * bio["immunogenicity"] / 100.0
        + WEIGHTS["risc"]    * bio["risc_loading"] / 100.0
        + WEIGHTS["thermo"]  * bio["thermo_stability"] / 100.0
    )
    composite_score = round(composite * 100, 2)

    return {
        "composite_score": composite_score,
        "efficacy_score": round(efficacy_raw, 2),
        "parameters": {
            "nuclease_resistance": round(bio["nuclease_resistance"], 1),
            "immunogenicity": round(bio["immunogenicity"], 1),
            "risc_loading": round(bio["risc_loading"], 1),
            "thermo_stability": round(bio["thermo_stability"], 1),
            "serum_stability": round(bio.get("serum_stability", 0), 1),
        }
    }
```

---

## 4. Integration Points

### 4.1 In `src/modification_engine.py` — Score Variants

Inside `multi_mod_scan.score_variants()` and `predict_modified()`:

```python
from .biophysics import compute_biophysics, compute_composite

# After computing efficacy score:
bio = compute_biophysics(v.sense, v.antisense,
                         v.parent_sense, v.parent_antisense)
composite = compute_composite(bio, float(s))

# Store on variant
v.composite = composite["composite_score"]
v.biophysics = composite["parameters"]
```

### 4.2 In `api/main.py` — Endpoint Response

Inside `/multi-mod-from-single` response formatting:

```python
results.append({
    ...
    "efficacy_score": round(v.efficacy_score, 2),
    "composite_score": getattr(v, 'composite', round(v.efficacy_score, 2)),
    "biophysics": getattr(v, 'biophysics', {}),
})
```

### 4.3 In `app.html` — UI Display

Multi-mod results table to include:

1. **Composite Score** column (primary sort, with colored bar)
2. **Expandable row detail** showing all 5 parameter breakdowns
3. **Radar view** (CSS-based) for the top candidate

UI sketch:
```
┌────┬──────────┬──────┬──────┬──────┬──────┬──────┬──────┬─────────────────────┐
│ #  │ Composite│ Eff. │ Nucl │ Imm  │ RISC │ Ther │ Serum│ Mod Pattern          │
├────┼──────────┼──────┼──────┼──────┼──────┼──────┼──────┼─────────────────────┤
│  1 │ ██ 88.5  │ 100  │  75  │  82  │  90  │  71  │  65  │ 8+M+F+S (5 mods)    │
│  2 │ ██ 84.2  │  98  │  80  │  78  │  75  │  80  │  70  │ M+F+L (3 mods)      │
└────┴──────────┴──────┴──────┴──────┴──────┴──────┴──────┴─────────────────────┘
```

---

## 5. Beam Search Change: Dynamic Early-Stop

Replace `max_mods` with dynamic early-stop in `multi_mod_scan()`:

```python
EARLY_STOP_THRESHOLD = 1.0  # points
NO_IMPROVEMENT_LIMIT = 2     # consecutive rounds

prev_best = -float('inf')
stale_rounds = 0

for n_mods in itertools.count(2):
    # ... beam expansion ...

    # Early stop check
    if len(scored) > 0:
        current_best = max(v.composite for v in scored if hasattr(v, 'composite'))
        if current_best - prev_best < EARLY_STOP_THRESHOLD:
            stale_rounds += 1
            if stale_rounds >= NO_IMPROVEMENT_LIMIT:
                break
        else:
            stale_rounds = 0
        prev_best = current_best

    if n_mods >= 21:  # absolute max
        break
```

---

## 6. Implementation Order

1. **Create** `src/biophysics.py` with all 5 scoring functions + composite
2. **Integrate** into `modification_engine.py` — add composite to beam search variants
3. **Integrate** into `predictor.py` — add composite to `predict_modified` results
4. **Update** `api/main.py` — surface composite + biophysics in API responses
5. **Update** `app.html` — display parameter breakdown in multi-mod results
6. **Update** `cli/run.py` — show biophysics breakdown in CLI output
7. **Test** — unit tests for each scoring function + integration tests
