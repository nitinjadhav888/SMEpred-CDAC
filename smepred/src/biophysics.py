"""
biophysics.py — Biophysical penalties that adjust the main efficacy score.

Each function returns a penalty (0 to max_value) reflecting how a candidate
violates a biophysical design principle. The total penalty is subtracted
from the raw LightGBM efficacy score to produce a realistic final score.

This captures trade-offs: a modification that improves nuclease resistance
may also impair RISC loading. The adjusted score naturally ranks well-balanced
multi-mod designs above over-modified or under-protected ones.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple

__all__ = [
    "adjusted_efficacy_score",
    "nuclease_penalty",
    "immuno_penalty",
    "risc_penalty",
    "thermo_penalty",
    "serum_penalty",
]

# ─── helpers ──────────────────────────────────────────────────────────────

_MOD_2PRIME = frozenset("FMLEBD")


def _gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.upper().count("G") + seq.upper().count("C")) / len(seq) * 100.0


def _has_palindrome(seq: str, half: int = 4) -> bool:
    trans = str.maketrans("AUGC", "UACG")
    for i in range(len(seq) - 2 * half + 1):
        rc = seq[i:i + half][::-1].translate(trans)
        if rc in seq[i + half:]:
            return True
    return False


def _has_homopolymer(seq: str, n: int = 5) -> bool:
    for base in ("A", "U", "G", "C"):
        if base * n in seq.upper():
            return True
    return False


def _find_motifs(seq: str, motifs: List[str]) -> List[int]:
    """Return start positions of all unmodified motif occurrences."""
    positions = []
    for motif in motifs:
        idx = seq.find(motif)
        while idx != -1:
            positions.append(idx)
            idx = seq.find(motif, idx + 1)
    return positions


# ═════════════════════════════════════════════════════════════════════════
# Penalty functions  (0 = best, higher = worse)
# ═════════════════════════════════════════════════════════════════════════


def nuclease_penalty(sense: str, antisense: str,
                      base_sense: str, base_antisense: str) -> float:
    """
    Penalty for inadequate nuclease resistance (0–16).

    Confined to endonuclease defence (2'-mod density) and backbone PS
    coverage. Exonuclease (termini) protection is handled exclusively by
    serum_penalty to avoid cross-module double-counting.
    """
    total = 0.0

    # PS backbone coverage — too few linkages → endonuclease vulnerability
    ps = (sense + antisense).count("S")
    if ps == 0:
        total += 5
    elif ps < 3:
        total += 3

    # Low 2'-mod density → vulnerable to endonucleases
    combined = sense + antisense
    mod_2prime = sum(1 for c in combined if c in _MOD_2PRIME)
    density = mod_2prime / 42.0
    if density < 0.2:
        total += 4
    elif density < 0.4:
        total += 2

    return min(total, 16.0)


def immuno_penalty(sense: str, antisense: str,
                    base_sense: str, base_antisense: str) -> float:
    """
    Penalty for immunostimulatory features (0–28).
    Unmodified Uridine in seed = strongest signal.
    GU-rich motifs also trigger TLR-mediated response.

    Penalties are calibrated to avoid stacking: a single TLR-binding epitope
    (e.g. GUUGU) is penalized once, not once for each contained sub-motif.
    """
    total = 0.0

    # ── Unmodified U in antisense seed (positions 2–8) → TLR7/8 signal ──
    # Calibrated per Sioud & Sørensen (2004): single U drives signal,
    # but +4 per U was too aggressive — lowered to +2.0 (C-DAC review 2026).
    for i in range(1, min(8, len(antisense))):
        if base_antisense[i] == "U" and antisense[i] == base_antisense[i]:
            total += 2.0

    # ── Unmodified U in antisense tail (positions 9–21) ──
    # Endosomal exposure of the tail is secondary to seed recognition.
    for i in range(8, len(antisense)):
        if base_antisense[i] == "U" and antisense[i] == base_antisense[i]:
            total += 0.5

    # ── Unmodified U in sense strand ──
    # Passenger strand is rapidly degraded; lower immune weight justified.
    for i in range(len(sense)):
        if base_sense[i] == "U" and sense[i] == base_sense[i]:
            total += 1.0

    # ── GU-rich motifs — non-stacking hierarchical search ──
    # GUUGU is the true high-affinity ligand for human TLR8 (Hornung et al. 2005).
    # GUGU and UGU are weaker fragments bound by the same pocket. We check
    # longest-first and mask covered positions with a sentinel so the same
    # bases cannot trigger multiple motif penalties.
    base_combined = list(base_sense + base_antisense)
    combined = list(sense + antisense)
    covered = [False] * len(combined)       # mask: True = already penalized

    for motif in ["GUUGU", "GUGU", "UGU"]:
        mlen = len(motif)
        # Build a search string where covered positions are replaced with '.'
        # (never matches A/U/G/C so sub-motifs inside a prior hit are invisible)
        search_str = "".join(
            base_combined[i] if not covered[i] else "."
            for i in range(len(base_combined))
        )
        idx = 0
        while True:
            idx = search_str.find(motif, idx)
            if idx == -1:
                break
            # Check the window is still unmodified
            region_mod = combined[idx:idx + mlen]
            region_base = base_combined[idx:idx + mlen]
            if all(r == region_base[j] for j, r in enumerate(region_mod)):
                total += 3.0
                for j in range(idx, idx + mlen):
                    covered[j] = True
            idx += 1

    # ── Over-methylation advisory ──
    # Extreme 2'-OMe saturation (>24 of 42 nt) can engage alternative pathways
    # (Robbins et al. 2007). Clinical ESC designs operate at 25–27 safely.
    if (sense + antisense).count("M") > 24:
        total += 4

    return min(total, 28.0)


def risc_penalty(sense: str, antisense: str,
                  base_sense: str, base_antisense: str) -> float:
    """
    Penalty for impaired RISC loading / Ago2 activity (min: −10, max: 60).
    Seed-region modifications and over-modified antisense are most harmful.
    GNA at pos 6–8 and UNA at pos 7 can reduce penalty (beneficial).
    """
    total = 0.0

    # 5'-phosphate on antisense is essential for Ago2 loading
    if antisense[0] != "1":
        total += 5

    # PS at antisense position 1 reduces Ago2 affinity
    if antisense[0] == "S":
        total += 2

    # Modifications in the seed region (positions 2–8) impair target recognition
    # UNA at position 7 is exempt — it improves off-target profile (Bramsen 2010)
    seed_mods = 0
    for i in range(1, min(8, len(antisense))):
        if antisense[i] != base_antisense[i]:
            if not (antisense[i] == "6" and i == 6):  # UNA@7 exempt
                seed_mods += 1
    total += seed_mods * 2  # max +14

    # LNA in early seed (positions 2–4) blocks RISC loading
    for i in range(1, min(4, len(antisense))):
        if antisense[i] == "L":
            total += 5

    # MOE in guide strand positions 2–14 impairs Ago2 loading (bulky 2'-MOE)
    for i in range(1, min(14, len(antisense))):
        if antisense[i] == "E":
            total += 3

    # GNA — position-dependent (Schlegel 2022 ESC+: pos 6-8 beneficial)
    # GNA in positions 2–5 is disruptive (seed core)
    for i in range(1, min(5, len(antisense))):
        if antisense[i] == "8":
            total += 4
    # GNA in positions 6–8 is beneficial (ESC+, 6-8x therapeutic window)
    for i in range(5, min(8, len(antisense))):
        if antisense[i] == "8":
            total -= 2  # bonus

    # ENA (Y) — bicyclic nucleotide, stiffer than LNA
    # ENA in seed positions 2–8 → steric clash (analogous to LNA)
    for i in range(1, min(8, len(antisense))):
        if antisense[i] == "Y":
            total += 4
    # ENA in central body 9–14 → over-stabilizes guide-target duplex
    for i in range(8, min(14, len(antisense))):
        if antisense[i] == "Y":
            total += 2

    # TNA (9) — 4-carbon backbone, position-dependent (Mori 2025)
    # TNA in seed positions 2–6: backbone shift disrupts Ago2 register
    for i in range(1, min(6, len(antisense))):
        if antisense[i] == "9":
            total += 3
    # TNA at position 7 → 0 (beneficial, ESC+-like; exempt from penalty)
    # TNA in positions 8–14: mild body disruption
    for i in range(7, min(14, len(antisense))):
        if antisense[i] == "9":
            total += 1

    # Missing 2'-F on pyrimidines reduces Ago2 compatibility
    f_on_pyrs = sum(1 for i in range(len(antisense))
                    if antisense[i] == "F" and base_antisense[i] in "UC")
    pyrimidines = sum(1 for b in base_antisense if b in "UC")
    if pyrimidines > 0 and f_on_pyrs / pyrimidines < 0.2:
        total += 6
    elif pyrimidines > 0 and f_on_pyrs / pyrimidines < 0.4:
        total += 3

    # Exotic mod micro-penalties — differentiate rare chemistries in guide strand
    # Bulky aromatic (Benzyl), non-canonical bases (Inosine), and other sparse-training
    # modifications get +1 each to break ties and reflect greater biological uncertainty.
    exotic_mods = frozenset("BJVINOPRHKZQWX7")
    exotic_count = sum(1 for c in antisense if c in exotic_mods)
    if exotic_count > 0:
        total += exotic_count * 1.0
    # Extra penalty for bulkiest aryl modifications
    if "B" in antisense:
        total += 1
    if "J" in antisense:
        total += 1

    return min(max(total, -10.0), 60.0)


def thermo_penalty(sense: str, antisense: str,
                    base_sense: str, base_antisense: str) -> float:
    """
    Penalty for thermodynamically unfavorable sequences (0–20).
    Extreme GC, homopolymer runs, and palindromes destabilize.
    """
    total = 0.0
    base = base_sense.upper()

    gc = _gc_pct(base)
    if gc < 30.0 or gc > 55.0:
        total += 8
    elif gc < 35.0 or gc > 50.0:
        total += 3

    if _has_palindrome(base):
        total += 5

    if _has_homopolymer(base):
        total += 5

    if re.search(r"[GC]{6}", base):
        total += 3

    return min(total, 20.0)


def serum_penalty(sense: str, antisense: str,
                   base_sense: str, base_antisense: str) -> float:
    """
    Penalty for poor serum stability (0–17).

    Exclusively checks exonuclease vulnerability at termini.
    Modification density and backbone PS are confined to nuclease_penalty
    to maintain strict orthogonal separation of concerns.
    """
    total = 0.0

    # Unprotected antisense termini → exonuclease degradation
    # Antisense 5'-PO4 ("1") bound by Ago2 MID domain provides equivalent protection;
    # sense 3'-conjugate ("4") provides steric nuclease shield.
    if antisense[0] not in ("S", "1"):
        total += 4
    if antisense[20] not in ("S", "1"):
        total += 3

    # Unprotected sense termini
    if sense[0] not in ("S", "4"):
        total += 3
    if sense[20] not in ("S", "4"):
        total += 2

    return min(total, 17.0)


# ═════════════════════════════════════════════════════════════════════════
# Public API
# ═════════════════════════════════════════════════════════════════════════

# Scale factor — tunes how aggressively penalties reduce the raw score.
# 0.7 means a total penalty of 50 reduces raw score by 35 points.
_ADJUSTMENT_FACTOR = 0.70


def adjusted_efficacy_score(
    raw_score: float,
    sense: str,
    antisense: str,
    base_sense: str,
    base_antisense: str,
) -> Tuple[float, Dict[str, float], float]:
    """
    Apply all five biophysical penalties to a raw efficacy score.

    Parameters
    ----------
    raw_score     : raw LightGBM prediction (0–100)
    sense         : modified sense strand
    antisense     : modified antisense strand
    base_sense    : unmodified parent sense strand
    base_antisense: unmodified parent antisense strand

    Returns
    -------
    (adjusted_score, penalty_breakdown, total_penalty)
      adjusted_score  : raw_score minus scaled penalties, clipped to [0, 100]
      penalty_breakdown: dict with per-parameter penalty values
      total_penalty   : sum of all penalties before scaling
    """
    penalties = {
        "nuclease": nuclease_penalty(sense, antisense, base_sense, base_antisense),
        "immuno": immuno_penalty(sense, antisense, base_sense, base_antisense),
        "risc": risc_penalty(sense, antisense, base_sense, base_antisense),
        "thermo": thermo_penalty(sense, antisense, base_sense, base_antisense),
        "serum": serum_penalty(sense, antisense, base_sense, base_antisense),
    }
    total_penalty = sum(penalties.values())
    adjusted = raw_score - _ADJUSTMENT_FACTOR * total_penalty
    adjusted = max(0.0, min(100.0, adjusted))
    return adjusted, penalties, total_penalty
