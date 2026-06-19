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
    Exposed termini and low 2'-mod density are penalized.
    """
    total = 0.0

    # PS at termini protects against exonucleases
    if sense[0] != "S":
        total += 3
    if sense[20] != "S":
        total += 2
    if antisense[0] != "S":
        total += 3
    if antisense[20] != "S":
        total += 2

    # Too few PS linkages overall
    ps = (sense + antisense).count("S")
    if ps < 3:
        total += 3
    elif ps == 0:
        total += 2  # extra penalty for zero PS

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
    """
    total = 0.0

    # Unmodified U in antisense seed (positions 2–8) → strongest TLR7/8 signal
    for i in range(1, min(8, len(antisense))):
        if base_antisense[i] == "U" and antisense[i] == base_antisense[i]:
            total += 4

    # Unmodified U in antisense tail (positions 9–21)
    for i in range(8, len(antisense)):
        if base_antisense[i] == "U" and antisense[i] == base_antisense[i]:
            total += 1

    # Unmodified U in sense strand
    for i in range(len(sense)):
        if base_sense[i] == "U" and sense[i] == base_sense[i]:
            total += 1.5

    # GU-rich motifs (GUUGU, GUGU, UGU) — only penalize if completely unmodified
    combined = sense + antisense
    base_combined = base_sense + base_antisense
    for motif in ["GUUGU", "GUGU", "UGU"]:
        idx = base_combined.find(motif)
        while idx != -1:
            region = combined[idx:idx + len(motif)]
            if all(c == base_combined[idx + j] for j, c in enumerate(region)):
                total += 3
            idx = base_combined.find(motif, idx + 1)

    # Over-methylation penalty: >16 2'-OMe can trigger alternative immune pathways
    if combined.count("M") > 16:
        total += 4

    return min(total, 28.0)


def risc_penalty(sense: str, antisense: str,
                  base_sense: str, base_antisense: str) -> float:
    """
    Penalty for impaired RISC loading / Ago2 activity (0–31).
    Seed-region modifications and over-modified antisense are most harmful.
    """
    total = 0.0

    # 5'-phosphate on antisense is essential for Ago2 loading
    if antisense[0] != "1":
        total += 5

    # PS at antisense position 1 reduces Ago2 affinity
    if antisense[0] == "S":
        total += 2

    # Modifications in the seed region (positions 2–8) impair target recognition
    seed_mods = 0
    for i in range(1, min(8, len(antisense))):
        if antisense[i] != base_antisense[i]:
            seed_mods += 1
    total += seed_mods * 2  # max +14

    # LNA in early seed (positions 2–4) blocks RISC loading
    for i in range(1, min(4, len(antisense))):
        if antisense[i] == "L":
            total += 5

    # Over-modified antisense (>60%) reduces RISC loading efficiency
    as_mods = sum(1 for i in range(len(antisense)) if antisense[i] != base_antisense[i])
    if as_mods > 12:
        total += 5

    return min(total, 31.0)


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
    Exposed termini and low modification density reduce half-life in serum.
    """
    total = 0.0

    # Unprotected antisense termini → exonuclease degradation
    if antisense[0] != "S":
        total += 4
    if antisense[20] != "S":
        total += 3

    # Unprotected sense termini
    if sense[0] != "S":
        total += 3
    if sense[20] != "S":
        total += 2

    # Low overall modification density
    mod_count = sum(1 for a, b in zip(sense, base_sense) if a != b) + \
                sum(1 for a, b in zip(antisense, base_antisense) if a != b)
    mod_frac = mod_count / 42.0
    if mod_frac < 0.2:
        total += 4
    elif mod_frac < 0.35:
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
