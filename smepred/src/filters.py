"""
filters.py — Safety, toxicity, and functionality filters for siRNA candidates.

  • toxicity_score  — predicted cell viability (%) for the candidate's 6-mer seed
                      region (antisense positions 2–7). LOWER = more toxic.
                      Source: 4,097-entry siRNA seed → cell viability table
                      (Janas et al., Mol Cell 2018).
  • func_ok         — boolean. True if the candidate passes standard functional
                      criteria (GC in 30–65%, no 5-base homopolymer, no GC₆
                      run, no internal 4-base palindrome). False otherwise.
  • func_reason     — short reason string when func_ok is False.
"""
from __future__ import annotations

import itertools
import re
from functools import lru_cache
from pathlib import Path
from typing import List, Optional

from .utils import _gc_pct, _has_palindrome

import pandas as pd

_TOX_PATH = Path(__file__).parent.parent / "data" / "oligoformer" / "cell_viability.tsv"


# ─── toxicity lookup ─────────────────────────────────────────────────────────

@lru_cache(maxsize=1)
def _tox_table() -> dict:
    """Load the seed → cell-viability table once, cache as a dict."""
    df = pd.read_csv(_TOX_PATH, sep="\t")
    # the file uses RNA letters (A/U/G/C) for the Seed column
    return dict(zip(df["Seed"].str.upper(), df["cell_viability"].astype(float)))


def seed_of_antisense(antisense: str) -> str:
    """Return the 6-mer seed region of an antisense strand (positions 2–7, 1-based)."""
    s = antisense.upper().replace("T", "U")
    return s[1:7]


def toxicity_score(antisense: str) -> Optional[float]:
    """
    Predicted cell viability (%) for the candidate's seed.
    Returns None when the seed isn't in the lookup table.
    Lower = more toxic. >70% is generally considered safe.
    """
    return _tox_table().get(seed_of_antisense(antisense))


def toxicity_label(viability: Optional[float], safe_threshold: float = 70.0) -> str:
    """Plain-English bucket for the UI."""
    if viability is None:
        return "Unknown"
    if viability >= safe_threshold:
        return "Safe"
    if viability >= 50.0:
        return "Caution"
    return "Toxic"


# ─── modification-aware toxicity (for cm-siRNA from Single/Multi-Mod tabs) ────
#
# Seed-rescuing chemical modifications (well-documented in the siRNA-drug literature):
#   • M (2'-OMe)        — Jackson et al., RNA 2006: 2'-OMe at antisense pos 2 suppresses
#                         miRNA-like seed-mediated off-target silencing.
#   • F (2'-Fluoro)     — Bramsen & Kjems 2012: reduces off-target activity.
#   • L (LNA)           — disrupts seed pairing → reduces miRNA-like toxicity.
#   • E (2'-MOE)        — 2'-O-methoxyethyl, same family as 2'-OMe.
#
# Position 2 is the most impactful, but any seed-region (positions 2–7) placement of
# these modifications is generally beneficial. We flag the strongest effect at pos 2.
_SEED_RESCUING_MODS = {"M", "F", "L", "E"}
_MOD_NAMES = {"M": "2'-OMe", "F": "2'-Fluoro", "L": "LNA", "E": "2'-MOE"}


def seed_rescue_check(modified_antisense: str) -> tuple[list[tuple[int, str]], str]:
    """
    Look for seed-rescuing modifications in positions 2–7 of the (modified) antisense.

    Returns
    -------
    (rescue_mods, note)
        rescue_mods : list of (1-based position, symbol) e.g. [(2, "M")]
        note        : human-readable summary used as a tooltip
    """
    s = modified_antisense.upper()
    found = []
    for i in range(1, 7):  # positions 2..7 inclusive (0-based 1..6)
        if i < len(s) and s[i] in _SEED_RESCUING_MODS:
            found.append((i + 1, s[i]))
    if not found:
        return [], ""
    parts = [f"{_MOD_NAMES[sym]} @ pos {pos}" for pos, sym in found]
    note = "Seed off-target rescue: " + ", ".join(parts)
    return found, note


def toxicity_for_modified(modified_antisense: str, base_antisense: str
                          ) -> tuple[Optional[float], str, str]:
    """
    Toxicity for a chemically-modified siRNA.

    Strategy
    --------
     1. The seed-toxicity table is keyed on canonical bases, so we look up
       the seed using `base_antisense` (the unmodified strand) — this gives the BASELINE
       toxicity risk.
    2. We then scan the MODIFIED antisense seed region (pos 2–7) for known
       off-target-rescuing modifications. If present, we override the label to
       "Mitigated" — the underlying liability is the same, but a literature-backed
       rescue strategy is in place.

    Returns
    -------
    (viability_pct, label, note)
        viability_pct : the canonical-seed viability % (None if seed not in table)
        label         : Safe / Caution / Toxic / Mitigated / Unknown
        note          : empty string, or a tooltip explaining a Mitigated flag
    """
    base_viab = toxicity_score(base_antisense)
    base_label = toxicity_label(base_viab)
    rescues, note = seed_rescue_check(modified_antisense)
    if rescues and base_label in {"Toxic", "Caution"}:
        return base_viab, "Mitigated", note
    if rescues and base_label == "Safe":
        # already safe; still surface that a rescue mod is present for completeness
        return base_viab, "Safe", note
    return base_viab, base_label, ""


# ─── functional filter (Reynolds/Ui-Tei rules) ───────────────────────────────

_FIVE_RUN = re.compile(r"A{5}|U{5}|G{5}|C{5}")
_GC6 = [re.compile("".join(p)) for p in itertools.product("GC", repeat=6)]


def functional_check(siRNA_strand: str) -> tuple[bool, str]:
    """
    Returns (ok, reason). Implements standard siRNA design rules:
      1. GC content must be in [30%, 65%]
      2. No run of 5 identical bases (AAAAA / UUUUU / GGGGG / CCCCC)
      3. No 6-base GC-only run (any combination of G/C, length 6)
      4. No internal 4-base palindromic complementarity
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


# ─── batch helpers used by the predictor ─────────────────────────────────────

def annotate_candidates(senses: List[str], antisenses: List[str]) -> List[dict]:
    """Return a list of {toxicity_score, toxicity_label, func_ok, func_reason}
    aligned with the input lists. Used by predictor.rank_sirnas()."""
    out = []
    for sense, anti in zip(senses, antisenses):
        viab = toxicity_score(anti)
        ok, reason = functional_check(sense)
        out.append({
            "toxicity_score": None if viab is None else round(viab, 1),
            "toxicity_label": toxicity_label(viab),
            "func_ok": ok,
            "func_reason": reason,
        })
    return out
