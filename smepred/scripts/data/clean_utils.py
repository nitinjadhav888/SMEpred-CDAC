"""
clean_utils.py — Shared data-cleaning helpers for the raw siRNA datasets.

These functions implement the data-engineering rules used by the pipeline:

SEQUENCE CLEANING
  - Uppercase everything.
  - Convert DNA Thymine (T) to RNA Uracil (U) so every sequence is in RNA alphabet.
  - Keep only sequences made purely of A/U/G/C after conversion (drop junk like "1", "3").
  - Length filter: keep biologically sensible siRNA strand lengths (19-25 nt). The paper
    used 21-24; we allow 19-25 to retain the slightly shorter/longer real entries.

EFFICACY CLEANING
  - Parse a number out of free text like "69.5 percent target mRNA inhibition".
  - Drop rows whose efficacy cannot be parsed (null / non-numeric).
  - Clip negative inhibition values to 0 (a negative value means no silencing / assay
    noise; biologically efficacy cannot be below 0). Clip anything >100 down to 100.

MODIFICATION NAME -> SYMBOL MAPPING
  - Load the ordered alias_rules from modification_codes.json.
  - Normalize a raw modification name (lowercase, strip apostrophes/quotes) and test the
    alias rules in order; first substring match wins and returns its 1-letter symbol.
  - Names that match no rule return None (caller logs them as "unmapped").
"""

import json
import re
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd

# ─── load alphabet + alias rules once ─────────────────────────────────────────

_MOD_FILE = Path(__file__).parent.parent.parent / "data" / "modification_codes.json"
with _MOD_FILE.open(encoding="utf-8") as _f:
    _MOD_DATA = json.load(_f)

CANONICAL = set(_MOD_DATA["canonical_symbols"])                 # A U G C T
VALID_SYMBOLS = set(_MOD_DATA["canonical_symbols"]) | set(_MOD_DATA["modification_symbols"])
_ALIAS_RULES: List[dict] = _MOD_DATA["alias_rules"]            # ordered list


# ─── sequence cleaning ────────────────────────────────────────────────────────

_RNA_OK = re.compile(r"^[AUGC]+$")


def clean_sequence(seq: str) -> Optional[str]:
    """
    Clean a raw nucleotide sequence.
    Returns the cleaned RNA string, or None if it is invalid/empty.
    """
    if seq is None or (isinstance(seq, float) and np.isnan(seq)):
        return None
    s = str(seq).strip().upper()
    s = re.sub(r"[^A-Z]", "", s)        # remove spaces, digits, punctuation
    s = s.replace("T", "U")             # DNA -> RNA
    if not s or not _RNA_OK.match(s):
        return None
    return s


def valid_length(seq: str, lo: int = 19, hi: int = 25) -> bool:
    """True if the sequence length is within the accepted siRNA strand range."""
    return seq is not None and lo <= len(seq) <= hi


# ─── efficacy cleaning ────────────────────────────────────────────────────────

_NUM_RE = re.compile(r"[-+]?\d*\.?\d+")


def parse_efficacy(value) -> Optional[float]:
    """
    Extract a numeric efficacy from a cell that may be a number or free text
    like '69.5 percent target mRNA inhibition'.
    Returns a float clipped to [0, 100], or None if unparseable.
    """
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    m = _NUM_RE.search(str(value))
    if not m:
        return None
    try:
        v = float(m.group())
    except ValueError:
        return None
    # clip: negatives -> 0 (no silencing), cap at 100
    return float(min(100.0, max(0.0, v)))


# ─── modification name -> symbol ──────────────────────────────────────────────

def _normalize_mod_name(name: str) -> str:
    """Lowercase and strip apostrophes/primes so '2'-O-Methyl' matches '2-o-methyl'."""
    n = str(name).lower()
    n = n.replace("'", "").replace("’", "").replace("`", "")
    n = n.replace("′", "")  # prime symbol
    return n.strip()


def map_modification(name: str) -> Optional[str]:
    """
    Map a raw modification name to a 1-letter symbol using ordered alias rules.

    Returns
    -------
    str  : the matched symbol, OR
    ""   : if the name is a plain canonical base (A/U/G/C/T) or '0' (no modification), OR
    None : if no rule matches (caller should log as unmapped).
    """
    if name is None:
        return ""
    raw = str(name).strip()
    if raw in ("", "0", "nan", "NaN"):
        return ""                       # no modification at this position
    if raw.upper() in CANONICAL:
        return ""                       # canonical base -> keep underlying nucleotide
    norm = _normalize_mod_name(raw)
    for rule in _ALIAS_RULES:
        if rule["match"] in norm:
            return rule["symbol"]
    return None                         # unmapped


def write_unmapped_report(unmapped_counter, out_path: Path, title: str = "Unmapped modifications"):
    """Write a frequency-sorted report of modification names that matched no alias rule."""
    lines = [f"# {title}", f"# Total distinct unmapped names: {len(unmapped_counter)}", ""]
    for name, ct in unmapped_counter.most_common():
        lines.append(f"{ct:>8}  {name}")
    out_path.write_text("\n".join(lines), encoding="utf-8")
