"""
parser.py — Sequence input handler.

Accepts mRNA / gene input in three forms:
  1. FASTA file path  (.fa, .fasta, .fna)
  2. Raw text file containing a plain sequence
  3. A sequence string passed directly as a Python string

Always returns a clean RNA sequence (uppercase, U instead of T).
"""

import re
from pathlib import Path
from typing import Union


# ─── helpers ──────────────────────────────────────────────────────────────────

def _normalize(seq: str) -> str:
    """Strip whitespace/numbers, uppercase, convert DNA T → RNA U."""
    seq = re.sub(r"[^A-Za-z]", "", seq).upper()
    seq = seq.replace("T", "U")
    valid = set("AUGC")
    bad = set(seq) - valid
    if bad:
        raise ValueError(
            f"Sequence contains unexpected characters: {bad}. "
            "Only A, U, G, C (or T for DNA input) are allowed."
        )
    return seq


def _parse_fasta(text: str) -> str:
    """
    Pull the first sequence out of a FASTA block.
    FASTA format: first line starts with '>', rest is the sequence (may span lines).
    """
    lines = text.strip().splitlines()
    seq_lines = []
    in_seq = False
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            if in_seq:
                break          # stop at second record — we only want the first
            in_seq = True
        elif in_seq:
            seq_lines.append(line)
    if not seq_lines:
        raise ValueError("No sequence data found in FASTA input.")
    return "".join(seq_lines)


# ─── public API ───────────────────────────────────────────────────────────────

def load_sequence(source: Union[str, Path]) -> str:
    """
    Load an mRNA or gene sequence from a file path or inline string.

    Parameters
    ----------
    source : str or Path
        - A file path ending in .fa / .fasta / .fna / .txt
        - A raw multiline string (plain sequence or FASTA format)
        - A single-line sequence string

    Returns
    -------
    str
        Cleaned RNA sequence (uppercase, only A/U/G/C).
    """
    # ── if it looks like a file path, read the file ──
    path = Path(str(source))
    if path.suffix.lower() in (".fa", ".fasta", ".fna", ".txt") and path.exists():
        text = path.read_text(encoding="utf-8")
    else:
        text = str(source)

    # ── detect FASTA vs plain sequence ──
    if text.lstrip().startswith(">"):
        raw = _parse_fasta(text)
    else:
        raw = text

    seq = _normalize(raw)

    if len(seq) < 21:
        raise ValueError(
            f"Sequence is too short ({len(seq)} nt). "
            "Minimum length is 21 nt to generate at least one siRNA."
        )

    return seq
