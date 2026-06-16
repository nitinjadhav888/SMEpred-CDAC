"""
sirna_generator.py — 21-mer siRNA candidate generator.

How siRNA works (plain English):
  A siRNA is a 21-nucleotide double-stranded RNA molecule.
  One strand (the 'sense' or 'passenger' strand) matches the mRNA target.
  The other strand (the 'antisense' or 'guide' strand) is the reverse complement
  of the sense strand. The antisense strand is loaded into the RISC complex, which
  then finds and destroys matching mRNA, silencing the gene.

What this module does:
  Slides a window of width 21 across the full mRNA sequence one nucleotide at a time.
  Each window position gives one sense strand. We derive the antisense strand from it.
  Result: (length_of_mRNA − 20) candidate siRNA pairs.
"""

from dataclasses import dataclass
from typing import List


# ─── complement table (RNA) ───────────────────────────────────────────────────
_COMPLEMENT = str.maketrans("AUGC", "UACG")


def _reverse_complement(seq: str) -> str:
    """
    Derive the antisense strand.
    Step 1: complement each base  (A↔U, G↔C).
    Step 2: reverse the result    (5'→3' convention).
    """
    return seq.translate(_COMPLEMENT)[::-1]


# ─── data container ───────────────────────────────────────────────────────────

@dataclass
class SiRNACandidate:
    """One 21-mer siRNA candidate with both strands."""
    position: int        # 0-based start position in the mRNA
    sense: str           # 21-mer matching the mRNA (5' → 3')
    antisense: str       # reverse complement of sense (the guide strand, 5' → 3')

    def to_dict(self) -> dict:
        return {
            "position": self.position,
            "sense":    self.sense,
            "antisense": self.antisense,
        }


# ─── public API ───────────────────────────────────────────────────────────────

SIRNA_LEN = 21


def generate_candidates(mrna: str) -> List[SiRNACandidate]:
    """
    Slide a 21-nt window across the mRNA and return all siRNA candidates.

    Parameters
    ----------
    mrna : str
        Full mRNA/gene sequence (RNA, uppercase, only A/U/G/C).

    Returns
    -------
    List[SiRNACandidate]
        One entry per valid 21-mer position. Length = len(mrna) − 20.
    """
    if len(mrna) < SIRNA_LEN:
        raise ValueError(f"mRNA must be at least {SIRNA_LEN} nt long.")

    candidates = []
    for i in range(len(mrna) - SIRNA_LEN + 1):
        sense = mrna[i : i + SIRNA_LEN]
        antisense = _reverse_complement(sense)
        candidates.append(SiRNACandidate(position=i, sense=sense, antisense=antisense))

    return candidates
