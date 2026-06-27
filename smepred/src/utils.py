"""
utils.py — Shared utility functions for biophysical and functional checks.
"""

def _gc_pct(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.upper().count("G") + seq.upper().count("C")) / len(seq) * 100.0

def _has_palindrome(seq: str, half: int = 4) -> bool:
    """Internal palindrome: a 4-mer whose reverse-complement appears downstream."""
    trans = str.maketrans("AUGC", "UACG")
    for i in range(len(seq) - 2 * half + 1):
        rc = seq[i:i + half][::-1].translate(trans)
        if rc in seq[i + half:]:
            return True
    return False
