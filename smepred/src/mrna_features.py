"""
mRNA target-site feature extraction for HelixZero-CMS.

Computes features that help the model generalise across genes:
  - Normalised position in transcript (0–1)
  - GC% of target site, upstream 20 nt, downstream 20 nt
  - Target-site MFE via ViennaRNA RNAfold (local window)
"""

import json
import re
from pathlib import Path

import numpy as np

try:
    from Bio import Entrez, SeqIO
    _HAS_BIO = True
except ImportError:
    _HAS_BIO = False

try:
    import RNA
    _HAS_VIENNA = True
except ImportError:
    _HAS_VIENNA = False


_MRNA_CACHE = None  # lazy-loaded dict[gene_symbol] -> {accession, sequence}
_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
_MRNA_JSON = _DATA_DIR / "mrna_sequences.json"

# ─── public API ────────────────────────────────────────────────────────────


def mrna_features_from_genes(
    antisense_list: list[str],
    gene_list: list[str],
    mrna_cache: dict[str, dict] | None = None,
) -> np.ndarray:
    """Compute mRNA target features for a batch of siRNAs.

    Parameters
    ----------
    antisense_list
        List of antisense (guide) strand sequences (21-nt).
    gene_list
        List of target gene symbols (same length).
    mrna_cache
        Optional dict loaded from JSON cache; loaded automatically if None.

    Returns
    -------
    (N, M) array where M is the number of mRNA features.
    """
    if mrna_cache is None:
        mrna_cache = load_mrna_cache()

    rows = []
    for as_seq, gene in zip(antisense_list, gene_list):
        entry = mrna_cache.get(gene)
        if entry is None:
            rows.append(_nan_row())
            continue
        pos = find_target_position(as_seq, entry["sequence"])
        if pos is None:
            rows.append(_nan_row())
            continue
        feat = _compute(as_seq, entry["sequence"], pos)
        rows.append(feat)
    return np.array(rows, dtype=np.float32)


def load_mrna_cache(path: str | Path | None = None) -> dict[str, dict]:
    """Load the mRNA sequence cache from JSON."""
    global _MRNA_CACHE
    if _MRNA_CACHE is not None:
        return _MRNA_CACHE
    p = Path(path) if path else _MRNA_JSON
    if not p.exists():
        _MRNA_CACHE = {}
        return _MRNA_CACHE
    raw = json.loads(p.read_text(encoding="utf-8"))
    _MRNA_CACHE = {g: {"accession": d["accession"], "sequence": d["sequence"].upper().replace("T", "U")}
                   for g, d in raw.items()}
    return _MRNA_CACHE


def fetch_mrna_for_gene(
    gene: str,
    known_accession: str | None = None,
    email: str = "helixzero@research.local",
) -> dict | None:
    """Fetch RefSeq mRNA from NCBI Entrez for a human gene symbol.

    Returns {'accession': ..., 'sequence': ...} or None on failure.
    """
    if not _HAS_BIO:
        return None
    Entrez.email = email
    try:
        if known_accession:
            acc = known_accession
        else:
            query = (f"Homo sapiens[Organism] AND {gene}[Gene Name] "
                     f"AND refseq[Filter] AND mRNA[Filter] AND biomol_mrna[PROP]")
            handle = Entrez.esearch(db="nucleotide", term=query, retmax=3)
            record = Entrez.read(handle)
            handle.close()
            if not record["IdList"]:
                query = f"{gene} human RefSeq mRNA"
                handle = Entrez.esearch(db="nucleotide", term=query, retmax=3)
                record = Entrez.read(handle)
                handle.close()
            if not record["IdList"]:
                return None
            acc = record["IdList"][0]
        handle = Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        seq = str(seq_record.seq).upper().replace("T", "U")
        return {"accession": acc, "sequence": seq}
    except Exception:
        return None


# ─── target site alignment ─────────────────────────────────────────────────


def find_target_position(guide: str, mrna: str, max_mismatch: int = 2) -> int | None:
    """Find the 0-based start position of the siRNA target site in the mRNA.

    The antisense (guide) strand is complementary to the mRNA target site.
    Searches for the reverse complement of the guide in the mRNA.
    Allows up to `max_mismatch` mismatches (default 2).

    Returns None if no suitable match found.
    """
    rc = _revcomp(guide)
    L = len(rc)
    if len(mrna) < L:
        return None
    idx = mrna.find(rc)
    if idx != -1:
        return idx
    best_mm = max_mismatch + 1
    best_pos = None
    for i in range(len(mrna) - L + 1):
        mismatches = sum(1 for a, b in zip(rc, mrna[i:i + L]) if a != b)
        if mismatches < best_mm:
            best_mm = mismatches
            best_pos = i
    return best_pos if best_mm <= max_mismatch else None


# ─── feature computation ───────────────────────────────────────────────────


def _compute(guide: str, mrna: str, pos: int) -> np.ndarray:
    """Return a 6-d feature vector for one siRNA targeting the given mRNA."""
    n = len(mrna)
    L = len(guide)

    # 1 — normalised position
    pos_norm = pos / max(n, 1)

    # 2 — target-site GC%
    target = mrna[pos:pos + L]
    tgc = _gc_frac(target)

    # 3 — upstream GC%
    up = mrna[max(0, pos - 20):pos]
    ugc = _gc_frac(up)

    # 4 — downstream GC%
    dn = mrna[pos + L:min(n, pos + L + 20)]
    dgc = _gc_frac(dn)

    # 5 — target MFE via ViennaRNA
    mfe = _local_mfe(mrna, pos, L)

    return np.array([pos_norm, tgc, ugc, dgc, mfe], dtype=np.float32)


def _gc_frac(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _local_mfe(mrna: str, pos: int, sirna_len: int, flank: int = 50) -> float:
    """Compute MFE (kcal/mol) of the region around the target site.

    Extracts a window of 2*flank + sirna_len centred on the target site
    and folds it with ViennaRNA. Returns 0.0 if ViennaRNA is unavailable
    or the window is too short.
    """
    if not _HAS_VIENNA:
        return 0.0
    start = max(0, pos - flank)
    end = min(len(mrna), pos + sirna_len + flank)
    window = mrna[start:end]
    if len(window) < 21:
        return 0.0
    fc = RNA.fold_compound(window)
    mfe_str, mfe = fc.mfe()
    return mfe


# ─── helpers ───────────────────────────────────────────────────────────────


_COMP = str.maketrans("AUGC", "UACG")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _nan_row() -> np.ndarray:
    return np.full(5, np.nan, dtype=np.float32)


def feature_count() -> int:
    """Number of mRNA features (for dimension tracking)."""
    return 5
