"""
features.py — Feature extraction for LightGBM models.

Three feature types are extracted from siRNA sequences (Paper Table 1):

1. MNC — Mononucleotide Composition
   The fraction of each nucleotide type in a sequence.
   Example: for a 21-nt sequence with 5 A's → frequency(A) = 5/21 = 0.238
   We have 35 nucleotide types (5 canonical + 30 chemical modification symbols).
   Sense strand gives 35 values, antisense strand gives 35 values → 70-d vector.
   This is the BEST single feature (PCC 0.80, Model-A).

2. DNC — Dinucleotide Composition
   The fraction of each pair of consecutive nucleotides.
   35 × 35 = 1225 possible pairs per strand → 2450-d vector total.
   Less predictive than MNC alone for this task.

3. BIN — Binary Pattern (Position Encoding)
   Tells the model *where* in the sequence each nucleotide type appears.
   For each of the 35 nucleotide types, we write a binary string of length equal
   to the siRNA strand length (24 nt max as used in paper): 1 at positions where
   that type is present, 0 elsewhere. Each position uses a 6-bit code.
   Used in Model-B (antisense 5'-end 13 nt) and Model-C (antisense 3'-end last 8 nt).

Why composition works so well:
   The nucleotide composition captures GC content, AU-richness, and the presence/absence
   of chemical modifications — all of which directly affect how well RISC loads the guide
   strand and cleaves the target mRNA.
"""

import json
from pathlib import Path
from typing import List, Optional

import numpy as np


# ─── load modification codes ──────────────────────────────────────────────────

_MOD_FILE = Path(__file__).parent.parent / "data" / "modification_codes.json"
with _MOD_FILE.open() as _f:
    _MOD_DATA = json.load(_f)

ALL_SYMBOLS: List[str] = [m["symbol"] for m in _MOD_DATA["modifications"]]  # 35 total
SYMBOL_INDEX: dict = {s: i for i, s in enumerate(ALL_SYMBOLS)}               # symbol → index
N_SYMBOLS: int = len(ALL_SYMBOLS)   # 35
SIRNA_LEN: int = 21                 # standard siRNA length (21 nt duplex)


# ─── MNC ──────────────────────────────────────────────────────────────────────

def mnc_vector(seq: str) -> np.ndarray:
    """
    Mononucleotide Composition for one strand.

    Returns a 35-dimensional vector of nucleotide frequencies (0.0 – 1.0).
    Each element = count(symbol_i) / len(seq).
    """
    counts = np.zeros(N_SYMBOLS, dtype=np.float32)
    n = len(seq)
    if n == 0:
        return counts
    for ch in seq:
        idx = SYMBOL_INDEX.get(ch)
        if idx is not None:
            counts[idx] += 1
    return counts / n


def mnc_full(sense: str, antisense: str) -> np.ndarray:
    """
    70-dimensional MNC vector = MNC(sense) + MNC(antisense) concatenated.
    This is the feature vector for Model-A.
    """
    return np.concatenate([mnc_vector(sense), mnc_vector(antisense)])


# ─── DNC ──────────────────────────────────────────────────────────────────────

def dnc_vector(seq: str) -> np.ndarray:
    """
    Dinucleotide Composition for one strand.

    Returns a (35×35) = 1225-dimensional vector of dinucleotide frequencies.
    Each element = count(symbol_i + symbol_j consecutive) / (len(seq)-1).
    """
    counts = np.zeros(N_SYMBOLS * N_SYMBOLS, dtype=np.float32)
    n = len(seq)
    if n < 2:
        return counts
    for k in range(n - 1):
        i = SYMBOL_INDEX.get(seq[k])
        j = SYMBOL_INDEX.get(seq[k + 1])
        if i is not None and j is not None:
            counts[i * N_SYMBOLS + j] += 1
    denom = n - 1
    return counts / denom


def dnc_full(sense: str, antisense: str) -> np.ndarray:
    """
    2450-dimensional DNC vector = DNC(sense) + DNC(antisense) concatenated.
    """
    return np.concatenate([dnc_vector(sense), dnc_vector(antisense)])


# ─── BIN ──────────────────────────────────────────────────────────────────────

def bin_vector(seq: str) -> np.ndarray:
    """
    Binary Pattern for one strand.

    For each nucleotide type (35 types) × each position in the sequence,
    we write 1 if that type is at that position, else 0.
    This gives a 35 × len(seq) flat binary vector.

    The paper uses 6-bit encoding per nucleotide position (35 types need 6 bits).
    We flatten this as: for each symbol type, a 1/0 at each position.
    Vector length = 35 × len(seq).
    """
    n = len(seq)
    vec = np.zeros(N_SYMBOLS * n, dtype=np.float32)
    for pos, ch in enumerate(seq):
        idx = SYMBOL_INDEX.get(ch)
        if idx is not None:
            vec[idx * n + pos] = 1.0
    return vec


def bin_antisense_seed(antisense: str, seed_len: int = 13) -> np.ndarray:
    """
    BIN of antisense seed region: first `seed_len` nucleotides from the 5'-end.
    Default seed_len=13 → used in Model-B (PCC 0.86 on independent set).
    The seed region is where RISC first contacts the mRNA target.
    """
    return bin_vector(antisense[:seed_len])


def bin_antisense_tail(antisense: str, tail_len: int = 8) -> np.ndarray:
    """
    BIN of antisense tail: last `tail_len` nucleotides from the 3'-end.
    Default tail_len=8 → used in Model-C (PCC 0.78 on independent set).
    """
    return bin_vector(antisense[-tail_len:])


# ─── composite feature builders for each model ────────────────────────────────

def _resolve_base(sense, antisense, base_sense, base_antisense):
    """If base strands aren't supplied, the modified strands ARE the base (unmodified)."""
    if base_sense is None:
        base_sense = sense
    if base_antisense is None:
        base_antisense = antisense
    return base_sense, base_antisense


def features_model_a(sense, antisense, base_sense=None, base_antisense=None) -> np.ndarray:
    """
    Model-A feature vector: MNC of the UNMODIFIED base strands + MNC of the MODIFIED
    strands, concatenated.
    Shape: (140,)  (70 base + 70 modified)

    Why both: near-fully-modified strands collapse to mostly M/F under the 35-symbol
    composition, losing the underlying sequence. Adding the base composition restores
    that sequence signal — empirically PCC 0.48 vs 0.37 for modified-only on this data.
    For unmodified siRNAs base == modified, so the two halves coincide.
    """
    base_sense, base_antisense = _resolve_base(sense, antisense, base_sense, base_antisense)
    return np.concatenate([mnc_full(base_sense, base_antisense), mnc_full(sense, antisense)])


def features_model_b(sense, antisense, base_sense=None, base_antisense=None) -> np.ndarray:
    """
    Model-B: Model-A features + BIN(modified antisense seed 13 nt from 5'-end).
    Shape: (140 + 35×13,) = (140 + 455,) = (595,)
    Adds the critical seed-region positional pattern of the antisense strand.
    """
    base = features_model_a(sense, antisense, base_sense, base_antisense)
    bin_seed = bin_antisense_seed(antisense, seed_len=13)
    return np.concatenate([base, bin_seed])


def features_model_c(sense, antisense, base_sense=None, base_antisense=None) -> np.ndarray:
    """
    Model-C: Model-A features + BIN(modified antisense last 8 nt from 3'-end).
    Shape: (140 + 35×8,) = (140 + 280,) = (420,)
    Adds 3'-end positional info which affects RISC clamp stability.
    """
    base = features_model_a(sense, antisense, base_sense, base_antisense)
    bin_tail = bin_antisense_tail(antisense, tail_len=8)
    return np.concatenate([base, bin_tail])


# ─── GBM feature builder (position-aware + experimental condition) ─────────────
#
# The LightGBM models use a richer, lower-dimensional feature set than the SVR MNC
# vectors. On top of base+modified MNC (140-d) we add:
#   • position-aware modification density (where mods sit, not just how many)
#   • base GC content per strand
#   • two experimental-condition features (dose, time) — the proven confound in the
#     HelixZero data. Fixing these at inference yields "predicted inhibition at a
#     standard assay condition", a well-defined quantity the app can report.
#
# Reference condition used at inference (the most common assay in the catalog):
REF_CONC_NM: float = 10.0
REF_TIME_H:  float = 24.0

# mRNA feature fallback means (from training data) — loaded lazily
_MRNA_MEAN_FILE = Path(__file__).parent.parent / "data" / "mrna_feature_means.json"
_MRNA_FALLBACK: Optional[np.ndarray] = None
_MRNA_FEAT_NAMES: list[str] = []
_MRNA_FEAT_CHECKED = False


def _load_feat_names() -> list[str]:
    """Load feature names from the means file (authoritative source)."""
    global _MRNA_FEAT_CHECKED
    if _MRNA_FEAT_CHECKED:
        return _MRNA_FEAT_NAMES
    try:
        import json
        raw = json.loads(_MRNA_MEAN_FILE.read_text())
        _MRNA_FEAT_NAMES.clear()
        _MRNA_FEAT_NAMES.extend(raw.keys())
    except Exception:
        pass
    _MRNA_FEAT_CHECKED = True
    return _MRNA_FEAT_NAMES


def _get_mrna_fallback() -> np.ndarray:
    global _MRNA_FALLBACK
    if _MRNA_FALLBACK is not None:
        return _MRNA_FALLBACK
    try:
        import json
        raw = json.loads(_MRNA_MEAN_FILE.read_text())
        _MRNA_FALLBACK = np.array(list(raw.values()), dtype=np.float32)
    except Exception:
        _MRNA_FALLBACK = np.array([], dtype=np.float32)
    return _MRNA_FALLBACK


def _gc_fraction(seq: str) -> float:
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)


def _mod_positions(mod_seq: str, base_seq: str) -> List[int]:
    """Indices (0-based) where the modified-symbol strand differs from its base."""
    return [i for i, (m, b) in enumerate(zip(mod_seq, base_seq)) if m != b]


def _mod_density_features(mod_seq: str, base_seq: str) -> List[float]:
    """
    Position-aware modification features for one strand:
      [overall mod fraction, seed-region (pos 1-8) mod fraction,
       3'-tail (last 3) mod fraction, modified-position count (normalized)].
    """
    n = max(len(mod_seq), 1)
    pos = _mod_positions(mod_seq, base_seq)
    if not pos:
        return [0.0, 0.0, 0.0, 0.0]
    seed = sum(1 for p in pos if p < 8) / 8.0
    tail = sum(1 for p in pos if p >= len(mod_seq) - 3) / 3.0
    return [len(pos) / n, seed, tail, len(pos) / 21.0]


def features_gbm(
    sense, antisense, base_sense=None, base_antisense=None,
    conc_nM: Optional[float] = None, time_h: Optional[float] = None,
    mrna_feat: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    Feature vector for the LightGBM models.
    Shape: 152 (standard) or 152 + N (with N extra features).

    When mrna_feat is provided (or available from fallback means JSON),
    it is appended after the condition features.
    """
    base_sense, base_antisense = _resolve_base(sense, antisense, base_sense, base_antisense)
    comp = features_model_a(sense, antisense, base_sense, base_antisense)  # 140-d

    extra = (
        _mod_density_features(sense, base_sense)
        + _mod_density_features(antisense, base_antisense)
        + [_gc_fraction(base_sense), _gc_fraction(base_antisense)]
    )

    c = REF_CONC_NM if conc_nM is None or not np.isfinite(conc_nM) else conc_nM
    t = REF_TIME_H if time_h is None or not np.isfinite(time_h) else time_h
    cond = [float(np.log10(c + 0.01)), t / 24.0]

    vec = np.concatenate([comp, np.array(extra + cond, dtype=np.float32)])
    mf = mrna_feat if mrna_feat is not None else _get_mrna_fallback()
    if len(mf):
        vec = np.concatenate([vec, mf.astype(np.float32)])
    return vec


def extract_batch_gbm(
    sense_list: List[str],
    antisense_list: List[str],
    base_sense_list: List[str] = None,
    base_antisense_list: List[str] = None,
    conc_list: List[float] = None,
    time_list: List[float] = None,
    mrna_feat_list: List[np.ndarray] = None,
) -> np.ndarray:
    """Batch version of features_gbm. Missing condition/base lists default to reference."""
    n = len(sense_list)
    bs = base_sense_list if base_sense_list is not None else [None] * n
    ba = base_antisense_list if base_antisense_list is not None else [None] * n
    cc = conc_list if conc_list is not None else [None] * n
    tt = time_list if time_list is not None else [None] * n
    ff = mrna_feat_list if mrna_feat_list is not None else [None] * n
    return np.array(
        [features_gbm(s, a, bs_i, ba_i, c_i, t_i, f_i)
         for s, a, bs_i, ba_i, c_i, t_i, f_i in zip(sense_list, antisense_list, bs, ba, cc, tt, ff)],
        dtype=np.float32,
    )


# ─── batch extractor ──────────────────────────────────────────────────────────

MODEL_FEATURE_FN = {
    "A": features_model_a,
    "B": features_model_b,
    "C": features_model_c,
}


def extract_batch(
    sense_list: List[str],
    antisense_list: List[str],
    model: str = "A",
    base_sense_list: List[str] = None,
    base_antisense_list: List[str] = None,
) -> np.ndarray:
    """
    Extract features for a batch of siRNA pairs.

    Parameters
    ----------
    sense_list          : modified-symbol sense strands
    antisense_list      : modified-symbol antisense strands
    model               : "A", "B", or "C"
    base_sense_list     : optional unmodified base sense strands (defaults to sense_list)
    base_antisense_list : optional unmodified base antisense strands

    Returns
    -------
    np.ndarray of shape (n_samples, n_features)
    """
    fn = MODEL_FEATURE_FN[model]
    n = len(sense_list)
    bs = base_sense_list if base_sense_list is not None else [None] * n
    ba = base_antisense_list if base_antisense_list is not None else [None] * n
    return np.array(
        [fn(s, a, bs_i, ba_i) for s, a, bs_i, ba_i in zip(sense_list, antisense_list, bs, ba)],
        dtype=np.float32,
    )
