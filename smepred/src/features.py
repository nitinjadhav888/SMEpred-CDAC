"""
features.py — Feature extraction for two models:
  - Model B v4 (HelixZero unified): positional binary features, 31 mod symbols, 1467-d
  - Naked V4: position one-hot + tri-nucleotide composition + GC, 214-d
"""

from typing import List
import numpy as np


# ══════════════════════════════════════════════════════════════════════
# Positional binary features — Model B v4
# ══════════════════════════════════════════════════════════════════════

_MOD_CHAR_MAP = {
    'F': 'is_2F', 'M': 'is_2OMe', 'L': 'is_LNA',
    'D': 'is_DNA', 'E': 'is_MOE',
    'B': 'is_Benzyl', 'N': 'is_4thio', 'I': 'is_FANA',
    'Z': 'is_ZOMe', 'Y': 'is_ENA',
    'S': 'is_PS', 'P': 'is_Borano',
    'R': 'is_MePhos', 'H': 'is_PhosAmid',
    'V': 'is_m5C', 'W': 'is_PseudoU',
    'J': 'is_Inosine', 'K': 'is_2thioU', 'O': 'is_DihydroU',
    '1': 'is_5Phos', '2': 'is_3P',
    '3': 'is_5OMe', '5': 'is_PEG',
    '6': 'is_UNA', '7': 'is_ANA',
    '8': 'is_GNA', '9': 'is_TNA',
    '4': 'is_Conj', 'Q': 'is_Abasic',
    'U': 'is_ModU', 'X': 'is_ModX',
}
_MOD_TYPES = sorted({v.replace('is_', '') for v in _MOD_CHAR_MAP.values()})
_POS_RANGE = range(1, 22)


def extract_positional_features_batch(
    sense_list: List[str],
    antisense_list: List[str],
    base_sense_list: List[str] = None,
    base_antisense_list: List[str] = None,
    conc_list: List[float] = None,
) -> np.ndarray:
    n = len(sense_list)
    bs = base_sense_list if base_sense_list is not None else [None] * n
    ba = base_antisense_list if base_antisense_list is not None else [None] * n
    cc = conc_list if conc_list is not None else [None] * n
    rows = []
    for s, a, b_s, b_a, c in zip(sense_list, antisense_list, bs, ba, cc):
        rows.append(_positional_features_one(s, a, b_s or s, b_a or a, c))
    return np.array(rows, dtype=np.float32)


def _positional_features_one(sense: str, antisense: str,
                              base_sense: str, base_antisense: str,
                              conc_nM: float = None) -> List[float]:
    feats = []

    for strand_key, seq, base_seq in [
        ("ss", sense, base_sense), ("as", antisense, base_antisense)
    ]:
        for pos in _POS_RANGE:
            i = pos - 1
            nt = seq[i] if i < len(seq) else ''
            base_nt = base_seq[i] if i < len(base_seq) else ''
            mod_char = nt if nt != base_nt else ''
            is_mod = int(mod_char != '')
            is_canon = 0 if is_mod else 1

            t = _MOD_CHAR_MAP.get(mod_char, '')
            for typename in _MOD_CHAR_MAP.values():
                feats.append(float(t == typename))
            feats.append(float(is_canon))
            feats.append(float(is_mod))

    for strand_key, seq, base_seq in [
        ("ss", sense, base_sense), ("as", antisense, base_antisense)
    ]:
        n = len(seq)
        counts = {t: 0 for t in _MOD_TYPES}
        n_mod = 0
        for i in range(min(n, 21)):
            nt = seq[i] if i < len(seq) else ''
            base_nt = base_seq[i] if i < len(base_seq) else ''
            if nt != base_nt:
                n_mod += 1
                type_name = _MOD_CHAR_MAP.get(nt, '').replace('is_', '')
                if type_name in counts:
                    counts[type_name] += 1

        frac_mod = n_mod / 21.0

        seed_positions = list(range(1, 8))
        seed_2f = sum(1 for p in seed_positions if p < len(seq) and seq[p] == 'F')
        seed_2ome = sum(1 for p in seed_positions if p < len(seq) and seq[p] == 'M')
        seed_frac_2f = seed_2f / 7.0
        seed_frac_2ome = seed_2ome / 7.0

        cleave_positions = list(range(8, 11))
        cleave_2f = sum(1 for p in cleave_positions if p < len(seq) and seq[p] == 'F')
        cleave_2ome = sum(1 for p in cleave_positions if p < len(seq) and seq[p] == 'M')
        cleave_lna = sum(1 for p in cleave_positions if p < len(seq) and seq[p] == 'L')

        gc_count = sum(1 for ch in base_seq[:21].upper() if ch in ('G', 'C'))
        gc_content = round(gc_count / min(len(base_seq), 21), 6) if base_seq else 0.5

        term_5_ps = 1.0 if (len(seq) > 0 and seq[0] == 'S') else 0.0
        term_3_ps = 1.0 if (len(seq) > 20 and seq[20] == 'S') else 0.0

        for mt in _MOD_TYPES:
            feats.append(float(counts[mt]))
        feats.extend([
            frac_mod, seed_frac_2f, seed_frac_2ome,
            float(cleave_2f), float(cleave_2ome), float(cleave_lna),
            gc_content, term_5_ps, term_3_ps,
        ])

    if conc_nM is not None and conc_nM > 0:
        log_conc = float(np.log1p(conc_nM))
    else:
        log_conc = float(np.log1p(10))
    feats.append(log_conc)

    return feats


# ══════════════════════════════════════════════════════════════════════
# V4 features — naked siRNA model (position one-hot + TNC + GC)
# ══════════════════════════════════════════════════════════════════════

_BASE_MAP_V4 = {"A": 0, "C": 1, "G": 2, "U": 3}


def _pad21_v4(seq: str) -> str:
    if len(seq) >= 21:
        return seq[:21]
    return seq + "A" * (21 - len(seq))


def extract_batch_v4(sense_list, antisense_list):
    n = len(sense_list)
    X = np.zeros((n, 214), dtype=np.float32)
    bm = _BASE_MAP_V4
    for i, (s, a) in enumerate(zip(sense_list, antisense_list)):
        sp = _pad21_v4(s)
        ap = _pad21_v4(a)
        for pos in range(21):
            idx = bm.get(sp[pos], 0)
            X[i, pos * 4 + idx] = 1.0
        for k in range(19):
            a1 = bm.get(sp[k], 0); a2 = bm.get(sp[k+1], 0); a3 = bm.get(sp[k+2], 0)
            X[i, 84 + a1*16 + a2*4 + a3] += 1.0
        X[i, 84:148] /= 19.0
        for k in range(19):
            a1 = bm.get(ap[k], 0); a2 = bm.get(ap[k+1], 0); a3 = bm.get(ap[k+2], 0)
            X[i, 148 + a1*16 + a2*4 + a3] += 1.0
        X[i, 148:212] /= 19.0
        X[i, 212] = (sp.count("G") + sp.count("C")) / 21.0
        X[i, 213] = (ap.count("G") + ap.count("C")) / 21.0
    return X
