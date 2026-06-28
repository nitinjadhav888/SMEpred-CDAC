"""
features.py — Feature Extraction pipelines for Machine Learning Models

This module transforms variable-length string sequences (siRNAs with or without 
chemical modifications) into structured, fixed-dimensional numeric tensors for 
ingestion by the LightGBM models.

Pipelines:
1. Model B v4 (HelixZero Unified Model):
   - Extracts 1,467-dimensional positional binary features.
   - Captures exact modification types (30 unique symbols) per position.
   - Extracts aggregate structural statistics (GC%, cleavage site LNA count, etc.).

2. Naked V4 Model (Baseline Unmodified Model):
   - Extracts 214-dimensional sequence-composition features.
   - Includes positional one-hot encoding for A/U/G/C.
   - Computes Tri-Nucleotide Composition (TNC) frequencies.
"""

from typing import List, Optional, Dict
import numpy as np


# ─── Model B (Modified) Feature Extractor ─────────────────────────────────────

# Mapping of raw modification symbols to semantic feature names
_MODIFICATION_MAP: Dict[str, str] = {
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

_MOD_CATEGORIES: List[str] = sorted(
    {value.replace('is_', '') for value in _MODIFICATION_MAP.values()}
)
_POSITION_RANGE = range(1, 22)


def extract_positional_features_batch(
    sense_list: List[str],
    antisense_list: List[str],
    base_sense_list: Optional[List[str]] = None,
    base_antisense_list: Optional[List[str]] = None,
    conc_list: Optional[List[float]] = None,
) -> np.ndarray:
    """
    Batch extraction of 1,467-dimensional positional features for chemically modified siRNAs.
    """
    num_samples = len(sense_list)
    
    # Handle optional lists by filling with None
    base_senses = base_sense_list if base_sense_list is not None else [None] * num_samples
    base_antisenses = base_antisense_list if base_antisense_list is not None else [None] * num_samples
    concentrations = conc_list if conc_list is not None else [None] * num_samples
    
    feature_rows = []
    
    for sense_seq, anti_seq, base_sense, base_anti, conc in zip(
        sense_list, antisense_list, base_senses, base_antisenses, concentrations
    ):
        # Fallback to modified sequence if base sequence is absent
        effective_base_sense = base_sense if base_sense is not None else sense_seq
        effective_base_anti = base_anti if base_anti is not None else anti_seq
        
        row_features = _extract_single_modified_features(
            sense_seq, anti_seq, effective_base_sense, effective_base_anti, conc
        )
        feature_rows.append(row_features)
        
    return np.array(feature_rows, dtype=np.float32)


def _extract_single_modified_features(
    sense: str, 
    antisense: str,
    base_sense: str, 
    base_antisense: str,
    conc_nM: Optional[float] = None
) -> List[float]:
    """
    Internal helper to extract the full feature vector for a single cm-siRNA duplex.
    """
    features: List[float] = []

    # 1. Positional Binary Features (per strand, per position)
    for strand_key, seq, base_seq in [
        ("ss", sense, base_sense), ("as", antisense, base_antisense)
    ]:
        for pos in _POSITION_RANGE:
            zero_index = pos - 1
            nucleotide = seq[zero_index] if zero_index < len(seq) else ''
            base_nucleotide = base_seq[zero_index] if zero_index < len(base_seq) else ''
            
            # If nucleotide differs from parent, it's a modification symbol
            modification_char = nucleotide if nucleotide != base_nucleotide else ''
            is_modified = int(modification_char != '')
            is_canonical = 0 if is_modified else 1

            mapped_type = _MODIFICATION_MAP.get(modification_char, '')
            
            for typename in _MODIFICATION_MAP.values():
                features.append(float(mapped_type == typename))
                
            features.append(float(is_canonical))
            features.append(float(is_modified))

    # 2. Aggregate Sequence Statistics (per strand)
    for strand_key, seq, base_seq in [
        ("ss", sense, base_sense), ("as", antisense, base_antisense)
    ]:
        seq_length = len(seq)
        mod_counts = {mod_type: 0 for mod_type in _MOD_CATEGORIES}
        total_modifications = 0
        
        for i in range(min(seq_length, 21)):
            nucleotide = seq[i] if i < len(seq) else ''
            base_nucleotide = base_seq[i] if i < len(base_seq) else ''
            
            if nucleotide != base_nucleotide:
                total_modifications += 1
                type_name = _MODIFICATION_MAP.get(nucleotide, '').replace('is_', '')
                if type_name in mod_counts:
                    mod_counts[type_name] += 1

        fraction_modified = total_modifications / 21.0

        # Sub-region analysis: Seed (2-8) and Cleavage (9-11)
        seed_indices = list(range(1, 8))
        seed_2f = sum(1 for p in seed_indices if p < seq_length and seq[p] == 'F')
        seed_2ome = sum(1 for p in seed_indices if p < seq_length and seq[p] == 'M')
        
        cleavage_indices = list(range(8, 11))
        cleave_2f = sum(1 for p in cleavage_indices if p < seq_length and seq[p] == 'F')
        cleave_2ome = sum(1 for p in cleavage_indices if p < seq_length and seq[p] == 'M')
        cleave_lna = sum(1 for p in cleavage_indices if p < seq_length and seq[p] == 'L')

        gc_count = sum(1 for char in base_seq[:21].upper() if char in ('G', 'C'))
        gc_content = round(gc_count / min(len(base_seq), 21), 6) if base_seq else 0.5

        # Terminus protections
        term_5_ps = 1.0 if (seq_length > 0 and seq[0] == 'S') else 0.0
        term_3_ps = 1.0 if (seq_length > 20 and seq[20] == 'S') else 0.0

        # Append aggregated features
        for mod_type in _MOD_CATEGORIES:
            features.append(float(mod_counts[mod_type]))
            
        features.extend([
            fraction_modified, 
            seed_2f / 7.0, 
            seed_2ome / 7.0,
            float(cleave_2f), 
            float(cleave_2ome), 
            float(cleave_lna),
            gc_content, 
            term_5_ps, 
            term_3_ps,
        ])

    # 3. Experimental Parameters
    if conc_nM is not None and conc_nM > 0:
        log_concentration = float(np.log1p(conc_nM))
    else:
        # Default proxy for standard 10nM transfection experiments
        log_concentration = float(np.log1p(10.0))
        
    features.append(log_concentration)

    return features


# ─── Naked V4 (Unmodified) Feature Extractor ──────────────────────────────────

_CANONICAL_MAP = {"A": 0, "C": 1, "G": 2, "U": 3}


def _pad_sequence_to_21(sequence: str) -> str:
    """Ensures sequences are strictly 21 nucleotides via 3' Poly-A padding."""
    if len(sequence) >= 21:
        return sequence[:21]
    return sequence + "A" * (21 - len(sequence))


def extract_batch_v4(sense_list: List[str], antisense_list: List[str]) -> np.ndarray:
    """
    Batch extraction of 214-dimensional features for unmodified siRNAs.
    Includes explicit A/U/G/C positional one-hot encoding, and Tri-Nucleotide 
    Composition (TNC) normalized frequencies.
    """
    num_samples = len(sense_list)
    feature_matrix = np.zeros((num_samples, 214), dtype=np.float32)
    base_map = _CANONICAL_MAP
    
    for row_idx, (sense_seq, anti_seq) in enumerate(zip(sense_list, antisense_list)):
        padded_sense = _pad_sequence_to_21(sense_seq)
        padded_anti = _pad_sequence_to_21(anti_seq)
        
        # Positional One-Hot Encoding (4 bases * 21 pos = 84 features per strand)
        for pos in range(21):
            base_idx = base_map.get(padded_sense[pos], 0)
            feature_matrix[row_idx, (pos * 4) + base_idx] = 1.0
            
        # Tri-Nucleotide Composition (Sense) -> 64 features (4^3)
        for k in range(19):
            base_1 = base_map.get(padded_sense[k], 0)
            base_2 = base_map.get(padded_sense[k+1], 0)
            base_3 = base_map.get(padded_sense[k+2], 0)
            # Index calculation: (b1 * 16) + (b2 * 4) + b3
            feature_matrix[row_idx, 84 + (base_1 * 16) + (base_2 * 4) + base_3] += 1.0
            
        feature_matrix[row_idx, 84:148] /= 19.0  # Normalize TNC counts to frequencies
        
        # Tri-Nucleotide Composition (Antisense) -> 64 features
        for k in range(19):
            base_1 = base_map.get(padded_anti[k], 0)
            base_2 = base_map.get(padded_anti[k+1], 0)
            base_3 = base_map.get(padded_anti[k+2], 0)
            feature_matrix[row_idx, 148 + (base_1 * 16) + (base_2 * 4) + base_3] += 1.0
            
        feature_matrix[row_idx, 148:212] /= 19.0

        # Global GC content (Sense and Antisense) -> 2 features
        feature_matrix[row_idx, 212] = (padded_sense.count("G") + padded_sense.count("C")) / 21.0
        feature_matrix[row_idx, 213] = (padded_anti.count("G") + padded_anti.count("C")) / 21.0
        
    return feature_matrix
