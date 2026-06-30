"""
filters.py — Candidate Safety, Toxicity, and Functionality Filters

Provides critical biological filtering mechanisms:
1. Seed Toxicity Prediction: Cross-references the candidate's 6-mer seed against 
   the Janas et al. (2018) empirical cell viability database (4,097 entries) to 
   predict off-target induced cytotoxicity.
2. Modification-Aware Mitigation: Detects if the user applied seed-rescuing chemical 
   modifications (e.g., 2'-OMe at position 2) that suppress innate miRNA-like toxicity.
3. Functional Rules: Enforces standard Reynolds/Ui-Tei biophysical design rules 
   (GC content limits, prevention of homopolymer runs and palindromes).
"""

from __future__ import annotations

import itertools
import re
import logging
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple, Dict

import pandas as pd

from .utils import calculate_gc_percentage, has_internal_palindrome

logger = logging.getLogger(__name__)

_TOX_PATH = Path(__file__).parent.parent / "data" / "oligoformer" / "cell_viability.tsv"


# ─── Seed Toxicity Lookup ─────────────────────────────────────────────────────

@lru_cache(maxsize=1)
def _load_toxicity_table() -> Dict[str, float]:
    """
    Loads the seed -> cell-viability mapping into a cached dictionary.
    
    Why: High-throughput screening (Janas et al., Mol Cell 2018) demonstrated that 
    specific 6-mer seeds are inherently toxic to human cells regardless of the target. 
    We cache this 4,000+ row table in memory for microsecond lookups during generation.
    """
    try:
        df = pd.read_csv(_TOX_PATH, sep="\t")
        # Enforce uppercase RNA formatting for reliable lookup keys
        return dict(zip(df["Seed"].str.upper(), df["cell_viability"].astype(float)))
    except Exception as e:
        logger.error(f"Failed to load toxicity table from {_TOX_PATH}: {e}")
        return {}


def _extract_seed(antisense: str) -> str:
    """
    Extracts the critical 6-mer seed region (positions 2-7, 1-indexed).
    """
    normalized_strand = antisense.upper().replace("T", "U")
    return normalized_strand[1:7]


def get_toxicity_score(antisense: str) -> Optional[float]:
    """
    Retrieves the predicted cell viability percentage for the candidate's seed.
    
    Args:
        antisense (str): The antisense strand sequence.
        
    Returns:
        Optional[float]: Cell viability percentage. Lower means more toxic. 
                         Returns None if the seed is undocumented.
    """
    seed_region = _extract_seed(antisense)
    score = _load_toxicity_table().get(seed_region)
    if score is None:
        logger.debug(f"Seed {seed_region} not found in empirical toxicity database.")
    return score


def get_toxicity_label(viability: Optional[float], safe_threshold: float = 70.0) -> str:
    """
    Translates raw cell viability percentages into human-readable clinical labels.
    """
    if viability is None:
        return "Unknown"
    if viability >= safe_threshold:
        return "Safe"
    if viability >= 50.0:
        return "Caution"
    return "Toxic"


# ─── Modification-Aware Toxicity Mitigation ───────────────────────────────────

# Modifications established in literature to suppress seed-mediated off-target binding:
# M (2'-OMe), F (2'-Fluoro), L (LNA), E (2'-MOE).
_SEED_RESCUING_MODS = frozenset({"M", "F", "L", "E"})
_MOD_NOMENCLATURE = {"M": "2'-OMe", "F": "2'-Fluoro", "L": "LNA", "E": "2'-MOE"}


def check_seed_rescue(modified_antisense: str) -> Tuple[List[Tuple[int, str]], str]:
    """
    Detects if seed-rescuing chemical modifications are present in the critical region.
    
    Why: A biologically toxic sequence can be "rescued" (rendered safe) if specific 
    steric modifications are placed in the seed region (positions 2-7), which disrupts 
    off-target miRNA-like binding (Jackson et al., RNA 2006).
    """
    upper_mod_strand = modified_antisense.upper()
    rescue_modifications = []
    
    # Scan positions 2 through 7 (indices 1 through 6)
    for i in range(1, min(7, len(upper_mod_strand))):
        if upper_mod_strand[i] in _SEED_RESCUING_MODS:
            rescue_modifications.append((i + 1, upper_mod_strand[i]))
            
    if not rescue_modifications:
        return [], ""
        
    mitigation_notes = [f"{_MOD_NOMENCLATURE[symbol]} @ pos {pos}" for pos, symbol in rescue_modifications]
    tooltip_note = "Seed off-target rescue: " + ", ".join(mitigation_notes)
    return rescue_modifications, tooltip_note


def toxicity_for_modified(
    modified_antisense: str, base_antisense: str
) -> Tuple[Optional[float], str, str]:
    """
    Evaluates toxicity for a chemically modified siRNA, applying mitigation overrides.
    
    Strategy: We first evaluate the unmodified (parent) baseline toxicity. Then we 
    scan the modified strand for rescuing chemistry. If a rescue is found in a Toxic 
    seed, we override the clinical label to "Mitigated".
    
    Returns:
        Tuple: (viability_percentage, clinical_label, mitigation_tooltip)
    """
    baseline_viability = get_toxicity_score(base_antisense)
    baseline_label = get_toxicity_label(baseline_viability)
    
    rescue_mods, mitigation_note = check_seed_rescue(modified_antisense)
    
    if rescue_mods:
        if baseline_label in {"Toxic", "Caution"}:
            logger.info("Toxic seed successfully mitigated via chemical modification.")
            return baseline_viability, "Mitigated", mitigation_note
        elif baseline_label == "Safe":
            # Pass the note forward even if already safe, for clinical completeness
            return baseline_viability, "Safe", mitigation_note
            
    return baseline_viability, baseline_label, ""


# ─── Functional Baseline Filters ──────────────────────────────────────────────

_HOMOPOLYMER_REGEX = re.compile(r"A{5}|U{5}|G{5}|C{5}")
_GC6_REGEXES = [re.compile("".join(p)) for p in itertools.product("GC", repeat=6)]


def check_functionality(sirna_strand: str) -> Tuple[bool, str]:
    """
    Evaluates whether the candidate violates baseline structural siRNA design rules.
    
    Why: A sequence may be non-toxic, but if it violates these thermodynamic boundaries, 
    it will fail to unwind or load into the RISC complex entirely, rendering it dead.
    """
    normalized_strand = sirna_strand.upper().replace("T", "U")
    
    gc_content = calculate_gc_percentage(normalized_strand)
    if not (30.0 <= gc_content <= 65.0):
        return False, f"GC {gc_content:.0f}% out of optimal 30-65% range"
        
    if _HOMOPOLYMER_REGEX.search(normalized_strand):
        return False, "5-base homopolymer run detected (prevents unwinding)"
        
    for gc_pattern in _GC6_REGEXES:
        if gc_pattern.search(normalized_strand):
            return False, "6-base contiguous GC run detected"
            
    if has_internal_palindrome(normalized_strand):
        return False, "Internal palindrome detected (forms stable hairpins)"
        
    return True, ""


# ─── Batch Annotation Helpers ─────────────────────────────────────────────────

def annotate_candidates(senses: List[str], antisenses: List[str]) -> List[Dict[str, Any]]:
    """
    Batch-annotates candidates with their toxicity scores and functional compliance flags.
    Used heavily by the `predictor` during sliding-window evaluation.
    """
    annotations = []
    for sense_strand, anti_strand in zip(senses, antisenses):
        viability = get_toxicity_score(anti_strand)
        is_functional_sense, reason_sense = check_functionality(sense_strand)
        is_functional_anti, reason_anti = check_functionality(anti_strand)
        is_functional = is_functional_sense and is_functional_anti
        failure_reason = reason_sense or reason_anti
        
        annotations.append({
            "toxicity_score": None if viability is None else round(viability, 1),
            "toxicity_label": get_toxicity_label(viability),
            "func_ok": is_functional,
            "func_reason": failure_reason,
        })
        
    return annotations
