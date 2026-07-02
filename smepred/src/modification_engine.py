"""
modification_engine.py — Chemical Modification Generator

This module applies chemical modifications to siRNA candidates. It supports
three distinct operation modes:

1. Single-Modification Scan
   Systematically applies each of the 30 chemical modification symbols to every 
   position (1-21) on both strands of a parent siRNA. This generates an exhaustive 
   1260-variant library to identify the single most effective modification point.

2. MultiModGen (Targeted Custom Modifications)
   Allows the user or downstream algorithms to apply specific modifications to 
   targeted positions across both strands simultaneously.

3. Beam Search Scan
   An intelligent, iterative search algorithm that combines top-performing single 
   modifications into multi-mod combinations, scoring them in rounds to find the 
   global biophysical optimum without brute-forcing millions of combinations.
"""

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Any, Set, Dict

logger = logging.getLogger(__name__)

# ─── Load Modification Definitions ──────────────────────────────────────────────

_MOD_FILE = Path(__file__).parent.parent / "data" / "modification_codes.json"
try:
    with _MOD_FILE.open("r", encoding="utf-8") as _f:
        _MOD_DATA = json.load(_f)
    CANONICAL_SYMBOLS: Set[str] = set(_MOD_DATA["canonical_symbols"])
    MODIFICATION_SYMBOLS: Set[str] = set(_MOD_DATA["modification_symbols"])
except Exception as e:
    logger.error(f"Failed to load modification codes: {e}")
    raise RuntimeError(f"Could not initialize modification engine: {e}")


# ─── Data Transfer Objects ────────────────────────────────────────────────────

@dataclass
class CmSiRNA:
    """
    Represents a Chemically Modified siRNA (cm-siRNA) variant.
    
    Attributes:
        sense (str): The chemically modified sense strand.
        antisense (str): The chemically modified antisense strand.
        mod_symbol (str): The symbol(s) representing the applied chemistry.
        mod_position (int): The 1-based index of the primary modification.
        mod_strand (str): The strand on which the modification occurs.
        parent_sense (str): The unmodified biological sense strand.
        parent_antisense (str): The unmodified biological antisense strand.
        mod_positions (str): Comma-separated list of all modified positions (for multi-mod).
        efficacy_score (float): The final biophysically adjusted efficacy score.
        delta_score (float): Efficacy improvement/loss relative to the parent.
        penalties (dict): Breakdown of biophysical penalties applied.
    """
    sense: str
    antisense: str
    mod_symbol: str
    mod_position: int
    mod_strand: str
    parent_sense: str
    parent_antisense: str
    mod_positions: str = ""
    efficacy_score: float = 0.0
    delta_score: float = 0.0
    penalties: Optional[Dict[str, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "sense": self.sense,
            "antisense": self.antisense,
            "mod_symbol": self.mod_symbol,
            "mod_position": self.mod_position,
            "mod_strand": self.mod_strand,
            "parent_sense": self.parent_sense,
            "parent_antisense": self.parent_antisense,
            "mod_positions": self.mod_positions,
        }
        if self.efficacy_score:
            result["efficacy_score"] = self.efficacy_score
        if self.delta_score:
            result["delta_score"] = self.delta_score
        if self.penalties:
            result["penalties"] = self.penalties
        return result


# ─── Helper Functions ─────────────────────────────────────────────────────────

def _apply_mod(sequence: str, position_1based: int, symbol: str) -> str:
    """
    Replaces a specific nucleotide with a chemical modification symbol.
    """
    if not (1 <= position_1based <= len(sequence)):
        logger.error(f"Position {position_1based} out of bounds for sequence length {len(sequence)}")
        raise ValueError(
            f"Position {position_1based} is out of range for sequence of length {len(sequence)}."
        )
    zero_indexed = position_1based - 1
    return sequence[:zero_indexed] + symbol + sequence[zero_indexed + 1:]


def _parse_multimod_input(
    mod_symbols_str: str, positions_str: str
) -> List[Tuple[str, List[int]]]:
    """
    Parses comma-delimited strings defining targeted modifications.
    Format: "F,,M" and "2,5,7,,10,12" -> Apply F at 2,5,7 and M at 10,12.
    """
    mod_groups = [m.strip() for m in mod_symbols_str.split(",,")]
    pos_groups = [p.strip() for p in positions_str.split(",,")]

    if len(mod_groups) != len(pos_groups):
        raise ValueError(
            "Mismatched group count: Number of modification groups must equal number of position groups."
        )

    parsed_instructions = []
    for symbol, pos_string in zip(mod_groups, pos_groups):
        clean_symbol = symbol.strip()
        if clean_symbol not in MODIFICATION_SYMBOLS | CANONICAL_SYMBOLS:
            logger.error(f"Unknown modification symbol detected: {clean_symbol}")
            raise ValueError(f"Unknown modification symbol: '{clean_symbol}'")
            
        parsed_positions = [int(p.strip()) for p in pos_string.split(",") if p.strip()]
        parsed_instructions.append((clean_symbol, parsed_positions))
        
    return parsed_instructions


# ─── Mode 1: Single-Modification Scan ─────────────────────────────────────────

def single_mod_scan(
    sense: str,
    antisense: str,
    target_symbols: Optional[List[str]] = None,
) -> List[CmSiRNA]:
    """
    Generates an exhaustive single-modification combinatorial library.
    """
    if target_symbols is None:
        target_symbols = list(MODIFICATION_SYMBOLS)

    generated_variants: List[CmSiRNA] = []

    for symbol in target_symbols:
        # Scan sense strand
        for pos in range(1, len(sense) + 1):
            modified_sense = _apply_mod(sense, pos, symbol)
            generated_variants.append(CmSiRNA(
                sense=modified_sense,
                antisense=antisense,
                mod_symbol=symbol,
                mod_position=pos,
                mod_strand="sense",
                parent_sense=sense,
                parent_antisense=antisense,
            ))
            
        # Scan antisense strand
        for pos in range(1, len(antisense) + 1):
            modified_antisense = _apply_mod(antisense, pos, symbol)
            generated_variants.append(CmSiRNA(
                sense=sense,
                antisense=modified_antisense,
                mod_symbol=symbol,
                mod_position=pos,
                mod_strand="antisense",
                parent_sense=sense,
                parent_antisense=antisense,
            ))

    return generated_variants


# ─── Mode 2: Targeted MultiModGen ─────────────────────────────────────────────

def multimod_gen(
    sense: str,
    antisense: str,
    sense_mods: str = "",
    sense_positions: str = "",
    antisense_mods: str = "",
    antisense_positions: str = "",
) -> CmSiRNA:
    """
    Applies precise, targeted modifications simultaneously across both strands.
    """
    mutable_sense = list(sense)
    mutable_antisense = list(antisense)

    if sense_mods and sense_positions:
        sense_instructions = _parse_multimod_input(sense_mods, sense_positions)
        for symbol, positions in sense_instructions:
            for pos in positions:
                if not (1 <= pos <= len(mutable_sense)):
                    raise ValueError(f"Sense position {pos} out of range.")
                mutable_sense[pos - 1] = symbol

    if antisense_mods and antisense_positions:
        antisense_instructions = _parse_multimod_input(antisense_mods, antisense_positions)
        for symbol, positions in antisense_instructions:
            for pos in positions:
                if not (1 <= pos <= len(mutable_antisense)):
                    raise ValueError(f"Antisense position {pos} out of range.")
                mutable_antisense[pos - 1] = symbol

    return CmSiRNA(
        sense="".join(mutable_sense),
        antisense="".join(mutable_antisense),
        mod_symbol="multi",
        mod_position=0,
        mod_strand="both",
        parent_sense=sense,
        parent_antisense=antisense,
    )


# ─── Mode 3: Combinatorial Beam Search Scan ───────────────────────────────────

def _is_sterically_viable(modified_strand: str, parent_strand: str) -> bool:
    """
    Rejects modification patterns with >= 3 consecutive bulky modifications
    (LNA, ENA, MOE) that create extreme backbone rigidity incompatible with
    Ago2 accommodation (Obad et al. 2011; ESC+ clinical guidelines).
    """
    consecutive_bulky = 0
    for i in range(len(modified_strand)):
        char = modified_strand[i]
        parent_char = parent_strand[i] if i < len(parent_strand) else char
        if char != parent_char and char in ('L', 'Y', 'E'):
            consecutive_bulky += 1
            if consecutive_bulky >= 3:
                return False
        else:
            consecutive_bulky = 0
    return True


def multi_mod_scan(
    sense: str,
    antisense: str,
    max_mods: int = 2,
    beam_width: int = 20,
    model_key: str = "B",
    full_scan: bool = False,
    single_results: Optional[List[Any]] = None,
    parent_score: Optional[float] = None,
    seed_variant: Optional[Any] = None,
    calibrator_key: Optional[str] = None,
    normalize_mode: str = "clip",
) -> List[CmSiRNA]:
    """
    Heuristically explores the vast combinatoric space of multi-modified siRNAs.
    Uses an iterative beam search to stack highly effective modifications while 
    pruning sub-optimal branches to avoid computational explosion.
    """
    # Lazy imports required to prevent circular dependency with predictor.py
    from .predictor import predict_modified, _get_model, _normalize_scores
    from .features import extract_positional_features_batch
    from .biophysics import calculate_adjusted_efficacy
    from collections import defaultdict

    logger.info("Starting combinatorial beam search.")

    if single_results is None:
        prediction_output = predict_modified(
            sense, antisense, mode="scan", model_key=model_key, full_scan=full_scan
        )
        parent_score = prediction_output.get("parent_score_raw", prediction_output["parent_score"])
        single_results = prediction_output["results"]
    elif parent_score is None:
        raise ValueError("parent_score must be provided when single_results is pre-calculated.")

    # Calculate baseline for delta comparisons
    parent_adjusted_score, _, _ = calculate_adjusted_efficacy(
        parent_score, sense, antisense, sense, antisense
    )

    def _score_variants_batch(variants: List[CmSiRNA], chunk_size: int = 200) -> List[CmSiRNA]:
        """Internal helper to batch-score variants using Model B, in chunks to limit memory."""
        if not variants:
            return []

        model = _get_model("B")
        scored_variants = []

        for i in range(0, len(variants), chunk_size):
            chunk = variants[i:i + chunk_size]
            s_list = [v.sense for v in chunk]
            a_list = [v.antisense for v in chunk]
            ps_list = [v.parent_sense for v in chunk]
            pa_list = [v.parent_antisense for v in chunk]

            feature_matrix = extract_positional_features_batch(s_list, a_list, ps_list, pa_list)
            raw_predictions = model.predict(feature_matrix)
            normalized_scores = _normalize_scores(raw_predictions, mode="rescale")

            for variant, raw_score in zip(chunk, normalized_scores):
                adj_score, penalties, _ = calculate_adjusted_efficacy(
                    float(raw_score), variant.sense, variant.antisense,
                    variant.parent_sense, variant.parent_antisense
                )
                variant.efficacy_score = round(adj_score, 2)
                variant.delta_score = round(adj_score - parent_adjusted_score, 2)
                variant.penalties = penalties
                scored_variants.append(variant)

        return scored_variants

    # Initialize the beam with diverse, high-performing single modifications
    mod_groups: Dict[str, List[Any]] = defaultdict(list)
    for result in single_results:
        mod_groups[result.mod_symbol].append(result)

    for symbol in mod_groups:
        mod_groups[symbol].sort(key=lambda r: r.efficacy_score, reverse=True)

    diversified_beam = []
    max_entries = max(len(lst) for lst in mod_groups.values())
    
    # Round-robin selection ensures chemical diversity in the starting beam
    for rank in range(max_entries):
        for symbol in sorted(mod_groups.keys()):
            if rank < len(mod_groups[symbol]):
                diversified_beam.append(mod_groups[symbol][rank])
            if len(diversified_beam) >= beam_width:
                break
        if len(diversified_beam) >= beam_width:
            break

    initial_beam: List[CmSiRNA] = []
    if seed_variant is not None:
        initial_beam.append(seed_variant)
        
    for result in diversified_beam:
        if len(initial_beam) >= beam_width:
            break
        variant = CmSiRNA(
            sense=result.sense,
            antisense=result.antisense,
            mod_symbol=result.mod_symbol,
            mod_position=result.mod_position,
            mod_strand=result.mod_strand,
            parent_sense=sense,
            parent_antisense=antisense,
        )
        variant.efficacy_score = result.efficacy_score
        variant.delta_score = result.delta_score
        initial_beam.append(variant)

    # Begin Expansion Rounds
    current_beam = _score_variants_batch(initial_beam)
    current_beam.sort(key=lambda x: x.efficacy_score, reverse=True)
    all_evaluated_variants = list(current_beam)

    pairing_pool = sorted(single_results, key=lambda r: r.efficacy_score, reverse=True)[:beam_width * 2]
    history_best_scores = [current_beam[0].efficacy_score if current_beam else 0.0]

    for iteration in range(2, max_mods + 1):
        round_best_score = current_beam[0].efficacy_score if current_beam else 0.0
        
        # Plateau detection: Stop searching if the optimum hasn't improved meaningfully
        if iteration >= 4 and (round_best_score - history_best_scores[-3]) < 0.5:
            logger.info("Beam search plateau detected. Stopping early.")
            break
            
        history_best_scores.append(round_best_score)
        round_candidates = []
        explored_pairs = set()

        def _generate_signature(v: Any) -> tuple:
            return (
                getattr(v, 'mod_symbol', ''), 
                getattr(v, 'mod_position', 0), 
                getattr(v, 'mod_strand', ''), 
                getattr(v, 'mod_positions', '')
            )

        for base_variant in current_beam:
            for addon_variant in pairing_pool:
                sig_1 = _generate_signature(base_variant)
                sig_2 = _generate_signature(addon_variant)
                pair_signature = tuple(sorted([sig_1, sig_2]))
                
                if pair_signature in explored_pairs:
                    continue
                    
                explored_pairs.add(pair_signature)

                # Merge modifications
                mutable_sense = list(base_variant.parent_sense)
                mutable_antisense = list(base_variant.parent_antisense)
                tracking_symbols = []
                tracking_positions = []
                tracking_strands = []

                # Restore base variant modifications
                for i in range(len(sense)):
                    if base_variant.sense[i] != base_variant.parent_sense[i]:
                        mutable_sense[i] = base_variant.sense[i]
                        tracking_symbols.append(base_variant.sense[i])
                        tracking_positions.append(i + 1)
                        tracking_strands.append("sense")
                        
                for i in range(len(antisense)):
                    if base_variant.antisense[i] != base_variant.parent_antisense[i]:
                        mutable_antisense[i] = base_variant.antisense[i]
                        tracking_symbols.append(base_variant.antisense[i])
                        tracking_positions.append(i + 1)
                        tracking_strands.append("antisense")

                # Apply new addon modification
                if addon_variant.mod_strand == "sense":
                    if mutable_sense[addon_variant.mod_position - 1] != sense[addon_variant.mod_position - 1]:
                        continue  # Position already modified, skip clash
                    mutable_sense[addon_variant.mod_position - 1] = addon_variant.mod_symbol
                else:
                    if mutable_antisense[addon_variant.mod_position - 1] != antisense[addon_variant.mod_position - 1]:
                        continue
                    mutable_antisense[addon_variant.mod_position - 1] = addon_variant.mod_symbol
                    
                tracking_symbols.append(addon_variant.mod_symbol)
                tracking_positions.append(addon_variant.mod_position)
                tracking_strands.append(addon_variant.mod_strand)

                # Check steric viability — reject consecutive bulky mods
                if not _is_sterically_viable("".join(mutable_antisense), antisense):
                    continue
                if not _is_sterically_viable("".join(mutable_sense), sense):
                    continue

                round_candidates.append(CmSiRNA(
                    sense="".join(mutable_sense),
                    antisense="".join(mutable_antisense),
                    mod_symbol="+".join(tracking_symbols),
                    mod_position=tracking_positions[0],
                    mod_positions=",".join(str(p) for p in tracking_positions),
                    mod_strand="+".join(tracking_strands),
                    parent_sense=sense,
                    parent_antisense=antisense,
                ))

        scored_candidates = _score_variants_batch(round_candidates)
        scored_candidates.sort(key=lambda v: v.efficacy_score, reverse=True)
        
        current_beam = scored_candidates[:beam_width]
        all_evaluated_variants.extend(scored_candidates)

    # Deduplicate based on exact sequence string to prevent permutations clogging the top 100
    unique_variants = {}
    for v in all_evaluated_variants:
        seq_key = v.sense + "|" + v.antisense
        # If we somehow have identical sequences with different scores, keep the highest
        if seq_key not in unique_variants or v.efficacy_score > unique_variants[seq_key].efficacy_score:
            unique_variants[seq_key] = v
            
    final_variants = list(unique_variants.values())
    final_variants.sort(key=lambda v: v.efficacy_score, reverse=True)
    
    logger.info(f"Beam search complete. Evaluated {len(all_evaluated_variants)} total permutations. Returning {len(final_variants)} unique sequences.")
    return final_variants
