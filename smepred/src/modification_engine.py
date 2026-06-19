"""
modification_engine.py — Chemical modification generator.

Two modes (inspired by the SMEpred approach, Dar et al. 2016):

MODE 1 — Single-Modification Scan
  For one chosen siRNA, systematically apply each of the 30 chemical
  modifications at every position (1–21) on both strands.
  Total variants = 30 modifications × 21 positions × 2 strands = 1260 cm-siRNAs.
  This lets you see which modification at which position maximally improves efficacy.

MODE 2 — MultiModGen (custom multiple modifications)
  User specifies one or more modification types and a list of positions for each.
  Multiple modification types are separated by ',,'.
  Example:
    siRNA:     GCAGCACGACUUCUUCAAGUU
    mods:      F,,M
    positions: 2,5,7,,10,12
  This means: apply 2'-Fluoro (F) at positions 2,5,7 AND 2'-OMe (M) at positions 10,12
  on the sense strand. Same for antisense if specified.
  Generates one cm-siRNA per combination submitted.

Sequence format:
  A cm-siRNA is stored as: <modified_sense>-<modified_antisense>
  Each position in the sequence is the nucleotide symbol (A/U/G/C) replaced by
  the modification symbol (e.g. F for 2'-Fluoro) at modified positions.
  Example: GFAGFACGACUUCUUCAAGUU-CUUGAAGAAGUCGUGCUGCUU
  (F at positions 2 and 4 on the sense strand).
"""

import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

# ─── load modification symbols ────────────────────────────────────────────────

_MOD_FILE = Path(__file__).parent.parent / "data" / "modification_codes.json"
with _MOD_FILE.open() as _f:
    _MOD_DATA = json.load(_f)

CANONICAL_SYMBOLS = set(_MOD_DATA["canonical_symbols"])        # A U G C T
MODIFICATION_SYMBOLS = set(_MOD_DATA["modification_symbols"])  # 30 mod symbols


# ─── data container ───────────────────────────────────────────────────────────

@dataclass
class CmSiRNA:
    """One chemically modified siRNA variant."""
    sense: str                  # modified sense strand
    antisense: str              # modified antisense strand (may also be modified)
    mod_symbol: str             # which modification was applied (e.g. "E" or "E+E")
    mod_position: int           # 1-based position of modification (first position for multi-mod)
    mod_strand: str             # "sense" or "antisense" (or "sense+sense" for multi-mod)
    parent_sense: str           # original unmodified sense strand
    parent_antisense: str       # original unmodified antisense strand
    mod_positions: str = ""     # all positions as comma-separated string for multi-mod (e.g. "4,6")

    def to_dict(self) -> dict:
        return {
            "sense":           self.sense,
            "antisense":       self.antisense,
            "mod_symbol":      self.mod_symbol,
            "mod_position":    self.mod_position,
            "mod_strand":      self.mod_strand,
            "parent_sense":    self.parent_sense,
            "parent_antisense": self.parent_antisense,
        }


# ─── helpers ──────────────────────────────────────────────────────────────────

def _apply_mod(seq: str, position_1based: int, symbol: str) -> str:
    """
    Replace the nucleotide at `position_1based` (1-indexed) with `symbol`.
    Returns the modified sequence string.
    """
    if not (1 <= position_1based <= len(seq)):
        raise ValueError(
            f"Position {position_1based} is out of range for sequence of length {len(seq)}."
        )
    idx = position_1based - 1
    return seq[:idx] + symbol + seq[idx + 1:]


def _parse_multimod_input(
    mods_str: str, positions_str: str
) -> List[Tuple[str, List[int]]]:
    """
    Parse MultiModGen input.

    mods_str      : modification symbols separated by ',,' for multiple types
                    e.g. "F,,M"
    positions_str : corresponding positions separated by ',' within a group,
                    groups separated by ',,'.
                    e.g. "2,5,7,,10,12"

    Returns list of (symbol, [positions]) tuples.
    """
    mod_groups = [m.strip() for m in mods_str.split(",,")]
    pos_groups = [p.strip() for p in positions_str.split(",,")]

    if len(mod_groups) != len(pos_groups):
        raise ValueError(
            "Number of modification groups must equal number of position groups. "
            "Separate groups with ',,'."
        )

    result = []
    for sym, pos_str in zip(mod_groups, pos_groups):
        sym = sym.strip()
        if sym not in MODIFICATION_SYMBOLS | CANONICAL_SYMBOLS:
            raise ValueError(f"Unknown modification symbol: '{sym}'")
        positions = [int(p.strip()) for p in pos_str.split(",") if p.strip()]
        result.append((sym, positions))
    return result


# ─── Mode 1: single-modification scan ─────────────────────────────────────────

def single_mod_scan(
    sense: str,
    antisense: str,
    mod_symbols: Optional[List[str]] = None,
) -> List[CmSiRNA]:
    """
    Apply each modification type at each position on both strands.

    Parameters
    ----------
    sense        : 21-nt sense strand sequence
    antisense    : 21-nt antisense strand sequence
    mod_symbols  : list of modification symbols to apply.
                   Defaults to all 30 modification symbols.

    Returns
    -------
    List[CmSiRNA] — up to 30 × 21 × 2 = 1260 variants.
    """
    if mod_symbols is None:
        mod_symbols = list(MODIFICATION_SYMBOLS)

    results: List[CmSiRNA] = []

    for sym in mod_symbols:
        # modify each position on sense strand
        for pos in range(1, len(sense) + 1):
            mod_sense = _apply_mod(sense, pos, sym)
            results.append(CmSiRNA(
                sense=mod_sense,
                antisense=antisense,
                mod_symbol=sym,
                mod_position=pos,
                mod_strand="sense",
                parent_sense=sense,
                parent_antisense=antisense,
            ))
        # modify each position on antisense strand (independent length)
        for pos in range(1, len(antisense) + 1):
            mod_antisense = _apply_mod(antisense, pos, sym)
            results.append(CmSiRNA(
                sense=sense,
                antisense=mod_antisense,
                mod_symbol=sym,
                mod_position=pos,
                mod_strand="antisense",
                parent_sense=sense,
                parent_antisense=antisense,
            ))

    return results


# ─── Mode 2: MultiModGen ──────────────────────────────────────────────────────

def multimod_gen(
    sense: str,
    antisense: str,
    sense_mods: str = "",
    sense_positions: str = "",
    antisense_mods: str = "",
    antisense_positions: str = "",
) -> CmSiRNA:
    """
    Apply multiple custom modifications to a siRNA at user-specified positions.

    Parameters
    ----------
    sense              : 21-nt sense strand
    antisense          : 21-nt antisense strand
    sense_mods         : modification symbols for sense strand e.g. "F,,M"
    sense_positions    : positions for sense strand mods e.g. "2,5,7,,10,12"
    antisense_mods     : modification symbols for antisense strand (same format)
    antisense_positions: positions for antisense strand mods

    Returns
    -------
    CmSiRNA with both strands modified as specified.
    """
    mod_sense = list(sense)
    mod_antisense = list(antisense)

    # apply sense strand modifications
    if sense_mods and sense_positions:
        ss_groups = _parse_multimod_input(sense_mods, sense_positions)
        for sym, positions in ss_groups:
            for pos in positions:
                if not (1 <= pos <= len(mod_sense)):
                    raise ValueError(f"Sense position {pos} out of range.")
                mod_sense[pos - 1] = sym

    # apply antisense strand modifications
    if antisense_mods and antisense_positions:
        as_groups = _parse_multimod_input(antisense_mods, antisense_positions)
        for sym, positions in as_groups:
            for pos in positions:
                if not (1 <= pos <= len(mod_antisense)):
                    raise ValueError(f"Antisense position {pos} out of range.")
                mod_antisense[pos - 1] = sym

    return CmSiRNA(
        sense="".join(mod_sense),
        antisense="".join(mod_antisense),
        mod_symbol="multi",
        mod_position=0,
        mod_strand="both",
        parent_sense=sense,
        parent_antisense=antisense,
    )


# ─── Mode 3: Multi-Modification Beam Search Scan ────────────────────────────────
# Builds on single-mod scan: take top-K single hits, combine, re-score, iterate.

def multi_mod_scan(
    sense: str,
    antisense: str,
    max_mods: int = 2,
    beam_width: int = 20,
    model_key: str = "B",
    full_scan: bool = False,
    single_results: Optional[List] = None,
    parent_score: Optional[float] = None,
    seed_variant: Optional[object] = None,
    calibrator_key: Optional[str] = None,
    normalize_mode: str = "clip",
) -> List[CmSiRNA]:
    """
    Beam-search multi-modification scan.

    1. Run single-mod scan (or mini-scan) to get top-K single-mod hits
    2. Combine top hits into 2-mod candidates (deduped, order-independent)
    3. Score 2-mod candidates, keep top-K
    4. Repeat for 3-mod, 4-mod, … up to max_mods
    5. Return all scored variants sorted by efficacy

    When seed_variant is provided, it is placed at the top of the initial beam
    so the search explores multi-mod combinations built on top of it.

    Parameters
    ----------
    sense, antisense : parent siRNA
    max_mods         : maximum modifications per variant (2+)
    beam_width       : keep top-K at each expansion step
    model_key        : "A", "B", or "C" for scoring
    full_scan        : if True, use full 1260 single-mod scan; else 40-variant mini-scan
    single_results   : optional pre-computed single-mod results (list of RankedCmSiRNA)
    parent_score     : optional pre-computed parent score (required if single_results provided)
    seed_variant     : optional single-mod variant to seed the beam
    calibrator_key   : calibrator key for _normalize_scores (used when normalize_mode="calibrate")
    normalize_mode   : "clip" — raw clip to [0,100]
                       "calibrate" — isotonic calibration
                       "rescale" — linear rescale to known model max (113.8)

    Returns
    -------
    List[CmSiRNA] sorted best -> worst by efficacy_score
    """
    # Lazy import to avoid circular deps
    from .predictor import predict_modified
    import numpy as np

    # Step 1: get single-mod results (includes parent score)
    if single_results is None:
        single_out = predict_modified(sense, antisense, mode="scan", model_key=model_key, full_scan=full_scan)
        parent_score = single_out.get("parent_score_raw", single_out["parent_score"])
        single_results = single_out["results"]
    elif parent_score is None:
        raise ValueError("parent_score must be provided when single_results is given")

    # Diversify beam: guarantee representation from every mod type
    from collections import defaultdict
    per_mod: dict = defaultdict(list)
    for r in single_results:
        per_mod[r.mod_symbol].append(r)

    for sym, entries in per_mod.items():
        entries.sort(key=lambda r: r.efficacy_score, reverse=True)

    # Round-robin pick: 1st from every mod type, then 2nd, etc. up to beam_width
    diversified: list = []
    max_per = max(len(lst) for lst in per_mod.values())
    for rank in range(max_per):
        for sym in sorted(per_mod.keys()):
            if rank < len(per_mod[sym]):
                diversified.append(per_mod[sym][rank])
            if len(diversified) >= beam_width:
                break
        if len(diversified) >= beam_width:
            break

    # Build initial beam — seed_variant first, then fill with diversified results
    beam_raw: list = []
    if seed_variant is not None:
        beam_raw.append(seed_variant)
    for r in diversified:
        if len(beam_raw) >= beam_width:
            break
        beam_raw.append(r)

    # Convert to CmSiRNA with scored efficacy (re-score with the correct normalize_mode)
    beam: List[CmSiRNA] = []
    for r in beam_raw:
        v = CmSiRNA(
            sense=r.sense,
            antisense=r.antisense,
            mod_symbol=r.mod_symbol,
            mod_position=r.mod_position,
            mod_strand=r.mod_strand,
            parent_sense=sense,
            parent_antisense=antisense,
        )
        v.efficacy_score = r.efficacy_score
        v.delta_score = r.delta_score
        beam.append(v)

    # Helper: predict efficacy for a batch of variants
    def score_variants(variants: List[CmSiRNA]) -> List[CmSiRNA]:
        if not variants:
            return []
        from .features import extract_positional_features_batch
        from .predictor import _get_model, _normalize_scores
        from .biophysics import adjusted_efficacy_score

        s_list = [v.sense for v in variants]
        a_list = [v.antisense for v in variants]
        bs_list = [v.parent_sense for v in variants]
        ba_list = [v.parent_antisense for v in variants]

        X = extract_positional_features_batch(s_list, a_list, bs_list, ba_list)
        model = _get_model("B")
        raw = model.predict(X)
        raw_scores = _normalize_scores(raw, mode="identity")

        out = []
        for v, rs in zip(variants, raw_scores):
            adj_s, penalties, _ = adjusted_efficacy_score(
                float(rs), v.sense, v.antisense,
                v.parent_sense, v.parent_antisense
            )
            v.efficacy_score = round(adj_s, 2)
            v.delta_score = round(adj_s - parent_adjusted, 2)
            v.penalties = penalties
            out.append(v)
        return out

    # Compute parent adjusted score once
    from .biophysics import adjusted_efficacy_score
    parent_adjusted, _, _ = adjusted_efficacy_score(
        parent_score, sense, antisense, sense, antisense
    )

    # Re-score initial beam — use adjusted scores from now on
    beam = score_variants(beam)
    beam.sort(key=lambda x: x.efficacy_score, reverse=True)
    all_scored = list(beam)

    # Beam expansion for multi-mod
    for n_mods in range(2, max_mods + 1):
        candidates = []

        def mod_key(v) -> tuple:
            sym = getattr(v, 'mod_symbol', '')
            pos = getattr(v, 'mod_position', 0)
            strand = getattr(v, 'mod_strand', '')
            positions = getattr(v, 'mod_positions', '')
            return (sym, pos, strand, positions)

        seen = set()
        for v1 in beam:
            for v2 in single_results:
                k1 = mod_key(v1)
                k2 = mod_key(v2)
                pair = tuple(sorted([k1, k2]))
                if pair in seen:
                    continue
                seen.add(pair)

                mod_sense = list(v1.parent_sense)
                mod_antisense = list(v1.parent_antisense)
                mod_symbols = []
                mod_positions = []
                mod_strands = []

                for i in range(len(sense)):
                    if v1.sense[i] != v1.parent_sense[i]:
                        mod_sense[i] = v1.sense[i]
                        mod_symbols.append(v1.sense[i])
                        mod_positions.append(i + 1)
                        mod_strands.append("sense")
                for i in range(len(antisense)):
                    if v1.antisense[i] != v1.parent_antisense[i]:
                        mod_antisense[i] = v1.antisense[i]
                        mod_symbols.append(v1.antisense[i])
                        mod_positions.append(i + 1)
                        mod_strands.append("antisense")

                if v2.mod_strand == "sense":
                    if mod_sense[v2.mod_position - 1] != sense[v2.mod_position - 1]:
                        continue
                    mod_sense[v2.mod_position - 1] = v2.mod_symbol
                else:
                    if mod_antisense[v2.mod_position - 1] != antisense[v2.mod_position - 1]:
                        continue
                    mod_antisense[v2.mod_position - 1] = v2.mod_symbol
                mod_symbols.append(v2.mod_symbol)
                mod_positions.append(v2.mod_position)
                mod_strands.append(v2.mod_strand)

                candidates.append(CmSiRNA(
                    sense="".join(mod_sense),
                    antisense="".join(mod_antisense),
                    mod_symbol="+".join(mod_symbols),
                    mod_position=mod_positions[0],
                    mod_positions=",".join(str(p) for p in mod_positions),
                    mod_strand="+".join(mod_strands),
                    parent_sense=sense,
                    parent_antisense=antisense,
                ))

        scored = score_variants(candidates)
        scored.sort(key=lambda v: v.efficacy_score, reverse=True)
        beam = scored[:beam_width]
        all_scored.extend(scored)

    all_scored.sort(key=lambda v: v.efficacy_score, reverse=True)
    return all_scored
