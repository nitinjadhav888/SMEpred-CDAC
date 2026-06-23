"""
predictor.py — Unified prediction interface.

This module ties together the parser, siRNA generator, feature extractor,
LightGBM models, and modification engine into two high-level workflows:

WORKFLOW 1 — rank_sirnas()
  Input  : mRNA/gene sequence or file path
  Output : list of 21-mer siRNA candidates sorted by predicted efficacy (high → low)
  Steps:
    1. Parse input sequence
    2. Generate all (N-20) 21-mer candidates with sliding window
    3. Extract GBM features for each
    4. Run through the normal siRNA LightGBM model
    5. Score 0-100, classify:
       - Score ≥ 90     : "Very High Efficacy"
       - Score 80-89    : "High Efficacy"
       - Score 70-79    : "Moderate Efficacy"
       - Score < 70     : "Low Efficacy"

WORKFLOW 2 — predict_modified()
  Input  : one siRNA (sense + antisense strings), modification mode
  Output : list of cm-siRNA variants sorted by predicted efficacy
  Steps:
    1. Generate 1260 cm-siRNAs (single-mod scan) OR user-defined cm-siRNA (MultiModGen)
    2. Extract features per chosen model (A, B, or C)
    3. Run through cm-siRNA LightGBM model
    4. Return sorted list with efficacy scores and delta vs parent

Score normalization:
  LightGBM output is already on the 0-100 inhibition scale. We simply clip to
  [0, 100] rather than min-max rescaling per batch.
"""

import warnings
warnings.filterwarnings('ignore', message='X does not have valid feature names')
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union

import numpy as np
import joblib

from .parser import load_sequence
from .sirna_generator import generate_candidates, SiRNACandidate
from .features import extract_batch_v4
from .modification_engine import single_mod_scan, multimod_gen, CmSiRNA
from .filters import annotate_candidates, toxicity_for_modified

# ─── model paths ──────────────────────────────────────────────────────────────

MODELS_DIR = Path(__file__).parent.parent / "models"

_MODEL_FILES = {
    "normal": MODELS_DIR / "model_normal.pkl",
    "B":      MODELS_DIR / "model_b.pkl",
}

_CALIBRATOR_FILES = {
    "normal": MODELS_DIR / "calibrator_naked.pkl",
}

_loaded_models      = {}   # lazy-loaded model cache
_loaded_calibrators = {}   # lazy-loaded calibrator cache


def _get_model(key: str):
    if key not in _loaded_models:
        path = _MODEL_FILES[key]
        if not path.exists():
            raise FileNotFoundError(
                f"Model file not found: {path}\n"
                "Run  python models/train_gbm_v3.py  to train and save models first."
            )
        loaded = joblib.load(path)
        # the naked model is saved as {"model": ..., "sources": [...]} so the inference
        # path can pad the source one-hot with the reference (largest) source's bit set.
        _loaded_models[key] = loaded
    return _loaded_models[key]


def _predict_naked(X_seq: np.ndarray) -> np.ndarray:
    """Run the naked model. Pads the source one-hot with the reference source bit set."""
    bundle = _get_model("normal")
    if isinstance(bundle, dict):
        model = bundle["model"]
        sources = bundle.get("sources", [])
        if sources:
            src_onehot = np.zeros((X_seq.shape[0], len(sources)), dtype=np.float32)
            ref = next((i for i, s in enumerate(sources) if "Hu" in s), 0)
            src_onehot[:, ref] = 1.0
            X = np.concatenate([X_seq, src_onehot], axis=1)
        else:
            X = X_seq
        return model.predict(X)
    return bundle.predict(X_seq)


# ─── score normalization ──────────────────────────────────────────────────────

def _get_calibrator(key: str):
    """Lazy-load an isotonic calibrator. Returns None if file doesn't exist."""
    if key not in _loaded_calibrators:
        path = _CALIBRATOR_FILES.get(key)
        if path is not None and path.exists():
            _loaded_calibrators[key] = joblib.load(path)
        else:
            _loaded_calibrators[key] = None
    return _loaded_calibrators[key]


def _normalize_scores(raw: np.ndarray, calibrator_key: str = None, mode: str = "clip") -> np.ndarray:
    """
    Convert raw model outputs to 0–100 scores.

    Parameters
    ----------
    raw             : raw model predictions
    calibrator_key  : isotonic calibrator key ("cm", "normal"), used when mode="calibrate"
    mode            : "clip"      — clip raw to [0, 100] (default)
                      "calibrate" — isotonic calibration then clip
                      "rescale"   — linear rescale: raw / 113.8 * 100, clipped to [0, 100]
                      "identity"  — clip to [0, 100] with no scaling (model already predicts 0-100)
    """
    if mode == "identity":
        return np.clip(raw, 0.0, 100.0)
    if mode == "rescale":
        RAW_MAX = 113.8  # known max raw output of the cm-siRNA model
        return np.clip(raw / RAW_MAX * 100.0, 0.0, 100.0)
    if mode == "calibrate" or calibrator_key is not None:
        cal = _get_calibrator(calibrator_key)
        if cal is not None:
            return np.clip(cal.transform(raw), 0.0, 100.0)
    return np.clip(raw, 0.0, 100.0)


def _efficacy_label(score: float) -> str:
    # Training-data percentiles: P50=48, P75=72, P84=80, P94=90.
    # Thresholds are calibrated to the prediction distribution (compressed vs labels).
    if score >= 80:
        return "Very High"
    elif score >= 70:
        return "High"
    elif score >= 55:
        return "Moderate"
    else:
        return "Low"


# ─── result containers ────────────────────────────────────────────────────────

@dataclass
class RankedSiRNA:
    rank:          int
    position:      int      # 0-based start position in mRNA
    sense:         str
    antisense:     str
    efficacy_score: float   # 0–100
    efficacy_label: str     # Very High / High / Moderate / Low
    # Safety / functionality annotations
    toxicity_score: Optional[float] = None    # predicted cell viability % (None = unknown seed)
    toxicity_label: str = "Unknown"            # Safe / Caution / Toxic / Unknown
    func_ok:        bool = True
    func_reason:    str = ""

    def to_dict(self) -> dict:
        return {
            "rank":              self.rank,
            "position":          self.position,
            "sense":             self.sense,
            "antisense":         self.antisense,
            "efficacy_score":    round(self.efficacy_score, 2),
            "efficacy_label":    self.efficacy_label,
            "toxicity_score":    self.toxicity_score,
            "toxicity_label":    self.toxicity_label,
            "func_ok":           self.func_ok,
            "func_reason":       self.func_reason,
        }


@dataclass
class RankedCmSiRNA:
    rank:            int
    sense:           str
    antisense:       str
    mod_symbol:      str
    mod_position:    int
    mod_strand:      str
    efficacy_score:  float   # 0–100 (biophysically adjusted)
    delta_score:     float   # adjusted_score − parent_adjusted_score
    efficacy_label:  str
    mod_positions:   str = ""    # all positions as comma-separated for multi-mod (e.g. "4,6")
    # Seed toxicity (canonical-base lookup + modification-aware mitigation flag)
    toxicity_score:  Optional[float] = None
    toxicity_label:  str = "Unknown"      # Safe / Caution / Toxic / Mitigated / Unknown
    toxicity_note:   str = ""             # tooltip explaining a Mitigated flag
    # Biophysical penalty breakdown (debug / tooltip)
    biophysics:      Optional[dict] = None

    def to_dict(self) -> dict:
        d = {
            "rank":            self.rank,
            "sense":           self.sense,
            "antisense":       self.antisense,
            "mod_symbol":      self.mod_symbol,
            "mod_position":    self.mod_position,
            "mod_strand":      self.mod_strand,
            "mod_positions":   self.mod_positions or str(self.mod_position),
            "efficacy_score":  round(self.efficacy_score, 2),
            "delta_score":     round(self.delta_score, 2),
            "efficacy_label":  self.efficacy_label,
            "toxicity_score":  self.toxicity_score,
            "toxicity_label":  self.toxicity_label,
            "toxicity_note":   self.toxicity_note,
        }
        if self.biophysics is not None:
            d["biophysics"] = self.biophysics
        return d


# ─── Workflow 1: rank unmodified siRNA candidates ─────────────────────────────

def rank_sirnas(
    source: Union[str, Path],
    top_n: Optional[int] = None,
    input_type: str = "gene",
) -> List[RankedSiRNA]:
    """
    From an mRNA/gene input, generate and rank all 21-mer siRNA candidates.

    Parameters
    ----------
    source : str or Path
        mRNA sequence, FASTA file path, or inline FASTA text.
    top_n  : optional int
        If set, return only the top N candidates.
    input_type : str
        "gene" (default) → sliding window across transcript.
        "dsirna"        → single 21-mer via Dicer cleavage rule.

    Returns
    -------
    List[RankedSiRNA] sorted best → worst by efficacy score.
    """
    # Step 1: parse
    seq = load_sequence(source)

    # Step 2: generate candidates (gene: sliding window; dsirna: Dicer rule)
    from .sirna_generator import generate_candidates, generate_dsirna_candidate
    if input_type == "dsirna":
        candidates: List[SiRNACandidate] = generate_dsirna_candidate(seq)
    else:
        candidates: List[SiRNACandidate] = generate_candidates(seq)

    # Step 3: extract features (V4 for naked model)
    sense_list     = [c.sense     for c in candidates]
    antisense_list = [c.antisense for c in candidates]
    X = extract_batch_v4(sense_list, antisense_list)

    # Step 4: predict (naked model is wrapped to handle the source one-hot)
    raw_scores = _predict_naked(X)
    scores = _normalize_scores(raw_scores, calibrator_key="normal")

    # Step 4b: safety + functionality annotations
    annotations = annotate_candidates(sense_list, antisense_list)

    # Step 5: rank
    order = np.argsort(scores)[::-1]  # highest first
    results = []
    for rank_i, idx in enumerate(order):
        c = candidates[idx]
        s = float(scores[idx])
        a = annotations[idx]
        results.append(RankedSiRNA(
            rank=rank_i + 1,
            position=c.position,
            sense=c.sense,
            antisense=c.antisense,
            efficacy_score=s,
            efficacy_label=_efficacy_label(s),
            toxicity_score=a["toxicity_score"],
            toxicity_label=a["toxicity_label"],
            func_ok=a["func_ok"],
            func_reason=a["func_reason"],
        ))

    if top_n is not None:
        results = results[:top_n]

    return results


# ─── Workflow 1b: rank by modification potential ──────────────────────────────

def _mini_mod_scan(sense: str, antisense: str) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Quick proxy for modification potential: test the 4 most effective modification
    symbols (E, D, Q, L) on antisense positions 1-10 — 40 variants per candidate.

    These 4 symbols span different chemistries and consistently rank highest
    in max modified score across diverse siRNA sequences.

    Returns (mod_senses, mod_antisenses, base_senses, base_antisenses).
    """
    from .modification_engine import _apply_mod
    s_list, a_list, bs_list, ba_list = [], [], [], []
    for sym in ('E', 'D', 'Q', 'L'):
        for pos in range(1, 11):
            mod_a = _apply_mod(antisense, pos, sym)
            s_list.append(sense)
            a_list.append(mod_a)
            bs_list.append(sense)
            ba_list.append(antisense)
    return s_list, a_list, bs_list, ba_list


def rank_by_naked_score(
    source: Union[str, Path],
    top_n: Optional[int] = None,
    input_type: str = "gene",
) -> List[RankedSiRNA]:
    """
    Rank siRNA candidates by their *naked (unmodified) silencing score*.

    This is the baseline efficacy of each siRNA backbone before any chemical
    modifications are applied. The naked model (PCC=0.55) scores each candidate
    using sequence-composition features only.

    Parameters
    ----------
    source : str or Path
        mRNA sequence, FASTA file path, or inline FASTA text.
    top_n  : optional int
        Number of results to return (default: all).
    input_type : str
        "gene" (default) → sliding window across transcript.
        "dsirna"        → single 21-mer via Dicer cleavage rule.

    Returns
    -------
    List[RankedSiRNA] sorted by naked score (best → worst).
    """
    seq = load_sequence(source)
    from .sirna_generator import generate_candidates, generate_dsirna_candidate
    if input_type == "dsirna":
        candidates: List[SiRNACandidate] = generate_dsirna_candidate(seq)
    else:
        candidates: List[SiRNACandidate] = generate_candidates(seq)

    if len(candidates) == 0:
        return []

    sense_list = [c.sense for c in candidates]
    antisense_list = [c.antisense for c in candidates]
    X_seq = extract_batch_v4(sense_list, antisense_list)
    raw = _predict_naked(X_seq)
    scores = _normalize_scores(raw, calibrator_key="normal")
    order = np.argsort(scores)[::-1]

    annotations = annotate_candidates(sense_list, antisense_list)

    results = []
    for rank_i, pos in enumerate(order):
        c = candidates[pos]
        s = float(scores[pos])
        a = annotations[pos]
        results.append(RankedSiRNA(
            rank=rank_i + 1,
            position=c.position,
            sense=c.sense,
            antisense=c.antisense,
            efficacy_score=s,
            efficacy_label=_efficacy_label(s),
            toxicity_score=a["toxicity_score"],
            toxicity_label=a["toxicity_label"],
            func_ok=a["func_ok"],
            func_reason=a["func_reason"],
        ))

    if top_n is not None:
        results = results[:top_n]

    return results


# ─── Workflow 2: predict modified siRNA efficacy ──────────────────────────────

def predict_modified(
    sense: str,
    antisense: str,
    mode: str = "scan",
    model_key: str = "B",
    full_scan: bool = False,
    # MultiModGen parameters (used when mode="multimod")
    sense_mods: str = "",
    sense_positions: str = "",
    antisense_mods: str = "",
    antisense_positions: str = "",
) -> List[RankedCmSiRNA]:
    """
    Predict efficacy of chemically modified variants of a given siRNA.
    Unified model (Model B / HelixZero) is always used regardless of model_key.

    Parameters
    ----------
    sense / antisense : the 21-nt parent siRNA strands.
    mode              : "scan"     → generate single-mod variants
                        "multimod" → apply user-specified multiple modifications
    model_key         : kept for backward compatibility, always uses unified B
    full_scan         : True  → generate all 1260 single-mod variants (30 mods × 21 pos × 2 strands)
                        False → generate 40-variant mini-scan (E/D/Q/L on antisense pos 1-10)
    sense_mods / sense_positions            : MultiModGen input for sense strand
    antisense_mods / antisense_positions    : MultiModGen input for antisense strand

    Returns
    -------
    List[RankedCmSiRNA] sorted best → worst.
    """
    # Step 1: get naked model parent score (matches what Rank tab shows)
    X_parent_v4 = extract_batch_v4([sense], [antisense])
    raw_naked = _predict_naked(X_parent_v4)
    parent_score = float(_normalize_scores(raw_naked, calibrator_key="normal")[0])
    # Also compute Model B baseline (different feature space — 1242 vs 614 dims)
    # so the UI can show the recalibrated baseline explicitly and avoid confusion
    from .features import extract_positional_features_batch
    X_parent_b = extract_positional_features_batch([sense], [antisense], [sense], [antisense])
    model_b = _get_model("B")
    raw_b = float(model_b.predict(X_parent_b)[0])
    model_b_parent_raw = float(_normalize_scores(np.array([raw_b]), mode="identity")[0])

    # Step 2: generate cm-siRNA variants
    if mode == "scan":
        if full_scan:
            variants: List[CmSiRNA] = single_mod_scan(sense, antisense)
        else:
            s_list, a_list, bs_list, ba_list = _mini_mod_scan(sense, antisense)
            from .modification_engine import CmSiRNA
            variants = []
            for i in range(len(s_list)):
                variants.append(CmSiRNA(
                    sense=s_list[i], antisense=a_list[i],
                    mod_symbol='E' if i < 10 else 'D' if i < 20 else 'Q' if i < 30 else 'L',
                    mod_position=(i % 10) + 1,
                    mod_strand='antisense',
                    parent_sense=bs_list[i], parent_antisense=ba_list[i],
                ))
    elif mode == "multimod":
        variants = [multimod_gen(
            sense, antisense,
            sense_mods=sense_mods,
            sense_positions=sense_positions,
            antisense_mods=antisense_mods,
            antisense_positions=antisense_positions,
        )]
    else:
        raise ValueError(f"mode must be 'scan' or 'multimod', got '{mode}'")

    if not variants:
        return []

    # Step 3: extract positional features (unified B path)
    s_list = [v.sense     for v in variants]
    a_list = [v.antisense for v in variants]
    bs_list = [v.parent_sense     for v in variants]
    ba_list = [v.parent_antisense for v in variants]
    from .features import extract_positional_features_batch
    X = extract_positional_features_batch(s_list, a_list, bs_list, ba_list)

    # Step 4: predict (unified model)
    model = _get_model("B")
    raw = model.predict(X)
    scores = _normalize_scores(raw, mode="identity")

    # Step 5: apply biophysical penalties + rank
    from .biophysics import adjusted_efficacy_score
    # Parent adjusted using Model B baseline (same model as variants) for fair delta
    parent_adjusted, parent_penalties, _ = adjusted_efficacy_score(
        model_b_parent_raw, sense, antisense, sense, antisense
    )
    # Also compute raw naked parent for the frontend to display both baselines
    raw_parent_adjusted, _, _ = adjusted_efficacy_score(
        parent_score, sense, antisense, sense, antisense
    )
    order = np.argsort(scores)[::-1]
    results = []
    for rank_i, idx in enumerate(order):
        v = variants[idx]
        raw_s = float(scores[idx])
        # Biophysical penalties adjust the main efficacy score
        adj_s, penalties, total_p = adjusted_efficacy_score(
            raw_s, v.sense, v.antisense, v.parent_sense, v.parent_antisense
        )
        # Modify score label based on adjusted score
        lab = _efficacy_label(adj_s)
        # Modification-aware seed toxicity
        viab, tox_label, tox_note = toxicity_for_modified(v.antisense, v.parent_antisense)
        results.append(RankedCmSiRNA(
            rank=rank_i + 1,
            sense=v.sense,
            antisense=v.antisense,
            mod_symbol=v.mod_symbol,
            mod_position=v.mod_position,
            mod_strand=v.mod_strand,
            mod_positions=v.mod_positions,
            efficacy_score=round(adj_s, 2),
            delta_score=round(adj_s - parent_adjusted, 2),
            efficacy_label=lab,
            toxicity_score=None if viab is None else round(viab, 1),
            toxicity_label=tox_label,
            toxicity_note=tox_note,
            biophysics=penalties,
        ))

    return {"results": results, "parent_score": round(parent_adjusted, 2),
            "parent_score_raw": round(parent_score, 2),
            "model_b_baseline": round(parent_adjusted, 2),
            "naked_baseline": round(raw_parent_adjusted, 2)}
