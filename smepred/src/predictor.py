"""
predictor.py — Unified Machine Learning Prediction Interface

This module acts as the central orchestration layer for the HelixZero-CMS pipeline. 
It ties together the sequence parser, candidate generator, feature extractor, 
LightGBM models, modification engine, and biophysical penalty algorithms.

Workflows:
1. rank_sirnas():
   Takes a raw mRNA/gene sequence, generates all possible unmodified 21-mer siRNA 
   candidates, extracts combinatorial features, and scores them using the baseline 
   LightGBM model (Model A). 

2. predict_modified():
   Takes a specific siRNA candidate and systematically applies chemical modifications 
   (either a single-mod scan or a specific multi-mod configuration). Features are 
   extracted using the positional-aware Model B, and final scores are heavily 
   penalized by the biophysics engine to enforce clinical realism.
"""

import warnings
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple, Union, Dict, Any

import numpy as np
import joblib

# Suppress sklearn feature name warnings when predicting from raw numpy arrays
warnings.filterwarnings('ignore', message='X does not have valid feature names')

from .parser import load_sequence
from .sirna_generator import generate_candidates, generate_dsirna_candidate, SiRNACandidate
from .features import extract_batch_v4, extract_positional_features_batch
from .modification_engine import single_mod_scan, multimod_gen, CmSiRNA, _apply_mod
from .filters import annotate_candidates, toxicity_for_modified
from .biophysics import calculate_adjusted_efficacy

logger = logging.getLogger(__name__)

# ─── Model Paths and Caching ──────────────────────────────────────────────────

MODELS_DIR = Path(__file__).parent.parent / "models"

_MODEL_FILES = {
    "normal": MODELS_DIR / "model_normal.pkl",
    "B":      MODELS_DIR / "model_b.pkl",
}

_CALIBRATOR_FILES = {
    "normal": MODELS_DIR / "calibrator_naked.pkl",
}

_loaded_models: Dict[str, Any] = {}
_loaded_calibrators: Dict[str, Any] = {}


def _get_model(key: str) -> Any:
    """Lazy-loads and caches LightGBM models from disk."""
    if key not in _loaded_models:
        path = _MODEL_FILES.get(key)
        if not path or not path.exists():
            logger.error(f"Model file not found: {path}")
            raise FileNotFoundError(
                f"Model file not found: {path}\n"
                "Run `python models/train_gbm_v3.py` to train and save models first."
            )
        _loaded_models[key] = joblib.load(path)
        logger.info(f"Successfully loaded model: {key}")
    return _loaded_models[key]


def _predict_naked(feature_matrix: np.ndarray) -> np.ndarray:
    """
    Executes inference using the baseline (unmodified) LightGBM model.
    Pads the source one-hot encoding array to match training structure.
    """
    model_bundle = _get_model("normal")
    
    if isinstance(model_bundle, dict):
        model = model_bundle["model"]
        sources = model_bundle.get("sources", [])
        if sources:
            source_onehot = np.zeros((feature_matrix.shape[0], len(sources)), dtype=np.float32)
            # Find the reference human source and set its bit to 1.0
            ref_idx = next((i for i, s in enumerate(sources) if "Hu" in s), 0)
            source_onehot[:, ref_idx] = 1.0
            input_matrix = np.concatenate([feature_matrix, source_onehot], axis=1)
        else:
            input_matrix = feature_matrix
        return model.predict(input_matrix)
        
    return model_bundle.predict(feature_matrix)


def _get_calibrator(key: str) -> Any:
    """Lazy-loads an isotonic calibrator. Returns None if file does not exist."""
    if key not in _loaded_calibrators:
        path = _CALIBRATOR_FILES.get(key)
        if path is not None and path.exists():
            _loaded_calibrators[key] = joblib.load(path)
            logger.info(f"Loaded isotonic calibrator for: {key}")
        else:
            _loaded_calibrators[key] = None
    return _loaded_calibrators[key]


def _normalize_scores(
    raw_predictions: np.ndarray, 
    calibrator_key: Optional[str] = None, 
    mode: str = "clip"
) -> np.ndarray:
    """
    Normalizes raw LightGBM output scores to a strict 0.0 - 100.0 scale.
    """
    if mode == "identity":
        return np.clip(raw_predictions, 0.0, 100.0)
        
    if mode == "rescale":
        # Known max raw output of the cm-siRNA model
        max_raw = 113.8  
        return np.clip((raw_predictions / max_raw) * 100.0, 0.0, 100.0)
        
    if mode == "calibrate" or calibrator_key is not None:
        calibrator = _get_calibrator(calibrator_key)
        if calibrator is not None:
            return np.clip(calibrator.transform(raw_predictions), 0.0, 100.0)
            
    return np.clip(raw_predictions, 0.0, 100.0)


def _get_efficacy_label(score: float) -> str:
    """
    Classifies a numerical efficacy score into human-readable categorical labels.
    """
    if score >= 80.0:
        return "Very High"
    elif score >= 70.0:
        return "High"
    elif score >= 55.0:
        return "Moderate"
    else:
        return "Low"


# ─── Data Transfer Objects ────────────────────────────────────────────────────

@dataclass
class RankedSiRNA:
    """DTO for a ranked, unmodified siRNA candidate."""
    rank: int
    position: int
    sense: str
    antisense: str
    efficacy_score: float
    efficacy_label: str
    toxicity_score: Optional[float] = None
    toxicity_label: str = "Unknown"
    func_ok: bool = True
    func_reason: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return {
            "rank": self.rank,
            "position": self.position,
            "sense": self.sense,
            "antisense": self.antisense,
            "efficacy_score": round(self.efficacy_score, 2),
            "efficacy_label": self.efficacy_label,
            "toxicity_score": self.toxicity_score,
            "toxicity_label": self.toxicity_label,
            "func_ok": self.func_ok,
            "func_reason": self.func_reason,
        }


@dataclass
class RankedCmSiRNA:
    """DTO for a ranked, chemically modified siRNA candidate."""
    rank: int
    sense: str
    antisense: str
    mod_symbol: str
    mod_position: int
    mod_strand: str
    efficacy_score: float
    delta_score: float
    efficacy_label: str
    mod_positions: str = ""
    toxicity_score: Optional[float] = None
    toxicity_label: str = "Unknown"
    toxicity_note: str = ""
    biophysics: Optional[Dict[str, float]] = None

    def to_dict(self) -> Dict[str, Any]:
        result = {
            "rank": self.rank,
            "sense": self.sense,
            "antisense": self.antisense,
            "mod_symbol": self.mod_symbol,
            "mod_position": self.mod_position,
            "mod_strand": self.mod_strand,
            "mod_positions": self.mod_positions or str(self.mod_position),
            "efficacy_score": round(self.efficacy_score, 2),
            "delta_score": round(self.delta_score, 2),
            "efficacy_label": self.efficacy_label,
            "toxicity_score": self.toxicity_score,
            "toxicity_label": self.toxicity_label,
            "toxicity_note": self.toxicity_note,
        }
        if self.biophysics is not None:
            result["biophysics"] = self.biophysics
        return result


# ─── Workflow 1: Unmodified siRNA Ranking ─────────────────────────────────────

def rank_sirnas(
    source: Union[str, Path],
    top_n: Optional[int] = None,
    input_type: str = "gene",
) -> List[RankedSiRNA]:
    """
    Parses an mRNA transcript, generates all combinatorial 21-mer candidates, 
    and ranks them by predicted naked efficacy.
    """
    logger.info("Starting rank_sirnas workflow.")
    sequence = load_sequence(source)

    if input_type == "dsirna":
        candidates = generate_dsirna_candidate(sequence)
    else:
        candidates = generate_candidates(sequence)

    if not candidates:
        logger.warning("No candidates generated.")
        return []

    sense_list = [c.sense for c in candidates]
    antisense_list = [c.antisense for c in candidates]
    
    # Extract structural features for the ML model
    feature_matrix = extract_batch_v4(sense_list, antisense_list)

    # Predict and normalize
    raw_scores = _predict_naked(feature_matrix)
    normalized_scores = _normalize_scores(raw_scores, calibrator_key="normal")

    # Annotate seed toxicity
    annotations = annotate_candidates(sense_list, antisense_list)

    # Rank by score (descending)
    sort_order = np.argsort(normalized_scores)[::-1]
    
    ranked_results = []
    for rank_idx, original_idx in enumerate(sort_order):
        cand = candidates[original_idx]
        score = float(normalized_scores[original_idx])
        annotation = annotations[original_idx]
        
        ranked_results.append(RankedSiRNA(
            rank=rank_idx + 1,
            position=cand.position,
            sense=cand.sense,
            antisense=cand.antisense,
            efficacy_score=score,
            efficacy_label=_get_efficacy_label(score),
            toxicity_score=annotation["toxicity_score"],
            toxicity_label=annotation["toxicity_label"],
            func_ok=annotation["func_ok"],
            func_reason=annotation["func_reason"],
        ))

    if top_n is not None:
        ranked_results = ranked_results[:top_n]

    logger.info(f"Successfully ranked {len(ranked_results)} siRNA candidates.")
    return ranked_results


def rank_by_naked_score(
    source: Union[str, Path],
    top_n: Optional[int] = None,
    input_type: str = "gene",
) -> List[RankedSiRNA]:
    """Alias for rank_sirnas."""
    return rank_sirnas(source, top_n, input_type)


# ─── Workflow 2: Modified siRNA Prediction ────────────────────────────────────

def _perform_mini_scan(
    sense: str, antisense: str
) -> Tuple[List[str], List[str], List[str], List[str]]:
    """
    Generates a rapid, low-compute 40-variant scan using the 4 most 
    effective modifications (E, D, Q, L) on antisense positions 1-10.
    """
    sense_variants, antisense_variants = [], []
    parent_sense, parent_antisense = [], []
    
    for symbol in ('E', 'D', 'Q', 'L'):
        for position in range(1, 11):
            modified_antisense = _apply_mod(antisense, position, symbol)
            
            sense_variants.append(sense)
            antisense_variants.append(modified_antisense)
            parent_sense.append(sense)
            parent_antisense.append(antisense)
            
    return sense_variants, antisense_variants, parent_sense, parent_antisense


def predict_modified(
    sense: str,
    antisense: str,
    mode: str = "scan",
    model_key: str = "B",
    full_scan: bool = False,
    sense_mods: str = "",
    sense_positions: str = "",
    antisense_mods: str = "",
    antisense_positions: str = "",
) -> Dict[str, Any]:
    """
    Predicts the efficacy of chemically modified siRNA variants.
    Applies the biophysical engine to penalize clinically non-viable modifications.
    """
    logger.info(f"Starting predict_modified workflow (mode: {mode}).")
    
    # 1. Establish parent baselines
    parent_v4_matrix = extract_batch_v4([sense], [antisense])
    raw_parent_score = float(_normalize_scores(_predict_naked(parent_v4_matrix), calibrator_key="normal")[0])
    
    parent_b_matrix = extract_positional_features_batch([sense], [antisense], [sense], [antisense])
    model_b = _get_model("B")
    raw_model_b_score = float(_normalize_scores(np.array([model_b.predict(parent_b_matrix)[0]]), mode="identity")[0])

    # 2. Generate variants
    if mode == "scan":
        if full_scan:
            variants = single_mod_scan(sense, antisense)
        else:
            s_vars, a_vars, ps_vars, pa_vars = _perform_mini_scan(sense, antisense)
            variants = []
            for i in range(len(s_vars)):
                variants.append(CmSiRNA(
                    sense=s_vars[i], 
                    antisense=a_vars[i],
                    mod_symbol='E' if i < 10 else 'D' if i < 20 else 'Q' if i < 30 else 'L',
                    mod_position=(i % 10) + 1,
                    mod_strand='antisense',
                    parent_sense=ps_vars[i], 
                    parent_antisense=pa_vars[i],
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
        raise ValueError(f"Invalid mode provided: {mode}")

    if not variants:
        return {"results": [], "parent_score": 0.0, "parent_score_raw": 0.0, "model_b_baseline": 0.0, "naked_baseline": 0.0}

    # 3. Extract features for variants
    s_list = [v.sense for v in variants]
    a_list = [v.antisense for v in variants]
    ps_list = [v.parent_sense for v in variants]
    pa_list = [v.parent_antisense for v in variants]
    
    feature_matrix = extract_positional_features_batch(s_list, a_list, ps_list, pa_list)

    # 4. Predict
    raw_variant_scores = model_b.predict(feature_matrix)
    normalized_scores = _normalize_scores(raw_variant_scores, mode="identity")

    # 5. Apply biophysical constraints and rank
    parent_adjusted_score, _, _ = calculate_adjusted_efficacy(
        raw_model_b_score, sense, antisense, sense, antisense
    )
    raw_parent_adjusted_score, _, _ = calculate_adjusted_efficacy(
        raw_parent_score, sense, antisense, sense, antisense
    )
    
    sort_order = np.argsort(normalized_scores)[::-1]
    ranked_results = []
    
    for rank_idx, original_idx in enumerate(sort_order):
        variant = variants[original_idx]
        score = float(normalized_scores[original_idx])
        
        # Apply physical penalties
        adj_score, penalties, _ = calculate_adjusted_efficacy(
            score, variant.sense, variant.antisense, variant.parent_sense, variant.parent_antisense
        )
        
        viability, tox_label, tox_note = toxicity_for_modified(variant.antisense, variant.parent_antisense)
        
        ranked_results.append(RankedCmSiRNA(
            rank=rank_idx + 1,
            sense=variant.sense,
            antisense=variant.antisense,
            mod_symbol=variant.mod_symbol,
            mod_position=variant.mod_position,
            mod_strand=variant.mod_strand,
            mod_positions=variant.mod_positions,
            efficacy_score=adj_score,
            delta_score=adj_score - parent_adjusted_score,
            efficacy_label=_get_efficacy_label(adj_score),
            toxicity_score=viability,
            toxicity_label=tox_label,
            toxicity_note=tox_note,
            biophysics=penalties,
        ))

    logger.info(f"Successfully evaluated {len(ranked_results)} modified siRNA variants.")
    return {
        "results": ranked_results, 
        "parent_score": round(parent_adjusted_score, 2),
        "parent_score_raw": round(raw_parent_score, 2),
        "model_b_baseline": round(parent_adjusted_score, 2),
        "naked_baseline": round(raw_parent_adjusted_score, 2)
    }
