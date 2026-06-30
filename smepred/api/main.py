"""
api/main.py — FastAPI REST API for HelixZero-CMS

This module serves as the primary gateway between the frontend UI and the 
HelixZero computational backend. It defines standard endpoints for siRNA 
ranking, chemical modification scanning, and transcriptomic safety validation.

Endpoints:
    POST /rank              : Rank unmodified siRNA candidates from a gene sequence.
    POST /rank/upload       : Same as /rank, but processes a raw FASTA file upload.
    POST /single-mod        : Generate 1,260 single-modification variants for a candidate.
    POST /multi-mod         : Evaluate a specific custom multi-modified cm-siRNA.
    POST /multi-mod-scan    : Combinatorial beam search for optimal multi-mod stacking.
    POST /offtarget-scan    : Run biological safety heuristics against human transcriptome.
    POST /generate-certificate : Generate a Markdown clinical safety dossier.
    GET  /modifications     : Retrieve supported chemical modification nomenclature.

Start Server:
    uvicorn api.main:app --reload --port 8000
"""

import logging
import json
from pathlib import Path
from typing import Optional, List, Dict, Any

from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field, ConfigDict

# Local internal imports
from src.predictor import (
    rank_by_naked_score, 
    predict_modified, 
    _get_efficacy_label,
    _predict_naked, 
    _normalize_scores, 
    _get_model
)
from src.biophysics import adjusted_efficacy_score
from src.filters import toxicity_score, toxicity_label, seed_of_antisense
from src.offtarget import get_offtarget_engine
from src.features import extract_batch_v4, extract_positional_features_batch
from src.modification_engine import multi_mod_scan

# Configure module-level logger
logging.basicConfig(
    level=logging.INFO, 
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
ROOT_DIR = Path(__file__).parent.parent
APP_HTML = ROOT_DIR / "app.html"


# ─── App Initialization ───────────────────────────────────────────────────────

app = FastAPI(
    title="HelixZero-CMS API",
    description=(
        "Production-grade REST API for Machine Learning-driven siRNA discovery, "
        "chemical modification optimization, and transcriptome-wide safety validation."
    ),
    version="2.1.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ─── Pydantic Data Models ─────────────────────────────────────────────────────

class RankRequest(BaseModel):
    sequence: str = Field(..., description="Target gene sequence (raw text or FASTA)")
    top_n: int = Field(20, ge=0, description="Limit results to Top-N (0 = return all)")
    input_type: str = Field("gene", description="Mode: 'gene' (sliding window) or 'dsirna' (Dicer)")

class SingleModRequest(BaseModel):
    sense: str = Field(..., description="21-nt sense strand")
    antisense: str = Field(..., description="21-nt antisense strand")
    model: str = Field("B", description="Model version key")
    top_n: int = Field(50, ge=0, description="Limit returned variants")
    full_scan: bool = Field(False, description="True=1260 variants, False=40-variant targeted scan")

class MultiModRequest(BaseModel):
    sense: str = Field(..., description="21-nt sense strand")
    antisense: str = Field(..., description="21-nt antisense strand")
    sense_mods: str = Field("", description="Modification symbols for sense strand (e.g. 'F,,M')")
    sense_positions: str = Field("", description="Positions for sense mods (e.g. '2,5,,10')")
    antisense_mods: str = Field("", description="Modification symbols for antisense strand")
    antisense_positions: str = Field("", description="Positions for antisense mods")
    model: str = Field("B", description="Model version key")

class MultiModScanRequest(BaseModel):
    sense: str
    antisense: str
    model: str = "B"
    max_mods: int = Field(2, ge=2, le=21)
    beam_width: int = Field(20, ge=5, le=50)
    full_scan: bool = False

class MultiModFromSingleRequest(BaseModel):
    sense: str
    antisense: str
    model: str = "B"
    max_mods: int = Field(5, ge=2, le=21)
    beam_width: int = Field(20, ge=5, le=100)
    full_scan: bool = True
    single_results: Optional[List[Dict[str, Any]]] = None
    parent_score: Optional[float] = None
    seed_variant: Optional[Dict[str, Any]] = None
    calibrator_key: Optional[str] = None
    normalize_mode: str = "rescale"

class OffTargetRequest(BaseModel):
    sense: str = Field(..., description="21-nt sense strand")
    antisense: str = Field(..., description="21-nt antisense strand")
    antisense_mods: str = Field("", description="Modification mask for antisense strand")


# ─── API Endpoints ────────────────────────────────────────────────────────────

@app.get("/")
def serve_frontend():
    """Serves the primary Single-Page Application (SPA) HTML."""
    return FileResponse(APP_HTML, headers={"Cache-Control": "no-cache, no-store, must-revalidate"})


@app.get("/app.html")
def serve_app_html():
    return serve_frontend()


@app.post("/offtarget-scan")
def offtarget_scan_endpoint(req: OffTargetRequest):
    """
    Executes a biological safety heuristic scan against the human transcriptome.
    """
    try:
        engine = get_offtarget_engine()
        result = engine.validate_safety(req.sense, req.antisense, req.antisense_mods)
        return result
    except Exception as e:
        logger.error(f"Transcriptome safety scan failed: {e}")
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")


@app.post("/generate-certificate")
def generate_certificate_endpoint(req: OffTargetRequest):
    """
    Generates a Markdown Clinical Safety Dossier based on the transcriptome scan.
    """
    try:
        engine = get_offtarget_engine()
        report = engine.validate_safety(req.sense, req.antisense, req.antisense_mods)
        cert_path = engine.generate_markdown_certificate(
            report, req.sense, req.antisense, req.antisense_mods
        )
        return {"success": True, "certificate_path": cert_path}
    except Exception as e:
        logger.error(f"Certificate generation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/rank")
def rank_endpoint(req: RankRequest):
    """
    Scores and ranks un-modified (naked) siRNA candidates utilizing Model A.
    """
    try:
        limit = req.top_n if req.top_n > 0 else None
        results = rank_by_naked_score(req.sequence, top_n=limit, input_type=req.input_type)
        return {
            "total_candidates": len(results),
            "input_type": req.input_type,
            "results": [r.to_dict() for r in results],
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail="Model file not found. Ensure models are compiled.")
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/rank/upload")
async def rank_upload_endpoint(file: UploadFile = File(...), top_n: int = 20):
    """
    Scores and ranks candidates ingested directly from a FASTA file upload.
    """
    try:
        content = (await file.read()).decode("utf-8")
        limit = top_n if top_n > 0 else None
        results = rank_by_naked_score(content, top_n=limit)
        return {
            "filename": file.filename,
            "total_candidates": len(results),
            "results": [r.to_dict() for r in results],
        }
    except Exception as e:
        raise HTTPException(status_code=422, detail=f"File processing failed: {str(e)}")


@app.post("/single-mod")
def single_mod_endpoint(req: SingleModRequest):
    """
    Exhaustively scans and evaluates all 1,260 single-point chemical modifications 
    across both strands of a parent siRNA candidate utilizing Model B.
    """
    try:
        output = predict_modified(
            req.sense, req.antisense,
            mode="scan",
            full_scan=req.full_scan,
            model_key=req.model
        )
        
        results = output["results"]
        parent_score = output["parent_score"]
        
        top_results = results[:req.top_n] if req.top_n > 0 else results
        
        # Calculate parent baseline toxicity
        parent_viability = toxicity_score(req.antisense)
        parent_seed = seed_of_antisense(req.antisense)
        
        # Calculate parent baseline transcriptome safety
        engine = get_offtarget_engine()
        parent_safety = engine.validate_safety(req.sense, req.antisense, "")

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": parent_score,
            "model_b_baseline": output.get("model_b_baseline", parent_score),
            "naked_baseline": output.get("naked_baseline", parent_score),
            "model": req.model,
            "total_variants": len(results),
            "full_scan": req.full_scan,
            "parent_toxicity": {
                "seed": parent_seed,
                "viability": round(parent_viability, 1) if parent_viability is not None else None,
                "label": toxicity_label(parent_viability),
            },
            "parent_safety": parent_safety,
            "results": [r.to_dict() for r in top_results],
        }
    except Exception as e:
        logger.error(f"Single-mod scan failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/multi-mod")
def multi_mod_endpoint(req: MultiModRequest):
    """
    Evaluates a highly specific, user-defined combinatorial modification pattern.
    """
    try:
        output = predict_modified(
            req.sense, req.antisense,
            mode="multimod",
            sense_mods=req.sense_mods,
            sense_positions=req.sense_positions,
            antisense_mods=req.antisense_mods,
            antisense_positions=req.antisense_positions,
        )
        results = output["results"]
        if not results:
            raise HTTPException(status_code=500, detail="Modification engine yielded no valid variants.")
            
        variant = results[0]
        variant_dict = variant.to_dict()
        
        # Apply safety scan and penalize efficacy if hazardous
        engine = get_offtarget_engine()
        safety_report = engine.validate_safety(
            variant.sense, req.antisense, variant.antisense, variant.sense
        )
        
        if safety_report["overallSafetyScore"] < 100:
            off_target_penalty_weight = (100 - safety_report["overallSafetyScore"]) * 0.2
            penalties = variant_dict.get("penalties", {})
            penalties["offtarget"] = round(off_target_penalty_weight, 1)
            
            variant_dict["penalties"] = penalties
            variant_dict["total_penalty"] = variant_dict.get("total_penalty", 0.0) + round(off_target_penalty_weight, 1)
            
            current_efficacy = variant_dict.get("efficacy_score", 0.0)
            adjusted_score = max(0.0, current_efficacy - off_target_penalty_weight)
            
            if not safety_report["isSafe"]:
                adjusted_score = 0.0
                
            variant_dict["efficacy_score"] = round(adjusted_score, 1)
            variant_dict["efficacy_label"] = _get_efficacy_label(variant_dict["efficacy_score"])

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": output["parent_score"],
            "model_b_baseline": output.get("model_b_baseline", output["parent_score"]),
            "naked_baseline": output.get("naked_baseline", output["parent_score"]),
            "model": req.model,
            "safety_report": safety_report,
            "result": variant_dict,
        }
    except Exception as e:
        logger.error(f"Multi-mod evaluation failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/multi-mod-scan")
def multi_mod_scan_endpoint(req: MultiModScanRequest):
    """
    Executes an autonomous beam search to stack synergistic modifications, 
    generating a highly optimized multi-mod library.
    """
    try:
        # Establish accurate baselines
        parent_features = extract_batch_v4([req.sense], [req.antisense])
        raw_naked = _predict_naked(parent_features)
        raw_parent_score = float(_normalize_scores(raw_naked, calibrator_key="normal")[0])
        naked_baseline_adj, _, _ = adjusted_efficacy_score(
            raw_parent_score, req.sense, req.antisense, req.sense, req.antisense
        )
        
        parent_features_b = extract_positional_features_batch(
            [req.sense], [req.antisense], [req.sense], [req.antisense]
        )
        model_b = _get_model(req.model)
        raw_b_score = float(model_b.predict(parent_features_b)[0])
        model_b_adj, _, _ = adjusted_efficacy_score(
            raw_b_score, req.sense, req.antisense, req.sense, req.antisense
        )

        variants = multi_mod_scan(
            req.sense, req.antisense,
            max_mods=req.max_mods,
            beam_width=req.beam_width,
            model_key=req.model,
            full_scan=req.full_scan,
        )

        formatted_results = []
        for idx, variant in enumerate(variants):
            penalties = getattr(variant, 'penalties', None) or {}
            total_penalty = sum(penalties.values())
            raw_efficacy = round(variant.efficacy_score + 0.70 * total_penalty, 2)
            
            formatted_results.append({
                "rank": idx + 1,
                "sense": variant.sense,
                "antisense": variant.antisense,
                "mod_symbol": variant.mod_symbol,
                "mod_position": variant.mod_position,
                "mod_strand": variant.mod_strand,
                "mod_positions": variant.mod_positions or str(variant.mod_position),
                "raw_efficacy_score": raw_efficacy,
                "efficacy_score": round(variant.efficacy_score, 2),
                "total_penalty": round(total_penalty, 1),
                "delta_score": round(variant.delta_score, 2),
                "efficacy_label": _get_efficacy_label(variant.efficacy_score),
                "penalties": {k: round(v, 1) for k, v in penalties.items()},
            })

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": round(model_b_adj, 2),
            "model_b_baseline": round(model_b_adj, 2),
            "naked_baseline": round(naked_baseline_adj, 2),
            "model": req.model,
            "total_variants": len(formatted_results),
            "results": formatted_results,
        }
    except Exception as e:
        logger.error(f"Multi-mod beam search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/multi-mod-from-single")
def multi_mod_from_single_endpoint(req: MultiModFromSingleRequest):
    """
    Advanced multi-mod search initialized from pre-computed single-modification results.
    Integrates safety scans inside the beam evaluation phase.
    """
    try:
        class ProxyVariant:
            """Translates raw JSON dictionaries back to Python objects for the engine."""
            def __init__(self, data: Dict[str, Any]):
                self.mod_symbol = data.get("mod_symbol") or data.get("modification", "")
                self.mod_position = data.get("mod_position") or data.get("position", 0)
                self.mod_strand = data.get("mod_strand") or data.get("strand", "")
                self.mod_positions = data.get("mod_positions") or str(self.mod_position)
                self.efficacy_score = data.get("efficacy_score") or data.get("score", 0.0)
                self.sense = data.get("sense", "")
                self.antisense = data.get("antisense", "")
                self.delta_score = data.get("delta_score", 0.0)

        single_proxies = [ProxyVariant(sr) for sr in req.single_results] if req.single_results else None
        seed_proxy = ProxyVariant(req.seed_variant) if req.seed_variant else None

        # Reconstruct Baseline
        if req.parent_score is None:
            features = extract_batch_v4([req.sense], [req.antisense])
            raw = float(_normalize_scores(_predict_naked(features), mode=req.normalize_mode)[0])
            naked_adj, _, _ = adjusted_efficacy_score(raw, req.sense, req.antisense, req.sense, req.antisense)
            parent_baseline = round(naked_adj, 2)
        else:
            parent_baseline = req.parent_score

        features_b = extract_positional_features_batch([req.sense], [req.antisense], [req.sense], [req.antisense])
        mb = _get_model("B")
        raw_b = float(mb.predict(features_b)[0])
        mb_adj, _, _ = adjusted_efficacy_score(raw_b, req.sense, req.antisense, req.sense, req.antisense)
        model_b_baseline = round(mb_adj, 2)

        variants = multi_mod_scan(
            req.sense, req.antisense,
            max_mods=req.max_mods,
            beam_width=req.beam_width,
            model_key=req.model,
            full_scan=req.full_scan,
            single_results=single_proxies,
            parent_score=parent_baseline,
            seed_variant=seed_proxy,
            calibrator_key=req.calibrator_key,
            normalize_mode=req.normalize_mode,
        )

        engine = get_offtarget_engine()
        formatted_results = []
        
        for idx, var in enumerate(variants):
            penalties = getattr(var, 'penalties', None) or {}
            
            # Integrate Transcriptome Safety Check
            safety = engine.validate_safety(var.sense, req.antisense, var.antisense, var.sense)
            if safety["overallSafetyScore"] < 100:
                offtarget_pen = (100 - safety["overallSafetyScore"]) * 0.2
                penalties["offtarget"] = round(offtarget_pen, 1)
                
            total_penalty = sum(penalties.values())
            
            # Recalculate raw score
            old_penalty = sum((getattr(var, 'penalties', None) or {}).values())
            if "offtarget" in penalties:
                old_penalty -= penalties["offtarget"]
                
            raw_score = round(var.efficacy_score + 0.70 * old_penalty, 2)
            adjusted_score = round(raw_score - 0.70 * total_penalty, 2)
            
            if not safety["isSafe"]:
                adjusted_score = 0.0
            adjusted_score = max(0.0, adjusted_score)

            formatted_results.append({
                "rank": 0,
                "sense": var.sense,
                "antisense": var.antisense,
                "mod_symbol": var.mod_symbol,
                "mod_position": var.mod_position,
                "mod_strand": var.mod_strand,
                "mod_positions": var.mod_positions or str(var.mod_position),
                "raw_efficacy_score": raw_score,
                "efficacy_score": adjusted_score,
                "total_penalty": round(total_penalty, 1),
                "delta_score": round(adjusted_score - model_b_baseline, 2),
                "efficacy_label": _get_efficacy_label(adjusted_score),
                "penalties": {k: round(v, 1) for k, v in penalties.items()},
                "offtarget_score": safety["overallSafetyScore"],
                "offtarget_status": safety["status"],
            })

        # Re-sort due to new safety penalties
        formatted_results.sort(key=lambda x: x["efficacy_score"], reverse=True)
        for idx, res in enumerate(formatted_results):
            res["rank"] = idx + 1

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": parent_baseline,
            "model_b_baseline": model_b_baseline,
            "naked_baseline": parent_baseline, # Simplified
            "model": req.model,
            "total_variants": len(formatted_results),
            "parent_safety": engine.validate_safety(req.sense, req.antisense, ""),
            "results": formatted_results,
        }
    except Exception as e:
        logger.error(f"Seeded multi-mod search failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/modifications")
def get_supported_modifications():
    """
    Returns the comprehensive dictionary of 30 supported chemical modifications.
    """
    try:
        mod_file = ROOT_DIR / "data" / "modification_codes.json"
        with mod_file.open("r", encoding="utf-8") as file:
            data = json.load(file)
            
        return {
            "canonical": [m for m in data["modifications"] if m["type"] == "canonical"],
            "modifications": [m for m in data["modifications"] if m["type"] != "canonical"],
        }
    except Exception as e:
        logger.error(f"Failed to load modification taxonomy: {e}")
        raise HTTPException(status_code=500, detail="Modification taxonomy file missing or corrupted.")
