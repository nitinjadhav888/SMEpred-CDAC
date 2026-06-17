"""
api/main.py — FastAPI REST API for HelixZero-CMS.

Exposes three endpoints:

  POST /rank
    Body: { "sequence": "AUGCAUG..." }  or  upload a FASTA file
    Returns: ranked list of siRNA candidates with efficacy scores

  POST /single-mod
    Body: { "sense": "...", "antisense": "...", "model": "A" }
    Returns: top 1260 cm-siRNA variants sorted by efficacy

  POST /multi-mod
    Body: { "sense": "...", "antisense": "...", "sense_mods": "F,,M",
            "sense_positions": "2,5,,10", "antisense_mods": "", "antisense_positions": "",
            "model": "A" }
    Returns: predicted efficacy of the custom cm-siRNA

  GET /modifications
    Returns: list of all 30 supported chemical modification symbols + names

Start the server:
  uvicorn api.main:app --reload --port 8000

Then visit http://localhost:8000/docs for the auto-generated Swagger UI.
"""

import sys
from pathlib import Path
from typing import Optional

from fastapi import FastAPI, UploadFile, File, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
import json

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.predictor import rank_sirnas, rank_by_naked_score, predict_modified, _efficacy_label

ROOT_DIR = Path(__file__).parent.parent
APP_HTML = ROOT_DIR / "app.html"

# ─── app setup ────────────────────────────────────────────────────────────────

app = FastAPI(
    title="HelixZero-CMS API",
    description=(
        "Rank siRNA candidates by predicted efficacy, scan 1,260 chemical modifications, "
        "and flag seed-based toxicity. Inspired by the SMEpred approach "
        "(Dar et al., RNA Biology, 2016)."
    ),
    version="2.0.0",
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)


# ─── request / response models ────────────────────────────────────────────────

class RankRequest(BaseModel):
    sequence: str = Field(..., description="mRNA or gene sequence (plain or FASTA format)")
    top_n: int    = Field(20, description="Number of top candidates to return (0 = all)")


class SingleModRequest(BaseModel):
    sense:     str = Field(..., description="21-nt sense strand")
    antisense: str = Field(..., description="21-nt antisense strand")
    model:     str = Field("A", description="Model: A (default), B, or C")
    top_n:     int = Field(50, description="Number of top cm-siRNA variants to return")
    full_scan: bool = Field(False, description="If true, run all 1260 variants. If false, run 40-variant mini-scan (E/D/Q/L on antisense pos 1-10).")


class MultiModRequest(BaseModel):
    sense:               str = Field(..., description="21-nt sense strand")
    antisense:           str = Field(..., description="21-nt antisense strand")
    sense_mods:          str = Field("", description="Mod symbols for sense strand, e.g. 'F,,M'")
    sense_positions:     str = Field("", description="Positions for sense mods, e.g. '2,5,,10,12'")
    antisense_mods:      str = Field("", description="Mod symbols for antisense strand")
    antisense_positions: str = Field("", description="Positions for antisense mods")
    model:               str = Field("A", description="Model: A, B, or C")


# ─── endpoints ────────────────────────────────────────────────────────────────

@app.get("/")
def root():
    return FileResponse(APP_HTML)


@app.get("/app.html")
def app_html():
    return FileResponse(APP_HTML)


@app.post("/rank")
def rank_endpoint(req: RankRequest):
    """
    Rank all 21-mer siRNA candidates by *naked (unmodified) silencing score*.

    Scores range 0–100. Very High ≥80, High 70–79, Moderate 55–69, Low <55.

    This is the baseline efficacy of each siRNA backbone; use the Single-Mod
    and Multi-Mod tabs to explore how chemical modifications improve silencing.
    """
    try:
        top = req.top_n if req.top_n > 0 else None
        results = rank_by_naked_score(req.sequence, top_n=top)

        return {
            "total_candidates": len(results),
            "results": [r.to_dict() for r in results],
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


@app.post("/rank/upload")
async def rank_upload(file: UploadFile = File(...), top_n: int = 20):
    """
    Same as /rank but accepts a FASTA file upload.
    """
    try:
        text = (await file.read()).decode("utf-8")
        top = top_n if top_n > 0 else None
        results = rank_by_naked_score(text, top_n=top)
        return {
            "filename": file.filename,
            "total_candidates": len(results),
            "results": [r.to_dict() for r in results],
        }
    except (FileNotFoundError, ValueError) as e:
        raise HTTPException(status_code=422, detail=str(e))


@app.post("/single-mod")
def single_mod_endpoint(req: SingleModRequest):
    """
    For one siRNA, generate and rank all 1260 single-modification cm-siRNA variants.

    Returns variants sorted by efficacy score (best first).
    delta_score = variant score − parent score (positive means improvement).

    Includes a `parent_toxicity` object so the UI can show the seed-baseline once at
    the top of the table instead of repeating the same % on every row (the seed
    canonical sequence is shared across all variants of one parent).
    """
    if req.model not in ("A", "B", "C"):
        raise HTTPException(status_code=422, detail="model must be A, B, or C")
    try:
        out = predict_modified(
            req.sense, req.antisense,
            mode="scan",
            model_key=req.model,
            full_scan=req.full_scan,
        )
        results = out["results"]
        parent_score = out["parent_score"]
        top = results[:req.top_n] if req.top_n > 0 else results
        # Parent seed toxicity (shared by all variants of this parent)
        from src.filters import toxicity_score, toxicity_label, seed_of_antisense
        parent_viab = toxicity_score(req.antisense)
        parent_seed = seed_of_antisense(req.antisense)
        return {
            "parent_sense":     req.sense,
            "parent_antisense": req.antisense,
            "parent_score":     parent_score,
            "model":            req.model,
            "total_variants":   len(results),
            "full_scan":        req.full_scan,
            "parent_toxicity": {
                "seed":            parent_seed,
                "viability":       None if parent_viab is None else round(parent_viab, 1),
                "label":           toxicity_label(parent_viab),
            },
            "results":          [r.to_dict() for r in top],
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


@app.post("/multi-mod")
def multi_mod_endpoint(req: MultiModRequest):
    """
    Predict the efficacy of a custom multi-modification cm-siRNA design.
    """
    if req.model not in ("A", "B", "C"):
        raise HTTPException(status_code=422, detail="model must be A, B, or C")
    try:
        out = predict_modified(
            req.sense, req.antisense,
            mode="multimod",
            model_key=req.model,
            sense_mods=req.sense_mods,
            sense_positions=req.sense_positions,
            antisense_mods=req.antisense_mods,
            antisense_positions=req.antisense_positions,
        )
        results = out["results"]
        parent_score = out["parent_score"]
        if not results:
            raise HTTPException(status_code=500, detail="No result generated.")
        r = results[0]
        return {
            "parent_sense":     req.sense,
            "parent_antisense": req.antisense,
            "parent_score":     parent_score,
            "model":            req.model,
            "result":           r.to_dict(),
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


# ─── Auto Multi-Mod Scan ───────────────────────────────────────────────────────

class MultiModScanRequest(BaseModel):
    sense: str
    antisense: str
    model: str = "A"
    max_mods: int = 2
    beam_width: int = 20
    full_scan: bool = False


@app.post("/multi-mod-scan")
def multi_mod_scan_endpoint(req: MultiModScanRequest):
    """
    Automatic multi-modification beam search.

    1. Runs single-mod scan to find top single-mod hits
    2. Combines top hits into multi-mod candidates (beam search)
    3. Returns ranked multi-mod variants with delta vs parent
    """
    if req.model not in ("A", "B", "C"):
        raise HTTPException(status_code=422, detail="model must be A, B, or C")
    if req.max_mods < 2 or req.max_mods > 3:
        raise HTTPException(status_code=422, detail="max_mods must be 2 or 3")
    if req.beam_width < 5 or req.beam_width > 50:
        raise HTTPException(status_code=422, detail="beam_width must be 5-50")

    try:
        from src.modification_engine import multi_mod_scan
        from src.predictor import _predict_naked, _normalize_scores
        from src.features import extract_batch_v4

        # Get parent score (naked model for consistency)
        X_parent = extract_batch_v4([req.sense], [req.antisense])
        raw_naked = _predict_naked(X_parent)
        parent_score = float(_normalize_scores(raw_naked, calibrator_key="normal")[0])

        # Run beam search
        variants = multi_mod_scan(
            req.sense,
            req.antisense,
            max_mods=req.max_mods,
            beam_width=req.beam_width,
            model_key=req.model,
            full_scan=req.full_scan,
        )

        # Format results
        results = []
        for i, v in enumerate(variants):
            results.append({
                "rank": i + 1,
                "sense": v.sense,
                "antisense": v.antisense,
                "mod_symbol": v.mod_symbol,
                "mod_position": v.mod_position,
                "mod_strand": v.mod_strand,
                "mod_positions": v.mod_positions or str(v.mod_position),
                "efficacy_score": round(v.efficacy_score, 2),
                "delta_score": round(v.delta_score, 2),
                "efficacy_label": _efficacy_label(v.efficacy_score),
            })

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": parent_score,
            "model": req.model,
            "total_variants": len(results),
            "results": results,
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


@app.get("/modifications")
def modifications_endpoint():
    """
    Return the list of 30 supported chemical modification symbols and their names.
    These symbols are used in the --sense-mods, --antisense-mods parameters.
    """
    mod_file = Path(__file__).parent.parent / "data" / "modification_codes.json"
    with mod_file.open() as f:
        data = json.load(f)
    return {
        "canonical": [m for m in data["modifications"] if m["type"] == "canonical"],
        "modifications": [m for m in data["modifications"] if m["type"] != "canonical"],
    }
