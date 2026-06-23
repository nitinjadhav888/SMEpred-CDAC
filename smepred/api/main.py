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
from src.predictor import rank_by_naked_score, predict_modified, _efficacy_label
from src.biophysics import adjusted_efficacy_score

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
    input_type: str = Field("gene", description="'gene' (sliding window) or 'dsirna' (Dicer cleavage)")


class SingleModRequest(BaseModel):
    sense:     str = Field(..., description="21-nt sense strand")
    antisense: str = Field(..., description="21-nt antisense strand")
    model:     str = Field("B", description="Model: B (default, unified HelixZero model)")
    top_n:     int = Field(50, description="Number of top cm-siRNA variants to return")
    full_scan: bool = Field(False, description="If true, run all 1260 variants. If false, run 40-variant mini-scan (E/D/Q/L on antisense pos 1-10).")


class MultiModRequest(BaseModel):
    sense:               str = Field(..., description="21-nt sense strand")
    antisense:           str = Field(..., description="21-nt antisense strand")
    sense_mods:          str = Field("", description="Mod symbols for sense strand, e.g. 'F,,M'")
    sense_positions:     str = Field("", description="Positions for sense mods, e.g. '2,5,,10,12'")
    antisense_mods:      str = Field("", description="Mod symbols for antisense strand")
    antisense_positions: str = Field("", description="Positions for antisense mods")
    model:               str = Field("B", description="Model: B (default, unified HelixZero model)")


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
        results = rank_by_naked_score(req.sequence, top_n=top, input_type=req.input_type)

        return {
            "total_candidates": len(results),
            "input_type": req.input_type,
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
    try:
        out = predict_modified(
            req.sense, req.antisense,
            mode="scan",
            full_scan=req.full_scan,
        )
        results = out["results"]
        parent_score = out["parent_score"]
        model_b_baseline = out.get("model_b_baseline", parent_score)
        naked_baseline = out.get("naked_baseline", parent_score)
        top = results[:req.top_n] if req.top_n > 0 else results
        # Parent seed toxicity (shared by all variants of this parent)
        from src.filters import toxicity_score, toxicity_label, seed_of_antisense
        parent_viab = toxicity_score(req.antisense)
        parent_seed = seed_of_antisense(req.antisense)
        return {
            "parent_sense":     req.sense,
            "parent_antisense": req.antisense,
            "parent_score":     parent_score,
            "model_b_baseline": model_b_baseline,
            "naked_baseline":   naked_baseline,
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
    try:
        out = predict_modified(
            req.sense, req.antisense,
            mode="multimod",
            sense_mods=req.sense_mods,
            sense_positions=req.sense_positions,
            antisense_mods=req.antisense_mods,
            antisense_positions=req.antisense_positions,
        )
        results = out["results"]
        parent_score = out["parent_score"]
        model_b_baseline = out.get("model_b_baseline", parent_score)
        naked_baseline = out.get("naked_baseline", parent_score)
        if not results:
            raise HTTPException(status_code=500, detail="No result generated.")
        r = results[0]
        return {
            "parent_sense":     req.sense,
            "parent_antisense": req.antisense,
            "parent_score":     parent_score,
            "model_b_baseline": model_b_baseline,
            "naked_baseline":   naked_baseline,
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
    model: str = "B"
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
    if req.max_mods < 2 or req.max_mods > 21:
        raise HTTPException(status_code=422, detail="max_mods must be 2–21")
    if req.beam_width < 5 or req.beam_width > 50:
        raise HTTPException(status_code=422, detail="beam_width must be 5-50")

    try:
        from src.modification_engine import multi_mod_scan
        from src.predictor import _predict_naked, _normalize_scores
        from src.features import extract_batch_v4

        # Get parent score (naked model for display)
        from src.biophysics import adjusted_efficacy_score
        X_parent = extract_batch_v4([req.sense], [req.antisense])
        raw_naked = _predict_naked(X_parent)
        raw_parent = float(_normalize_scores(raw_naked, calibrator_key="normal")[0])
        raw_naked_adj, _, _ = adjusted_efficacy_score(raw_parent, req.sense, req.antisense, req.sense, req.antisense)
        naked_baseline = round(raw_naked_adj, 2)
        # Also compute Model B baseline (different feature space - fair comparison)
        from src.features import extract_positional_features_batch
        from src.predictor import _get_model
        X_parent_b = extract_positional_features_batch([req.sense], [req.antisense], [req.sense], [req.antisense])
        model_b = _get_model("B")
        raw_b = float(model_b.predict(X_parent_b)[0])
        model_b_adj, _, _ = adjusted_efficacy_score(raw_b, req.sense, req.antisense, req.sense, req.antisense)
        parent_score = round(model_b_adj, 2)
        model_b_baseline = round(model_b_adj, 2)

        # Run beam search
        variants = multi_mod_scan(
            req.sense,
            req.antisense,
            max_mods=req.max_mods,
            beam_width=req.beam_width,
            model_key=req.model,
            full_scan=req.full_scan,
        )

        # Format results (already biophysically adjusted in score_variants)
        results = []
        for i, v in enumerate(variants):
            p = getattr(v, 'penalties', None) or {}
            total_pen = sum(p.values())
            raw_score = round(v.efficacy_score + 0.70 * total_pen, 2)
            entry = {
                "rank": i + 1,
                "sense": v.sense,
                "antisense": v.antisense,
                "mod_symbol": v.mod_symbol,
                "mod_position": v.mod_position,
                "mod_strand": v.mod_strand,
                "mod_positions": v.mod_positions or str(v.mod_position),
                "raw_efficacy_score": raw_score,
                "efficacy_score": round(v.efficacy_score, 2),
                "total_penalty": round(total_pen, 1),
                "delta_score": round(v.delta_score, 2),
                "efficacy_label": _efficacy_label(v.efficacy_score),
                "penalties": {k: round(v, 1) for k, v in p.items()},
            }
            results.append(entry)

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": parent_score,
            "model_b_baseline": model_b_baseline,
            "naked_baseline": naked_baseline,
            "model": req.model,
            "total_variants": len(results),
            "results": results,
        }
    except FileNotFoundError as e:
        raise HTTPException(status_code=503, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))


# ─── Multi-Mod from Single-Mod Results ──────────────────────────────────────────

class MultiModFromSingleRequest(BaseModel):
    sense: str
    antisense: str
    model: str = "B"
    max_mods: int = 5
    beam_width: int = 20
    full_scan: bool = True
    single_results: Optional[list] = None
    parent_score: Optional[float] = None
    seed_variant: Optional[dict] = None
    calibrator_key: Optional[str] = None
    normalize_mode: str = "rescale"


@app.post("/multi-mod-from-single")
def multi_mod_from_single_endpoint(req: MultiModFromSingleRequest):
    """
    Generate multi-modification candidates.

    If single_results are provided, they are used as building blocks for beam search.
    If not, the beam search runs its own single-mod scan internally — use this for
    a standalone multi-mod call directly from a parent siRNA (no single-mod prerequisite).
    """
    if req.max_mods < 2 or req.max_mods > 21:
        raise HTTPException(status_code=422, detail="max_mods must be 2–21")
    if req.beam_width < 5 or req.beam_width > 100:
        raise HTTPException(status_code=422, detail="beam_width must be 5–100")

    try:
        from src.modification_engine import multi_mod_scan, CmSiRNA
        from src.predictor import _efficacy_label

        # Convert raw dicts to lightweight objects matching RankedCmSiRNA interface
        class SingleModProxy:
            def __init__(self, d):
                self.mod_symbol = d.get("mod_symbol") or d.get("modification", "")
                self.mod_position = d.get("mod_position") or d.get("position", 0)
                self.mod_strand = d.get("mod_strand") or d.get("strand", "")
                self.mod_positions = d.get("mod_positions") or str(self.mod_position)
                self.efficacy_score = d.get("efficacy_score") or d.get("score", 0)
                self.sense = d.get("sense", "")
                self.antisense = d.get("antisense", "")
                self.delta_score = d.get("delta_score", 0)

        single_results = [SingleModProxy(sr) for sr in req.single_results] if req.single_results is not None else None
        seed = SingleModProxy(req.seed_variant) if req.seed_variant else None

        # When single_results not provided, compute parent_score from the naked model
        calibrator_key = req.calibrator_key
        if req.parent_score is None:
            from src.predictor import _predict_naked, _normalize_scores
            from src.features import extract_batch_v4
            from src.biophysics import adjusted_efficacy_score
            X_parent = extract_batch_v4([req.sense], [req.antisense])
            raw = float(_normalize_scores(_predict_naked(X_parent), calibrator_key=calibrator_key, mode=req.normalize_mode)[0])
            parent_adj, _, _ = adjusted_efficacy_score(raw, req.sense, req.antisense, req.sense, req.antisense)
            parent_score = round(parent_adj, 2)
            naked_baseline = round(parent_adj, 2)
        else:
            parent_score = req.parent_score
            from src.biophysics import adjusted_efficacy_score
            naked_adj, _, _ = adjusted_efficacy_score(parent_score, req.sense, req.antisense, req.sense, req.antisense)
            naked_baseline = round(naked_adj, 2)

        # Compute Model B baseline for fair multi-mod comparison
        from src.features import extract_positional_features_batch
        from src.predictor import _get_model
        X_parent_b = extract_positional_features_batch([req.sense], [req.antisense], [req.sense], [req.antisense])
        mb = _get_model("B")
        raw_b = float(mb.predict(X_parent_b)[0])
        # Model B output is already normalized (mode=identity)
        mb_adj, _, _ = adjusted_efficacy_score(raw_b, req.sense, req.antisense, req.sense, req.antisense)
        model_b_baseline = round(mb_adj, 2)

        variants = multi_mod_scan(
            req.sense,
            req.antisense,
            max_mods=req.max_mods,
            beam_width=req.beam_width,
            model_key=req.model,
            full_scan=req.full_scan,
            single_results=single_results,
            parent_score=parent_score,
            seed_variant=seed,
            calibrator_key=calibrator_key,
            normalize_mode=req.normalize_mode,
        )

        results = []
        for i, v in enumerate(variants):
            p = getattr(v, 'penalties', None) or {}
            total_pen = sum(p.values())
            raw_score = round(v.efficacy_score + 0.70 * total_pen, 2)
            entry = {
                "rank": i + 1,
                "sense": v.sense,
                "antisense": v.antisense,
                "mod_symbol": v.mod_symbol,
                "mod_position": v.mod_position,
                "mod_strand": v.mod_strand,
                "mod_positions": v.mod_positions or str(v.mod_position),
                "raw_efficacy_score": raw_score,
                "efficacy_score": round(v.efficacy_score, 2),
                "total_penalty": round(total_pen, 1),
                "delta_score": round(v.delta_score, 2),
                "efficacy_label": _efficacy_label(v.efficacy_score),
                "penalties": {k: round(v, 1) for k, v in p.items()},
            }
            results.append(entry)

        return {
            "parent_sense": req.sense,
            "parent_antisense": req.antisense,
            "parent_score": parent_score,
            "model_b_baseline": model_b_baseline,
            "naked_baseline": naked_baseline,
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
