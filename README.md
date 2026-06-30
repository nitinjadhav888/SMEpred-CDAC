# HelixZero-CMS: Chemical Modification Scanner & Safety Engine

HelixZero-CMS is a production-grade, state-of-the-art computational engine and API designed to optimize Small Interfering RNA (siRNA) therapeutics. Unlike traditional models that only evaluate raw sequences, **its primary novelty lies in predicting the true clinical efficacy of *chemically modified* siRNAs** (via Model B). It provides end-to-end Machine Learning-driven prediction of both unmodified (naked) and heavily modified siRNA silencing efficacy, paired with exhaustive chemical modification scanning and transcriptome-wide biological safety validation.

Designed for clinical development teams and RNAi scientists, this framework drastically reduces the search space for potent, safe, and stable siRNA drugs.

## Key Capabilities

*   **Baseline Efficacy Ranking:** Evaluates every 21-mer candidate derived from a target transcript, applying thermodynamic boundaries (Reynolds/Ui-Tei rules) and a trained LightGBM model to predict raw silencing efficacy.
*   **Chemical Modification Scanning:** Rapidly scans the combinatoric chemical space of an siRNA candidate. Evaluates all 1,260 single-modification variants across 30 clinically relevant chemical moieties (e.g., 2'-OMe, 2'-F, LNA, GalNAc, Phosphorothioate linkages).
*   **Beam Search Combinatorial Design:** Autonomously stacks the most potent single modifications using a heuristic beam search to yield synergistic multi-modified siRNA designs.
*   **Biophysical & Pharmacokinetic Penalties:** Predicts sequence stability and loading potential, applying penalizations for lack of nuclease shielding, missing delivery conjugates, and steric RISC hindrance.
*   **Transcriptome-Wide Safety Validation:** Scans the siRNA seed and full length against the entire human transcriptome to mitigate catastrophic off-target slicing and innate immunostimulatory (TLR7/8) activation.

## System Architecture

The repository is modularized following strict multi-national corporation (MNC) production standards:

```text
HelixZero-CMS/
│
├── api/
│   └── main.py                 # FastAPI REST application exposing core functionality
│
├── src/
│   ├── predictor.py            # ML orchestration: prediction, feature building, and ranking
│   ├── modification_engine.py  # Combinatorial generation and beam search logic
│   ├── offtarget.py            # Transcriptome-wide heuristic safety validation engine
│   ├── biophysics.py           # Biophysical penalty calculation (PK, RISC, Immuno, Nuclease)
│   ├── features.py             # Feature extraction for Models (Naked & Unified)
│   ├── filters.py              # Tox-seed mitigation and sequence baseline functional filtering
│   ├── parser.py               # Robust sequence and FASTA parsing utilities
│   ├── utils.py                # Core constants, utilities, and helper functions
│   └── download_transcriptome.py # Utility to acquire human reference sequences
│
├── models/
│   └── lgb_..._model.txt       # Pre-trained LightGBM gradient boosting artifacts
│
├── data/
│   ├── modification_codes.json # Nomenclature library mapping chemical moieties
│   └── human_transcriptome.fasta # Transcriptomic baseline for safety scans
│
├── app.html                    # Single-Page Application (SPA) frontend interface
└── docs/                       # Auto-generated Clinical Safety Certificates
```

## Scientific Logic & Validation

This framework enforces rigorous scientific heuristics based on modern siRNA literature:

1.  **Thermodynamic Asymmetry:** Calculates 5' end ΔG differentials to predict Strand Bias, ensuring the antisense strand is preferentially loaded into the Ago2 RISC complex.
2.  **Slicer-Mediated (15-mer) Off-Targets:** Employs optimized string matching to absolutely reject any candidate sharing a 15-nucleotide contiguous identical match with an unintended transcript.
3.  **Seed-Mediated (miRNA-like) Mitigation:** Cross-references the seed region against the transcriptome to quantify partial matches. The penalty is drastically mitigated if specific chemical modifications (e.g., 2'-OMe) are placed at precise seed locations (pos 2 or 7) to structurally inhibit miRNA-like binding (Jackson et al., 2006).
4.  **Innate Immunogenicity Masking:** Scans for recognized Toll-Like Receptor (TLR7/TLR8) sequence motifs (e.g., `UGGC`, `GUUC`) and enforces 2'-O-methyl masking at the precise motif indices to evade interferon response cascades.
5.  **Toxicity Lookup:** Uses empirical cell-viability metrics (Janas et al., 2018) for 6-mer seed regions to predict base cytotoxicity, overriding warnings when clinical modifications are applied.

## Setup & Installation

**Prerequisites:** Python 3.10+

1.  **Install Requirements:**
    ```bash
    pip install fastapi uvicorn pydantic pandas numpy lightgbm joblib
    ```

2.  **Acquire Transcriptome Baseline:**
    Before utilizing the off-target module, download the reference transcriptome:
    ```bash
    python src/download_transcriptome.py
    ```

3.  **Start the REST API:**
    ```bash
    uvicorn api.main:app --reload --port 8000
    ```
    *The Single Page Application is served automatically at `http://localhost:8000/`.*

## License & Usage

This project is intended for research and discovery optimization within computational biology and RNAi therapeutics. All models and generated candidates must be empirically validated in *in vitro* and *in vivo* clinical settings.
