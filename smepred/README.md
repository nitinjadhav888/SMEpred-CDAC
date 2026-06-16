# SMEpred — siRNA Modification Efficacy Predictor

Predict which chemical modification pattern will silence your target gene best — before you spend money on synthesis. Given an mRNA target, SMEpred ranks all 1,260 modification variants per siRNA candidate by predicted efficacy (0–100), flags seed-based toxicity, and checks off-target complementarity.

---

## Quick Start

```bash
# Install
pip install -r requirements.txt

# Start the web app
uvicorn api.main:app --reload --port 8000
# → Open http://localhost:8000/docs

# Or use the command line
python cli/run.py rank --sequence AUGGAGGAGCCGCAGUCAGAUCCUAG --top 10
python cli/run.py single-mod --sense GCAGCACGACUUCUUCAAGUU --antisense CUUGAAGAAGUCGUGCUGCUU --top 20
```

On Windows, double-click `start_server.bat`.

---

## What It Solves

Designing an siRNA therapeutic means choosing a modification pattern from **30 modification types × 21 positions × 2 strands = 1,260 variants per siRNA**. Testing even 5 in the lab costs $1,000–$2,500 and takes 2–4 weeks.

**SMEpred eliminates the guesswork.** It ranks all 1,260 variants computationally using LightGBM:

| Model | Scope | Algorithm | Accuracy |
|-------|-------|-----------|----------|
| **LightGBM** | Modified siRNAs (all 1,260 variants) | Gradient-boosted trees, 152-d features | Within-gene PCC **0.68**, MAE **16.4%** |

---

## How It Works (End-to-End)

```
mRNA/gene sequence
        │
        ▼
[1] Generate all 21-mer siRNA candidates (sliding window)
        │
        ▼
[2] Extract 152-d feature vector per candidate:
    ├── Mononucleotide composition (base + modified)     140-d
    ├── Position-aware modification density               8-d
    ├── GC content                                        2-d
    └── Assay condition (dose, time)                      2-d
        │
        ▼
[3] LightGBM predicts efficacy (0–100)
        │
        ▼
[4] For your chosen siRNA → generate 1,260 modification variants
        │
        ▼
[5] Rank all variants by predicted efficacy + toxicity flags
```

---

## Model Training Details

### Training Environment

| Component | Value |
|-----------|-------|
| **Hardware** | CPU (no GPU required) |
| **Algorithm** | LightGBM LGBMRegressor (gradient-boosted trees) |
| **Training paradigm** | Boosting rounds (not epochs like neural networks) |
| **Boosting rounds** | 800 trees, early stopping at 50 rounds of no improvement |
| **Typical stopped at** | ~799 rounds (random holdout), ~26 rounds (gene-grouped holdout) |
| **Training time** | ~2 minutes on CPU for 25k rows × 152 features |
### What "Boosting Rounds" Means

Gradient-boosted trees train in **boosting rounds** — each round adds one new decision tree that corrects the errors of all previous trees. A model with 800 trees has iteratively added 800 trees, each one fixing what the previous ensemble got wrong.

The model uses **early stopping**: if the validation score doesn't improve for 50 consecutive rounds, training halts to prevent overfitting. The final model's iteration count (799 out of 800 max) shows it used nearly all trees, indicating the problem benefits from the full ensemble capacity.

### LightGBM Hyperparameters

| Parameter | Value | Purpose |
|-----------|-------|---------|
| `n_estimators` | 800 | Maximum boosting rounds |
| `learning_rate` | 0.03 | Step size per tree |
| `num_leaves` | 63 | Tree complexity cap |
| `subsample` | 0.8 | Row sampling (prevents overfit) |
| `colsample_bytree` | 0.7 | Feature sampling (prevents overfit) |
| `reg_lambda` | 1.0 | L2 regularization |
| `min_child_samples` | 30 | Minimum data per leaf |

---

## Current Performance

### cm-siRNA Model (25,763 rows, 13 genes)

| Metric | Value | What It Measures |
|--------|-------|-----------------|
| **Within-gene PCC** (random 5-fold CV) | **0.6789** | How well the model ranks modifications for a known gene |
| **Within-gene Spearman** | **0.6736** | Rank-order agreement for modification ranking |
| **Within-gene MAE** | **16.42%** | Average prediction error in inhibition percentage points |
| **Random-split PCC** | **0.6789** | Modification ranking for known siRNA (use case) |
| **Best single gene** (PCSK9) | **0.8465** | Some genes are very predictable |
| **Worst single gene** (PLN, n=95) | **0.4553** | Small-data genes are harder |

### Naked siRNA Model (4,060 rows, 4 sources)

| Source | n | PCC | Notes |
|--------|---|-----|-------|
| All (source-aware) | 4,060 | **0.5543** | Uses source one-hot encoding |
| **Taka** (independent) | 699 | **0.6905** | Best generalization — peer-reviewed data |
| Huesken (largest set) | 2,361 | 0.4179 | Public benchmark |
| Mix | 462 | 0.4622 | Mixed public data |
| **smepred_existing** | **538** | **0.0486** | Poor quality — labels may be unreliable |

### Homo Model (4,716 rows)
PCC **0.7370** — retrained as LightGBM (was SVR with unknown metrics).

---

## Key Limitations (Honest)

| Limitation | Impact | Workaround |
|------------|--------|------------|
| **smepred_existing labels** (PCC=0.05) | Our internal catalog data is noisy | We train with source one-hot; Taka and Huesken data are reliable |
| **smepred_existing labels** (PCC=0.05) | Our internal catalog data is noisy | We train with source one-hot; Taka and Huesken data are reliable |
| **No prospective wet-lab validation** | Scientific credibility requires experimental confirmation | Choose top-5 predictions for your gene, test in your lab, report hit rate |
| — | — | — |

---

## Improving Feature Quality

Our current 152-d features are **composition-based** — they encode *which* nucleotides and modifications are present, not *what the sequence means biologically*. Planned improvements include target mRNA features (secondary structure accessibility, GC content of target region) and learned embeddings (per RN.Ai-Predict findings).

---

## What Each File Does

| Path | Purpose |
|------|---------|
| `src/features.py` | Converts siRNA sequences to 152-d numerical feature vectors |
| `src/predictor.py` | Orchestrates ranking and modification prediction workflows |
| `src/parser.py` | Handles FASTA files, inline sequences, DNA→RNA conversion |
| `src/sirna_generator.py` | Sliding window to generate all 21-mer candidates |
| `src/modification_engine.py` | Generates 1,260 single-mod variants + custom multi-mod designs |
| `src/filters.py` | Seed toxicity lookup + functional filter rules |
| `cli/run.py` | Command-line interface (rank / single-mod / multi-mod) |
| `api/main.py` | FastAPI REST server with /rank, /single-mod, /multi-mod endpoints |
| `models/train_gbm_v3.py` | Training script (corrected: 152-d features, proper early stopping) |
| `models/model_a.pkl` | Deployed LightGBM (cm-siRNA) — 152-d, 799 trees |
| `models/model_normal.pkl` | Deployed LightGBM (naked siRNA) — 156-d with source one-hot |
| `models/model_homo.pkl` | Deployed LightGBM (dose-controlled cm-siRNA) |
| (removed) | — |
| `data/hetero_train_2728.csv` | 23,187 cm-siRNA training rows |
| `data/hetero_val_303.csv` | 2,576 cm-siRNA validation rows |
| `data/normal_siRNA_extended.csv` | 4,060 naked siRNA rows with source labels |

---

## Future Roadmap

| Priority | Item | Expected Gain |
|----------|------|---------------|
| 1 | Learned embeddings (per RN.Ai-Predict) for feature improvement | Within-gene MAE 16.4 → ~12 |
| 2 | More training genes (50+ via CRO partnership) | Broader gene coverage for cm model |
| 3 | Prospective wet-lab validation (20 siRNAs, ~$5k) | Scientific credibility, bioRxiv paper |
| 4 | Per-condition models (separate GBM per dose/cell line) | Within-gene MAE 16.4 → ~12 |
| 5 | Multi-target optimization (pick siRNA that hits multiple isoforms) | Clinical utility |

---

## Citation

If you use SMEpred in your research:

```bibtex
@software{smepred2026,
  title = {SMEpred: siRNA Modification Efficacy Predictor},
  year = {2026},
  description = {LightGBM for ranking chemical modification patterns in siRNA therapeutics}
}
```

---

## License

MIT
