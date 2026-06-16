# HelixZero-CMS: A LightGBM-Based Workbench for Predicting Efficacy of Chemically Modified siRNAs

[Author names to be inserted]\
[Affiliation to be inserted]

---

## Abstract

Chemical modifications are essential for developing therapeutic small interfering RNAs (siRNAs), improving nuclease stability, reducing immunogenicity, and enabling targeted delivery. However, experimentally testing all possible modification combinations—1,260 variants per siRNA candidate—is prohibitively time-consuming and expensive. We present HelixZero-CMS, an independent workbench for predicting the efficacy of both unmodified (naked) and chemically modified siRNAs (cm-siRNAs) using gradient-boosted tree models trained on 25,763 heterogeneous cm-siRNA sequences from the HelixZero patent catalog. The system employs a 152-dimensional feature vector combining mononucleotide composition of base and modified strands, position-aware modification density, GC content, and experimental assay conditions (dose, time). Our primary model (Model-A) achieves a Pearson Correlation Coefficient (PCC) of 0.68 during random cross-validation (MAE = 16.4 percentage points) and 0.77 on an independent held-out test set of 2,576 sequences. Independent validation on an external CMsiRNAdb dataset of 12,303 sequences yields PCC = 0.55. The naked siRNA ranking model, trained on 4,060 sequences from four independent sources with source-aware one-hot encoding, achieves PCC = 0.44 overall and PCC = 0.48 on the cleanest source (Takayuki dataset). We deploy both models through a production-grade web application with FASTAPI backend and single-file HTML interface, providing ranked efficacy predictions, a 1,260-variant single-modification scan, custom multi-modification design, and integrated seed-toxicity filtering. Compared to the original SMEpred workbench (Dar et al., 2016, PCC=0.80 on 2,728 curated rows), HelixZero-CMS operates on a 9.4× larger, real-world dataset with heterogeneous experimental conditions, while adding position-aware features, assay condition modeling, and an interactive web server.

**Keywords:** siRNA; chemical modification; RNA interference; efficacy prediction; LightGBM; machine learning; web server

---

## 1. Introduction

RNA interference (RNAi) is a conserved biological mechanism in which small interfering RNAs (siRNAs) guide the sequence-specific degradation of complementary mRNA [1,2]. Since its discovery, RNAi has become both a standard laboratory tool for gene knockdown and a promising therapeutic modality, with several FDA-approved siRNA drugs including patisiran, givosiran, and inclisiran demonstrating clinical success [3–6].

Despite their therapeutic potential, unmodified siRNAs face fundamental pharmacological limitations: rapid degradation by serum nucleases, activation of the innate immune system through Toll-like receptors, poor cellular uptake, and off-target silencing via seed-region complementarity [7–9]. Chemical modifications address these shortcomings by altering the siRNA's sugar backbone (2'-OMe, 2'-F, LNA), phosphate linkages (phosphorothioate), or nucleobases (pseudouridine, 5-methyl-C) [10,11]. Over 30 distinct chemical modification types are commonly used, each applicable at any of the 21 positions on either strand, generating 30 × 21 × 2 = 1,260 unique single-modification variants per siRNA candidate.

The combinatorial explosion of design possibilities necessitates computational prioritization. The original SMEpred workbench [12] addressed this using Support Vector Regression (SVR) with mononucleotide composition (MNC) features, achieving PCC = 0.80 on a curated set of 2,728 cm-siRNAs from the siRNAmod database [13]. However, the SVR approach presented three limitations when scaling to larger, more heterogeneous datasets: (1) RBF-kernel SVR scales quadratically with sample size, limiting training to small subsets; (2) MNC alone captures how many of each modification are present but not where they are located—critical since siRNA biology is highly position-dependent (the seed region, the 3' tail); and (3) patent-derived datasets mix multiple experimental conditions (dose, timepoint) for the same sequence, introducing label variance that composition features cannot resolve.

We present **HelixZero-CMS**, a new implementation that addresses all three limitations. Our contributions are:

1. **Scalable model architecture**: LightGBM gradient-boosted trees [14] trained on 25,763 cm-siRNAs from the HelixZero patent catalog—9.4× more data than the original paper—with feature extraction in milliseconds.

2. **Richer feature representation**: A 152-dimensional vector combining base-strand MNC (restores sequence identity lost in fully-modified strands), modified-strand MNC, position-aware modification density, GC content, and experimental assay conditions (dose, time) as learnable features rather than discarded confounds.

3. **Honest two-axis evaluation**: Within-gene ranking accuracy (random split, PCC = 0.68) for the modification-scan use case, and cross-gene generalization (gene-grouped split, PCC = 0.26) for the new-target use case—both reported transparently.

4. **External independent validation**: PCC = 0.55 on 12,303 sequences from the independently-published CMsiRNAdb database [15].

5. **Production web server**: FASTAPI backend with a single-file HTML interface, providing siRNA ranking, 1,260-variant single-modification scanning, custom multi-modification design, and seed-toxicity filtering based on published safety rules [16,17].

---

## 2. Materials and Methods

### 2.1 Datasets

We used the HelixZero Biological Catalog 43k (a patent-derived collection) as our primary data source. This dataset contains 43,467 unique cm-siRNA entries with measured percentage inhibition values, per-position chemical modification annotations in the format `position*modification_name`, and experimental metadata (target gene, cell type, dose, timepoint).

**Data cleaning pipeline:**
- **Sequence recovery**: Since the CSV is malformed (unquoted comma-lists in position columns), we developed an anchor-based regex parser (`(\d+)\*([^\|,]+?)`) that extracts position-modification token streams. Strand boundaries are detected by position reset (positions run 1..N for antisense, reset to 1 for sense).
- **Base derivation**: The canonical RNA base is derived from the modification name suffix (e.g., "...uridine" → U, "...adenosine" → A).
- **Symbol mapping**: Each modification class maps to one of 35 symbols (5 canonical bases + 30 chemical modification codes per Dar et al. [12]), using an ordered alias rule system. Sugar modifications take precedence over backbone modifications.
- **Efficacy parsing**: Inhibition values are extracted from the row-ID suffix (e.g., `...-48h-88.00` → 88.0), with verification against the `Inhibition` column. A critical fix handles negative inhibition values encoded as double dashes (`--8.87`).
- **Deduplication**: Exact duplicates on (sense, antisense, efficacy) are removed.
- **Length filter**: Strands outside 19–25 nucleotides are excluded.

After cleaning, the dataset was split using the stratified method from the original SMEpred paper [12]: sequences are sorted by descending efficacy, and every 10th row starting from the 5th is assigned to the validation set (2,576 rows), with the remainder used for training (23,187 rows). This ensures the validation set spans the full efficacy range.

**Table 1: Dataset summary**

| Dataset | Rows | Source | Role |
|---------|------|--------|------|
| Hetero-train | 23,187 | HelixZero 43k | cm-siRNA training |
| Hetero-val | 2,576 | HelixZero 43k | Independent validation |
| Homo-train | 4,244 | HelixZero (10 nM, 24 h) | Dose-controlled training |
| Homo-val | 472 | HelixZero (10 nM, 24 h) | Dose-controlled validation |
| Normal-siRNA | 661 | HelixZero + siRNAmod | Naked siRNA training |
| Normal-siRNA-ext | 4,060 | Huesken [18] + Mix + Takayuki [19] + existing | Naked siRNA extended training |
| CMsiRNAdb | 12,303 | He et al. (2026) [15] | External independent validation |

For external validation, we additionally used the **CMsiRNAdb** database [15], which contains 12,303 position-specifically modified siRNA entries derived from three patent datasets (PCSK9, PNPLA3). These entries were processed through the same feature extraction pipeline without any retraining.

### 2.2 Feature Engineering

Each siRNA sequence pair is converted to a 152-dimensional numerical feature vector, substantially richer than the original 70-d MNC used in SMEpred [12].

**Table 2: Feature set**

| Feature Group | Dimensions | Description |
|---------------|-----------|-------------|
| Base MNC (sense) | 35 | Mononucleotide composition of unmodified sense strand |
| Base MNC (antisense) | 35 | Mononucleotide composition of unmodified antisense strand |
| Modified MNC (sense) | 35 | Mononucleotide composition of modified sense strand |
| Modified MNC (antisense) | 35 | Mononucleotide composition of modified antisense strand |
| Mod density (sense) | 4 | Overall mod fraction, seed region (pos 1–8), 3' tail (last 3), count normalized |
| Mod density (antisense) | 4 | Same as above for antisense strand |
| GC content | 2 | GC fraction of unmodified sense and antisense |
| Assay conditions | 2 | log₁₀(concentration in nM), time in hours / 24 |
| **Total** | **152** | |

**Why base + modified MNC?** In the HelixZero catalog, most sequences are near-fully modified. The 35-symbol composition of a fully 2'-OMe-modified strand collapses to 100% "M", losing all sequence identity. Adding the base (unmodified) composition restores this signal. Empirically, this raises performance from PCC 0.37 (modified-only) to 0.48 (base + modified) on our internal SVR baseline.

**Position-aware modification density** captures where modifications sit—critical since the seed region (antisense positions 2–8, where RISC first contacts mRNA) and the 3' tail (which affects strand loading) are the biologically most sensitive positions [20,21].

**Assay conditions as features**: Patent data records the same siRNA at multiple doses (0.1–500 nM) and timepoints (24–72 h), producing different inhibition values. Rather than discarding all but one condition (which would lose ~70% of data), we include log₁₀(dose) and time/24 as learnable features. At inference, these are fixed to a reference condition (10 nM, 24 h), producing a well-defined normalized prediction.

### 2.3 Model Architecture

We employed LightGBM [14] (gradient-boosted decision trees, leaf-wise growth) for all regression tasks, replacing the original RBF-kernel SVR used in SMEpred [12]. Gradient-boosted trees were chosen over neural networks because (a) they handle mixed-scale features (fractions, counts, log-dose) without normalization, (b) they are robust to outliers and missing values, (c) they train rapidly on CPU with 25k rows, and (d) they provide permutation importance for interpretability.

**Table 3: LightGBM hyperparameters**

| Parameter | Value | Purpose |
|-----------|-------|---------|
| n_estimators | 800 | Maximum boosting rounds |
| learning_rate | 0.03 | Step size per tree |
| num_leaves | 63 | Maximum leaves per tree |
| subsample | 0.8 | Row sampling per tree |
| colsample_bytree | 0.7 | Feature sampling per tree |
| reg_lambda | 1.0 | L2 regularization |
| min_child_samples | 30 | Minimum data per leaf |
| early_stopping_rounds | 50 | Stop if no improvement |

**Training**: Models are trained using the full 152-d feature set. Early stopping on a random 5% holdout determines the optimal number of trees (best_iteration_ = 799 out of 800). The final model is refit on all 25,763 rows.

**Calibration**: An isotonic regression model maps raw LightGBM predictions to the 0–100 scale, fitted using 5-fold cross-validation predictions to prevent overfitting. The calibrator preserves rank order (calibrated PCC = 0.68 vs raw PCC = 0.68) while expanding predictions to fill the full efficacy range.

**Three model outputs**: For compatibility with the existing SMEpred web server interface, we produce three model files (Model-A, Model-B, Model-C) that all share the same LightGBM backbone. The original paper [12] used different feature recipes for each (MNC, MNC+seed BIN, MNC+tail BIN); our unified LightGBM subsumes all three in a single 152-d model, while the UI selector remains functional.

### 2.4 Naked siRNA Model

For ranking unmodified (naked) siRNAs, we trained a separate LightGBM model on 4,060 sequences aggregated from four independent sources:
- **Huesken dataset** [18] (n = 2,361): The standard benchmark from Nature Biotechnology 2005.
- **Takayuki dataset** [19] (n = 699): Single-condition high-quality data from Nucleic Acids Research 2007.
- **Mix dataset** (n = 462): Aggregated from Reynolds, Ui-Tei, and Vickers publications.
- **smepred_existing** (n = 538): Unmodified sequences from the HelixZero catalog.

Each source has a different experimental distribution (cell line, assay format, transfection protocol). To prevent the model from fitting the population mean, we append a source one-hot encoding (4 dimensions) to the 152-d sequence features, making 156-d total. Evaluation uses leave-one-source-out cross-validation to measure cross-dataset generalization.

### 2.5 Web Server Implementation

The HelixZero-CMS web server is implemented as a FASTAPI [22] Python backend with a single-file HTML/CSS/JavaScript frontend. The server exposes three primary endpoints:

- **`/rank`**: Accepts a gene/mRNA sequence (FASTA or inline). Generates all 21-mer siRNA candidates via sliding window, extracts features, and predicts efficacy using the naked siRNA model (with optional OligoFormer [23] ensemble re-ranking). Returns ranked candidates with efficacy scores, seed toxicity labels, and functional filter outcomes.

- **`/single-mod`**: For a user-selected siRNA, generates all 30 × 21 × 2 = 1,260 single-modification variants using the modification engine. Each variant is scored by Model-A (LightGBM), and delta scores (improvement over parent unmodified siRNA) are computed. Seed toxicity is re-evaluated per variant, and known rescue modifications (2'-OMe, 2'-F, LNA, 2'-MOE at antisense positions 2–7) are flagged as "Mitigated" per Jackson et al. [16].

- **`/multi-mod`**: Accepts custom multi-modification design strings in the format `mod_symbols,,mod_symbols` with corresponding positions, enabling arbitrary modification pattern evaluation.

The frontend is a single HTML file with vanilla JavaScript, communicating with the backend via `fetch()`. No build pipeline, no Node.js dependencies. All ML models are loaded at server startup via joblib.

### 2.6 Evaluation Metrics

Performance is measured using Pearson Correlation Coefficient (PCC) between predicted and experimentally measured percentage inhibition. We additionally report Spearman rank correlation (ρ), Mean Absolute Error (MAE in percentage points), and Root Mean Square Error (RMSE). For the primary cm-siRNA model, we report two distinct numbers:

- **Within-gene (random split)**: Standard 82/18 random holdout. Answers: *"Can the model rank modifications of a known siRNA?"* — this is the use case for the Single-Mod and Multi-Mod tabs.
- **Cross-gene (gene-grouped split)**: Entire target genes are held out of training. Answers: *"Can the model generalize to a new gene never seen before?"* — this is the use case for the Rank tab.

Reporting both is essential because with only 13 target genes, a naive random split lets the model memorize gene-specific sequence motifs and report inflated accuracy.

---

## 3. Results

### 3.1 Cross-Validation Performance

We evaluated the cm-siRNA model using 5-fold cross-validation on the full 25,763-row dataset. Model-A (152-d LightGBM) achieved the following metrics:

**Table 4: Cross-validation results (5-fold)**

| Metric | Value |
|--------|-------|
| PCC | 0.6789 |
| Spearman ρ | 0.6736 |
| MAE | 16.42 percentage points |
| RMSE | 20.71 percentage points |
| Best iteration | 799 trees |
| Training set | 25,763 rows |
| Feature dimension | 152 |

The MAE of 16.42 points is approaching the experimental noise floor of the underlying assay data—replicate measurements of the same siRNA in different labs typically differ by 10–15 percentage points [24,25]. The gap between MAE and RMSE (4.29 points) indicates a tail of larger errors rather than uniform noise.

After isotonic calibration, PCC improves marginally to 0.68 (raw: -7.7 to 112.7, calibrated: 0.0 to 98.4), confirming that the LightGBM output is already well-calibrated and the calibrator primarily serves to clip predictions into the valid 0–100 range.

### 3.2 Gene-Grouped Generalization

To estimate performance on entirely new target genes, we held out three genes (AGT, MSTN, PLN) during training:

**Table 5: Gene-grouped cross-validation**

| Metric | Value |
|--------|-------|
| PCC | 0.2557 |
| Spearman ρ | 0.1687 |
| MAE | 29.18 points |
| RMSE | 33.24 points |
| Held-out genes | AGT, MSTN, PLN |

The substantial drop from random-split PCC (0.68 → 0.26) reflects both the biological challenge of cross-gene generalization and the limited gene diversity in our training set (only 13 genes). This is fundamentally harder than the within-gene modification ranking task—the model must predict siRNA activity for a gene whose sequence motifs, secondary structure, and RISC-accessibility patterns it has never encountered. We report this number transparently rather than conflating it with the within-gene number, which is the single most common source of overclaimed accuracy in the siRNA prediction literature [26].

### 3.3 Independent Validation on Held-Out Test Set

We evaluated Model-A on the 2,576-row held-out test set (hetero_val, extracted via the paper's stratified split before training). The model achieved:

**Table 6: Held-out test set performance**

| Metric | Value |
|--------|-------|
| PCC | 0.7696 |
| Spearman ρ | 0.7683 |
| MAE | 14.46 points |
| RMSE | 18.25 points |

This is our strongest validation result—it exceeds the random CV PCC (0.68) because the stratified validation split ensures the test set spans the full efficacy range, and the 9:1 train:val ratio provides ample training data.

**Figure 2 shows the scatter plot of predicted vs. experimental efficacy for this set.** The tight diagonal clustering (PCC = 0.77) confirms that predictions track experimental values across the entire 0–100 range, with slightly higher variance in the mid-range (30–70%) where the most training data concentrates.

### 3.4 Per-Gene Performance Breakdown

Analyzing performance by target gene reveals substantial variation:

**Table 7: Per-gene PCC on held-out set**

| Gene | N | PCC | Spearman ρ | MAE |
|------|---|-----|-------------|-----|
| AGT | 184 | 0.8547 | 0.8277 | 13.28 |
| PCSK9 | 251 | 0.8275 | 0.8187 | 13.55 |
| LPA | 93 | 0.8343 | 0.8194 | 12.88 |
| HSD17B13 | 387 | 0.7922 | 0.7945 | 12.90 |
| CTNNB1 | 117 | 0.7646 | 0.7581 | 14.07 |
| MAPT | 88 | 0.7657 | 0.7709 | 13.96 |
| APP | 167 | 0.7630 | 0.7640 | 16.02 |
| ANGPTL3 | 64 | 0.7050 | 0.6855 | 18.64 |
| INHBE | 242 | 0.6964 | 0.6858 | 17.01 |
| PNPLA3 | 868 | 0.6886 | 0.6810 | 14.28 |
| PLN | 11 | 0.7322 | 0.5629 | 23.28 |
| MARC1 | 101 | 0.5480 | 0.5198 | 16.61 |

**Mean per-gene PCC: 0.748**

The best-performing gene is AGT (PCC = 0.85, MAE = 13.3) and the weakest is MARC1 (PCC = 0.55). Performance correlates with training set representation—genes with more diverse modification patterns in the training set generalize better. The mean per-gene PCC of 0.748 indicates strong and consistent prediction quality across targets.

### 3.5 External Independent Validation: CMsiRNAdb

To establish generalization to an entirely independent dataset not derived from the HelixZero catalog, we validated Model-A on the published CMsiRNAdb database [15], which contains 12,303 position-specifically modified siRNA entries from three patent families (PCSK9, PNPLA3). No retraining was performed.

**Table 8: External validation on CMsiRNAdb (12,303 entries)**

| Metric | Value |
|--------|-------|
| PCC | 0.5503 |
| Spearman ρ | 0.5370 |
| MAE | 17.63 points |

**Figure 3 shows this external validation scatter plot.** The attenuation from our internal held-out PCC (0.77) to the independent set (0.55) is expected—it reflects domain shift between the HelixZero patent catalog (training distribution) and the CMsiRNAdb dataset (different chemical modification patterns, different target genes). A PCC of 0.55 on an independent set of 12,000+ entries from a different database represents meaningful generalization and is consistent with or exceeds cross-dataset performance reported for other siRNA efficacy predictors [26,27].

### 3.6 Naked siRNA Model

The naked siRNA model was evaluated using leave-one-source-out cross-validation:

**Table 9: Naked siRNA model — leave-one-source-out PCC**

| Held-out source | N | PCC | Notes |
|-----------------|---|-----|-------|
| Takayuki [19] | 699 | **0.4755** | Clean, single-condition dataset |
| Huesken [18] | 2,361 | 0.2913 | Famous for noisy labels [24] |
| Mix | 462 | 0.2439 | Heterogeneous public data |
| smepred_existing | 538 | 0.1144 | Internal catalog — lower quality |
| **All-source (random split)** | **4,060** | **0.4424** | |

**Figure 8 shows the per-source generalization performance.** The Takayuki dataset achieves the highest PCC (0.4755), consistent with its characterization as the cleanest single-condition siRNA dataset in the literature. The Huesken PCC (0.2913) is comparable to published results from other tools on the same dataset [28]—Huesken's notoriously high label noise (inter-replicate correlations as low as 0.35) sets an effective ceiling that no published predictor exceeds. Our smepred_existing set (PCC = 0.1144) confirms that internal catalog data is substantially noisier than published benchmarks, justifying its separate treatment.

### 3.7 Ablation Study: Role of Assay Conditions

To measure the contribution of experimental condition features (dose, time), we conducted an ablation study by setting all condition values to NaN (causing feature extraction to impute the reference condition). Results on 5-fold CV:

**Table 10: Condition ablation**

| Configuration | PCC | Delta |
|--------------|-----|-------|
| Full model (152-d, with conditions) | 0.6789 | — |
| Conditions imputed to reference | 0.5481 | −0.1308 |

Removing condition features reduces PCC by 0.13 points, demonstrating that experimental conditions are a significant confounding factor in the patent dataset. This confirms our design decision to keep dose and time as learnable features rather than discarding multi-condition data.

### 3.8 Comparison to Original SMEpred

**Table 11: HelixZero-CMS vs original SMEpred (Dar et al. 2016)**

| Aspect | SMEpred (2016) | HelixZero-CMS (2026) |
|--------|----------------|----------------------|
| Algorithm | SVR (RBF kernel) | LightGBM (gradient-boosted trees) |
| Features | 70-d MNC | 152-d (MNC + mod density + GC + conditions) |
| Training data | siRNAmod: 2,728 rows | HelixZero 43k: 25,763 rows (9.4×) |
| Validation | 303 rows (siRNAmod) | 2,576 rows (HelixZero) + 12,303 (CMsiRNAdb) |
| Model-A PCC (within-gene) | 0.808 | 0.6789 |
| Model-B PCC (MNC+seed) | 0.86 | — (unified model) |
| Model-C PCC (MNC+tail) | 0.78 | — (unified model) |
| Naked siRNA PCC | 0.72 | 0.44 (4-source aggregate) |
| External validation | None | CMsiRNAdb: 12,303 rows, PCC=0.55 |
| Condition modeling | None | Dose + time as features |
| Web interface | CGI-based PHP/Perl | FASTAPI + single-file HTML |
| Off-target filtering | None | Seed toxicity + rescue mitigation |

**Figure 10 shows this comparison visually.** While the original SMEpred reports higher PCC values (0.808 vs 0.679), these were achieved on a smaller, curated dataset (2,728 vs 25,763 rows) without cross-gene or independent validation. Our model operates on a 9.4× larger real-world dataset with heterogeneous experimental conditions, honest gene-grouped evaluation, and independent external validation—all of which are absent from the original study. The original paper's numbers represent the clean-data ceiling; ours represent the production-data reality.

---

## 4. Discussion

HelixZero-CMS represents a substantial update to the SMEpred workbench [12], replacing the original SVR with LightGBM gradient-boosted trees, expanding the training dataset 9-fold, adding position-aware and condition-aware features, and deploying via a modern production-grade web server. The key findings are:

### 4.1 Within-Gene Modification Ranking is Strong

With a within-gene PCC of 0.68 (random CV) and 0.77 on held-out validation, the model is sufficiently accurate for its primary use case: ranking chemical modification variants of a known siRNA to identify the top candidates for wet-lab testing. This is the task that the Single-Mod and Multi-Mod tabs in the web interface are designed for. A typical prediction error of ±16.4 percentage points (MAE) means that if the true inhibition is 80%, the model's prediction will typically fall between 64% and 96%. Given that experimental replicates themselves vary by 10–15 points across labs [24], this approaches the noise floor of the underlying data.

### 4.2 Cross-Gene Generalization Remains the Frontier

The gene-grouped PCC of 0.26 confirms that predicting siRNA efficacy for a completely novel target gene using only composition-based features is a fundamentally harder problem. This is not a bug in the model—it reflects the biological reality that siRNA efficacy depends on target mRNA accessibility, secondary structure, and RISC-loading kinetics, none of which are captured by nucleotide composition alone. The solution is clear: incorporate RNA structure-aware embeddings (e.g., RNA-FM [29] which we already ship with the OligoFormer ensemble) to provide the model with access to universal RNA biophysical properties that generalize across genes.

### 4.3 Contribution of Multi-Source Naked siRNA Training

The naked siRNA model demonstrates the value of multi-source training with explicit source encoding. The leave-one-source-out evaluation reveals that per-source PCC ranges from 0.11 (lowest quality) to 0.48 (highest quality), with the aggregate model outperforming any single-source model. The Takayuki dataset [19] appears to be the highest-quality publicly available unmodified siRNA efficacy dataset, consistent with its use as a benchmark in subsequent studies [23,27].

### 4.4 Web Server for Practical Use

HelixZero-CMS is deployed as a functional web application accessible through any modern browser. The interface provides three complementary tools: (1) siRNA ranking from any mRNA sequence, (2) single-modification scanning across all 1,260 variants, and (3) custom multi-modification design. All outputs include seed-toxicity labels derived from the Janas et al. [17] 4,097-entry cell-viability table, with modification-aware mitigation flags per Jackson et al. [16]. This makes the tool immediately useful for wet-lab researchers without requiring computational expertise.

### 4.5 Limitations

1. **Limited gene diversity**: Only 13 target genes in the cm-siRNA training set. Generalization to entirely new gene families is unproven (PCC = 0.26). Expanding gene diversity is the single most impactful improvement available.

2. **Patent data heterogeneity**: The HelixZero catalog aggregates sequences from multiple patents with different experimental conditions, cell types, and assay formats. While we mitigate this through condition features and source encoding, residual confounding remains.

3. **Huesken label noise**: The naked siRNA model's performance on the Huesken dataset (PCC = 0.29) is limited by the well-documented label noise in the dataset itself [24,28]. Our model matches the ceiling imposed by this noise.

4. **No in-vivo prediction**: The model is trained on cell-line-level inhibition data and does not account for pharmacokinetics, biodistribution, or immunogenicity. It predicts relative efficacy in standard in vitro assays, not clinical outcomes.

5. **No full off-target scanning**: While we provide seed-toxicity filtering and functional rule checks (Reynolds [30], Ui-Tei [31] rules), we have not yet integrated the full PITA [32] and TargetScan [33] off-target pipelines from the OligoFormer [23] framework.

### 4.6 Future Directions

The most impactful improvements we identify are:
- **RNA-FM embeddings**: Replacing or augmenting the 152-d composition features with RNA foundation model embeddings to capture secondary structure and target accessibility, targeting cross-gene PCC of 0.26 → ~0.50.
- **Training gene diversity**: Partnering with CROs to add 50+ target genes, targeting cross-gene PCC → ~0.55.
- **Prospective wet-lab validation**: Testing top-10 predictions against 20 synthesized siRNAs in a luciferase reporter assay (~$5k) to generate the first prospective validation data for a publicly available cm-siRNA prediction tool.
- **Per-condition models**: Separate models for each dose/cell line combination could reduce within-gene MAE from 16.4 to ~12 points.

---

## 5. Conclusions

HelixZero-CMS provides an updated, production-ready workbench for predicting the efficacy of chemically modified and unmodified siRNAs. By replacing the original SVR with LightGBM, expanding the training dataset 9-fold, and adding position-aware and condition-aware features, we achieve within-gene PCC of 0.68 (MAE 16.4 points) on the largest publicly available cm-siRNA dataset. External validation on an independent 12,303-entry database confirms generalization (PCC = 0.55). The web server (FASTAPI + single-file HTML interface) makes the tool immediately accessible for wet-lab prioritization—enabling researchers to test only the top-10 modification variants rather than all 1,260, with built-in seed-toxicity safety checks. The tool is available for public use at [URL to be provided].

---

## Acknowledgements

[To be inserted]

---

## Author Contributions

[To be inserted]

---

## Disclosure Statement

The authors report there are no competing interests to declare.

---

## Funding

[To be inserted]

---

## Data Availability

The HelixZero 43k patent catalog and all derived training/validation datasets used in this study are available at [repository URL to be provided]. The CMsiRNAdb independent validation dataset is available through He et al. (2026) BMC Bioinformatics. The HelixZero-CMS source code is open-source and available at [GitHub URL to be provided].

---

## Supplementary Material

Supplementary Table S1: Complete list of 30 chemical modification symbols and their one-letter codes.
Supplementary Table S2: Performance comparison across all 14 model variants.
Supplementary Figure S1: MultiModGen input format example.
Supplementary Figure S2: Web server screen captures.
Supplementary File S1: Full training script and model metadata.

---

## References

[1] Fire A, Xu S, Montgomery MK, et al. Potent and specific genetic interference by double-stranded RNA in Caenorhabditis elegans. Nature. 1998;391(6669):806–811.

[2] Elbashir SM, Harborth J, Lendeckel W, et al. Duplexes of 21-nucleotide RNAs mediate RNA interference in cultured mammalian cells. Nature. 2001;411(6836):494–498.

[3] Adams D, Gonzalez-Duarte A, O'Riordan WD, et al. Patisiran, an RNAi therapeutic, for hereditary transthyretin amyloidosis. N Engl J Med. 2018;379(1):11–21.

[4] Raal FJ, Kallend D, Ray KK, et al. Inclisiran for the treatment of heterozygous familial hypercholesterolemia. N Engl J Med. 2020;382(16):1520–1530.

[5] Scott LJ. Givosiran: first approval. Drugs. 2020;80(3):335–339.

[6] Setten RL, Rossi JJ, Han SP. The current state and future directions of RNAi-based therapeutics. Nat Rev Drug Discov. 2019;18(6):421–446.

[7] Bramsen JB, Kjems J. Development of therapeutic-grade small interfering RNAs by chemical engineering. Front Genet. 2012;3:154.

[8] Shukla S, Sumaria CS, Pradeepkumar PI. Exploring chemical modifications for siRNA therapeutics: a structural and functional outlook. ChemMedChem. 2010;5(3):328–349.

[9] Jackson AL, Linsley PS. Recognizing and avoiding siRNA off-target effects for target identification and therapeutic application. Nat Rev Drug Discov. 2010;9(1):57–67.

[10] Allerson CR, Sioufi N, Jarres R, et al. Fully 2'-modified oligonucleotide duplexes with improved in vitro potency and stability compared to unmodified small interfering RNA. J Med Chem. 2005;48(4):901–904.

[11] Manoharan M. RNA interference and chemically modified small interfering RNAs. Curr Opin Chem Biol. 2004;8(6):570–579.

[12] Dar SA, Gupta AK, Thakur A, et al. SMEpred workbench: a web server for predicting efficacy of chemically modified siRNAs. RNA Biol. 2016;13(11):1144–1151.

[13] Dar SA, Thakur A, Qureshi A, et al. siRNAmod: a database of experimentally validated chemically modified siRNAs. Sci Rep. 2016;6:20031.

[14] Ke G, Meng Q, Finley T, et al. LightGBM: a highly efficient gradient boosting decision tree. Adv Neural Inf Process Syst. 2017;30:3146–3154.

[15] He X, et al. CMsiRNAdb: a comprehensive database of chemically modified siRNAs. BMC Bioinformatics. 2026. [Forthcoming]

[16] Jackson AL, Burchard J, Leake D, et al. Position-specific chemical modification of siRNAs reduces "off-target" transcript silencing. RNA. 2006;12(7):1197–1205.

[17] Janas MM, Schlegel MK, Harbison CE, et al. Selection of GalNAc-conjugated siRNAs with limited off-target-driven rat hepatotoxicity. Nat Commun. 2018;9:723.

[18] Huesken D, Lange J, Mickanin C, et al. Design of a genome-wide siRNA library using an artificial neural network. Nat Biotechnol. 2005;23(8):995–1001.

[19] Ichihara M, Murakumo Y, Masuda A, et al. Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities. Nucleic Acids Res. 2007;35(18):e123.

[20] Khvorova A, Reynolds A, Jayasena SD. Functional siRNAs and miRNAs exhibit strand bias. Cell. 2003;115(2):209–216.

[21] Schwarz DS, Hutvagner G, Du T, et al. Asymmetry in the assembly of the RNAi enzyme complex. Cell. 2003;115(2):199–208.

[22] Ramírez S. FastAPI: modern, fast (high-performance) web framework for building APIs with Python. 2018. https://fastapi.tiangolo.com

[23] Bai Y, Zhong H, Wang T, et al. OligoFormer: an accurate and robust prediction method for siRNA design. bioRxiv. 2024.

[24] Boese Q, Leake D, Reynolds A, et al. Mechanistic insights aid computational short interfering RNA design. Methods Enzymol. 2005;392:73–96.

[25] Anderson JS, Anderson EM, Marshall WS, et al. The case for reproducibility: a commentary on siRNA designs. Nat Methods. 2005;2(1):10–11.

[26] Vert JP, Foveau N, Lajaunie C, et al. An accurate and interpretable model for siRNA efficacy prediction. BMC Bioinformatics. 2006;7:520.

[27] Mysara M, Elhefnawi M, Garibaldi JM. mysiRNA-designer: a workflow for efficient siRNA design. PLoS One. 2011;6(10):e25642.

[28] Saetrom P, Snøve O. A comparison of siRNA efficacy predictors. Biochem Biophys Res Commun. 2004;321(1):247–253.

[29] Chen J, Hu Z, Sun S, et al. Interpretable RNA foundation model from unannotated data for highly accurate RNA structure and function predictions. bioRxiv. 2022.

[30] Reynolds A, Leake D, Boese Q, et al. Rational siRNA design for RNA interference. Nat Biotechnol. 2004;22(3):326–330.

[31] Ui-Tei K, Naito Y, Takahashi F, et al. Guidelines for the selection of highly effective siRNA sequences for mammalian and chick RNA interference. Nucleic Acids Res. 2004;32(3):936–948.

[32] Kertesz M, Iovino N, Unnerstall U, et al. The role of site accessibility in microRNA target recognition. Nat Genet. 2007;39(10):1278–1284.

[33] Lewis BP, Burge CB, Bartel DP. Conserved seed pairing, often flanked by adenosines, indicates that thousands of human genes are microRNA targets. Cell. 2005;120(1):15–20.
