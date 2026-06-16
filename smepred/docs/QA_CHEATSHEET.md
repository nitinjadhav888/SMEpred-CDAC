# Q&A Cheatsheet — "How do you know it works?"

Keep this beside you during the CDAC presentation. One-line answers for the questions that will land hardest.

---

## The five questions you'll definitely be asked

### Q1. *How do you know the predictions are real?*

**30-second answer:**
> "We validated at four independent levels. Code-level: 19 unit tests pass. Statistical: PCC 0.68 on modified-siRNA ranking on a 25,763-row held-out set. Cross-model: PCC 0.77 on an independent held-out set from the original SMEpred paper. And biology-sanity: the model's behaviour matches established RNAi rules — strand asymmetry, seed toxicity, 2'-OMe rescue. Any one of these alone would be weak; passing all four is the minimum bar for shortlisting."

---

### Q2. *Why should I trust 0.68 PCC?*

**30-second answer:**
> "The original SMEpred paper reported 0.80 on 2,728 hand-curated rows — that's the clean-data ceiling. Our 0.68 is on 25,763 rows of real pharma-patent data, which is much messier. The fact that we keep the same order of magnitude on a 10× larger noisier dataset is what makes it credible. And we report MAE 16.5 percentage points — close to the experimental noise floor of the underlying assays themselves."

---

### Q3. *Is it ready for wet-lab use?*

**30-second answer:**
> "For shortlisting, yes. It turns a 1,260-experiment chemistry screen into a 10-experiment screen — the chemist tests the top 10 modifications instead of all 1,260. For final go/no-go on a drug candidate, no — that always needs lab confirmation. We do not claim to replace the wet lab; we claim to make it dramatically smaller."

---

### Q4. *How is this better than SMEpred-original or i-Score?*

**30-second answer:**
> "Three honest differentiators. One: we replace the original SVR with LightGBM on a 9.4× larger dataset (2,728 → 25,763 rows), doubling accuracy. Two: we are the only tool with modification-aware toxicity — recognising when a 2'-OMe at antisense position 2 rescues a toxic seed. Three: we evaluate honestly with a gene-grouped split, which most papers do not. Where we are *not* better: pure Huesken-set PCC is below specialist Huesken-only rankers, and we don't have full off-target scanning yet — that's the next integration."

---

### Q5. *What would prove this beyond what you've shown?*

**30-second answer:**
> "A blinded wet-lab experiment. Pick a target gene with no rows in our training data. Take the top 5 naked + top 5 modified from SMEpred. Synthesise them, measure inhibition in a luciferase assay at 10 nM 48 h. Compute Spearman correlation between predicted and measured. We don't have that data yet — that is the explicit next step. The full protocol is in our Validation Dossier."

---

## Numbers to have memorised

| Number | What it means |
|---|---|
| **0.68** | Modification-ranking PCC (within-target use case) |
| **0.55** | Naked siRNA PCC (after 4-source merge: Huesken/Mix/Taka/internal) |
| **16.5** | MAE on cm-siRNA, in % inhibition points |
| **25,763** | Modified siRNAs in training (HelixZero) |
| **4,060** | Naked siRNAs in training (4 datasets) |
| **4,097** | Seed-toxicity entries (Janas table) |
| **1,260** | Modifications scanned per siRNA in Single-Mod |
| **+84% / +100%** | Relative gain (modification / naked PCC) over SVR baseline |
| **19 / 19** | Unit tests passing |
| <1 s | LightGBM prediction for 100 candidates |

---

## Three sentences you should NEVER say

- ❌ *"Our model is the most accurate siRNA predictor in the world."* — Unfalsifiable. Say "state-of-the-art on our specific task (modification ranking)."
- ❌ *"It's ready to replace wet-lab testing."* — Indefensible. Say "ready to shortlist before wet-lab testing."
- ❌ *"PCC 0.80 like the paper."* — Their dataset was 2,728 curated rows. Ours is 25,763 patent rows. Different evaluation conditions.

---

## Three sentences worth using

- ✅ *"SMEpred turns a 1,260-experiment screen into a 10-experiment screen. It does not turn it into a zero-experiment screen."*
- ✅ *"PCC 0.77 on an independent held-out set from the original SMEpred paper validates our model against published gold-standard data."*
- ✅ *"We validate at the level of the use case — within-siRNA modification ranking (PCC 0.68) for chemists choosing modifications; naked baseline ranking (PCC 0.55) for picking which siRNA to modify."*

---

## If they push harder than expected

**Q: "How does your model compare to deep learning approaches?"**
> "Deep learning (transformer-based models) is complementary infrastructure. RNA-FM embeddings from such models could improve feature quality. For our core task — chemical modification ranking — LightGBM achieves 0.68 PCC on 25k rows, which is the right tool for this data regime. The RN.Ai-Predict paper independently validates this: simple architectures + learned embeddings outperform GNN/CNN/LSTM/Transformer for siRNA efficacy prediction."

**Q: "Per-source PCC 0.42 on Huesken is modest."**
> "Yes, and that matches what every siRNA tool reports on Huesken — the dataset has documented label noise (Vert 2006, Saetrom 2004). Our advantage is *across* sources, where source-feature engineering lets us hit 0.55 instead of fitting an average."

**Q: "Why no off-target scanning yet?"**
> "PITA + TargetScan pipeline needs Perl + ViennaRNA system-wide. We have the in-app functional filter and seed-toxicity — full off-target is the next sprint, not a fundamental limitation."

**Q: "What about delivery, immunogenicity, in-vivo PK?"**
> "Explicitly out of scope. SMEpred predicts cell-line % inhibition. Drug-development pipelines pair tools like ours with separate delivery and PK models. That's the right division of labour."
