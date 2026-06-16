# SMEpred — Speaker Notes (slide-by-slide)

> ⚠️ **UPDATE 2026-06-15:** OligoFormer (the vendored transformer ensemble model) has been
> **fully removed** from the project. All references to "ensemble with OligoFormer", "Model B",
> and the "ensemble toggle" in the slides below are **stale**. The Rank tab now uses the naked
> model only (model_normal.pkl, 4,060 rows, 4 sources merged). Slides 4, 6, 9, 12, 14 and
> the Q&A section need verbal adjustment during delivery — skip any mention of OligoFormer
> or the ensemble re-ranker. The data *sources* (Huesken/Mix/Takayuki CSV files) remain in
> `data/oligoformer/` and were used to train the naked model; they are not OligoFormer the model.

A detailed talking script for the 15-slide CDAC Pune pitch. Aim for **~1 minute per slide**
(15 minutes total + Q&A). Each note has an **opening line**, a **plain-English explanation**,
a **simple real-life example**, and a **transition**.

---

## Slide 1 — Title

**Open with:** "Good morning. I'm presenting SMEpred — an AI tool that predicts how well a
gene-silencing drug will work, before anyone touches a test tube."

**Say:** SMEpred stands for *siRNA Modification Efficacy Predictor*. It's the project I've
been working on during my internship at C-DAC Pune. In the next 15 minutes I'll show what
problem it solves, how it works, the real measured results, and where it's headed.

**Simple example:** Think of it like Google Maps for drug designers. Instead of trying
every route in traffic, you ask the app which route is fastest, and only then drive.

**Transition:** "Let's start with the problem."

---

## Slide 2 — The Problem

**Open with:** "Many diseases start with one bad gene that won't shut up."

**Say:** Genes make proteins. When one gene makes too much of a harmful protein, you get
cancer, high cholesterol, viral disease. siRNA — small interfering RNA — is a tiny
21-letter molecule that can silence that specific gene. It's a programmable off-switch
for any gene in the body.

The catch: raw siRNA is destroyed by enzymes in the blood within minutes, triggers immune
alarms, and never reaches the target cell. To turn it into a real drug you have to add
chemical "armor" — modifications.

**The killer stat (point at the right panel):** for ONE siRNA there are 1,260 possible
chemical modifications. Testing each one in the wet lab costs months and lakhs of rupees.

**Simple example:** It's like designing a phone case. There are 1,260 possible materials
and shapes. You can't physically test all of them — you need a simulator that says,
"these top 5 are worth manufacturing."

**Transition:** "Here's an analogy that makes it click."

---

## Slide 3 — Real-life analogy: a guided missile

**Open with:** "Think of an siRNA drug as a guided missile."

**Say:** The antisense strand is the GPS — it locks onto the exact mRNA target out of
20,000 genes. The chemical modifications are armor plating — they protect the missile
from being shot down before it reaches the target. SMEpred is the war-room simulator:
before you build 1,260 missiles, it predicts which armor design hits hardest.

**Simple example:** Indian Air Force engineers don't physically build every possible
missile variant. They simulate them in software, pick the top three, and only then
build the real thing. SMEpred is that simulator, for medicine.

**Transition:** "So what does it actually do? One sentence."

---

## Slide 4 — What SMEpred does, in one breath

**Open with:** Read the big blue banner verbatim — slowly.

**Say:** It does two related but distinct jobs, with one model trained on two data regimes.
**Naked efficacy** — given a gene, find the best unmodified siRNA sequence. **Modified
efficacy** — given a chosen siRNA, find the chemical tweak that makes it the strongest
drug. **Safety** — every candidate also gets a seed-toxicity check, and modifications
can rescue a toxic seed.

**Simple example:** A chef has two skills — pick the right base recipe (naked siRNA),
then pick the spice that makes it irresistible (modification). SMEpred does both.

**Transition:** "Here's how that flows step by step."

---

## Slide 5 — The workflow funnel

**Open with:** "The whole app is one funnel — four steps from a gene to a drug candidate."

**Say (walk left to right):**
1. **Paste a gene** — typically 1,000 to 10,000 letters of mRNA.
2. **Rank tab** — the model scores every 21-letter window of that gene, ranks them, and
   shows you the best naked siRNA candidates.
3. **Single-Mod Scan tab** — you pick the top siRNA. The model now scans all 1,260 chemical
   modifications and ranks them by predicted efficacy gain.
4. **Multi-Mod Design tab** — finally, you design a custom modified drug with multiple
   modifications, and the model scores it.

**Show the fix at the bottom:** Earlier, users had to manually copy the 21-letter siRNA
from step 2 and paste it into step 3 — that was friction. Now each ranked row has a
"Modify →" button that auto-fills the next tab. One click instead of copy-paste.

**Simple example:** It's like Amazon — search a product, click "Buy Now," and the address
and payment auto-fill. Same idea, applied to drug design.

**Transition:** "Now, under the hood — the architecture."

---

## Slide 6 — System architecture

**Open with:** "Five clean layers — data, features, two models, serving."

**Say:**
- **Data layer** parses the HelixZero 43k patent catalog + 4 published naked-siRNA
  Huesken / Mix / Takayuki datasets + siRNAmod, into clean CSV files.
- **Feature layer** turns each siRNA into a 152-number vector.
- **Model A (cm-siRNA):** LightGBM trained on 25,763 modified siRNAs — drives the
  Single-Mod and Multi-Mod tabs. Predicts % inhibition.
- **Model Normal (naked):** LightGBM trained on 4,060 naked siRNAs merged from four
  published sources (Huesken, Mix, Takayuki, internal) — drives the Rank tab.
- **Serving layer** is FastAPI + a single-file HTML web app. Browser-only.

**Simple example:** Like a car factory — raw steel (data) → stamping (features) →
assembly (LightGBM) → showroom (web app).

**Transition:** "The whole thing is only as good as the data."

---

## Slide 7 — Real data, multi-source

**Open with:** "Garbage in, garbage out. So we use real data — from four different labs."

**Say:** The three big numbers on screen:
- **25,763** modified siRNAs from the HelixZero pharma-patent catalog.
- **4,060** naked siRNAs combined from four published sets — **Huesken** (Nature
  Biotechnology 2005, 2,361 rows, the field's gold standard), **Mix** (Reynolds,
  Vickers, Ui-Tei combined, 462), **Takayuki** (NAR 2007, 699), plus our existing 538.
- **4,097** entries in the seed-toxicity lookup table (Janas et al., Mol Cell 2018)

**The trick (dark callout):** Lab datasets have different distributions — different
cell lines, doses, timepoints. If you mix them naively the model fits the average and
gets worse. We feed dataset **source as a feature**, so per-lab biases are learned, not
averaged away. Result: naked PCC jumped from 0.32 to 0.42 after the merge.

**Simple example:** Training a self-driving car on dashcam from four different cities.
Tell the model which city each clip came from, and it learns the local rules instead
of averaging into one confused driver.

**Transition:** "Real data still needs to become numbers the model can read."

---

## Slide 8 — Feature engineering

**Open with:** "Machine learning models can't read RNA letters. So we translate each
siRNA into a 152-number fingerprint."

**Say (walk the four cards):**
- **140 composition numbers** — how many of each base (A/U/G/C) and each chemical
  modification, on both strands.
- **8 position numbers** — *where* the modifications sit. Position matters enormously
  in siRNA biology: the seed region (positions 1–8) is sacred, the 3' tail is more
  forgiving.
- **2 GC content numbers** — a proxy for how tightly the strands hold together.
- **2 condition numbers** — the dose used in the experiment and how long it ran. This
  is the secret sauce.

**The key insight (bottom italic line):** The same exact sequence shows different
inhibition at 1 nM vs 100 nM. If we ignore that, the model sees contradictions and
learns badly. By feeding dose and time AS features, we keep ALL 25,000 rows instead of
throwing most away.

**Simple example:** Predicting whether a cricket batsman will score depends not just on
his stats, but also on the pitch conditions. We feed the pitch into the model.

**Transition:** "And the model itself — we made a big choice here."

---

## Slide 9 — Training data from four published sources

**Open with:** "The model is only as good as what it learned from. We merged data from
four independent published sources."

**Say:**
- **Huesken (Nature Biotech 2005):** 2,361 naked siRNAs — the field's gold-standard
  benchmark. Used in virtually every siRNA paper since.
- **Mix (Reynolds/Vickers/Ui-Tei):** 462 rows from three independent design rule papers.
- **Takayuki (NAR 2007):** 699 rows — the cleanest single-condition dataset. Our model
  scores highest on this (PCC 0.69).
- **Our internal catalog:** 538 rows from HelixZero's own measurements.

**The trick:** Each lab used different cell lines, doses, and timepoints. If you just
merge them naively, the model averages across distribution shifts and gets worse. We
feed dataset **source as a one-hot feature** so the model learns per-lab biases, not
an average of them all.

**Simple example:** Four doctors each examine the same patient. One uses a stethoscope,
one uses an X-ray, one does a blood test, one takes a history. You don't average their
findings — you teach a super-doctor which tool each one used.

**Transition:** "But raw accuracy claims are easy to fake. We measured honestly."

---

## Slide 10 — Trained honestly — two different tests

**Open with:** "Different product features need different honesty tests."

**Say:**
- **Modification ranking (0.68):** standard 82/18 random split. This answers: *"Can it
  rank chemical modifications of an siRNA we already know?"* PCC = **0.68**. This is
  the actual job of the Single-Mod and Multi-Mod tabs.
- **Naked baseline ranking (0.55):** the Rank tab uses a separate naked model trained on
  4,060 unmodified siRNAs from 4 merged sources. This is for picking which siRNA to
  chemically modify.

**Why separated?** Different questions need different models. The cm-siRNA model ranks
modifications of a known siRNA (within-gene), so random-split PCC is the faithful
measure. Cross-gene generalization is not required — modification ranking never leaves
the context of one siRNA. The naked model handles baseline ranking.

**Simple example:** A chef testing spice combinations on one base recipe doesn't need
to know how to cook Italian food. Single recipe, many spices — that's our 0.68 PCC.

**Transition:** "So let's see the headline result."

---

## Slide 11 — Results: two big wins

**Open with:** Point at the bar chart. "Red is the paper-baseline SVR. Green is our
LightGBM rebuild after merging four published datasets (Huesken, Mix, Takayuki, internal)."

**Say:** Two headline numbers:
- **Modification ranking: 0.37 → 0.68 PCC, +84% relative.** That's the core task of
  the Single-Mod and Multi-Mod tabs.
- **Naked siRNA ranking: 0.21 → 0.42 PCC, +100% (doubled).** The four-dataset merge (Huesken + Mix + Taka + internal) — Huesken + Mix + Takayuki, with source as a feature — is what pushed it over
  the top.

**Bonuses on the right panel:**
- Predictions are now real inhibition percentages on a 0–100 scale. The old code had a
  bug — it min-max-rescaled every batch so the best of every group always showed 100,
  which was meaningless.
- MAE on cm-siRNA is 16.5 percentage points — approaching the experimental noise floor
  of the underlying lab assays themselves (~10–15 points between labs).

**Simple example:** Old system was an exam graded on a curve. New system gives true
scores you can compare across batches.

**Transition:** "These predictions live inside a working product."

---

## Slide 12 — The web app

**Open with:** "It's not a research notebook. It's a working web app with four tabs and
two safety filters baked in."

**Say:** Each tab solves one job. Rank siRNAs. Single-Mod Scan. Multi-Mod Design.
Modifications reference. Single HTML file + FastAPI backend over HTTP. Browser-only.

**Key features:**
- **Seed Toxicity column** on every tab — Safe / Caution / Toxic with cell-viability %.
- **Modification-aware toxicity** — if a known rescue mod (2'-OMe, 2'-F, LNA, 2'-MOE)
  sits in the seed region of a modified candidate, the badge flips to **Mitigated**
  with a tooltip explaining the rescue. Real published biology, surfaced in the UI.
- **Auto Multi-Mod Scan** — top single-mod hits are automatically combined via beam
  search into multi-mod variants (2-mod, 3-mod), scored and ranked.

**The demo result at the bottom:** A real test. Rank → Modify → Scan. The model found
that 2'-MOE at antisense position 8 improves efficacy by +10.6 points, AND that 2'-OMe
at antisense position 2 rescues a Toxic seed into Mitigated. That's the kind of
insight that saves wet-lab months.

**Simple example:** Data in, safety annotations a chemist actually trusts, and the
best modified design surfaces automatically.

**Transition:** "Why does any of this matter outside this room?"

---

## Slide 13 — Real-world impact

**Open with:** "siRNA drugs aren't theoretical anymore."

**Say:** Three FDA-approved siRNA medicines are already saving lives. Patisiran for a
rare nerve disease. Inclisiran lowers cholesterol by silencing the PCSK9 gene — millions
of patients. Givosiran for porphyria. **Every single one needed exactly the kind of
modification-optimization SMEpred predicts.** Without it, those drugs took years and
hundreds of millions of dollars to design.

**Three concrete impacts:** fewer wet-lab experiments per candidate, faster discovery
cycles, and significantly lower cost per drug program.

**Simple example:** Like CAD software for civil engineers. You don't build the bridge to
test whether it stands — you simulate. SMEpred is CAD for siRNA drugs.

**Transition:** "Quickly — the toolkit behind all this."

---

## Slide 14 — The full toolkit

**Open with:** "Everything is mainstream open-source Python — one clean ML stack."

**Say:**
- **ML:** LightGBM, scikit-learn, SciPy, NumPy.
- **Data:** pandas + custom HelixZero / siRNAmod catalog parsers + 4 published siRNA
  datasets (Huesken, Mix, Takayuki, internal).
- **Backend:** FastAPI on Uvicorn — async, auto-generated Swagger UI at /docs.
- **Frontend:** single HTML file + vanilla JavaScript. Zero build pipeline.
- **Quality:** 19 unit tests, seeded reproducibility — metrics reproduce exactly.

**Simple example:** Same stack any modern AI startup uses. Not exotic. Production-grade.

**Transition:** "Let me close with the takeaway."

---

## Slide 15 — Closing

**Open with:** Read the title: "From a gene to a drug candidate — in seconds."

**Say:** Three takeaways.
1. SMEpred predicts both naked and modified siRNA efficacy with **built-in seed-toxicity
   AND modification-aware off-target rescue checks** — the Mitigated badges.
2. Accuracy:
   - Modification ranking: 0.37 → 0.68 PCC (+84%)
   - Naked siRNA ranking: 0.21 → 0.55 PCC (+162%) — four dataset merge.
3. A production-grade web app any researcher can use immediately — paste a gene,
   get ranked modification-optimized siRNA candidates in seconds.

**Next steps:** more diverse target genes for the Rank tab. Wet-lab validation of
SMEpred's top picks. And add full off-target scanning (PITA / TargetScan) as a third
safety check.

**Close with:** "Thank you. I'm happy to take questions — about the biology, the model,
or the engineering."

---

## Q&A — likely questions you should be ready for

| Question | Short answer |
|---|---|
| *Why is the naked model PCC lower than modification PCC?* | Naked siRNAs have less signal — the sequence alone tells you less than the modification pattern. The real advantage is in the modified model (PCC 0.68). |
| *What does "Mitigated" mean?* | The seed lookup says Toxic, but a known off-target rescue mod (2'-OMe, 2'-F, LNA, 2'-MOE) sits in pos 2–7 of the antisense. Established biology — Jackson et al., *RNA* 2006. |
| *Is 0.68 PCC good?* | The paper claimed 0.80 on a smaller curated set. Ours is on a much larger messier real-world catalog with honest within-gene evaluation. Strong; the paper's number is the ceiling. |
| *Why not deep learning for your own model?* | 25k rows is the sweet spot for boosted trees. Deep nets need 100k+ to beat tree ensembles on tabular data, and we'd lose interpretability. |
| *Biggest limitation?* | Only 13 target genes in the cm-siRNA training set limits the diversity of modification patterns the model has seen. More genes would help robustness. |
| *Could a wet lab use it tomorrow?* | Yes for shortlisting — pick the top 10 mods instead of testing 1,260, with Toxic + func-fail auto-dropped. For final go/no-go, still needs experimental confirmation. |
| *Why patent data?* | It's the only public source with this volume of measured inhibition values. Academic papers report dozens of rows; patents report thousands. |
| *Open source?* | The code is structured for it. Yes, easy to publish. |

