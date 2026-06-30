# HelixZero-CMS: Architecture & Workflow Overview

## What is HelixZero-CMS?
HelixZero-CMS is a computational pipeline designed to predict the efficacy and safety of siRNA (small interfering RNA) therapeutics. Unlike older models that only predict how well a "naked" (unmodified) sequence works, our model specifically simulates **clinically modified** siRNAs—the exact type of chemically enhanced drugs that actually enter the human body and get approved by the FDA.

## The Core Workflow

Our model processes data through a strict series of scientific filters to mimic human biology:

### 1. Sequence Parsing & Generation
- **What it does:** The system takes a target human gene (mRNA) and slices it into every possible 21-nucleotide segment. 
- **Why it's needed:** An siRNA drug must be exactly 21-23 nucleotides long to fit perfectly into the human cell's machinery (Argonaute/RISC). 
- **The Science:** Dicer is a human enzyme that naturally chops RNA into ~21 base-pair fragments. Our generator mimics Dicer cleavage, focusing specifically on the 3' end.

### 2. Machine Learning Efficacy Scoring (Model A vs Model B)
- **Model A (Naked Baseline):** First, we score all 21-mer sequences to see how well their raw, unmodified genetic code silences the target gene. This creates our `Naked Model` baseline score.
- **Model B (Chemical Awareness):** This is our advanced AI. It doesn't just look at the genetic sequence; it looks at exactly which chemical modifications (like 2'-OMe, 2'-F) are placed at which positions. It recalculates the score based on how these chemicals interact with the RNA.
- **Why we need both:** Comparing Model A vs Model B allows us to calculate the **Delta**—the exact improvement or degradation in efficacy caused by the applied chemistry. 

### 3. The Biophysics Engine
- **What it does:** The AI score is passed through a strict physics engine that penalizes sequences if they violate known laws of human biology. 
- **Why it's needed:** An AI might predict that a certain sequence is 100% effective, but if that sequence would trigger a fatal immune response in a human patient, the drug is useless. The Biophysics Engine grounds the AI in clinical reality.

### 4. Output Generation
The final output is an **Adjusted Score**. This score represents the true clinical viability of the drug after surviving simulated immune system attacks, blood enzymes, and cellular machinery checks.
