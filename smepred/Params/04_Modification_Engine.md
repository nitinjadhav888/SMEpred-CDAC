# The Modification Engine and Beam Search

This document explains how the model explores billions of chemical combinations to find the perfect drug design.

## 1. Experimental Chemistry Application
- **The Concept:** Natural RNA breaks down instantly in the body. To create a drug, we must replace the natural OH-groups on the RNA backbone with synthetic chemicals like 2'-Fluoro (F), 2'-O-Methyl (M), or Phosphorothioates (S).
- **The Challenge:** An siRNA drug is 42 nucleotides long (21 on each strand). There are over 30 different chemical modifications available. Trying every possible combination would take billions of years of computing time.

## 2. Multi-Mod Beam Search (The AI Workflow)
- **What it is:** "Beam Search" is an advanced AI search algorithm. Instead of trying every combination randomly, it acts like a chess grandmaster planning moves ahead.
- **How it works:** 
  1. It first tests single modifications one at a time across all 42 positions (generating a library of 1,260 variants).
  2. It takes the top 20 best-performing modifications (the "Beam").
  3. It starts stacking them. It takes the best single modification and pairs it with the next best, continually testing combinations.
  4. It throws away combinations that cause steric clashes (e.g., $\ge$ 3 consecutive bulky modifications) and keeps the ones that increase the Model B score while keeping Biophysical Penalties low.
- **The Result:** We get the `Top 30 of 100 beam variants`—the absolute most potent and safest chemical architectures discovered by the AI.

## 3. Highlighting and Terminology
When you look at the model's visual outputs, you will see specific color coding or highlighting:

- **Highlighted = modified:** Any position that has a synthetic chemical applied to it is highlighted. 
- **Seed Region (pos 2-8):** This specific stretch of 7 nucleotides is highlighted because it is the "danger zone". If we apply the wrong chemistry here, the drug becomes toxic (Seed Toxicity). The model explicitly verifies that protective chemistry (like GNA at position 7) is applied here.
- **Terminus:** The very ends of the drug (Position 1 and Position 21). These are highlighted because they are the attack points for blood enzymes. They must be protected with specific "cap" chemistries.

## 4. Why this matters to non-scientists
In the past, scientists had to guess where to put these chemicals and then spend 6 months and $100,000 testing a single configuration in a lab. Our AI tests thousands of highly probable, physically valid combinations in seconds, outputting a drug that is already computationally proven to survive the human bloodstream, hide from the immune system, and perfectly silence the disease gene.
