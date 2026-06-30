# Organ Targeting and Delivery: Is the Model Liver-Specific?

A critical question from stakeholders is often: *"Is this model only useful for liver diseases?"*

The answer is **No. The core intelligence of the model applies universally to all human cells.** However, the model currently enforces a specific *delivery* requirement for the liver because that is the most clinically validated route in modern medicine.

Here is the breakdown of what is universal versus what is liver-specific:

## 1. Universal Biology (Organ-Agnostic)
The vast majority of the model's physics engine evaluates biological laws that are identical in every cell in the human body, whether in the brain, lungs, heart, or liver:

* **Efficacy & Gene Silencing:** The physical laws of RNA interference (RISC loading, Argonaute 2 slicing, nearest-neighbor thermodynamics) are universal. A sequence that successfully loads into RISC in a liver cell will also load into RISC in a lung cell.
* **Toxicity & Off-Target Safety:** The engine scans against the *entire* human transcriptome, checking for off-target slicing across all human genes. 
* **Immune Evasion:** The system masks the drug from Toll-Like Receptors (TLR7/8), which are part of the innate immune system present throughout the bloodstream and entire body.
* **Seed Toxicity:** The 4,000+ empirical seed-toxicity experiments the model relies on were tested on generic human cell lines, not strictly liver cells.

## 2. Liver-Specific Delivery (GalNAc Pharmacokinetics)
While the drug's *function* is universal, the physical challenge of *delivering* the drug into the target organ is immense. 

Currently, the model's Pharmacokinetic (PK) check assumes a **hepatic (liver) delivery route**.
* **The Conjugate Requirement:** Because it defaults to the liver, the model strictly checks if you have attached the **GalNAc** conjugate (represented by modification symbol `4`) to the drug.
* **What is GalNAc?** GalNAc (N-Acetylgalactosamine) is a specific sugar molecule that acts like a tracking beacon and a biological key. It binds exclusively to ASGPR receptors, which are found almost entirely on the surface of liver cells (hepatocytes).
* **The Impact:** If you design a brilliant, 100% effective drug but forget to attach GalNAc, the model will apply a heavy pharmacokinetic penalty. Why? Because without that conjugate, the drug will simply float around in the human bloodstream until it is flushed out by the kidneys. It will never enter the liver cells to do its job.

### Summary
The model designs drugs that are universally potent and safe for human biology. It currently enforces GalNAc delivery because liver-targeting (used in FDA-approved drugs like *Inclisiran* and *Givosiran*) is the most successful, clinically proven method for delivering RNA therapeutics today.
