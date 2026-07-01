import fitz, json, os

results = {}

paper_files = [
    ("1.pdf", "Weingartner 2020 - GalNAc conjugates"),
    ("2.pdf", "Sakamuri 2020 - PS stereochemistry"),
    ("3.pdf", "Weingartner 2020 - GalNAc (duplicate check)"),
    ("ML.pdf", "ML paper"),
    ("OligoFormer.pdf", "OligoFormer"),
    ("Chemical and structural modifications of RNAi therapeutics.pdf", "Chem modifications"),
    ("Design of siRNA Therapeutics from the Molecular Scale.pdf", "Design from molecular scale"),
    ("Thermodynamics.pdf", "Thermodynamics"),
    ("Challenges In SIRNA Design.pdf", "Challenges"),
    ("Chemical modification resolves the asymmetry of siRNA.pdf", "Asymmetry"),
    ("main.pdf", "main"),
]

for fname, desc in paper_files:
    path = f"papers/{fname}"
    if not os.path.exists(path):
        continue
    doc = fitz.open(path)
    text = ""
    for page in doc:
        text += page.get_text()
        if len(text) > 12000:
            break
    doc.close()
    results[desc] = text[:12000]

with open("papers_text.json", "w", encoding="utf-8") as f:
    json.dump(results, f, ensure_ascii=False)
print("Done!", list(results.keys()))
