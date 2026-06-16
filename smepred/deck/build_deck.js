/* SMEpred pitch deck — 15 slides, CDAC Pune. Run: node build_deck.js */
const pptxgen = require("pptxgenjs");
const React = require("react");
const ReactDOMServer = require("react-dom/server");
const sharp = require("sharp");
const Fa = require("react-icons/fa");

// ─── palette (scientific teal / midnight) ───────────────────────────────────
const C = {
  navy:  "0B2239",   // dark slide bg
  deep:  "0E5A82",   // primary
  teal:  "1C7293",
  mint:  "02C39A",   // accent
  gold:  "F4B860",   // secondary accent (sparingly)
  light: "F2F7FB",   // content bg
  card:  "FFFFFF",
  ink:   "1E293B",   // body text on light
  muted: "64748B",
  white: "FFFFFF",
  redbad:"E2574C",
};
const HF = "Georgia";       // header font
const BF = "Calibri";       // body font

// ─── icon rasteriser ─────────────────────────────────────────────────────────
async function icon(Comp, color, size = 256) {
  const svg = ReactDOMServer.renderToStaticMarkup(
    React.createElement(Comp, { color, size: String(size) }));
  const png = await sharp(Buffer.from(svg)).png().toBuffer();
  return "image/png;base64," + png.toString("base64");
}
const shadow = () => ({ type: "outer", color: "0B2239", blur: 9, offset: 3, angle: 135, opacity: 0.18 });

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";          // 13.3 x 7.5
const W = 13.3, H = 7.5;
pres.author = "SMEpred";
pres.title = "SMEpred — AI siRNA Efficacy Predictor";

// ─── shared slide furniture ────────────────────────────────────────────────
function footer(slide, n) {
  slide.addShape(pres.shapes.RECTANGLE, { x: 0, y: H - 0.32, w: W, h: 0.32, fill: { color: C.navy } });
  slide.addText("SMEpred  ·  siRNA Efficacy Predictor", { x: 0.5, y: H - 0.32, w: 8, h: 0.32, fontFace: BF, fontSize: 9, color: "AEC6D6", valign: "middle", margin: 0 });
  slide.addText(`${n} / 15   ·   CDAC Pune`, { x: W - 4.5, y: H - 0.32, w: 4, h: 0.32, fontFace: BF, fontSize: 9, color: "AEC6D6", align: "right", valign: "middle", margin: 0 });
}
// content-slide header with icon-in-circle motif
function header(slide, iconImg, title, kicker) {
  slide.addShape(pres.shapes.OVAL, { x: 0.5, y: 0.45, w: 0.7, h: 0.7, fill: { color: C.deep } });
  slide.addImage({ data: iconImg, x: 0.63, y: 0.58, w: 0.44, h: 0.44 });
  if (kicker) slide.addText(kicker.toUpperCase(), { x: 1.4, y: 0.42, w: 11, h: 0.3, fontFace: BF, fontSize: 11, color: C.mint, bold: true, charSpacing: 2, margin: 0 });
  slide.addText(title, { x: 1.4, y: kicker ? 0.68 : 0.55, w: 11.4, h: 0.6, fontFace: HF, fontSize: 28, color: C.navy, bold: true, margin: 0 });
}
function card(slide, x, y, w, h, fill = C.card) {
  slide.addShape(pres.shapes.RECTANGLE, { x, y, w, h, fill: { color: fill }, line: { color: "DCE7F0", width: 1 }, shadow: shadow() });
}

(async () => {
  const I = {
    dna:    await icon(Fa.FaDna, "#FFFFFF"),
    bug:    await icon(Fa.FaBug, "#FFFFFF"),
    target: await icon(Fa.FaCrosshairs, "#FFFFFF"),
    flask:  await icon(Fa.FaFlask, "#FFFFFF"),
    flow:   await icon(Fa.FaProjectDiagram, "#FFFFFF"),
    db:     await icon(Fa.FaDatabase, "#FFFFFF"),
    cogs:   await icon(Fa.FaCogs, "#FFFFFF"),
    brain:  await icon(Fa.FaBrain, "#FFFFFF"),
    chart:  await icon(Fa.FaChartLine, "#FFFFFF"),
    laptop: await icon(Fa.FaLaptopCode, "#FFFFFF"),
    pills:  await icon(Fa.FaPills, "#FFFFFF"),
    layers: await icon(Fa.FaLayerGroup, "#FFFFFF"),
    rocket: await icon(Fa.FaRocket, "#FFFFFF"),
    check:  await icon(Fa.FaCheckCircle, "#02C39A"),
    cross:  await icon(Fa.FaTimesCircle, "#E2574C"),
    arrow:  await icon(Fa.FaArrowRight, "#0E5A82"),
    micro:  await icon(Fa.FaMicroscope, "#FFFFFF"),
    seed:   await icon(Fa.FaSeedling, "#FFFFFF"),
  };

  // ════════════════ SLIDE 1 — TITLE ════════════════
  let s = pres.addSlide();
  s.background = { color: C.navy };
  // decorative helix dots
  for (let i = 0; i < 14; i++) {
    const yy = 0.5 + i * 0.48;
    s.addShape(pres.shapes.OVAL, { x: 11.4 + Math.sin(i)*0.5, y: yy, w: 0.16, h: 0.16, fill: { color: i%2? C.mint : C.teal } });
    s.addShape(pres.shapes.OVAL, { x: 12.4 - Math.sin(i)*0.5, y: yy, w: 0.16, h: 0.16, fill: { color: i%2? C.teal : C.mint } });
  }
  s.addShape(pres.shapes.OVAL, { x: 0.9, y: 1.9, w: 1.5, h: 1.5, fill: { color: C.deep } });
  s.addImage({ data: I.dna, x: 1.27, y: 2.27, w: 0.76, h: 0.76 });
  s.addText("SMEpred", { x: 0.8, y: 3.6, w: 10, h: 1.1, fontFace: HF, fontSize: 60, color: C.white, bold: true, margin: 0 });
  s.addText("AI-Powered siRNA Efficacy & Chemical-Modification Predictor", { x: 0.85, y: 4.7, w: 11, h: 0.6, fontFace: BF, fontSize: 20, color: C.mint, margin: 0 });
  s.addText("Predicting which gene-silencing drugs will work — before the lab.", { x: 0.85, y: 5.3, w: 11, h: 0.5, fontFace: BF, fontSize: 14, color: "AEC6D6", italic: true, margin: 0 });
  s.addText([
    { text: "Internship Project  ·  ", options: { color: "AEC6D6" } },
    { text: "C-DAC Pune", options: { color: C.white, bold: true } },
  ], { x: 0.85, y: 6.3, w: 11, h: 0.4, fontFace: BF, fontSize: 13, margin: 0 });

  // ════════════════ SLIDE 2 — THE PROBLEM ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.bug, "The Problem: most diseases start with a faulty gene", "Why this matters");
  const probText = [
    { t: "Genes → proteins → disease", d: "A faulty gene makes a harmful protein. Cancer, high cholesterol, viral infections — many trace back to one over-active gene." },
    { t: "siRNA can switch it off", d: "Small interfering RNA is a 21-letter molecule that tells the cell to destroy a specific gene's message (mRNA). A programmable 'off-switch'." },
    { t: "But raw siRNA fails as a drug", d: "Enzymes in blood shred it in minutes, it triggers immune alarms, and never reaches the target. It needs chemical armor." },
  ];
  probText.forEach((p, i) => {
    const y = 1.6 + i * 1.65;
    card(s, 0.6, y, 7.3, 1.45);
    s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y, w: 0.1, h: 1.45, fill: { color: C.teal } });
    s.addText(p.t, { x: 0.95, y: y + 0.18, w: 6.7, h: 0.4, fontFace: HF, fontSize: 17, bold: true, color: C.deep, margin: 0 });
    s.addText(p.d, { x: 0.95, y: y + 0.62, w: 6.8, h: 0.7, fontFace: BF, fontSize: 12.5, color: C.ink, margin: 0 });
  });
  // right callout
  card(s, 8.3, 1.6, 4.4, 5.0, C.navy);
  s.addText("The expensive bottleneck", { x: 8.6, y: 1.9, w: 3.9, h: 0.5, fontFace: HF, fontSize: 18, bold: true, color: C.white, margin: 0 });
  s.addText("There are", { x: 8.6, y: 2.7, w: 3.9, h: 0.3, fontFace: BF, fontSize: 13, color: "AEC6D6", margin: 0 });
  s.addText("1,260", { x: 8.6, y: 3.0, w: 3.9, h: 0.9, fontFace: HF, fontSize: 54, bold: true, color: C.mint, margin: 0 });
  s.addText("possible chemical modifications for a single siRNA.", { x: 8.6, y: 3.95, w: 3.9, h: 0.7, fontFace: BF, fontSize: 13, color: C.white, margin: 0 });
  s.addText("Testing each one in a lab costs time & money. SMEpred predicts the winners first.", { x: 8.6, y: 4.8, w: 3.9, h: 1.3, fontFace: BF, fontSize: 13, color: "CFE0EC", italic: true, margin: 0 });
  footer(s, 2);

  // ════════════════ SLIDE 3 — REAL-LIFE ANALOGY ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.target, "Think of it like a guided missile", "Real-life analogy");
  const anal = [
    { ic: I.target, t: "The guide strand = GPS", d: "The antisense strand is the targeting system. It finds the exact mRNA to destroy, ignoring the other 20,000 genes." },
    { ic: I.flask, t: "Modifications = armor plating", d: "Chemical tweaks (2'-F, 2'-OMe...) protect the missile from being destroyed mid-flight by enzymes." },
    { ic: I.brain, t: "SMEpred = the war-room simulator", d: "Before you build 1,260 missiles, the simulator predicts which armor design hits the target hardest." },
  ];
  anal.forEach((a, i) => {
    const x = 0.6 + i * 4.15;
    card(s, x, 1.7, 3.85, 4.0);
    s.addShape(pres.shapes.OVAL, { x: x + 1.4, y: 2.0, w: 1.05, h: 1.05, fill: { color: i===2? C.mint : C.deep } });
    s.addImage({ data: a.ic, x: x + 1.68, y: 2.28, w: 0.5, h: 0.5 });
    s.addText(a.t, { x: x + 0.25, y: 3.25, w: 3.35, h: 0.7, fontFace: HF, fontSize: 16, bold: true, color: C.navy, align: "center", margin: 0 });
    s.addText(a.d, { x: x + 0.3, y: 4.0, w: 3.25, h: 1.5, fontFace: BF, fontSize: 12.5, color: C.ink, align: "center", margin: 0 });
  });
  s.addText("\"Don't build 1,260 prototypes. Simulate them, then build the best 3.\"", { x: 1, y: 6.1, w: 11.3, h: 0.5, fontFace: BF, fontSize: 15, italic: true, bold: true, color: C.deep, align: "center", margin: 0 });
  footer(s, 3);

  // ════════════════ SLIDE 4 — WHAT SMEPRED DOES ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.dna, "What SMEpred does, in one breath", "The solution");
  card(s, 0.6, 1.7, 12.1, 1.5, C.deep);
  s.addText("You give it a gene. It finds every possible siRNA, predicts how well each silences the gene, then for your chosen one it predicts which chemical modification makes the best drug. Two ML models (our LightGBM + OligoFormer transformer) ensemble for cross-validated picks, with built-in toxicity & off-target safety filters.", { x: 0.95, y: 1.85, w: 11.4, h: 1.2, fontFace: BF, fontSize: 14, color: C.white, valign: "middle", margin: 0 });
  const two = [
    { ic: I.seed, t: "Naked efficacy", d: "How well does the unmodified siRNA silence the gene? (Rank tab)", col: C.teal },
    { ic: I.flask, t: "Modified efficacy", d: "How much does a chemical modification improve or hurt it? (Single & Multi-Mod tabs)", col: C.mint },
  ];
  two.forEach((p, i) => {
    const x = 0.6 + i * 6.25;
    card(s, x, 3.5, 5.85, 3.1);
    s.addShape(pres.shapes.RECTANGLE, { x, y: 3.5, w: 5.85, h: 0.12, fill: { color: p.col } });
    s.addShape(pres.shapes.OVAL, { x: x + 0.35, y: 3.85, w: 0.95, h: 0.95, fill: { color: p.col } });
    s.addImage({ data: p.ic, x: x + 0.62, y: 4.12, w: 0.42, h: 0.42 });
    s.addText(p.t, { x: x + 1.5, y: 3.95, w: 4.1, h: 0.6, fontFace: HF, fontSize: 22, bold: true, color: C.navy, margin: 0 });
    s.addText(p.d, { x: x + 0.4, y: 5.0, w: 5.1, h: 1.4, fontFace: BF, fontSize: 14, color: C.ink, margin: 0 });
  });
  footer(s, 4);

  // ════════════════ SLIDE 5 — THE WORKFLOW FUNNEL ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.flow, "One workflow, start to finish (a funnel)", "How it flows");
  const steps = [
    { n: "1", t: "Paste a gene", d: "mRNA / gene sequence", col: C.teal },
    { n: "2", t: "RANK", d: "Score every 21-mer siRNA → pick the best naked candidate", col: C.deep },
    { n: "3", t: "SINGLE-MOD", d: "Scan all 1,260 modifications → find the best single tweak", col: C.teal },
    { n: "4", t: "MULTI-MOD", d: "Score your final custom drug design", col: C.mint },
  ];
  steps.forEach((st, i) => {
    const x = 0.55 + i * 3.18;
    card(s, x, 2.2, 2.85, 2.9);
    s.addShape(pres.shapes.OVAL, { x: x + 1.0, y: 2.45, w: 0.85, h: 0.85, fill: { color: st.col } });
    s.addText(st.n, { x: x + 1.0, y: 2.45, w: 0.85, h: 0.85, fontFace: HF, fontSize: 30, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    s.addText(st.t, { x: x + 0.15, y: 3.45, w: 2.55, h: 0.5, fontFace: HF, fontSize: 16, bold: true, color: C.navy, align: "center", margin: 0 });
    s.addText(st.d, { x: x + 0.2, y: 3.95, w: 2.45, h: 1.1, fontFace: BF, fontSize: 11.5, color: C.ink, align: "center", margin: 0 });
    if (i < 3) s.addImage({ data: I.arrow, x: x + 2.92, y: 3.45, w: 0.32, h: 0.32 });
  });
  card(s, 0.55, 5.5, 12.15, 1.05, C.navy);
  s.addText([
    { text: "The fix we shipped:  ", options: { bold: true, color: C.mint } },
    { text: "each ranked result now has a ", options: { color: C.white } },
    { text: "\"Modify →\"", options: { bold: true, color: C.gold } },
    { text: " button that auto-fills the next stage — no re-typing sequences.", options: { color: C.white } },
  ], { x: 0.9, y: 5.5, w: 11.5, h: 1.05, fontFace: BF, fontSize: 14, valign: "middle", margin: 0 });
  footer(s, 5);

  // ════════════════ SLIDE 6 — SYSTEM ARCHITECTURE ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.cogs, "System architecture (end-to-end)", "Under the hood");
  const arch = [
    { ic: I.db, t: "Data", d: "HelixZero 43k +\nHuesken/Mix/Taka\nparse · clean · dedup", col: C.teal },
    { ic: I.layers, t: "Features", d: "152-d vector:\ncomposition · position\n· assay condition", col: C.deep },
    { ic: I.brain, t: "Model A", d: "LightGBM\ngradient-boosted\n% inhibition", col: C.teal },
    { ic: I.micro, t: "Model B", d: "OligoFormer\ntransformer\n+RNA-FM (re-rank)", col: C.deep },
    { ic: I.laptop, t: "Serve", d: "FastAPI +\nHTML app\nensemble in seconds", col: C.mint },
  ];
  arch.forEach((a, i) => {
    const x = 0.45 + i * 2.55;
    card(s, x, 2.3, 2.35, 3.2);
    s.addShape(pres.shapes.OVAL, { x: x + 0.7, y: 2.6, w: 0.9, h: 0.9, fill: { color: a.col } });
    s.addImage({ data: a.ic, x: x + 0.95, y: 2.85, w: 0.4, h: 0.4 });
    s.addText(a.t, { x: x + 0.1, y: 3.65, w: 2.15, h: 0.45, fontFace: HF, fontSize: 16, bold: true, color: C.navy, align: "center", margin: 0 });
    s.addText(a.d, { x: x + 0.15, y: 4.1, w: 2.05, h: 1.3, fontFace: BF, fontSize: 11, color: C.ink, align: "center", margin: 0 });
    if (i < 4) s.addImage({ data: I.arrow, x: x + 2.40, y: 3.55, w: 0.22, h: 0.22 });
  });
  s.addText("Python 3.11  ·  LightGBM · PyTorch (OligoFormer + RNA-FM) · scikit-learn · pandas  ·  FastAPI", { x: 0.6, y: 6.0, w: 12, h: 0.5, fontFace: BF, fontSize: 13, italic: true, color: C.muted, align: "center", margin: 0 });
  footer(s, 6);

  // ════════════════ SLIDE 7 — THE DATA ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.db, "Real data — not synthetic guesses", "The fuel");
  const stats = [
    { big: "25,763", lab: "modified siRNAs\n(HelixZero patent catalog)" },
    { big: "4,060", lab: "naked siRNAs\n(Huesken + Mix + Taka + ours)" },
    { big: "4,097", lab: "seed-toxicity entries\n(OligoFormer Janas et al.)" },
  ];
  stats.forEach((p, i) => {
    const x = 0.6 + i * 4.13;
    card(s, x, 1.75, 3.85, 2.1);
    s.addText(p.big, { x: x + 0.2, y: 1.95, w: 3.45, h: 0.85, fontFace: HF, fontSize: 42, bold: true, color: C.deep, align: "center", margin: 0 });
    s.addText(p.lab, { x: x + 0.25, y: 2.85, w: 3.35, h: 0.9, fontFace: BF, fontSize: 12.5, color: C.ink, align: "center", margin: 0 });
  });
  card(s, 0.6, 4.15, 12.1, 2.45, C.navy);
  s.addText("Why real, multi-source data wins", { x: 0.95, y: 4.35, w: 11, h: 0.5, fontFace: HF, fontSize: 18, bold: true, color: C.mint, margin: 0 });
  s.addText([
    { text: "Modified siRNAs: 25,763 experimentally measured % inhibitions from pharma patents.  Naked siRNAs: 4,060 from four published lab datasets including Huesken (the field's gold standard).  ", options: { color: C.white, breakLine: true } },
    { text: "We feed dataset SOURCE as a model feature so per-lab biases are learned, not averaged away. Result: naked PCC jumped from 0.32 → 0.42 after the multi-source merge.", options: { color: "CFE0EC", italic: true } },
  ], { x: 0.95, y: 4.9, w: 11.4, h: 1.6, fontFace: BF, fontSize: 13, margin: 0 });
  footer(s, 7);

  // ════════════════ SLIDE 8 — FEATURE ENGINEERING ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.layers, "Turning RNA letters into numbers", "Feature engineering");
  s.addText("A machine learning model needs numbers. Each siRNA becomes one 152-number fingerprint:", { x: 0.6, y: 1.55, w: 12, h: 0.4, fontFace: BF, fontSize: 14, color: C.ink, margin: 0 });
  const feats = [
    { d: "140", t: "Composition", s: "How many of each base & modification on both strands (base + modified)" },
    { d: "8", t: "Position density", s: "WHERE mods sit: seed region (pos 1–8), 3' tail, overall" },
    { d: "2", t: "GC content", s: "Duplex stability proxy for each strand" },
    { d: "2", t: "Assay condition", s: "Dose (nM) & time (h) — the hidden confound, now a feature" },
  ];
  feats.forEach((f, i) => {
    const x = 0.6 + (i % 2) * 6.1;
    const y = 2.15 + Math.floor(i / 2) * 1.7;
    card(s, x, y, 5.7, 1.5);
    s.addShape(pres.shapes.OVAL, { x: x + 0.25, y: y + 0.32, w: 0.9, h: 0.9, fill: { color: i===3? C.mint : C.deep } });
    s.addText(f.d, { x: x + 0.25, y: y + 0.32, w: 0.9, h: 0.9, fontFace: HF, fontSize: 22, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
    s.addText(f.t, { x: x + 1.35, y: y + 0.22, w: 4.2, h: 0.4, fontFace: HF, fontSize: 16, bold: true, color: C.navy, margin: 0 });
    s.addText(f.s, { x: x + 1.35, y: y + 0.62, w: 4.2, h: 0.8, fontFace: BF, fontSize: 12, color: C.ink, margin: 0 });
  });
  s.addText([
    { text: "Key insight:  ", options: { bold: true, color: C.deep } },
    { text: "the same sequence shows different inhibition at 1 nM vs 100 nM. Feeding dose+time as features let us keep ALL the data instead of throwing it away.", options: { color: C.ink } },
  ], { x: 0.6, y: 5.7, w: 12, h: 0.7, fontFace: BF, fontSize: 13.5, italic: true, align: "center", margin: 0 });
  footer(s, 8);

  // ════════════════ SLIDE 9 — SMEPRED + OLIGOFORMER ENSEMBLE ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.brain, "Two models, one verdict — ensemble re-ranking", "Cross-validation in silico");
  // Two complementary model cards
  card(s, 0.6, 1.75, 5.85, 3.2);
  s.addShape(pres.shapes.RECTANGLE, { x: 0.6, y: 1.75, w: 5.85, h: 0.65, fill: { color: C.teal } });
  s.addText("Our LightGBM (this project)", { x: 0.6, y: 1.75, w: 5.85, h: 0.65, fontFace: HF, fontSize: 16, bold: true, color: C.white, align: "center", valign: "middle", margin: 0 });
  const mineRows = [
    "Trained on 25,763 modified + 4,060 naked siRNAs",
    "Handles chemical modifications (the cm-siRNA tabs)",
    "Predicts % inhibition directly (regression)",
    "Fast: scores 1,000 candidates in <100 ms",
  ];
  mineRows.forEach((t, i) => {
    s.addImage({ data: I.check, x: 0.85, y: 2.55 + i * 0.56, w: 0.26, h: 0.26 });
    s.addText(t, { x: 1.22, y: 2.48 + i * 0.56, w: 5.0, h: 0.55, fontFace: BF, fontSize: 12, color: C.ink, valign: "middle", margin: 0 });
  });

  card(s, 6.85, 1.75, 5.85, 3.2);
  s.addShape(pres.shapes.RECTANGLE, { x: 6.85, y: 1.75, w: 5.85, h: 0.65, fill: { color: C.mint } });
  s.addText("OligoFormer (Bai et al., Tsinghua 2024)", { x: 6.85, y: 1.75, w: 5.85, h: 0.65, fontFace: HF, fontSize: 16, bold: true, color: C.navy, align: "center", valign: "middle", margin: 0 });
  const ofRows = [
    "Transformer + RNA-FM foundation model (1.1 GB)",
    "Naked siRNA only (no modification support)",
    "Predicts activity probability (deep learning)",
    "Slower: re-ranks the top-50 from LightGBM (~11 s)",
  ];
  ofRows.forEach((t, i) => {
    s.addImage({ data: I.check, x: 7.1, y: 2.55 + i * 0.56, w: 0.26, h: 0.26 });
    s.addText(t, { x: 7.47, y: 2.48 + i * 0.56, w: 5.0, h: 0.55, fontFace: BF, fontSize: 12, color: C.ink, valign: "middle", margin: 0 });
  });

  // ensemble strip below
  card(s, 0.6, 5.15, 12.1, 1.4, C.navy);
  s.addText("Ensemble = within-batch percentile-rank average", { x: 0.95, y: 5.3, w: 11.5, h: 0.4, fontFace: HF, fontSize: 16, bold: true, color: C.mint, margin: 0 });
  s.addText([
    { text: "The two models score on different scales (% inhibition vs activity probability). Averaging raw scores would over-weight one; ", options: { color: C.white } },
    { text: "averaging each model's percentile within the batch is calibration-invariant — \"top 5% by both\" surfaces real cross-validated picks.", options: { color: "CFE0EC", italic: true } },
  ], { x: 0.95, y: 5.7, w: 11.5, h: 0.8, fontFace: BF, fontSize: 12.5, valign: "middle", margin: 0 });
  footer(s, 9);

  // ════════════════ SLIDE 10 — TRAINING & HONEST EVAL ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.micro, "Trained honestly — two different tests", "Training & validation");
  const ev = [
    { t: "Gene-grouped split", q: "\"Can it predict a brand-new gene?\"", d: "Whole target genes (AGT, MSTN, PLN) held out — the model never saw them in training. The honest hard test. Used for the Rank tab.", col: C.teal, sc: "0.26" },
    { t: "Random split", q: "\"Can it rank modifications of a known siRNA?\"", d: "The actual job of the Single/Multi-Mod tabs — compare chemical variants of a sequence you already have.", col: C.mint, sc: "0.68" },
  ];
  ev.forEach((e, i) => {
    const x = 0.6 + i * 6.25;
    card(s, x, 1.75, 5.85, 3.4);
    s.addShape(pres.shapes.RECTANGLE, { x, y: 1.75, w: 0.12, h: 3.4, fill: { color: e.col } });
    s.addText(e.t, { x: x + 0.35, y: 1.95, w: 5.3, h: 0.45, fontFace: HF, fontSize: 18, bold: true, color: C.navy, margin: 0 });
    s.addText(e.q, { x: x + 0.35, y: 2.42, w: 5.3, h: 0.45, fontFace: BF, fontSize: 13.5, italic: true, color: C.deep, margin: 0 });
    s.addText(e.d, { x: x + 0.35, y: 2.95, w: 5.3, h: 1.4, fontFace: BF, fontSize: 12.5, color: C.ink, margin: 0 });
    s.addText([{ text: "PCC ", options: { fontSize: 14, color: C.muted } }, { text: e.sc, options: { fontSize: 30, bold: true, color: e.col } }], { x: x + 0.35, y: 4.4, w: 5.3, h: 0.6, fontFace: HF, margin: 0 });
  });
  card(s, 0.6, 5.45, 12.1, 1.1, C.navy);
  s.addText([
    { text: "Why two tests?  ", options: { bold: true, color: C.mint } },
    { text: "With only 13 genes, a naive split lets the model memorise gene-specific motifs and report fake accuracy. Holding out whole genes proves real generalisation.", options: { color: C.white } },
  ], { x: 0.95, y: 5.45, w: 11.5, h: 1.1, fontFace: BF, fontSize: 13.5, valign: "middle", margin: 0 });
  footer(s, 10);

  // ════════════════ SLIDE 11 — RESULTS (chart) ════════════════
  s = pres.addSlide(); s.background = { color: C.navy };
  s.addShape(pres.shapes.OVAL, { x: 0.5, y: 0.45, w: 0.7, h: 0.7, fill: { color: C.mint } });
  s.addImage({ data: await icon(Fa.FaChartLine, "#0B2239"), x: 0.63, y: 0.58, w: 0.44, h: 0.44 });
  s.addText("RESULTS", { x: 1.4, y: 0.42, w: 11, h: 0.3, fontFace: BF, fontSize: 11, color: C.mint, bold: true, charSpacing: 2, margin: 0 });
  s.addText("We nearly doubled accuracy", { x: 1.4, y: 0.68, w: 11, h: 0.6, fontFace: HF, fontSize: 28, bold: true, color: C.white, margin: 0 });
  s.addChart(pres.charts.BAR, [
    { name: "Old SVR (paper baseline)", labels: ["Modification ranking", "Naked siRNA"], values: [0.37, 0.21] },
    { name: "New LightGBM + OligoFormer data", labels: ["Modification ranking", "Naked siRNA"], values: [0.68, 0.42] },
  ], {
    x: 0.7, y: 1.7, w: 7.6, h: 5.0, barDir: "col",
    chartColors: [C.redbad, C.mint],
    chartArea: { fill: { color: C.navy } }, plotArea: { fill: { color: C.navy } },
    catAxisLabelColor: "CFE0EC", valAxisLabelColor: "CFE0EC", catAxisLabelFontSize: 13, valAxisLabelFontSize: 11,
    valGridLine: { color: "1C3B52", size: 0.5 }, catGridLine: { style: "none" },
    valAxisMaxVal: 0.8, valAxisMinVal: 0, valAxisMajorUnit: 0.2, valAxisLabelFormatCode: "0.0",
    showValue: true, dataLabelPosition: "outEnd", dataLabelColor: "FFFFFF", dataLabelFontSize: 13, dataLabelFontBold: true, dataLabelFormatCode: "0.00",
    showLegend: true, legendPos: "b", legendColor: "CFE0EC", legendFontSize: 12,
    showTitle: false,
  });
  card(s, 8.7, 1.9, 4.0, 4.4, C.deep);
  s.addText("+84%", { x: 8.9, y: 2.05, w: 3.6, h: 0.9, fontFace: HF, fontSize: 44, bold: true, color: C.mint, margin: 0 });
  s.addText("modification ranking (0.37 → 0.68)", { x: 8.9, y: 2.95, w: 3.6, h: 0.5, fontFace: BF, fontSize: 12, color: C.white, margin: 0 });
  s.addText("+100%", { x: 8.9, y: 3.55, w: 3.6, h: 0.9, fontFace: HF, fontSize: 44, bold: true, color: C.gold, margin: 0 });
  s.addText("naked siRNA (0.21 → 0.42)\nafter merging Huesken + Mix + Taka", { x: 8.9, y: 4.45, w: 3.6, h: 0.9, fontFace: BF, fontSize: 12, color: C.white, margin: 0 });
  s.addText("MAE on cm-siRNA: 16.5 pts ·  approaching the noise floor.", { x: 8.9, y: 5.55, w: 3.6, h: 0.6, fontFace: BF, fontSize: 11, color: "CFE0EC", italic: true, margin: 0 });
  footer(s, 11);

  // ════════════════ SLIDE 12 — THE WEB APP ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.laptop, "A working web app, not just a notebook", "The product");
  const tabs = [
    { ic: I.micro, t: "Rank siRNAs", d: "Paste a gene → ranked table of the best naked siRNA candidates with efficacy scores." },
    { ic: I.flask, t: "Single-Mod Scan", d: "Pick an siRNA → all 1,260 modifications ranked, each with its gain/loss vs the original (Δ score)." },
    { ic: I.dna, t: "Multi-Mod Design", d: "Design a custom modified drug → get its predicted efficacy instantly." },
    { ic: I.layers, t: "Modifications", d: "Reference of all 30 chemical modification symbols — click to copy." },
  ];
  tabs.forEach((t, i) => {
    const x = 0.6 + (i % 2) * 6.1;
    const y = 1.7 + Math.floor(i / 2) * 1.78;
    card(s, x, y, 5.7, 1.58);
    s.addShape(pres.shapes.RECTANGLE, { x, y, w: 0.12, h: 1.58, fill: { color: i%2? C.mint : C.teal } });
    s.addShape(pres.shapes.OVAL, { x: x + 0.32, y: y + 0.4, w: 0.78, h: 0.78, fill: { color: i%2? C.mint : C.teal } });
    s.addImage({ data: t.ic, x: x + 0.53, y: y + 0.61, w: 0.36, h: 0.36 });
    s.addText(t.t, { x: x + 1.28, y: y + 0.2, w: 4.3, h: 0.45, fontFace: HF, fontSize: 16, bold: true, color: C.deep, margin: 0 });
    s.addText(t.d, { x: x + 1.28, y: y + 0.62, w: 4.3, h: 0.85, fontFace: BF, fontSize: 12, color: C.ink, margin: 0 });
  });
  card(s, 0.6, 5.4, 12.1, 1.15, C.navy);
  s.addText([
    { text: "Live demo result:  ", options: { bold: true, color: C.mint } },
    { text: "Rank → \"Modify →\" → Scan picks 2'-MOE @ antisense pos 8 → +10.6 efficacy. The same scan finds 2'-OMe @ pos 2 ", options: { color: C.white } },
    { text: "rescues a Toxic seed → Mitigated", options: { bold: true, color: C.gold } },
    { text: " — the kind of insight that saves wet-lab months.", options: { color: C.white } },
  ], { x: 0.95, y: 5.4, w: 11.5, h: 1.15, fontFace: BF, fontSize: 13, valign: "middle", margin: 0 });
  footer(s, 12);

  // ════════════════ SLIDE 13 — REAL-LIFE IMPACT ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.pills, "This is already changing medicine", "Real-world impact");
  card(s, 0.6, 1.7, 12.1, 1.6, C.deep);
  s.addText([
    { text: "siRNA drugs are real and FDA-approved.  ", options: { bold: true, color: C.white } },
    { text: "Patisiran (rare nerve disease), Inclisiran (lowers cholesterol via the PCSK9 gene), Givosiran — all are chemically-modified siRNAs. Every one needed exactly the modification-optimisation SMEpred predicts.", options: { color: "E8F1F8" } },
  ], { x: 0.95, y: 1.85, w: 11.4, h: 1.3, fontFace: BF, fontSize: 14.5, valign: "middle", margin: 0 });
  const impact = [
    { ic: I.flask, t: "Fewer lab experiments", d: "Screen 1,260 designs in silico; test only the top few in the wet lab." },
    { ic: I.rocket, t: "Faster drug discovery", d: "Cut months off the design-test-repeat cycle for new gene targets." },
    { ic: I.check, t: "Lower cost", d: "Each avoided synthesis & assay saves significant time and money." },
  ];
  impact.forEach((a, i) => {
    const x = 0.6 + i * 4.13;
    card(s, x, 3.6, 3.85, 2.9);
    s.addShape(pres.shapes.OVAL, { x: x + 1.42, y: 3.9, w: 1.0, h: 1.0, fill: { color: C.teal } });
    s.addImage({ data: a.ic, x: x + 1.69, y: 4.17, w: 0.46, h: 0.46 });
    s.addText(a.t, { x: x + 0.25, y: 5.05, w: 3.35, h: 0.5, fontFace: HF, fontSize: 15, bold: true, color: C.navy, align: "center", margin: 0 });
    s.addText(a.d, { x: x + 0.3, y: 5.55, w: 3.25, h: 0.9, fontFace: BF, fontSize: 12, color: C.ink, align: "center", margin: 0 });
  });
  footer(s, 13);

  // ════════════════ SLIDE 14 — TECH STACK SUMMARY ════════════════
  s = pres.addSlide(); s.background = { color: C.light };
  header(s, I.cogs, "The full toolkit", "Tech stack");
  const stack = [
    { ic: I.brain, t: "ML — own model", d: "LightGBM · scikit-learn · SciPy · NumPy" },
    { ic: I.micro, t: "ML — ensemble", d: "PyTorch · OligoFormer + RNA-FM (vendored)" },
    { ic: I.db, t: "Data", d: "pandas · custom catalog parsers · joblib" },
    { ic: I.laptop, t: "Backend / API", d: "FastAPI · Uvicorn · Pydantic" },
    { ic: I.flow, t: "Frontend", d: "Single HTML · vanilla JS (fetch)" },
    { ic: I.check, t: "Quality", d: "19 unit tests · seeded reproducibility" },
  ];
  stack.forEach((a, i) => {
    const x = 0.6 + (i % 3) * 4.13;
    const y = 1.95 + Math.floor(i / 3) * 2.15;
    card(s, x, y, 3.85, 1.9);
    s.addShape(pres.shapes.OVAL, { x: x + 0.3, y: y + 0.55, w: 0.8, h: 0.8, fill: { color: C.deep } });
    s.addImage({ data: a.ic, x: x + 0.52, y: y + 0.77, w: 0.36, h: 0.36 });
    s.addText(a.t, { x: x + 1.25, y: y + 0.45, w: 2.45, h: 0.5, fontFace: HF, fontSize: 15, bold: true, color: C.navy, margin: 0 });
    s.addText(a.d, { x: x + 1.25, y: y + 0.95, w: 2.5, h: 0.85, fontFace: BF, fontSize: 11.5, color: C.ink, margin: 0 });
  });
  footer(s, 14);

  // ════════════════ SLIDE 15 — CONCLUSION ════════════════
  s = pres.addSlide(); s.background = { color: C.navy };
  for (let i = 0; i < 12; i++) s.addShape(pres.shapes.OVAL, { x: 0.4 + i * 1.08, y: 0.4, w: 0.14, h: 0.14, fill: { color: i%2? C.mint : C.teal } });
  s.addShape(pres.shapes.OVAL, { x: 0.9, y: 1.5, w: 1.3, h: 1.3, fill: { color: C.mint } });
  s.addImage({ data: await icon(Fa.FaDna, "#0B2239"), x: 1.23, y: 1.83, w: 0.64, h: 0.64 });
  s.addText("From a gene to a drug candidate — in seconds.", { x: 0.85, y: 3.05, w: 11.5, h: 1.0, fontFace: HF, fontSize: 32, bold: true, color: C.white, margin: 0 });
  const takeaways = [
    "Predicts BOTH naked and modified siRNA efficacy with built-in toxicity & off-target safety filters",
    "Modification PCC 0.37 → 0.68 (+84%) · Naked PCC 0.21 → 0.42 (+100%) — both via real published data",
    "Two ML models ensemble: our LightGBM + the OligoFormer transformer, cross-validated picks in one click",
  ];
  takeaways.forEach((t, i) => {
    s.addImage({ data: I.check, x: 1.0, y: 4.25 + i * 0.62, w: 0.34, h: 0.34 });
    s.addText(t, { x: 1.5, y: 4.2 + i * 0.62, w: 10.5, h: 0.5, fontFace: BF, fontSize: 15, color: C.white, valign: "middle", margin: 0 });
  });
  s.addText("Next: more target genes to lift new-gene accuracy  ·  experimental wet-lab validation of top picks", { x: 1.0, y: 6.25, w: 11.3, h: 0.5, fontFace: BF, fontSize: 12.5, italic: true, color: "AEC6D6", margin: 0 });
  s.addText("Thank you  —  C-DAC Pune", { x: 0.85, y: 6.75, w: 11.5, h: 0.5, fontFace: HF, fontSize: 18, bold: true, color: C.mint, margin: 0 });

  await pres.writeFile({ fileName: "C:/Helixx/smepred/deck/SMEpred_Pitch_Deck.pptx" });
  console.log("DECK WRITTEN");
})();
