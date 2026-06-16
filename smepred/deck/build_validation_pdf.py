"""Render docs/VALIDATION.md to a clean A4 PDF for the C-DAC review panel."""
from pathlib import Path
import re
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor
from reportlab.lib.units import inch
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak,
)

SRC = Path(__file__).parent.parent / "docs" / "VALIDATION.md"
OUT = Path(__file__).parent / "SMEpred_Validation_Dossier.pdf"

NAVY  = HexColor("#0B2239")
TEAL  = HexColor("#1C7293")
MINT  = HexColor("#02C39A")
MUTED = HexColor("#64748B")
LIGHT = HexColor("#F2F7FB")

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontSize=18, leading=22,
                    textColor=NAVY, spaceBefore=12, spaceAfter=8)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontSize=13, leading=17,
                    textColor=TEAL, spaceBefore=14, spaceAfter=6)
H3 = ParagraphStyle("H3", parent=ss["Heading3"], fontSize=11, leading=14,
                    textColor=NAVY, spaceBefore=10, spaceAfter=4)
P  = ParagraphStyle("P", parent=ss["BodyText"], fontSize=9.5, leading=13,
                    textColor=NAVY, spaceAfter=5, allowOrphans=1)
BUL = ParagraphStyle("BUL", parent=P, leftIndent=14, bulletIndent=6)
QUOTE = ParagraphStyle("Q", parent=P, leftIndent=14, rightIndent=14,
                       textColor=TEAL, fontName="Helvetica-Oblique")
CODE = ParagraphStyle("C", parent=P, fontName="Courier", fontSize=9, textColor=NAVY,
                      backColor=LIGHT)
TBL_CELL = ParagraphStyle("TBL", parent=P, fontSize=8.5, leading=11,
                          textColor=NAVY, spaceAfter=0)
TBL_HDR  = ParagraphStyle("TBLH", parent=TBL_CELL, textColor=HexColor("#FFFFFF"),
                          fontName="Helvetica-Bold")


def inline_md(text: str) -> str:
    """Minimal markdown inline formatting → ReportLab HTML-ish."""
    text = text.replace("&", "&amp;")
    text = re.sub(r"`([^`]+)`", r'<font face="Courier" color="#0B2239">\1</font>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*\n]+)\*", r"<i>\1</i>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"<u>\1</u>", text)
    return text


def parse_table(block_lines):
    rows = []
    for ln in block_lines:
        if not ln.strip().startswith("|"):
            continue
        if set(ln.replace("|", "").replace(":", "").strip()) <= {"-", " "}:
            continue
        cells = [c.strip() for c in ln.strip().strip("|").split("|")]
        rows.append(cells)
    return rows


def build():
    story = []
    text = SRC.read_text(encoding="utf-8")
    lines = text.split("\n")

    i = 0
    in_code = False
    para_buf = []
    list_buf = []

    def flush_para():
        nonlocal para_buf
        if para_buf:
            joined = " ".join(s.strip() for s in para_buf if s.strip())
            if joined:
                story.append(Paragraph(inline_md(joined), P))
            para_buf = []

    def flush_list():
        nonlocal list_buf
        for item in list_buf:
            story.append(Paragraph("• " + inline_md(item), BUL))
        list_buf = []

    while i < len(lines):
        ln = lines[i]
        # fenced code skip
        if ln.startswith("```"):
            in_code = not in_code
            i += 1
            continue
        if in_code:
            story.append(Paragraph(ln.replace(" ", "&nbsp;"), CODE))
            i += 1
            continue

        s = ln.rstrip()
        # blockquote
        if s.startswith(">"):
            flush_para(); flush_list()
            story.append(Paragraph(inline_md(s.lstrip("> ").strip()), QUOTE))
            i += 1; continue

        # heading
        if s.startswith("# "):
            flush_para(); flush_list()
            story.append(Paragraph(inline_md(s[2:]), H1))
            i += 1; continue
        if s.startswith("## "):
            flush_para(); flush_list()
            story.append(Paragraph(inline_md(s[3:]), H2))
            i += 1; continue
        if s.startswith("### "):
            flush_para(); flush_list()
            story.append(Paragraph(inline_md(s[4:]), H3))
            i += 1; continue

        # table
        if s.startswith("|"):
            flush_para(); flush_list()
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i]); i += 1
            rows = parse_table(block)
            if rows:
                table_data = []
                # header
                table_data.append([Paragraph(inline_md(c), TBL_HDR) for c in rows[0]])
                for r in rows[1:]:
                    if len(r) != len(rows[0]):
                        r = (r + [""] * len(rows[0]))[:len(rows[0])]
                    table_data.append([Paragraph(inline_md(c), TBL_CELL) for c in r])
                avail = A4[0] - 1.4*inch
                col_w = avail / len(rows[0])
                t = Table(table_data, colWidths=[col_w]*len(rows[0]), repeatRows=1)
                t.setStyle(TableStyle([
                    ("BACKGROUND", (0,0), (-1,0), TEAL),
                    ("ROWBACKGROUNDS", (0,1), (-1,-1), [HexColor("#FFFFFF"), LIGHT]),
                    ("GRID", (0,0), (-1,-1), 0.4, MUTED),
                    ("VALIGN", (0,0), (-1,-1), "TOP"),
                    ("LEFTPADDING", (0,0), (-1,-1), 5),
                    ("RIGHTPADDING", (0,0), (-1,-1), 5),
                    ("TOPPADDING", (0,0), (-1,-1), 4),
                    ("BOTTOMPADDING", (0,0), (-1,-1), 4),
                ]))
                story.append(Spacer(1, 4))
                story.append(t)
                story.append(Spacer(1, 4))
            continue

        # list
        m = re.match(r"^\s*[-*]\s+(.*)$", s)
        m_num = re.match(r"^\s*\d+\.\s+(.*)$", s)
        if m or m_num:
            flush_para()
            list_buf.append((m or m_num).group(1))
            i += 1; continue

        # horizontal rule
        if re.match(r"^-{3,}$", s):
            flush_para(); flush_list()
            story.append(Spacer(1, 6))
            i += 1; continue

        # blank line = paragraph break
        if not s.strip():
            flush_para(); flush_list()
            i += 1; continue

        para_buf.append(s)
        i += 1

    flush_para(); flush_list()

    doc = SimpleDocTemplate(
        str(OUT), pagesize=A4,
        leftMargin=0.6*inch, rightMargin=0.6*inch,
        topMargin=0.6*inch, bottomMargin=0.6*inch,
        title="SMEpred — Validation Dossier",
    )
    doc.build(story)
    print(f"Wrote {OUT}  ({OUT.stat().st_size/1024:.0f} KB)")


if __name__ == "__main__":
    build()
