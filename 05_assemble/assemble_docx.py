# -*- coding: utf-8 -*-
"""y70: 拼装完整稿 v1（JTM 格式，docx）
输入：四稿 md ＋ 图注 md ＋ 表 1/表 2 csv ＋ 附表说明 md
输出：manuscript/output/ALMP_manuscript_v1_20261023.docx
规则：正文严格按 md 散文段落抽取；跳过内部批注块（>）、自查表、变更日志；
     参考文献从 Abstract&Intro 的引文表解析；图表置于正文之后。
"""
import io, os, re, csv
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

O = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output"
OUT_DOCX = os.path.join(O, "ALMP_manuscript_v1_20261023.docx")


def read(fn):
    return io.open(os.path.join(O, fn), encoding="utf-8").read()


def is_body(line):
    s = line.strip()
    if not s:
        return False
    if s[0] in "#>|":
        return False
    if set(s) <= set("-—–=*_ "):
        return False
    if s.startswith("**") and s.endswith("**") and len(s) < 60:
        return False   # 小标题/段标签
    if s.startswith("```"):
        return False
    return True


def sections(text, cut=None):
    """→ [(heading, [paras])]，按 # 标题切分；同一标题下多段合并"""
    if cut:
        i = text.find(cut)
        if i > 0:
            text = text[:i]
    out, cur, buf = [], "_head", []
    for line in text.splitlines():
        if line.startswith("#"):
            if buf:
                out.append((cur, " ".join(buf))); buf = []
            cur = line.lstrip("# ").strip(); continue
        if is_body(line):
            buf.append(line.strip())
        else:
            if buf:
                out.append((cur, " ".join(buf))); buf = []
    if buf:
        out.append((cur, " ".join(buf)))
    merged = []
    for h, t in out:
        if merged and merged[-1][0] == h:
            merged[-1] = (h, merged[-1][1] + "\n\n" + t)
        else:
            merged.append((h, t))
    return merged


AI = read("论文初稿_Abstract与Introduction_v1.md")
RS = read("论文初稿_Results_v1.md")
DS = read("论文初稿_Discussion_v1.md")
MS = read("论文初稿_Methods_v1.md")
LG = read("图注_七条_JTM_v1.md")


def between(text, start, end):
    a = text.find(start)
    b = text.find(end, a + 1) if end else len(text)
    return text[a:b]


# ---------------- 取各部分
def paras(text):
    return [l.strip() for l in text.splitlines() if is_body(l)]


title_txt = paras(between(AI, "## 一、Title", "## 二、Abstract"))[0]
title_txt = re.sub(r"^\*\*|\*\*$", "", title_txt).strip()

abs_block = paras(between(AI, "## 二、Abstract", "## 三、Introduction"))
intro_block = paras(between(AI, "## 三、Introduction", "## 四、"))

# 引文表 → References
refs = []
for line in between(AI, "## 四、引文清单", "## 五、").splitlines():
    m = re.match(r"\|\s*\**\[(\d+)\]\**\s*\|(.*?)\|(.*?)\|(.*?)\|", line)
    if m:
        num, cit, pmid, rest = m.groups()
        cit = cit.strip().strip("*").strip()
        rest = rest.strip().strip("*").strip().replace("·", ".")
        parts = [p.strip() for p in rest.split("|")]
        refs.append((int(num), cit, pmid.strip(), rest))

doc = Document()

# ---------------- 页面与样式
for sec in doc.sections:
    sec.page_height, sec.page_width = Cm(29.7), Cm(21.0)
    sec.top_margin = sec.bottom_margin = Cm(2.54)
    sec.left_margin = sec.right_margin = Cm(2.54)

st = doc.styles["Normal"]
st.font.name = "Times New Roman"
st.font.size = Pt(12)
st.element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
st.paragraph_format.line_spacing = 2.0
st.paragraph_format.space_after = Pt(0)


def P(text="", bold=False, size=None, align=None, space_before=0, space_after=0,
      italic=False, hang=False):
    p = doc.add_paragraph()
    for k, chunk in enumerate(re.split(r"\*\*", text)):
        if not chunk:
            continue
        r = p.add_run(chunk)
        r.bold = bold or (k % 2 == 1)
        r.italic = italic
    if size:
        r.font.size = Pt(size)
    if align == "c":
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    elif align == "j":
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    p.paragraph_format.space_before = Pt(space_before)
    p.paragraph_format.space_after = Pt(space_after)
    if hang:
        p.paragraph_format.left_indent = Cm(0.75)
        p.paragraph_format.first_line_indent = Cm(-0.75)
    return p


def H(text, level=1):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    r.font.size = Pt(14 if level == 1 else 12)
    p.paragraph_format.space_before = Pt(12 if level == 1 else 8)
    p.paragraph_format.space_after = Pt(4)
    return p


# ---------------- Title page
P("A decomposable adipokine–lipid programme in cutaneous melanoma: a stromal "
  "adipokine core and a stroma-independent tumour-intrinsic axis",
  bold=True, size=14, align="c", space_after=12)
P("[Authors — to be completed]", align="c", space_after=4)
P("[Affiliations — to be completed]", align="c", space_after=4)
P("Corresponding author: [name, email]", align="c", space_after=4)
P("Word count (Introduction + Results + Discussion): 5,943 · Abstract: 299 · "
  "Figures: 7 · Tables: 2 · Supplementary tables: S1, S3, S4 (S2 in preparation)",
  align="c", size=10, space_after=18)

# ---------------- Abstract
H("Abstract")
for para in abs_block:
    m = re.match(r"^\*\*(.+?)\*\*\s*(.*)$", para)
    if m:
        p = doc.add_paragraph()
        r1 = p.add_run(m.group(1) + ": ")
        r1.bold = True
        p.add_run(m.group(2))
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
    else:
        P(para, align="j")
P("Keywords: cutaneous melanoma; adipokine; lipid metabolism; tumour microenvironment; "
  "deconvolution; Mendelian randomization", space_before=6, align="j")

# ---------------- Body
def emit(section_list, startswith, level=1):
    for h, t in section_list:
        if not re.match(startswith, h):
            continue
        H(re.sub(r"^\d+(\.\d+)?\s*", "", h) if level > 1 else h, level)
        for para in [x.strip() for x in t.split("\n\n") if x.strip()]:
            P(para, align="j", space_after=6)


H("Background")
for para in intro_block:
    P(para, align="j", space_after=6)

H("Methods")
emit(sections(MS, "## 9."), r"^2\.\d", 2)
H("Results")
emit(sections(RS, "## 3.6"), r"^3\.[1-5] ", 2)
H("Discussion")
emit(sections(DS, "## 4.8"), r"^4\.[1-7] ", 2)

# ---------------- Declarations
H("Declarations")
for sub, body in [
    ("Ethics approval and consent to participate",
     "Not applicable. This study used only publicly available, de-identified datasets and "
     "published summary statistics; no human participants or animals were involved."),
    ("Consent for publication", "Not applicable."),
    ("Availability of data and materials",
     "All datasets are public and listed with their accessions in the Data Availability "
     "statement below. Derived score matrices, covariate-adjusted correlation tables and "
     "figure source data are provided as Supplementary Tables S1–S4. Analysis code is "
     "available from the corresponding author on reasonable request."),
    ("Competing interests", "[To be completed by the authors.]"),
    ("Funding", "[To be completed by the authors.]"),
    ("Authors' contributions", "[To be completed: CRediT roles.]"),
    ("Acknowledgements", "[To be completed.]"),
]:
    p = doc.add_paragraph()
    r = p.add_run(sub + ". ")
    r.bold = True
    p.add_run(body)
    p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY

H("Data availability")
da = between(read("附表S1S3S4与数据可得性_v1.md"), "## Data Availability", "> **待补两项")
for para in [x.strip("> ").strip() for x in da.splitlines() if x.strip().startswith(">")]:
    if para and not para.startswith("**Data Availability"):
        P(para.replace("**", ""), align="j", size=11, space_after=4)

# ---------------- References
doc.add_page_break()
H("References")
for num, cit, pmid, rest in refs:
    P(f"{num}. {cit}. {rest} PMID: {pmid}.", size=10, hang=True, align="j")

# ---------------- Figure legends（图注后附图，便于通读；投稿时另附单独图件）
FIGF = os.path.join(O, "figures")
FIGMAP = {
    1: ["fig1_design_hierarchy.png"],
    2: ["fig2_almp_construction.png"],
    3: ["fig3_axis_decomposition.png"],
    4: ["fig4_group.png"],
    5: ["fig5_singlecell.png"],
    6: ["fig6a_causal_forest.png", "fig6bcd_causal_audit.png"],
    7: ["fig7_exploratory.png"],
}

doc.add_page_break()
H("Figure legends")
for block in LG.split("\n## ")[1:]:
    head_line, _, body = block.partition("\n")
    P(head_line.strip(), bold=True, space_before=8, space_after=2)
    for para in [x for x in body.splitlines()
                 if x.strip() and not x.startswith(">")
                 and not re.fullmatch(r"[-–—\s|]+", x)]:
        P(para.strip().replace("**", ""), size=10, align="j", space_after=4)
    m = re.match(r"Figure\s+(\d+)\.", head_line.strip())
    if m:
        for fn in FIGMAP.get(int(m.group(1)), []):
            fp = os.path.join(FIGF, fn)
            if os.path.exists(fp):
                doc.add_picture(fp, width=Cm(15.5))
                doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
                P("", space_after=6)

# ---------------- Tables
def add_csv_table(fn, caption):
    P(caption, bold=True, space_before=12, space_after=4)
    rows = list(csv.reader(io.open(os.path.join(O, fn), encoding="utf-8-sig")))
    rows = [r for r in rows if any(c.strip() for c in r)]
    if not rows:
        return
    t = doc.add_table(rows=len(rows), cols=len(rows[0]))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, row in enumerate(rows):
        for j, cell in enumerate(row[:len(rows[0])]):
            c = t.cell(i, j)
            c.text = ""
            run = c.paragraphs[0].add_run(cell)
            run.font.size = Pt(7)
            if i == 0:
                run.bold = True
    P("", space_after=6)


doc.add_page_break()
H("Tables")
add_csv_table("表1_队列与数据源基线_v1.csv",
              "Table 1. Cohorts and data sources.")
add_csv_table("表2_因果层32组合判定矩阵_v1.csv",
              "Table 2. The 32 tissue × locus combinations of the causal audit.")

doc.add_page_break()
H("Supplementary material")
for s in [
    "Table S1. The 308-gene adipokine–lipid metabolic programme, axis assignment and "
    "adipose-specificity criteria (supp/S1_ALMP_308gene_list_and_axis_v1.csv).",
    "Table S2. Sensitivity analyses (deconvolution framework, scoring algorithm, "
    "preprocessing, instrument pruning). [To be compiled.]",
    "Table S3. All immune checkpoint blockade response and TMB-stratified tests, with "
    "pooled estimates (supp/S3_ICB_response_and_TMB_tests_v1.csv).",
    "Table S4. Reporting-standard crosswalk (STROBE, STROBE-MR, TRIPOD+AI, REMARK, SAGER, "
    "RECORD) and the SAP version chain.",
]:
    P(s, align="j", size=10, space_after=4)

doc.save(OUT_DOCX)
print("saved:", OUT_DOCX)
print("段落:", len(doc.paragraphs), "| 表:", len(doc.tables), "| 参考文献:", len(refs))
