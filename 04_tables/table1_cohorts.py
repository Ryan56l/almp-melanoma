# -*- coding: utf-8 -*-
"""表 1：队列与数据源基线（Panel A 研究队列 / Panel B GWAS 汇总统计）"""
import os, sys, csv, io
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
OUT = os.path.join(ROOT, "output")

PANEL_A = [
    # cohort, layer, n, sample, platform, endpoint, source
    ("TCGA-SKCM", "Discovery", "472 tumours (368 metastases, 103 primary, 1 unclassified)",
     "Bulk tumour RNA-seq", "Illumina HiSeq", "Niche coupling; TMB orthogonality", "GDC"),
    ("Gide 2019", "ICB cohort", "91 (73 with RECIST response)", "Pre-treatment biopsy",
     "RNA-seq", "Response; coupling", "cBioPortal"),
    ("Liu 2019", "ICB cohort", "122", "Pre-treatment biopsy", "RNA-seq", "Response; coupling", "cBioPortal"),
    ("Riaz 2017", "ICB cohort", "107 (49 with pre-treatment response)", "Pre-/on-treatment biopsy",
     "RNA-seq", "Response; coupling", "cBioPortal"),
    ("UCLA 2016 (Hugo)", "ICB cohort", "27", "Pre-treatment biopsy", "RNA-seq",
     "Coupling only", "cBioPortal"),
    ("Van Allen 2015", "ICB cohort, post hoc", "40", "Pre-treatment biopsy", "RNA-seq",
     "Coupling only; no response field", "cBioPortal"),
    ("GSE65904", "External prognostic", "214 (150 with survival annotation)", "Primary tumours",
     "Microarray", "DMFS; DSS", "GEO"),
    ("GSE72056", "Single cell", "4,645 cells / 19 patients", "Tumour and adjacent tissue",
     "scRNA-seq", "Cell-type localisation", "GEO"),
]

GWAS = [
    # id, trait, role, n, cases/controls, population, k, F range, PMID
    # k / F range are instrument-strength fields and are filled only for exposure,
    # negative-control and mediator rows. Outcome rows carry no instruments, so both
    # fields are em-dashed rather than reporting the outcome GWAS's own significant-locus
    # count, which is a different quantity (see manuscript Table 1 footnote).
    ("GCST90475578", "Cutaneous melanoma (PheCode)", "Primary outcome",
     "446,751", "10,468 / 436,283", "European (MVP)", "—", "—", "39024449"),
    ("ebi-a-GCST90007319", "Circulating leptin, BMI-adjusted", "Exposure (primary)",
     "56,708", "—", "Mixed", "5", "38.6–100.7", "32917775"),
    ("ieu-a-1", "Adiponectin", "Exposure (primary)",
     "39,883", "—", "Mixed", "10", "33.0–368.2", "22479202"),
    ("ebi-a-GCST90012034", "Resistin", "Exposure (primary)",
     "21,758", "—", "European", "13", "30.1–192.1", "33067605"),
    ("ieu-b-40", "Body mass index (Yengo)", "Exposure (negative-control branch)",
     "681,275", "—", "European", "501", "28.6–2197.3", "30124842"),
    ("ebi-a-GCST90025984", "Lymphocyte count (Barton, MVP-free)", "Mediator (primary)",
     "443,762", "—", "European", "379", "25.7–1825.7", "34226706"),
    ("ebi-a-GCST90025985", "White blood cell count (Barton, MVP-free)", "Mediator (primary)",
     "444,734", "—", "European", "—", "—", "34226706"),
    ("ieu-b-4969", "Cutaneous melanoma (UK Biobank)", "Outcome (replication)",
     "375,767", "—", "European", "—", "—", "—"),
    ("finn-b-C3_SKIN", "Malignant skin neoplasm (FinnGen)", "Outcome (sensitivity)",
     "—", "—", "Finnish", "—", "—", "36653562"),
    ("ieu-b-4959", "Non-melanoma skin cancer", "Outcome (sensitivity)",
     "395,710", "—", "European", "—", "—", "—"),
]

md = io.StringIO()
md.write("# 表 1　队列与数据源基线（Table 1. Cohorts and data sources）\n\n")
md.write("> **口径**：n 为进入分析的样本量；TCGA 的 368 + 103 = 471 与总数 472 差 1，"
         "该 1 例 `sample_type` 未归类（见 Methods §9.2 待核验第 1 项）。\n"
         "> **样本重叠**：主结局为纯 MVP 队列；三暴露（Yaghootkar 2020 外显子芯片、ADIPOGen、"
         "Folkersen 2020）与 MVP 无个体重叠；免疫中介端采用 Barton 2021（纯 UK Biobank，无 MVP），"
         "含 MVP 的 Vuckovic / Chen 仅作敏感性。\n\n")

md.write("## A. 研究队列\n\n")
md.write("| Cohort | Layer | n | Sample | Platform | Endpoint used | Source |\n")
md.write("|---|---|---|---|---|---|---|\n")
for r in PANEL_A:
    md.write("| " + " | ".join(r) + " |\n")

md.write("\n## B. GWAS 汇总统计\n\n")
md.write("| Dataset | Trait | Role | n | Cases / controls | Population | k | F range | PMID |\n")
md.write("|---|---|---|---|---|---|---|---|---|\n")
for r in GWAS:
    md.write("| " + " | ".join(r) + " |\n")

md.write("\n> **注**：`ieu-b-4969` 与 `ieu-b-4959` 的效应量为线性概率尺度（BOLT-LMM），"
         "已按 logOR = β / [p(1−p)] 转换（因子 101.19 与 17.77），转换经四锚点验证"
         "（MC1R rs1805007 OR 1.787 等）；正文不报告任何 naive exp(β)。\n"
         "> 完整 23 个数据集清单（含 7 个未进入主分析的扩展暴露与 4 个敏感性暴露）见附表 S1。\n")

path_md = os.path.join(OUT, "表1_队列与数据源基线_v1.md")
io.open(path_md, "w", encoding="utf-8").write(md.getvalue())

# CSV
rows = []
for r in PANEL_A:
    rows.append({"panel": "A", "cohort_or_dataset": r[0], "layer": r[1], "n": r[2],
                 "sample": r[3], "platform_or_cases": r[4], "endpoint": r[5], "source_or_population": r[6],
                 "k_instruments": "", "F_range": "", "pmid": ""})
for r in GWAS:
    rows.append({"panel": "B", "cohort_or_dataset": r[0], "layer": r[2], "n": r[3],
                 "sample": r[4], "platform_or_cases": r[4], "endpoint": r[1], "source_or_population": r[5],
                 "k_instruments": r[6], "F_range": r[7], "pmid": r[8]})
df = pd.DataFrame(rows)
path_csv = os.path.join(OUT, "表1_队列与数据源基线_v1.csv")
df.to_csv(path_csv, index=False, encoding="utf-8-sig")

print("written:", path_md)
print("written:", path_csv)
print(df.shape)
