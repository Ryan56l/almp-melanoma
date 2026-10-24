# -*- coding: utf-8 -*-
"""y90: docx ↔ md 词数对账（v3.3）"""
import os
from docx import Document

O = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output"
doc = Document(os.path.join(O, "ALMP_manuscript_v3_20261023.docx"))
paras = [p.text.strip() for p in doc.paragraphs]

marks = {}
for i, p in enumerate(paras):
    if p in ("Abstract", "Background", "Methods", "Results", "Discussion",
             "Declarations", "References", "Figure legends", "Tables",
             "Supplementary material"):
        marks.setdefault(p, i)
print("锚点:", marks)


def wc(a, b):
    return sum(len(p.split()) for p in paras[a:b] if p)


seg = {
    "Abstract": wc(marks["Abstract"] + 1, marks["Background"]),
    "Background": wc(marks["Background"] + 1, marks["Methods"]),
    "Methods": wc(marks["Methods"] + 1, marks["Results"]),
    "Results": wc(marks["Results"] + 1, marks["Discussion"]),
    "Discussion": wc(marks["Discussion"] + 1, marks["Declarations"]),
}

md = {"Abstract": 298, "Background": 880, "Methods": 2154, "Results": 3291,
      "Discussion": 1980}
# Δ 恒等于该节小节标题的"词数"（不是条数）。逐节标题词数：
#   Abstract  11 = 关键词行 "Keywords: ..."
#   Methods   48 = 8 条小节标题合计
#   Results   75 = 12 + 15 + 11 + 14 + 12 + 11（六条）
#   Discussion 35 = 7 条小节标题合计
heads = {"Abstract": 11, "Background": 0, "Methods": 48, "Results": 75, "Discussion": 35}

print("\n%-12s %8s %8s %8s %8s" % ("节", "docx", "md", "标题词数", "实际Δ"))
ok = True
for k in ["Abstract", "Background", "Methods", "Results", "Discussion"]:
    d = seg[k]
    exp = md[k] + heads[k]
    print("%-12s %8d %8d %8d %8d %s" % (k, d, md[k], heads[k], d - md[k],
                                        "OK" if d == exp else "MISMATCH"))
    ok &= (d == exp)

imrad_docx = seg["Background"] + seg["Results"] + seg["Discussion"]
print("\ndocx IMRaD =", imrad_docx, " | md IMRaD =", 880 + 3291 + 1980)
print("docx 全部含摘要方法 =", seg["Abstract"] + imrad_docx + seg["Methods"])
print("段落:", len(paras), "| 表:", len(doc.tables),
      "| 内嵌图:", len(doc.inline_shapes))
print("总判定:", "PASS" if ok else "FAIL")
