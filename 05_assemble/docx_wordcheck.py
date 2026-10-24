# -*- coding: utf-8 -*-
"""y90: docx ↔ md 词数对账（v3）"""
import io, re, os
from docx import Document

O = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output"
doc = Document(os.path.join(O, "ALMP_manuscript_v3_20261023.docx"))
paras = [p.text.strip() for p in doc.paragraphs]

# 定位各一级标题（H() 生成的是普通段落，加粗 14pt）
marks = {}
for i, p in enumerate(paras):
    if p in ("Abstract", "Background", "Methods", "Results", "Discussion",
             "Declarations", "References", "Figure legends", "Tables",
             "Supplementary material"):
        marks.setdefault(p, i)
print("锚点:", marks)

def wc(a, b):
    n = 0
    for p in paras[a:b]:
        if p:
            n += len(p.split())
    return n

order = ["Abstract", "Background", "Methods", "Results", "Discussion", "Declarations"]
seg = {}
for k in order:
    a = marks[k] + 1
    b = marks[order[order.index(k) + 1]]
    if k == "Abstract":
        b = marks["Background"]
    seg[k] = wc(a, b)

# Abstract 与 Background 相邻
seg["Abstract"] = wc(marks["Abstract"] + 1, marks["Background"])

md = {"Abstract": 298, "Background": 880, "Methods": 2056, "Results": 3272,
      "Discussion": 1957}
heads = {"Abstract": 0, "Background": 0, "Methods": 8, "Results": 6, "Discussion": 7}

print("\n%-12s %8s %8s %8s %8s" % ("节", "docx", "md", "预期Δ", "实际Δ"))
ok = True
for k in ["Abstract", "Background", "Methods", "Results", "Discussion"]:
    d = seg[k]
    exp = md[k] + heads[k]
    print("%-12s %8d %8d %8d %8d %s" % (k, d, md[k], exp - md[k], d - md[k],
                                        "OK" if d == exp else "MISMATCH"))
    ok &= (d == exp)

imrad_docx = seg["Background"] + seg["Results"] + seg["Discussion"]
print("\ndocx IMRaD =", imrad_docx, " | md IMRaD =", 880 + 3272 + 1957)
print("docx 全部含摘要方法 =", seg["Abstract"] + imrad_docx + seg["Methods"])
print("总判定:", "PASS" if ok else "FAIL")
