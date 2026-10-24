# -*- coding: utf-8 -*-
"""① GSE115978 第二队列复现（v2：全转录组秩，修正 v1 的子集宇宙错误）

v1 的错误：只抽了 295 个目标基因，再在这 295 个基因内部算分位 → 得 0.46–0.50 的
全体噪声底，属于 SAP v1.11-J 明令禁止的「子集矩阵下算集合内部相对位次」。

v2 做法：读入完整 23,686 × 7,186 计数矩阵，按细胞（列）在全转录组内算分位，
再取目标基因的分位均值。与 GSE72056 的 UCell 口径一致。
"""
import os, gzip, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
GD = os.path.join(os.path.expanduser("~"), "WorkBuddy", "2026-10-07-22-15-36",
                  ".workbuddy", "data", "gse115978")
OUT = os.path.join(ROOT, "P5-T3_GSE115978复现v1")
OD = os.path.join(OUT, "data")
os.makedirs(OD, exist_ok=True)

s1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v2.csv"))
axm = dict(zip(s1.gene_symbol, s1.axis))
sgm = dict(zip(s1.gene_symbol, s1.segment_primary))
core = [g for g in s1.gene_symbol if sgm.get(g) == "SEG1 adipokine core"]
intr = [g for g in s1.gene_symbol if axm.get(g) == "tumour-intrinsic"]
comp = [g for g in s1.gene_symbol if axm.get(g) == "tissue-compositional"]
hi = pd.read_csv(os.path.join(D, "P5T1_cell_autonomous_intrinsic_subset_v2.tsv"), sep="\t")
ca166 = hi.gene_symbol.tolist()  # 171 genes: denominator corrected to the 50 melanoma lines

ann = pd.read_csv(os.path.join(GD, "GSE115978_cell.annotations.csv.gz")).set_index("cells")
print("annotations:", ann.shape, flush=True)

print("reading full matrix (genes x cells, float32) ...", flush=True)
_h = pd.read_csv(os.path.join(GD, "GSE115978_counts.csv.gz"), nrows=0)
_dt = {c: np.float32 for c in _h.columns[1:]}
M = pd.read_csv(os.path.join(GD, "GSE115978_counts.csv.gz"), index_col=0, dtype=_dt)
print("matrix:", M.shape, flush=True)
M.columns = [c.strip('"') for c in M.columns]
M.index = [g.strip('"') for g in M.index]

cell = [c for c in M.columns if c in ann.index]
M = M[cell]
ann = ann.loc[cell]
print("analysis cells:", M.shape, flush=True)

want = [g for g in (core + intr + comp + ca166) if g in M.index]
want = list(dict.fromkeys(want))
print("target genes found:", len(want),
      "| core", len([g for g in core if g in want]),
      "| intrinsic", len([g for g in intr if g in want]),
      "| 171", len([g for g in ca166 if g in want]),
      "| composition", len([g for g in comp if g in want]), flush=True)

# 按列（细胞）分块，在每个细胞内部对全部 23,686 个基因排序
ARR = np.empty((M.shape[1], len(want)), dtype=np.float32)
CH = 600
for s in range(0, M.shape[1], CH):
    blk = M.iloc[:, s:s + CH]
    rk = blk.rank(axis=0, pct=True)          # axis=0 = 每个细胞内跨基因排序
    ARR[s:s + CH, :] = rk.loc[want, :].to_numpy(dtype=np.float32).T
    if s % 3000 == 0:
        print("  ranked", s, "/", M.shape[1], flush=True)
R = pd.DataFrame(ARR, index=M.columns, columns=want)
print("rank matrix:", R.shape, flush=True)
del M

sc = pd.DataFrame(index=R.index)
sc["core"] = R[[g for g in core if g in R.columns]].mean(axis=1)
sc["intrinsic"] = R[[g for g in intr if g in R.columns]].mean(axis=1)
sc["intrinsic171"] = R[[g for g in ca166 if g in R.columns]].mean(axis=1)
sc["composition"] = R[[g for g in comp if g in R.columns]].mean(axis=1)
sc["cell_type"] = ann["cell.types"]
sc["patient"] = ann["samples"]
sc["cohort"] = ann["Cohort"]
sc.to_csv(os.path.join(OD, "P5T3_cell_scores.tsv"), sep="\t")

rows = []
for coh in ["New", "Tirosh", "ALL"]:
    sub = sc if coh == "ALL" else sc[sc.cohort == coh]
    for ct in sorted(sub.cell_type.unique()):
        if ct == "?":
            continue
        for ax_ in ["core", "intrinsic", "intrinsic171", "composition"]:
            m = sub.cell_type == ct
            rows.append(dict(cohort=coh, cell_type=ct, axis=ax_, n_cells=int(m.sum()),
                             median=float(sub.loc[m, ax_].median()),
                             mean=float(sub.loc[m, ax_].mean())))
C = pd.DataFrame(rows)
C.to_csv(os.path.join(OD, "P5T3_celltype_medians.tsv"), sep="\t", index=False)

pair = []
CTS = ["T.CD8", "T.CD4", "T.cell", "B.cell", "NK", "Macrophage", "CAF", "Endo."]
for coh in ["New", "ALL"]:
    sub = sc if coh == "ALL" else sc[sc.cohort == coh]
    for ct in CTS:
        for ax_ in ["core", "intrinsic", "intrinsic171", "composition"]:
            pv = sub.pivot_table(index="patient", columns="cell_type", values=ax_, aggfunc="median")
            if "Mal" not in pv.columns or ct not in pv.columns:
                continue
            d = pv[["Mal", ct]].dropna()
            rec = dict(cohort=coh, contrast=f"Mal vs {ct}", axis=ax_, n_patients=len(d),
                       median_Mal=float(d["Mal"].median()), median_other=float(d[ct].median()),
                       p=np.nan, note="n<5" if len(d) < 5 else "")
            if len(d) >= 5:
                if (d["Mal"] == d[ct]).all():
                    rec["note"] = "all ties"
                else:
                    rec["p"] = float(stats.wilcoxon(d["Mal"], d[ct]).pvalue)
            pair.append(rec)
P = pd.DataFrame(pair)
P["q"] = np.nan
for coh in P.cohort.unique():
    m = (P.cohort == coh) & P["p"].notna()
    if m.sum():
        P.loc[m, "q"] = stats.false_discovery_control(P.loc[m, "p"].values, method="bh")
P.to_csv(os.path.join(OD, "P5T3_paired_wilcoxon.tsv"), sep="\t", index=False)

print("\n=== 细胞类型中位（Cohort=New） ===")
print(C[(C.cohort == "New") & C.axis.isin(["core", "intrinsic", "intrinsic171"])]
      .pivot(index="cell_type", columns="axis", values="median").round(4).to_string())
print("\n=== 细胞类型中位（ALL, n=32 例） ===")
print(C[(C.cohort == "ALL") & C.axis.isin(["core", "intrinsic", "intrinsic171"])]
      .pivot(index="cell_type", columns="axis", values="median").round(4).to_string())
for coh in ["New", "ALL"]:
    for ax_ in ["core", "intrinsic", "intrinsic171"]:
        print(f"\n=== 配对 Wilcoxon（{coh}, {ax_}） ===")
        print(P[(P.cohort == coh) & (P.axis == ax_)].round(4)
              .drop(columns=["cohort"]).to_string(index=False))
print("\nsaved ->", OD)
