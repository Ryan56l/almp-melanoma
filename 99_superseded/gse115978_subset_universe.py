# -*- coding: utf-8 -*-
"""① GSE115978（Jerby-Arnon 2018）第二队列复现

- 细胞级评分：UCell 式秩评分（每细胞内按表达排名，取基因集平均分位）
- 主检验：患者内配对 Wilcoxon（恶性细胞 vs 各非恶性类型），以患者为单位
- 细胞级 p 只作描述
- 独立性处理：Cohort=="New" 为主分析（17 例），"Tirosh" 子集仅作技术对照（与 GSE72056 有重叠风险）
"""
import os, gzip, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
GD = r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\data\gse115978"
OUT = os.path.join(ROOT, "P5-T3_GSE115978复现v1")
OD = os.path.join(OUT, "data")
os.makedirs(OD, exist_ok=True)

s1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v2.csv"))
ax = dict(zip(s1.gene_symbol, s1.axis))
seg = dict(zip(s1.gene_symbol, s1.segment_primary))
core = [g for g in s1.gene_symbol if seg.get(g) == "SEG1 adipokine core"]
intr = [g for g in s1.gene_symbol if ax.get(g) == "tumour-intrinsic"]
comp = [g for g in s1.gene_symbol if ax.get(g) == "tissue-compositional"]
hi = pd.read_csv(os.path.join(D, "P5T1_cell_autonomous_intrinsic_subset.tsv"), sep="\t")
ca166 = hi.gene_symbol.tolist()

ann = pd.read_csv(os.path.join(GD, "GSE115978_cell.annotations.csv.gz"))
ann = ann.set_index("cells")
print("annotations:", ann.shape)

# ---------- 读矩阵（判断方向）
with gzip.open(os.path.join(GD, "GSE115978_counts.csv.gz"), "rt") as f:
    head = f.readline().rstrip("\n")
cols = head.split(",")
print("header fields:", len(cols), "first3:", cols[:3])
index_is_cells = cols[1].strip('"') in set(ann.index)  # cols[0] 为行名前的空字段

usecols = None
if index_is_cells:          # 行=基因, 列=细胞
    genes_in_file = None
frames = []
with gzip.open(os.path.join(GD, "GSE115978_counts.csv.gz"), "rt") as f:
    head = f.readline()
    ncol = len(head.rstrip("\n").split(","))
    if index_is_cells:
        # 逐行流式，只保留需要的基因
        need = set(core + intr + ca166 + comp) | {"COL1A1", "DCN", "PTPRC", "MLANA", "S100B", "SOX10", "TYR"}
        keep = {}
        for line in f:
            p = line.rstrip("\n").split(",")
            g = p[0].strip('"')
            if g in need:
                keep[g] = p[1:]
        M = pd.DataFrame(keep, index=[c.strip('"') for c in head.rstrip("\n").split(",")[1:]])
    else:
        raise SystemExit("矩阵方向未知，需先检查")

print("matrix:", M.shape)
M.index = M.index.str.strip('"')
M = M.astype(float)
print("cells in matrix:", len(M), "mapped to annotation:", len(set(M.index) & set(ann.index)))

ann = ann.loc[[c for c in M.index if c in ann.index]]
M = M.loc[ann.index]
print("analysis cells:", M.shape)


def ucell_score(mat, genes):
    g = [x for x in genes if x in mat.columns]
    if not g:
        return pd.Series(np.nan, index=mat.index), 0
    r = mat.rank(axis=1, pct=True)
    return r[g].mean(axis=1), len(g)


sc = pd.DataFrame(index=M.index)
sc["core"], n1 = ucell_score(M, core)
sc["intrinsic"], n2 = ucell_score(M, intr)
sc["intrinsic166"], n3 = ucell_score(M, ca166)
sc["composition"], n4 = ucell_score(M, comp)
sc["cell_type"] = ann["cell.types"]
sc["patient"] = ann["samples"]
sc["cohort"] = ann["Cohort"]
print(f"genes used: core {n1} intrinsic {n2} 166 {n3} composition {n4}")

rows = []
for coh in ["New", "Tirosh", "ALL"]:
    sub = sc if coh == "ALL" else sc[sc.cohort == coh]
    for ct in sorted(sub.cell_type.unique()):
        if ct in ("?",):
            continue
        for axis in ["core", "intrinsic", "intrinsic166", "composition"]:
            rows.append(dict(cohort=coh, cell_type=ct, axis=axis,
                             n_cells=int((sub.cell_type == ct).sum()),
                             median=float(sub.loc[sub.cell_type == ct, axis].median())))
C = pd.DataFrame(rows)
C.to_csv(os.path.join(OD, "P5T3_celltype_medians.tsv"), sep="\t", index=False)

# ---------- 患者内配对
pair = []
for coh in ["New", "ALL"]:
    sub = sc if coh == "ALL" else sc[sc.cohort == coh]
    for ct in ["T.CD8", "T.CD4", "T.cell", "B.cell", "NK", "Macrophage", "CAF", "Endo."]:
        for axis in ["core", "intrinsic", "intrinsic166", "composition"]:
            pv = sub.pivot_table(index="patient", columns="cell_type", values=axis, aggfunc="median")
            if "Mal" not in pv or ct not in pv:
                continue
            d = pv[["Mal", ct]].dropna()
            if len(d) < 5:
                pair.append(dict(cohort=coh, contrast=f"Mal vs {ct}", axis=axis,
                                 n_patients=len(d), median_Mal=float(d.Mal.median()),
                                 median_other=float(d[ct].median()), p=np.nan, note="n<5"))
                continue
            w = stats.wilcoxon(d.Mal, d[ct])
            pair.append(dict(cohort=coh, contrast=f"Mal vs {ct}", axis=axis,
                             n_patients=len(d), median_Mal=float(d.Mal.median()),
                             median_other=float(d[ct].median()), p=float(w.pvalue), note=""))
P = pd.DataFrame(pair)
P["q"] = np.nan
for coh in P.cohort.unique():
    m = (P.cohort == coh) & P.p.notna()
    if m.sum():
        P.loc[m, "q"] = stats.false_discovery_control(P.loc[m, "p"].values, method="bh")
P.to_csv(os.path.join(OD, "P5T3_paired_wilcoxon.tsv"), sep="\t", index=False)

print("\n=== 细胞类型中位（Cohort=New） ===")
print(C[(C.cohort == "New") & (C.axis.isin(["core", "intrinsic", "intrinsic166"]))]
      .pivot(index="cell_type", columns="axis", values="median").round(4).to_string())
print("\n=== 患者内配对 Wilcoxon（Cohort=New，axis=intrinsic） ===")
print(P[(P.cohort == "New") & (P.axis == "intrinsic")].round(4).to_string(index=False))
print("\n=== 同（axis=core） ===")
print(P[(P.cohort == "New") & (P.axis == "core")].round(4).to_string(index=False))
print("\nsaved ->", OD)
