# -*- coding: utf-8 -*-
"""②c／新发现：细胞系面板谱系审计 + 成纤维系作为核心段阳性对照

发现：57 株「皮肤黑色素瘤细胞系」中有 7 株（HS600T/HS688AT/HS834T/HS839T/HS895T/
HS934T/HS940T）实为成纤维细胞系 —— DepMap 标注为 SKCM，但 MLANA/TYR/SOX10/S100B
≈0、COL1A1 577–4,869 RPKM。

处理：按标记把面板拆为「50 株黑素瘤 ＋ 7 株成纤维」，主分析用 50；
7 株成纤维改作核心段的**细胞类型阳性对照**（全文缺的正是这一层证据）。
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
OUT = os.path.join(ROOT, "P5-T2_跨平台秩校准v1")
OD = os.path.join(OUT, "data")
os.makedirs(OD, exist_ok=True)

gm = json.load(open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gm["entrez2symbol"].items()}
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)
P = cc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                             values="rpkm", aggfunc="mean")
gl = pd.read_csv(os.path.join(D, "P5T1_gene_level_detection.tsv"), sep="\t")
core = [g for g in gl[gl.segment_primary == "SEG1 adipokine core"].gene_symbol if g in P.columns]
intr = [g for g in gl[gl.axis == "tumour-intrinsic"].gene_symbol if g in P.columns]

# ---------- 谱系判定（只看标记，不看名字）
ctl = pd.read_csv(os.path.join(D, "ccle_melanoma_controlgens_rpkm.tsv"), sep="\t")
C = ctl.pivot_table(index="sampleId", columns="symbol", values="rpkm", aggfunc="mean")
MEL_MK = ["MLANA", "TYR", "SOX10", "S100B", "PMEL", "DCT"]  # MITF 在成纤维系可非零，已剔除
mel_mk = [g for g in MEL_MK if g in C.columns]
fib_mk = [g for g in ["COL1A1", "COL1A2", "DCN", "LUM", "THY1", "PDGFRB"] if g in C.columns]
aud = pd.DataFrame({
    "melanocytic_max": C[mel_mk].max(axis=1),
    "fibroblast_median": C[fib_mk].median(axis=1),
}).loc[P.index]
aud["call"] = np.where((aud.melanocytic_max < 1) & (aud.fibroblast_median >= 10),
                       "fibroblast", "melanoma")
print(aud.call.value_counts().to_string())
print(aud.sort_values("fibroblast_median", ascending=False).head(10).round(1).to_string())

fib = aud.index[aud.call == "fibroblast"].tolist()
mel = aud.index[aud.call == "melanoma"].tolist()
print("\nfibroblast lines:", fib)
print("melanoma lines:", len(mel))


def med_det(sub, genes):
    """与正文同口径：每个基因的检出率（≥1 RPKM），取基因间中位数"""
    return float((sub[genes] >= 1).mean(axis=0).median())


def rho_axes(sub):
    a = np.log1p(sub[core].clip(lower=0)).mean(axis=1)
    b = np.log1p(sub[intr].clip(lower=0)).mean(axis=1)
    r = stats.spearmanr(a, b)
    return float(r.correlation), float(r.pvalue)


rows = []
for tag, sub in [("all 57 (as analysed)", P), ("melanoma 50 (audited)", P.loc[mel]),
                 ("fibroblast 7 (control)", P.loc[fib])]:
    r, p = rho_axes(sub)
    rows.append(dict(group=tag, n_lines=len(sub),
                     core_detection=med_det(sub, core),
                     intrinsic_detection=med_det(sub, intr),
                     rho_core_intrinsic=r, p_rho=p))
S = pd.DataFrame(rows)
print("\n=== 分组主统计 ===")
print(S.round(4).to_string(index=False))

# ---------- 逐基因：成纤维 vs 黑素瘤
G = pd.DataFrame({
    "gene": core,
    "mel50_median_RPKM": [float(P.loc[mel, g].median()) for g in core],
    "fib7_median_RPKM": [float(P.loc[fib, g].median()) for g in core],
    "mel50_det": [float((P.loc[mel, g] >= 1).mean()) for g in core],
    "fib7_det": [float((P.loc[fib, g] >= 1).mean()) for g in core],
})
G["fib_over_mel"] = np.log2((G.fib7_median_RPKM + 0.1) / (G.mel50_median_RPKM + 0.1))
G = G.sort_values("fib_over_mel", ascending=False)
print("\n=== 核心 22 基因：成纤维/黑素瘤 ===")
print(G.round(2).to_string(index=False))
print("\n成纤维高出 ≥8 倍（log2 ≥3）的基因:",
      G[G.fib_over_mel >= 3].gene.tolist())
print("黑素瘤高出 ≥8 倍的基因:",
      G[G.fib_over_mel <= -3].gene.tolist())

# ---------- 核心—内在 ρ 的 2x2（全 288 vs 166 细胞自主子集）
tc = pd.read_csv(os.path.join(os.path.expanduser("~"), "WorkBuddy", "2026-10-07-22-15-36", ".workbuddy", "data", "cbio_skcm", "mrna_data.tsv"), sep="\t")
tc["symbol"] = tc.entrezGeneId.astype(int).map(e2s)
T = tc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                             values="value", aggfunc="mean")
hi = pd.read_csv(os.path.join(D, "P5T1_cell_autonomous_intrinsic_subset_v2.tsv"), sep="\t")
ca = [g for g in hi.gene_symbol if g in T.columns and g in P.columns]


def rr(sub, genes):
    genes = [g for g in genes if g in sub.columns]
    co = [g for g in core if g in sub.columns]
    a = np.log1p(sub[co].clip(lower=0)).mean(axis=1)
    b = np.log1p(sub[genes].clip(lower=0)).mean(axis=1)
    r = stats.spearmanr(a, b)
    return float(r.correlation), float(r.pvalue), len(genes)


q = []
for tag, sub in [("tumours (472)", T), ("melanoma lines (audited)", P.loc[mel]),
                 ("fibroblast lines (control)", P.loc[fib])]:
    for gname, genes in [("all intrinsic", intr), ("171 cell-autonomous", ca)]:
        r, p, k = rr(sub, genes)
        q.append(dict(dataset=tag, gene_set=gname, n_genes=k, rho=r, p=p))
Q = pd.DataFrame(q)
print("\n=== 核心—内在 ρ 2x2 ===")
print(Q.round(4).to_string(index=False))
Q.to_csv(os.path.join(OD, "P5T2_rho_2x2_audited.tsv"), sep="\t", index=False)

S.to_csv(os.path.join(OD, "P5T2_panel_audit_summary.tsv"), sep="\t", index=False)
G.to_csv(os.path.join(OD, "P5T2_core_genes_fibroblast_vs_melanoma.tsv"), sep="\t", index=False)
aud.to_csv(os.path.join(OD, "P5T2_lineage_audit.tsv"), sep="\t")
print("\nsaved ->", OD)
