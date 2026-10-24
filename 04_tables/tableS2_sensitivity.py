# -*- coding: utf-8 -*-
"""S2 表修订：按 50 株黑素瘤 ＋ 7 株成纤维审计口径替换 e/f 块，并新增 g 块"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
T2 = os.path.join(ROOT, "P5-T2_跨平台秩校准v1", "data")
S2 = os.path.join(ROOT, "supp", "S2_sensitivity_analyses_v1.csv")

gm = json.load(open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gm["entrez2symbol"].items()}
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)
P = cc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                             values="rpkm", aggfunc="mean")
gl = pd.read_csv(os.path.join(D, "P5T1_gene_level_detection.tsv"), sep="\t")
core = [g for g in gl[gl.segment_primary == "SEG1 adipokine core"].gene_symbol if g in P.columns]
intr = [g for g in gl[gl.axis == "tumour-intrinsic"].gene_symbol if g in P.columns]
aud = pd.read_csv(os.path.join(T2, "P5T2_lineage_audit.tsv"), sep="\t").set_index("sampleId")
mel = aud.index[aud["call"] == "melanoma"].tolist()
fib = aud.index[aud["call"] == "fibroblast"].tolist()
Q = pd.read_csv(os.path.join(T2, "P5T2_rho_2x2_audited.tsv"), sep="\t")
G = pd.read_csv(os.path.join(T2, "P5T2_core_genes_fibroblast_vs_melanoma.tsv"), sep="\t")

d = pd.read_csv(S2)
d = d[~d.analysis_block.str.startswith(("S2e", "S2f"))].copy()

rows = []
for th in [1.0, 2.0, 5.0]:
    cm = float((P.loc[mel, core] >= th).mean(axis=0).median())
    im = float((P.loc[mel, intr] >= th).mean(axis=0).median())
    rows.append(dict(analysis_block="S2e detection threshold (cell-line validation, audited panel)",
                     comparison=f"RPKM >= {th}", level="50 melanoma cell lines",
                     statistic="median detection fraction (core / intrinsic)",
                     value=f"{cm:.3f} / {im:.3f}", n=50,
                     note=f"difference {cm - im:+.3f}"))
cf = float((P.loc[fib, core] >= 1.0).mean(axis=0).median())
inf = float((P.loc[fib, intr] >= 1.0).mean(axis=0).median())
rows.append(dict(analysis_block="S2e detection threshold (cell-line validation, audited panel)",
                 comparison="RPKM >= 1.0", level="7 fibroblast lines (positive control)",
                 statistic="median detection fraction (core / intrinsic)",
                 value=f"{cf:.3f} / {inf:.3f}", n=7,
                 note=f"core detection 1.9x that of melanoma lines ({cf:.3f} vs {cm:.3f})"
                      if (cm := float((P.loc[mel, core] >= 1.0).mean(axis=0).median())) else ""))

lab = {"tumours (472)": "472 TCGA-SKCM tumours",
       "melanoma lines (audited)": "50 melanoma cell lines (marker audit)",
       "fibroblast lines (control)": "7 fibroblast lines (positive control)"}
sets = {"all intrinsic": "core vs intrinsic, all 288 assessable intrinsic genes",
        "166 cell-autonomous": "core vs intrinsic, 166 cell-autonomous genes"}
for _, r in Q.iterrows():
    rows.append(dict(analysis_block="S2f gene-set restriction (axis correlation, audited panel)",
                     comparison=sets[r.gene_set], level=lab[r.dataset],
                     statistic="Spearman rho", value=f"{r.rho:+.3f}", n=int(
                         {"tumours (472)": 472, "melanoma lines (audited)": 50,
                          "fibroblast lines (control)": 7}[r.dataset]),
                     note=f"p = {r.p:.3g}; {r.n_genes} genes"))
top = G.assign(absr=G.fib_over_mel.abs()).nlargest(3, "absr")
for _, r in top.iterrows():
    rows.append(dict(analysis_block="S2g lineage audit (deposited label vs marker expression)",
                     comparison=f"{r.gene} median RPKM, fibroblast / melanoma",
                     level="7 fibroblast vs 50 melanoma lines", statistic="fold difference",
                     value=f"{r.fib7_median_RPKM:.2f} / {r.mel50_median_RPKM:.2f}",
                     n=57, note=f"{r.fib7_median_RPKM / r.mel50_median_RPKM:.1f}x"
                     if r.mel50_median_RPKM > 0 else "melanoma median = 0"))
rows.append(dict(analysis_block="S2g lineage audit (deposited label vs marker expression)",
                 comparison="lines called melanocytic vs fibroblastic",
                 level="57 lines deposited as cutaneous melanoma", statistic="line counts",
                 value="50 / 7", n=57,
                 note="melanocytic = max(MLANA,TYR,SOX10,S100B,PMEL,DCT) >= 1 RPKM; "
                      "fibroblastic = all six < 1 and median(COL1A1,COL1A2,DCN,LUM,THY1,PDGFRB) >= 10"))
rows.append(dict(analysis_block="S2g lineage audit (deposited label vs marker expression)",
                 comparison="cross-platform position of the core after rank conversion",
                 level="cell lines vs tumours", statistic="median within-sample percentile",
                 value="0.276 / 0.289", n=57,
                 note="same 308-gene universe on both sides; absolute-threshold contrast is withdrawn"))
d = pd.concat([d, pd.DataFrame(rows)], ignore_index=True)
d.to_csv(S2, index=False, encoding="utf-8-sig")
print(d.shape)
print(d[d.analysis_block.str.startswith(("S2e", "S2f", "S2g"))].to_string())
