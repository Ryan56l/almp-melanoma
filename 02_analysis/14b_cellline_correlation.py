# -*- coding: utf-8 -*-
"""P5-T1 追加：细胞自主子集下的相关对比 + 低检出组在肿瘤中的表现"""
import json, io, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
gl = pd.read_csv(os.path.join(D, "P5T1_gene_level_detection.tsv"), sep="\t")

core = gl[gl.segment_primary == "SEG1 adipokine core"]
intr = gl[gl.axis == "tumour-intrinsic"]
low = intr[intr.detect_ccle < 0.2]
hi = intr[intr.detect_ccle >= 0.9]
print("低检出组 n=%d  TCGA 检出率中位 %.3f | 表达中位 %.2f"
      % (len(low), low.detect_tcga.median(), low.median_tcga.median()))
print("高检出组 n=%d  TCGA 检出率中位 %.3f | 表达中位 %.2f"
      % (len(hi), hi.detect_tcga.median(), hi.median_tcga.median()))

# 重算相关：仅用细胞自主内在子集
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
gmap = json.load(io.open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gmap["entrez2symbol"].items()}
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)
M = cc.pivot_table(index="sampleId", columns="symbol", values="rpkm", aggfunc="mean")

tc = pd.read_csv(r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\data\cbio_skcm\mrna_data.tsv",
                 sep="\t")
tc["symbol"] = tc.entrezGeneId.astype(int).map(e2s)
T = tc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                            values="value", aggfunc="mean")

def score(mat, genes):
    g = [x for x in genes if x in mat.columns]
    return np.log1p(mat[g].clip(lower=0)).mean(axis=1), len(g)

res = {}
for tag, genes in [("all", intr.gene_symbol.tolist()),
                   ("cell_autonomous", hi.gene_symbol.tolist())]:
    a_cc, n_cc = score(M, core.gene_symbol.tolist())
    b_cc, _ = score(M, genes)
    a_tc, _ = score(T, core.gene_symbol.tolist())
    b_tc, _ = score(T, genes)
    r1 = stats.spearmanr(a_tc, b_tc); r2 = stats.spearmanr(a_cc, b_cc)
    print(f"\n[{tag}] 内在基因 {n_cc} 个")
    print(f"  实体瘤  n={len(T):3d}  rho={r1.correlation:+.3f}  p={r1.pvalue:.2g}")
    print(f"  细胞系  n={len(M):3d}  rho={r2.correlation:+.3f}  p={r2.pvalue:.2g}")
    res[tag] = dict(rho_tcga=float(r1.correlation), p_tcga=float(r1.pvalue),
                    rho_ccle=float(r2.correlation), p_ccle=float(r2.pvalue),
                    n_genes=int(n_cc), n_lines=int(len(M)), n_tcga=int(len(T)))
    np.savetxt(os.path.join(D, f"P5T1_scatter_{tag}.tsv"),
               np.c_[a_tc.values, b_tc.values], delimiter="\t", header="core\tintrinsic", comments="")
    np.savetxt(os.path.join(D, f"P5T1_scatter_ccle_{tag}.tsv"),
               np.c_[a_cc.values, b_cc.values], delimiter="\t", header="core\tintrinsic", comments="")

json.dump(res, io.open(os.path.join(D, "P5T1_rho_compare.json"), "w"), indent=1)
print("\nsaved")
