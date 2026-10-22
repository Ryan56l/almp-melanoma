# -*- coding: utf-8 -*-
"""
P0-d：六集合 × 10 MCP 群体，单独控「成纤维细胞含量」——判定基质条件由哪个轴驱动
并检验髓系偶联是否由成纤维含量中介
"""
import os, io
import numpy as np
import pandas as pd
from scipy import stats

WS = r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36"
SC = r"C:\Rtmp\p2t12\results50"
MCP = os.path.join(WS, "output/P2-T5b_多去卷积一致性v1/data")
OUT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output\P2-T20_纯度去混杂v1"
os.makedirs(OUT, exist_ok=True)
ALMP6 = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
         "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]


def zr(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / (v.std(ddof=0) + 1e-12)


def partial_rho(x, y, covs):
    n, k = len(x), len(covs)
    C = np.column_stack([np.ones(n)] + list(covs)) if k else np.ones((n, 1))
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    r, _ = stats.pearsonr(ex, ey)
    dfree = n - 2 - k
    t = r * np.sqrt(dfree / max(1e-12, 1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dfree)


rows = []
for alg, fn in [("ssGSEA", "scores_ssgsea_fulltranscriptome.tsv"),
                ("GSVA", "scores_gsva_fulltranscriptome.tsv")]:
    S = pd.read_csv(os.path.join(SC, fn), sep="\t", index_col=0)
    S.columns = [c.strip() for c in S.columns]
    for mcp_tag, mfn in [("noQN", "P2T5b_mcpcounter_tcga_noQN.tsv"),
                         ("QN", "P2T5b_mcpcounter_tcga_QN.tsv")]:
        M = pd.read_csv(os.path.join(MCP, mfn), sep="\t", index_col=0)
        M.columns = [c.strip().rstrip("A") for c in M.columns]
        common = [c for c in S.columns if c in M.columns]
        fib_raw = M.loc["Fibroblasts", common].values
        fib = zr(fib_raw)
        for almp in ALMP6:
            x = S.loc[almp, common].values
            # 与成纤维含量的偶联（原始 / 控内皮，因内皮与成纤维高度共变）
            endo = zr(M.loc["Endothelial cells", common].values)
            r_f, p_f = partial_rho(x, fib_raw, [])
            r_f2, p_f2 = partial_rho(x, fib_raw, [endo])
            rows.append(dict(algorithm=alg, mcp=mcp_tag, almp=almp, target="Fibroblasts",
                             rho=r_f, p=p_f, rho_c_endo=r_f2, p_c_endo=p_f2, n=len(common)))
            # 各群体：原始 / 控成纤维
            for pop in M.index:
                if pop == "Fibroblasts":
                    continue
                y = M.loc[pop, common].values
                r0, p0 = partial_rho(x, y, [])
                r1, p1 = partial_rho(x, y, [fib])
                rows.append(dict(algorithm=alg, mcp=mcp_tag, almp=almp, target=pop,
                                 rho=r0, p=p0, rho_c_endo=np.nan, p_c_endo=np.nan,
                                 rho_c_fib=r1, p_c_fib=p1, n=len(common)))

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "P2T20_F_fibroblast_mediation.tsv"), sep="\t", index=False)

L = ["# P0-d　成纤维含量：是哪个轴驱动，以及是否中介髓系偶联\n",
     "> 目的：Seg1 与成纤维含量在全部校正下存活（ρ≈0.45）。若该条由 16 基因成分轴（ADIPOQ/CFD/FABP4/ITLN1 等脂肪组织来源基因）驱动，"
     "则提示评分部分在报告自身组织来源；若由 292 基因内在轴驱动，则更像肿瘤细胞自身程序。\n"]
for alg in ["ssGSEA", "GSVA"]:
    for tag in ["noQN", "QN"]:
        s = R[(R.algorithm == alg) & (R.mcp == tag) & (R.target == "Fibroblasts")]
        if s.empty:
            continue
        L.append(f"\n## 表 F1-{alg}-{tag}　各集合 × 成纤维含量\n")
        L.append("| 集合 | ρ | p | 控内皮后 ρ | p |")
        L.append("|---|---|---|---|---|")
        for _, r in s.iterrows():
            L.append(f"| {r.almp} | {r.rho:+.3f} | {r.p:.1e} | {r.rho_c_endo:+.3f} | {r.p_c_endo:.1e} |")

L.append("\n## 表 F2　控成纤维后各群体偶联（ssGSEA／noQN，SEG1）\n")
s = R[(R.algorithm == "ssGSEA") & (R.mcp == "noQN") & (R.almp == "SEG1_ADIPOKINE_CORE")
      & (R.target != "Fibroblasts")]
L.append("| 群体 | 原始 ρ | p | 控成纤维 ρ | p | 保留 |")
L.append("|---|---|---|---|---|---|")
for _, r in s.iterrows():
    keep = r.rho_c_fib / r.rho * 100 if abs(r.rho) > 1e-9 else np.nan
    L.append(f"| {r.target} | {r.rho:+.3f} | {r.p:.1e} | {r.rho_c_fib:+.3f} | {r.p_c_fib:.1e} | {keep:.0f}% |")

io.open(os.path.join(OUT, "P2T20_fibroblast_mediation_v1.md"), "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
