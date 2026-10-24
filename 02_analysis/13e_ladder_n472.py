# -*- coding: utf-8 -*-
"""
P0-c v2：纯度／含量去混杂「阶梯」
  c_imm  = 免疫限定含量（MCP 去 Fibro/Endo，**留一**：测谁就把它从协变量里剔除）
  c_str  = 基质含量（Fibro+Endo，**留一**）★测成纤维/内皮时必须剔靶标，否则自控
  c_pur  = 恶性含量（8 个黑色素谱系标记 z 分均值）
  c_all  = c_str + c_imm（同时控）
Part A: 10 个 MCP 群体 × 5 组协变量 × 2 评分算法 × 2 MCP 变体
Part B: 分段梯度（6 集合 × 两套读数：MCP 10 群体 / 27 免疫读数）中位 ρ
"""
import os
import io
import numpy as np
import pandas as pd
from scipy import stats

WS = r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36"
SC = r"C:\Rtmp\p2t12\results50"
TPM = r"C:\Rtmp\deconv\TCGA-SKCM.star_tpm.tsv.gz"
MCP = os.path.join(WS, "output/P2-T5b_多去卷积一致性v1/data")
OUT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output\P2-T20_纯度去混杂v1"
os.makedirs(OUT, exist_ok=True)

MARKERS = ["ENSG00000077498", "ENSG00000120215", "ENSG00000185664", "ENSG00000080166",
           "ENSG00000187098", "ENSG00000160307", "ENSG00000107165", "ENSG00000164199"]
ALMP6 = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
         "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]
STROMA = ["Fibroblasts", "Endothelial cells"]


def zr(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / (v.std(ddof=0) + 1e-12)


def partial_rho(x, y, covs, k):
    n = len(x)
    C = np.column_stack([np.ones(n)] + list(covs)) if k else np.ones((n, 1))
    rx = pd.Series(x).rank().values
    ry = pd.Series(y).rank().values
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    r, _ = stats.pearsonr(ex, ey)
    dfree = n - 2 - k
    t = r * np.sqrt(dfree / max(1e-12, 1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dfree)


print("[1] 恶性含量代理 …")
X = pd.read_csv(TPM, sep="\t", index_col=0)
X.index = [i.split(".")[0] for i in X.index]
X.columns = [str(c).strip()[:15] for c in X.columns]
sub = X.loc[[g for g in MARKERS if g in X.index]]
malig = pd.Series(np.mean([zr(sub.loc[g].values) for g in sub.index], axis=0), index=sub.columns)
print("    命中 %d/%d 标记" % (len(sub), len(MARKERS)))

rowsA, rowsB = [], []
for alg, fn in [("ssGSEA", "scores_ssgsea_fulltranscriptome.tsv"),
                ("GSVA", "scores_gsva_fulltranscriptome.tsv")]:
    S = pd.read_csv(os.path.join(SC, fn), sep="\t", index_col=0)
    S.columns = [c.strip() for c in S.columns]
    READS = [i for i in S.index if "::" in i]          # 27 条免疫读数
    for mcp_tag, mfn in [("noQN", "P2T5b_mcpcounter_tcga_noQN.tsv"),
                         ("QN", "P2T5b_mcpcounter_tcga_QN.tsv")]:
        M = pd.read_csv(os.path.join(MCP, mfn), sep="\t", index_col=0)
        M.columns = [str(c).strip()[:15] for c in M.columns]
        common = [c for c in S.columns if c in M.columns and c in malig.index]
        Zo = pd.DataFrame({p: zr(M.loc[p, common].values) for p in M.index}, index=common)
        pur = malig[common].values
        n = len(common)

        # ---- Part A ----
        x = S.loc["SEG1_ADIPOKINE_CORE", common].values
        for pop in M.index:
            y = M.loc[pop, common].values
            imm_loo = Zo[[p for p in M.index if p not in STROMA and p != pop]].mean(axis=1).values
            str_loo = Zo[[p for p in STROMA if p != pop]].mean(axis=1).values
            rec = dict(algorithm=alg, mcp=mcp_tag, target=pop, n=n)
            for tag, covs, k in [("raw", [], 0), ("c_imm", [imm_loo], 1), ("c_str", [str_loo], 1),
                                 ("c_pur", [pur], 1), ("c_all", [str_loo, imm_loo], 2),
                                 ("c_str_pur", [str_loo, pur], 2)]:
                r, p = partial_rho(x, y, covs, k)
                rec["rho_" + tag], rec["p_" + tag] = r, p
            rowsA.append(rec)

        # ---- Part B ----
        for almp in ALMP6:
            if almp not in S.index:
                continue
            xa = S.loc[almp, common].values
            for fam, targets in [("MCP10", list(M.index)), ("READ27", READS)]:
                for tag in ["raw", "c_imm", "c_str", "c_all"]:
                    rr = []
                    for t in targets:
                        if fam == "MCP10":
                            yv = M.loc[t, common].values
                            imm = Zo[[p for p in M.index if p not in STROMA and p != t]].mean(axis=1).values
                            st = Zo[[p for p in STROMA if p != t]].mean(axis=1).values
                        else:
                            yv = S.loc[t, common].values
                            imm = Zo[[p for p in M.index if p not in STROMA]].mean(axis=1).values
                            st = Zo[STROMA].mean(axis=1).values
                        cv = {"raw": [], "c_imm": [imm], "c_str": [st],
                              "c_all": [st, imm]}[tag]
                        rr.append(partial_rho(xa, yv, cv, len(cv))[0])
                    rowsB.append(dict(algorithm=alg, mcp=mcp_tag, almp=almp, family=fam, cov=tag,
                                      median_rho=float(np.median(rr)),
                                      n_sig=int(sum(1 for v in rr if v > 0)),
                                      n_eval=len(rr), n=n))

RA, RB = pd.DataFrame(rowsA), pd.DataFrame(rowsB)
RA.to_csv(os.path.join(OUT, "P2T20_D_covariate_ladder.tsv"), sep="\t", index=False)
RB.to_csv(os.path.join(OUT, "P2T20_E_segment_ladder.tsv"), sep="\t", index=False)

L = ["# P0-c　纯度／含量去混杂阶梯（TCGA-SKCM，n=472）\n"]
L.append("> **协变量定义（全部留一）**：`c_imm`=免疫限定含量（MCP 去 Fibro/Endo 的 z 分均值，"
         "**测谁就剔谁**）；`c_str`=基质含量（Fibro+Endo，**测成纤维/内皮时剔靶标**）；"
         "`c_pur`=恶性含量（TYR/MLANA/PMEL/DCT/MITF/S100B/TYRP1/SLC45A2 八标记 z 分均值）；"
         "`c_all`=基质+免疫同时控。\n")
L.append("> 评分与去卷积均为**全转录组输入**（SAP v1.11 条款 J）。偏 Spearman＝秩变换＋线性残差化，"
         "p 取自 t 分布，df=n−2−k。\n")

for alg in ["ssGSEA", "GSVA"]:
    for tag in ["noQN", "QN"]:
        s = RA[(RA.algorithm == alg) & (RA.mcp == tag)]
        L.append(f"\n## 表 D-{alg}-{tag}　SEG1 × MCP 群体　n={int(s.n.iloc[0])}\n")
        L.append("| 群体 | 原始 | 控免疫 | 控基质 | 控恶性 | 控基质+免疫 | p | 控基质+恶性 |")
        L.append("|---|---|---|---|---|---|---|---|")
        for _, r in s.iterrows():
            L.append(f"| {r.target} | {r.rho_raw:+.3f} | {r.rho_c_imm:+.3f} | {r.rho_c_str:+.3f} | "
                     f"{r.rho_c_pur:+.3f} | **{r.rho_c_all:+.3f}** | {r.p_c_all:.1e} | {r.rho_c_str_pur:+.3f} |")

for fam in ["MCP10", "READ27"]:
    L.append(f"\n## 表 E-{fam}　分段梯度在四组协变量下的中位 ρ\n")
    L.append("| 算法 | MCP | 集合 | 原始 | 控免疫 | 控基质 | 控基质+免疫 | 正号数 原始→控后 |")
    L.append("|---|---|---|---|---|---|---|---|")
    g0 = RB[RB.family == fam]
    for (alg, tag, almp), g in g0.groupby(["algorithm", "mcp", "almp"], sort=False):
        d = dict(zip(g["cov"], g["median_rho"]))
        np_ = dict(zip(g["cov"], g["n_sig"]))
        ne = int(g["n_eval"].iloc[0])
        L.append(f"| {alg} | {tag} | {almp} | {d.get('raw', np.nan):+.3f} | {d.get('c_imm', np.nan):+.3f} | "
                 f"{d.get('c_str', np.nan):+.3f} | **{d.get('c_all', np.nan):+.3f}** | "
                 f"{np_.get('raw', 0)}/{ne}→{np_.get('c_all', 0)}/{ne} |")

io.open(os.path.join(OUT, "P2T20_purity_ladder_v1.md"), "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
