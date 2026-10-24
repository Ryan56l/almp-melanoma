# -*- coding: utf-8 -*-
"""
P0：纯度／整体细胞含量去混杂（TCGA-SKCM 层）

两个部分：
  A. 生态位主张（SEG1 × MCP-counter 十群体）：原始 ρ vs 偏 ρ
     协变量：①成分轴 16 基因评分 ②整体细胞含量（其余九群体 z 分均值，leave-one-out，避免把目标群体从自身均值里扣掉）
  B. 段间梯度（六集合 × 27 免疫读出）：原始中位 ρ vs 偏中位 ρ
  C. 中性粒框架依赖性（两算法 × MCP noQN/QN）

输出 output/P2-T20_纯度去混杂v1/
"""
import io, os
import numpy as np
import pandas as pd
from scipy import stats

WS = r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36"
SC = r"C:\Rtmp\p2t12\results50"
MCP = os.path.join(WS, "output/P2-T5b_多去卷积一致性v1/data")
OUT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36/.workbuddy/tmp/n472"
os.makedirs(OUT, exist_ok=True)

ALMP6 = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
         "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]


def load_scores(fn):
    df = pd.read_csv(os.path.join(SC, fn), sep="\t", index_col=0)
    df.columns = [c.strip() for c in df.columns]
    return df


def load_mcp(fn):
    df = pd.read_csv(os.path.join(MCP, fn), sep="\t", index_col=0)
    df.columns = [str(c).strip()[:15] for c in df.columns]
    return df


def partial_rho(x, y, covars, n_cov):
    """秩变换后的偏相关；covars 为 n×k 数组"""
    n = len(x)
    C = np.column_stack([np.ones(n), covars]) if n_cov else np.ones((n, 1))
    rx = x.rank().values
    ry = y.rank().values
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    r, _ = stats.pearsonr(ex, ey)
    dfree = n - 2 - n_cov
    t = r * np.sqrt(dfree / max(1e-12, 1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), dfree)


rows_a, rows_b, rows_c = [], [], []
nsamp = {}

for alg, fn in [("ssGSEA", "scores_ssgsea_fulltranscriptome.tsv"),
                ("GSVA", "scores_gsva_fulltranscriptome.tsv")]:
    S = load_scores(fn)
    immune = [r for r in S.index if r not in ALMP6]
    for tag, mfn in [("noQN", "P2T5b_mcpcounter_tcga_noQN.tsv"),
                     ("QN", "P2T5b_mcpcounter_tcga_QN.tsv")]:
        M = load_mcp(mfn)
        common = [c for c in S.columns if c in M.columns]
        nsamp[(alg, tag)] = len(common)
        Mc = M.loc[:, common]
        Z = Mc.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) + 1e-12), axis=1)   # 每群体 z 分
        comp = S.loc["ALMP_COMPOSITIONAL_AXIS", common]

        for a in ALMP6:
            x = S.loc[a, common]

            # ---- A：与十群体（LOO 整体含量） ----
            STROMAL = [p for p in ["Fibroblasts", "Endothelial cells"] if p in Mc.index]
            immune_only = [p for p in Mc.index if p not in STROMAL]
            g_imm = Z.loc[immune_only].mean(axis=0)          # 免疫限定含量（不含基质）
            for pop in Mc.index:
                y = Mc.loc[pop, common]
                g_loo = Z.drop(index=pop).mean(axis=0)        # 全量 leave-one-out（过校正）
                r0, p0 = stats.spearmanr(x, y)
                r_imm, p_imm = partial_rho(x, y, np.column_stack([comp.values, g_imm.values]), 2)
                r1, p1 = partial_rho(x, y, np.column_stack([comp.values, g_loo.values]), 2)
                r_g, p_g = partial_rho(x, y, g_imm.values.reshape(-1, 1), 1)   # 仅控免疫限定含量（不与 SEG1 共享基因）
                rows_a.append(dict(algorithm=alg, mcp=tag, almp=a, target=pop, layer="MCP",
                                   n=len(common), rho=r0, p=p0,
                                   rho_ctrl_immuneonly=r_g, p_ctrl_immuneonly=p_g,
                                   rho_ctrl_immune=r_imm, p_ctrl_immune=p_imm,
                                   rho_ctrl_both=r1, p_ctrl_both=p1))

            # ---- B：与 27 免疫读出 ----
            G = Z.mean(axis=0)          # 此处整体含量相对于免疫读出是外部量，用全均值
            for imm in immune:
                y = S.loc[imm, common]
                r0, p0 = stats.spearmanr(x, y)
                r1, p1 = partial_rho(x, y, np.column_stack([comp.values, G.values]), 2)
                rows_b.append(dict(algorithm=alg, almp=a, immune=imm, n=len(common),
                                   rho=r0, p=p0, rho_ctrl_both=r1, p_ctrl_both=p1))

            # ---- C：中性粒 ----
            if a == "SEG1_ADIPOKINE_CORE" and "Neutrophils" in Mc.index:
                y = Mc.loc["Neutrophils", common]
                g_loo = Z.drop(index="Neutrophils").mean(axis=0)
                g_imm = Z.loc[[p for p in Mc.index if p not in ("Fibroblasts", "Endothelial cells")]].mean(axis=0)
                r0, p0 = stats.spearmanr(x, y)
                r1, p1 = partial_rho(x, y, np.column_stack([comp.values, g_loo.values]), 2)
                ri, pi = partial_rho(x, y, np.column_stack([comp.values, g_imm.values]), 2)
                rows_c.append(dict(algorithm=alg, mcp=tag, rho=r0, p=p0,
                                   rho_ctrl_immune=ri, p_ctrl_immune=pi,
                                   rho_ctrl_both=r1, p_ctrl_both=p1, n=len(common)))

A = pd.DataFrame(rows_a); B = pd.DataFrame(rows_b); C = pd.DataFrame(rows_c)
A.to_csv(os.path.join(OUT, "P2T20_A_niche_partial.tsv"), sep="\t", index=False)
B.to_csv(os.path.join(OUT, "P2T20_B_segment_gradient.tsv"), sep="\t", index=False)
C.to_csv(os.path.join(OUT, "P2T20_C_neutrophil.tsv"), sep="\t", index=False)

L = []
L.append("# P0　纯度／整体细胞含量去混杂（TCGA-SKCM 层）v1\n")
L.append("> **日期**：2026-10-22 ｜ **问题**：SEG1 与基质／中性粒的偶联，是不是只是肿瘤纯度或基质含量？\n")
L.append("> **数据**：全转录组评分（33 × 472，`C:/Rtmp/p2t12/results50/`）× MCP-counter 十群体（472 样本，`P2T5b`）\n")
L.append("> **方法**：秩变换后线性偏相关（＝偏 Spearman）。主协变量＝**免疫限定整体含量**（MCP 十群体去掉 Fibro/Endo 后各自 z 分的均值），"
         "该变量与 SEG1 **不共享任何基因**；次协变量加成分轴 16 基因评分（注意：成分轴与 SEG1 共享 ADIPOQ/CFD/FABP4/ITLN1 四基因，属部分自控，仅作次要参考）\n")
L.append("> **对齐样本**：" + "；".join(f"{k[0]}-{k[1]} n={v}" for k, v in nsamp.items()) + "\n")

L.append("\n## A　生态位主张（SEG1 × 十群体）\n")
for alg in ["ssGSEA", "GSVA"]:
    L.append(f"\n### {alg}｜MCP-counter noQN\n")
    L.append("| 群体 | 原始 ρ | p | **仅控免疫限定含量**（无共享基因） | p | 控成分轴＋免疫含量 | p |")
    L.append("|---|---|---|---|---|---|---|")
    s = A[(A.algorithm == alg) & (A.mcp == "noQN") & (A.almp == "SEG1_ADIPOKINE_CORE")]
    for pop in s.target:
        r = s[s.target == pop].iloc[0]
        L.append(f"| {pop} | {r.rho:+.3f} | {r.p:.2e} | {r.rho_ctrl_immuneonly:+.3f} | {r.p_ctrl_immuneonly:.2e} | "
                 f"{r.rho_ctrl_immune:+.3f} | {r.p_ctrl_immune:.2e} |")

L.append("\n## B　段间梯度（六集合 × 27 免疫读出，中位 ρ）\n")
L.append("| 算法 | 集合 | 原始中位 ρ | 控制后中位 ρ | 保留 | 27 读出中 p<0.05 个数（原始→控制后） |")
L.append("|---|---|---|---|---|---|")
for alg in ["ssGSEA", "GSVA"]:
    for a in ALMP6:
        s = B[(B.algorithm == alg) & (B.almp == a)]
        k0 = int((s.p < 0.05).sum()); k1 = int((s.p_ctrl_both < 0.05).sum())
        L.append(f"| {alg} | {a} | {s.rho.median():+.3f} | {s.rho_ctrl_both.median():+.3f} | "
                 f"{(s.rho_ctrl_both.median()/s.rho.median()*100):.0f}% | {k0} → {k1} |")

L.append("\n## C　中性粒框架依赖性\n")
L.append("| 算法 | 去卷积 | 原始 ρ | p | 控免疫限定含量 ρ | p | 控全量含量 ρ | p |")
L.append("|---|---|---|---|---|---|---|---|")
for _, r in C.iterrows():
    L.append(f"| {r.algorithm} | MCP-{r.mcp} | {r.rho:+.3f} | {r.p:.2e} | {r.rho_ctrl_immune:+.3f} | "
             f"{r.p_ctrl_immune:.2e} | {r.rho_ctrl_both:+.3f} | {r.p_ctrl_both:.2e} |")

io.open(os.path.join(OUT, "P2T20_purity_deconfound_v1.md"), "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
