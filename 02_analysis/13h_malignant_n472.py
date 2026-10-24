# -*- coding: utf-8 -*-
"""
P0-b：以「恶性细胞含量」为纯度代理，重做 MCP 层的偏相关
纯度代理 = 黑色素细胞谱系标记（TYR/MLANA/PMEL/DCT/MITF/S100B/TYRP1/SLC45A2）的 z 分均值
控制 = ①免疫限定含量 ②恶性含量  —— 刻意不控基质，因为基质正是被检验的信号
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
# 成稿期单写：新产出写回 manuscript/output/
OUT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36/.workbuddy/tmp/n472"
os.makedirs(OUT, exist_ok=True)

MARKERS = {
    "ENSG00000077498": "TYR", "ENSG00000120215": "MLANA", "ENSG00000185664": "PMEL",
    "ENSG00000080166": "DCT", "ENSG00000187098": "MITF", "ENSG00000160307": "S100B",
    "ENSG00000107165": "TYRP1", "ENSG00000164199": "SLC45A2",
}

print("[1] 读 TPM 矩阵 …")
X = pd.read_csv(TPM, sep="\t", index_col=0)
X.index = [i.split(".")[0] for i in X.index]
X.columns = [str(c).strip()[:15] for c in X.columns]
print("    维度:", X.shape)

sub = X.loc[[g for g in MARKERS if g in X.index]]
print("[2] 标记基因命中 %d/%d" % (len(sub), len(MARKERS)))
print(sub.mean(axis=1).round(1).rename(index=MARKERS).to_string())
C = sub.T.corr(method="spearman")
print("[3] 标记间 Spearman 中位:", round(np.median(C.values[np.triu_indices_from(C, 1)]), 3))

Z = sub.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) + 1e-12), axis=1)
malig = Z.mean(axis=0)          # 恶性含量（纯度代理）

rows = []
for alg, fn in [("ssGSEA", "scores_ssgsea_fulltranscriptome.tsv"),
                ("GSVA", "scores_gsva_fulltranscriptome.tsv")]:
    S = pd.read_csv(os.path.join(SC, fn), sep="\t", index_col=0)
    S.columns = [c.strip() for c in S.columns]
    for tag, mfn in [("noQN", "P2T5b_mcpcounter_tcga_noQN.tsv"),
                     ("QN", "P2T5b_mcpcounter_tcga_QN.tsv")]:
        M = pd.read_csv(os.path.join(MCP, mfn), sep="\t", index_col=0)
        M.columns = [str(c).strip()[:15] for c in M.columns]
        common = [c for c in S.columns if c in M.columns and c in malig.index]
        Mc = M.loc[:, common]
        Zc = Mc.apply(lambda r: (r - r.mean()) / (r.std(ddof=0) + 1e-12), axis=1)
        imm = Zc.loc[[p for p in Mc.index if p not in ("Fibroblasts", "Endothelial cells")]].mean(axis=0)
        pur = malig[common]
        x = S.loc["SEG1_ADIPOKINE_CORE", common]

        def pcor(y, covs, k):
            n = len(common)
            Cm = np.column_stack([np.ones(n)] + list(covs))
            rx = x.rank().values; ry = y.rank().values
            ex = rx - Cm @ np.linalg.lstsq(Cm, rx, rcond=None)[0]
            ey = ry - Cm @ np.linalg.lstsq(Cm, ry, rcond=None)[0]
            r, _ = stats.pearsonr(ex, ey)
            t = r * np.sqrt((n - 2 - k) / max(1e-12, 1 - r ** 2))
            return r, 2 * stats.t.sf(abs(t), n - 2 - k)

        for pop in Mc.index:
            y = Mc.loc[pop, common]
            r0, p0 = stats.spearmanr(x, y)
            r1, p1 = pcor(y, [imm.values], 1)
            r2, p2 = pcor(y, [imm.values, pur.values], 2)
            r3, p3 = pcor(y, [pur.values], 1)
            rows.append(dict(algorithm=alg, mcp=tag, target=pop, n=len(common),
                             rho=r0, p=p0, rho_ctrl_imm=r1, p_ctrl_imm=p1,
                             rho_ctrl_imm_pur=r2, p_ctrl_imm_pur=p2,
                             rho_ctrl_pur=r3, p_ctrl_pur=p3))

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "P2T20_B_malignant_adjusted.tsv"), sep="\t", index=False)

L = []
L.append("# P0-b　以恶性细胞含量为纯度代理的偏相关（TCGA-SKCM，SEG1）\n")
L.append("> **纯度代理**：黑色素细胞谱系标记（TYR／MLANA／PMEL／DCT／MITF／S100B／TYRP1／SLC45A2）"
         "各自 z 分后取均值，越高＝肿瘤细胞比例越高。\n")
L.append("> **控制策略**：①免疫限定含量（MCP 去 Fibro/Endo）②恶性含量。**刻意不控基质**——基质正是被检验的信号，"
         "控了就是循环论证。\n")
for alg in ["ssGSEA", "GSVA"]:
    for tag in ["noQN", "QN"]:
        s = R[(R.algorithm == alg) & (R.mcp == tag)]
        L.append(f"\n### {alg}｜MCP-{tag}（n={int(s.n.iloc[0])}）\n")
        L.append("| 群体 | 原始 ρ | p | 控免疫含量 | 控恶性含量 | **控免疫＋恶性** | p | 保留 |")
        L.append("|---|---|---|---|---|---|---|---|")
        for _, r in s.iterrows():
            L.append(f"| {r.target} | {r.rho:+.3f} | {r.p:.1e} | {r.rho_ctrl_imm:+.3f} | "
                     f"{r.rho_ctrl_pur:+.3f} | **{r.rho_ctrl_imm_pur:+.3f}** | {r.p_ctrl_imm_pur:.1e} | "
                     f"{(r.rho_ctrl_imm_pur / r.rho * 100 if r.rho else float('nan')):.0f}% |")
io.open(os.path.join(OUT, "P2T20_B_malignant_adjusted.md"), "w", encoding="utf-8").write("\n".join(L))
print("\n".join(L))
