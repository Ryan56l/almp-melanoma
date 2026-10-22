# -*- coding: utf-8 -*-
"""
y51: TCGA-SKCM 全转录组口径下 SEG1 x MCP-counter 耦合（补 P2-T6b-R 的 TCGA 锚点）
目的：图 3a 五队列森林图中 TCGA 一行的 rho 必须与 ICB 四队列同宇宙（SAP v1.11 条款 J）。
输入：
  C:/Rtmp/p2t12/results50/scores_ssgsea_fulltranscriptome.tsv  (33 x 472, 全转录组 ssGSEA)
  C:/Rtmp/p2t12/results50/scores_gsva_fulltranscriptome.tsv    (33 x 472, 全转录组 GSVA)
  output/P2-T5b_多去卷积一致性v1/data/P2T5b_mcpcounter_tcga_noQN.tsv
输出：
  output/P2-T6b-R_五队列基因宇宙复查v1/data/P2T6bR_tcga_full_coupling.csv
"""
import os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
OUT = os.path.join(ROOT, "output/P2-T6b-R_五队列基因宇宙复查v1/data")
P50 = r"C:/Rtmp/p2t12/results50"
MCP = os.path.join(ROOT, "output/P2-T5b_多去卷积一致性v1/data/P2T5b_mcpcounter_tcga_noQN.tsv")


def read_scores(fn):
    df = pd.read_csv(fn, sep="\t", index_col=0)
    df.columns = [str(c)[:15] for c in df.columns]   # 截到 15 位，与 MCP 侧对齐
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def read_mcp(fn):
    df = pd.read_csv(fn, sep="\t", index_col=0)
    df.columns = [str(c)[:15] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return df


mcp = read_mcp(MCP)
mcp.columns = [str(c).replace("A", "") if False else c for c in mcp.columns]
# MCP 侧列名形如 TCGA-ER-A199-06A（16 位），截 15 位后与评分侧 TCGA-ER-A199-06 一致
rows = []
for alg, fn in (("ssGSEA", os.path.join(P50, "scores_ssgsea_fulltranscriptome.tsv")),
                ("GSVA", os.path.join(P50, "scores_gsva_fulltranscriptome.tsv"))):
    sc = read_scores(fn).loc["SEG1_ADIPOKINE_CORE"]
    sh = [s for s in mcp.columns if s in sc.index]
    print(f"[{alg}] 交集样本 {len(sh)} / MCP {mcp.shape[1]} / 评分 {sc.shape[0]}")
    a_all = sc.loc[sh].astype(float)
    for pop in mcp.index:
        v = mcp.loc[pop, sh].astype(float)
        ok = np.isfinite(v) & np.isfinite(a_all)
        r, p = stats.spearmanr(a_all[ok], v[ok])
        rows.append(dict(algorithm=alg, population=pop, n=int(ok.sum()), rho=round(r, 4), p=p))
out = pd.DataFrame(rows)
out["q"] = np.nan
for alg, g in out.groupby("algorithm"):
    out.loc[g.index, "q"] = np.round(stats.false_discovery_control(g["p"], method="bh"), 5) \
        if hasattr(stats, "false_discovery_control") else np.nan
os.makedirs(OUT, exist_ok=True)
out.to_csv(os.path.join(OUT, "P2T6bR_tcga_full_coupling.csv"), index=False, encoding="utf-8-sig")
pd.set_option("display.width", 200)
print(out.sort_values(["algorithm", "rho"], ascending=[True, False]).to_string(index=False))
