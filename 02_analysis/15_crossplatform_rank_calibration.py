# -*- coding: utf-8 -*-
"""② 跨平台秩校准 + 非黑色素瘤细胞系对照

设计要点
--------
1. RPKM 与 RSEM 的绝对阈值（≥1）不可比 -> 放弃绝对阈值，改用**秩**。
2. 秩在平台内计算，且**两侧使用同一个 308 基因宇宙**（符合 ALMP v1.11-J「两侧同一宇宙」），
   因此 RPKM/RSEM 的归一化差异被完全消掉（秩对单调变换不变）。
3. 对照：非黑色素瘤 CCLE 细胞系，检验 7 个「被检出」核心基因是否为谱系通用基因。
"""
import sys, os, json, time
import numpy as np
import pandas as pd

sys.path.insert(0, r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\tmp")
from netlib import get, post, API

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
OUT = os.path.join(ROOT, "P5-T2_跨平台秩校准v1")
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
OD = os.path.join(OUT, "data")

# ---------- 0. 基因与轴定义
s1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v2.csv"))
gm = json.load(open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gm["entrez2symbol"].items()}
sym2e = {v: k for k, v in e2s.items()}
prog = [g for g in s1.gene_symbol if g in sym2e]
ax = dict(zip(s1.gene_symbol, s1.axis))
seg = dict(zip(s1.gene_symbol, s1.segment_primary))
print("programme genes mapped:", len(prog))

core = [g for g in prog if seg[g] == "SEG1 adipokine core"]
intr = [g for g in prog if ax[g] == "tumour-intrinsic"]
comp = [g for g in prog if ax[g] == "tissue-compositional"]
print("core", len(core), "intrinsic", len(intr), "compositional", len(comp))

# ---------- 1. 已有矩阵
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)
CC = cc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                              values="rpkm", aggfunc="mean")
tc = pd.read_csv(r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\data\cbio_skcm\mrna_data.tsv",
                 sep="\t")
tc["symbol"] = tc.entrezGeneId.astype(int).map(e2s)
TC = tc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                              values="value", aggfunc="mean")
print("CCLE mel matrix", CC.shape, "TCGA", TC.shape)

# ---------- 2. 非黑色素瘤对照细胞系
rows = get(f"{API}/studies/ccle_broad_2019/clinical-data"
           "?clinicalDataType=SAMPLE&attributeId=CANCER_TYPE_DETAILED")
lin = {r["sampleId"]: (r["value"] or "").strip() for r in rows}
mel = set(json.load(open(os.path.join(D, "ccle_melanoma_lines.json"))))
cnt = pd.Series([v for k, v in lin.items() if k not in mel]).value_counts()
rng = np.random.default_rng(20261023)
ctrl = []
for lg, n in cnt.items():
    if n < 10:
        continue
    cand = sorted([k for k, v in lin.items() if v == lg and k not in mel])
    pick = list(rng.choice(cand, size=min(2, len(cand)), replace=False))
    ctrl += [(s, lg) for s in pick]
print("control lines:", len(ctrl), "lineages:", len({l for _, l in ctrl}))

# ---------- 3. 取对照细胞系 308 基因
sel = [s for s, _ in ctrl]
have = sorted(set(prog))
ent = [sym2e[g] for g in have]
records = []
B = 30
for i in range(0, len(ent), B):
    body = {"sampleIds": sel, "entrezGeneIds": ent[i:i + B]}
    records += post(f"{API}/molecular-profiles/ccle_broad_2019_rna_seq_mrna/"
                    f"molecular-data/fetch?projection=SUMMARY", body)
    time.sleep(0.5)
print("records fetched:", len(records))
df = pd.DataFrame(records)
df["symbol"] = df.entrezGeneId.astype(int).map(e2s)
CT = df.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                              values="value", aggfunc="mean")
CT = CT.loc[[s for s in sel if s in CT.index]]
print("control matrix", CT.shape)
CT.to_csv(os.path.join(OD, "ccle_control_308genes.tsv"), sep="\t")


def pct_within(m):
    """每个样本内，各基因在该样本 308 基因向量中的分位（0-1）。秩对 RPKM/RSEM 差异免疫。"""
    return m.rank(axis=1, pct=True, na_option="keep")


P_cc, P_tc, P_ct = pct_within(CC), pct_within(TC), pct_within(CT)


def summarise(P, tag, groups):
    out = []
    for name, genes in groups.items():
        gg = [g for g in genes if g in P.columns]
        v = P[gg].stack().values
        out.append(dict(dataset=tag, set=name, n_genes=len(gg),
                        median=float(np.median(v)),
                        q25=float(np.percentile(v, 25)), q75=float(np.percentile(v, 75)),
                        frac_upper_half=float(np.mean(v >= 0.5))))
    return out


groups = {"adipokine core (22)": core, "intrinsic axis (292)": intr,
          "compositional axis (16)": comp, "all 308": prog}
S = []
for P, tag in [(P_cc, "melanoma cell lines (57)"), (P_tc, "TCGA-SKCM tumours (472)"),
               (P_ct, f"non-melanoma cell lines ({len(CT)})")]:
    S += summarise(P, tag, groups)
S = pd.DataFrame(S)
S.to_csv(os.path.join(OD, "P5T2_axis_percentiles.tsv"), sep="\t", index=False)
print(S.to_string(index=False))

# ---------- 4. 7 个「被检出」核心基因 + 6 个零检出基因的谱系对照
seven = ["ADIPOR1", "ADIPOR2", "NAMPT", "VEGFA", "GPNMB", "SERPINE1", "ANGPTL4"]
zero = ["ADIPOQ", "ITLN1", "LEP", "RETN", "SERPINA12", "TNF"]
G = []
for g in seven + zero:
    if g not in P_ct.columns:
        continue
    G.append(dict(gene=g, group="detected" if g in seven else "undetected",
                  mel_median=float(P_cc[g].median()) if g in P_cc else np.nan,
                  ctrl_median=float(P_ct[g].median()),
                  ctrl_q25=float(P_ct[g].quantile(.25)),
                  ctrl_q75=float(P_ct[g].quantile(.75)),
                  ctrl_frac_upper=float((P_ct[g] >= .5).mean())))
G = pd.DataFrame(G)
G.to_csv(os.path.join(OD, "P5T2_core_genes_lineage_generality.tsv"), sep="\t", index=False)
print(G.to_string(index=False))
print("saved ->", OD)
