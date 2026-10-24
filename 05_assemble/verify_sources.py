# -*- coding: utf-8 -*-
"""③ 同源核查：Figure S1 六面板、S5 v2、S2 的 e/f/g 块是否同一套数字

同时重算秩分位——原 P5T2_axis_percentiles.tsv 用的是全部 57 株（含 7 株成纤维），
与正文 50 株口径不一致，须重算并核对。
"""
import os, json
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
T2 = os.path.join(ROOT, "P5-T2_跨平台秩校准v1", "data")

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
tc = pd.read_csv(os.path.join(os.path.expanduser("~"), "WorkBuddy", "2026-10-07-22-15-36",
                              ".workbuddy", "data", "cbio_skcm", "mrna_data.tsv"), sep="\t")
tc["symbol"] = tc.entrezGeneId.astype(int).map(e2s)
T = tc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                             values="value", aggfunc="mean")
U = [g for g in P.columns if g in T.columns]   # 两侧同一宇宙

print("宇宙交集基因数:", len(U), "| 核心在其中:", len([g for g in core if g in U]))


def pct(sub):
    """样本内 308 基因宇宙分位"""
    r = sub[U].rank(axis=1, pct=True)
    return r


res = {}
for tag, sub in [("mel50", P.loc[mel]), ("fib7", P.loc[fib]), ("all57", P), ("tcga", T)]:
    r = sub[U].rank(axis=1, pct=True, na_option="keep")
    cc_ = [g for g in core if g in U]; ii = [g for g in intr if g in U]
    res[tag] = dict(core_genemed=float(r[cc_].median(axis=0).median()),
                    core_flat=float(np.median(r[cc_].stack().values)),
                    intr_flat=float(np.median(r[ii].stack().values)),
                    n=len(sub))
    print(f"{tag:6s} n={len(sub):3d} core 基因间中位={res[tag]['core_genemed']:.4f} "
          f"扁平中位={res[tag]['core_flat']:.4f} intrinsic扁平={res[tag]['intr_flat']:.4f}")

CHK = []


def chk(name, value, expect, tol):
    ok = abs(value - expect) <= tol
    CHK.append((name, f"{value:.4g}", f"{expect:.4g}", "PASS" if ok else "**FAIL**"))


# ---- 1. S5 v2 逐基因
s5 = pd.read_csv(os.path.join(ROOT, "supp", "S5_ALMP_gene_detection_celllines_v2.csv"))
g = s5.set_index("gene_symbol")
core_det_mel = float((P.loc[mel, core] >= 1).mean(axis=0).median())
core_det_fib = float((P.loc[fib, core] >= 1).mean(axis=0).median())
intr_det_mel = float((P.loc[mel, intr] >= 1).mean(axis=0).median())
intr_det_fib = float((P.loc[fib, intr] >= 1).mean(axis=0).median())
chk("S5 core 检出中位 (50)", s5[s5.segment == "SEG1 adipokine core"].detection_fraction_50_melanoma.median(), core_det_mel, 1e-6)
chk("S5 core 检出中位 (7)", s5[s5.segment == "SEG1 adipokine core"].detection_fraction_7_fibroblast.median(), core_det_fib, 1e-6)
chk("S5 SERPINE1 黑素瘤中位", g.loc["SERPINE1", "median_RPKM_50_melanoma"], float(P.loc[mel, "SERPINE1"].median()), 1e-4)
chk("S5 SERPINE1 成纤维中位", g.loc["SERPINE1", "median_RPKM_7_fibroblast"], float(P.loc[fib, "SERPINE1"].median()), 1e-4)
chk("S5 IL6 倍数(原始)", g.loc["IL6", "fold_fibroblast_over_melanoma"], float(P.loc[fib, "IL6"].median() / P.loc[mel, "IL6"].median()), 1e-6)

# ---- 2. S2 e/f/g
s2 = pd.read_csv(os.path.join(ROOT, "supp", "S2_sensitivity_analyses_v1.csv"))
e = s2[s2.analysis_block.str.startswith("S2e")].set_index("level")
def eval_(lvl):
    v = e.loc[lvl, "value"]
    if hasattr(v, "iloc"):
        v = v.iloc[0]
    return float(str(v).split(" / ")[0])
chk("S2e 50 株 core/intrinsic", eval_("50 melanoma cell lines"), core_det_mel, 1e-3)
chk("S2e 7 株 core", eval_("7 fibroblast lines (positive control)"), core_det_fib, 1e-3)

Q = pd.read_csv(os.path.join(T2, "P5T2_rho_2x2_audited.tsv"), sep="\t")
f = s2[s2.analysis_block.str.startswith("S2f")]
m50 = f[(f.level.str.startswith("50")) & (f.comparison.str.contains("all 288"))].value.iloc[0]
chk("S2f 50 株 ρ(全288)", float(m50), Q[(Q.dataset == "melanoma lines (audited)") &
    (Q.gene_set == "all intrinsic")].rho.iloc[0], 5e-3)
m166 = f[(f.level.str.startswith("50")) & (f.comparison.str.contains("166"))].value.iloc[0]
chk("S2f 50 株 ρ(166)", float(m166), Q[(Q.dataset == "melanoma lines (audited)") &
    (Q.gene_set == "166 cell-autonomous")].rho.iloc[0], 5e-3)
gg = s2[s2.analysis_block.str.startswith("S2g")].set_index("comparison")
chk("S2g SERPINE1 倍数", float(gg.loc["SERPINE1 median RPKM, fibroblast / melanoma", "note"].rstrip("x")),
    float(P.loc[fib, "SERPINE1"].median() / P.loc[mel, "SERPINE1"].median()), 0.5)

# ---- 3. 秩分位：50 株 vs 旧 57 株
print("\n=== 秩分位核对 ===")
print(f"旧表（57 株，扁平）= 0.2763 ｜ 50 株扁平 = {res['mel50']['core_flat']:.4f} ｜ 7 株扁平 = {res['fib7']['core_flat']:.4f}")
print(f"肿瘤 core 扁平 = {res['tcga']['core_flat']:.4f}（旧报 0.2885）")

print("\n=== 同源核查 ===")
for r in CHK:
    print(f"{r[3]:9s} {r[0]:28s} 实测 {r[1]:>10s}  登记 {r[2]:>10s}")
bad = [r for r in CHK if "FAIL" in r[3]]
print("\n失败项:", len(bad))
pd.DataFrame(res).T.to_csv(os.path.join(T2, "P5T2_rank_percentiles_audited.tsv"), sep="\t")
print("saved -> P5T2_rank_percentiles_audited.tsv")
