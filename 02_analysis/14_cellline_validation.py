# -*- coding: utf-8 -*-
"""P5-T1 细胞系内在轴验证：CCLE 黑色素瘤细胞系 vs TCGA-SKCM 实体瘤"""
import json, io, os
import numpy as np
import pandas as pd
from scipy import stats

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
S1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v1.csv"))

# ---------- 1. 载入
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
gmap = json.load(io.open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gmap["entrez2symbol"].items()}
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)

tc = pd.read_csv(r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\data\cbio_skcm\mrna_data.tsv",
                 sep="\t")
tc["symbol"] = tc.entrezGeneId.astype(int).map(e2s)
tc = tc.dropna(subset=["symbol"])

print("CCLE 尺度: max=%.1f p99=%.1f" % (cc.rpkm.max(), cc.rpkm.quantile(.99)))
LOG_CC = cc.rpkm.quantile(.99) < 40
print("CCLE 视为 log2 值:", LOG_CC)

cc["val"] = np.where(LOG_CC, 2 ** cc.rpkm - 1, cc.rpkm)
tc["val"] = tc.value

# ---------- 2. 每条细胞系的可用基因数
M = cc.pivot_table(index="sampleId", columns="symbol", values="val", aggfunc="mean")
T = tc.pivot_table(index="sampleId", columns="symbol", values="val", aggfunc="mean")
print("CCLE 矩阵", M.shape, "| TCGA 矩阵", T.shape)
print("CCLE 每系基因数中位:", int(M.notna().sum(axis=1).median()))

# ---------- 3. 逐基因检出率
det_cc = (M >= 1).mean(axis=0)
det_tc = (T >= 1).mean(axis=0)
gl = pd.DataFrame({"gene_symbol": det_cc.index,
                   "detect_ccle": det_cc.values,
                   "detect_tcga": det_tc.reindex(det_cc.index).values,
                   "median_ccle": M.median(axis=0).values,
                   "median_tcga": T.median(axis=0).reindex(det_cc.index).values})
gl = gl.merge(S1[["gene_symbol", "axis", "segment_primary", "HPA_adipose_nTPM",
                  "adipose_to_nonadipose_ratio", "meets_adipose_criterion"]],
              on="gene_symbol", how="left")
gl["delta_detect"] = gl.detect_tcga - gl.detect_ccle
gl = gl.sort_values(["axis", "detect_ccle"], ascending=[True, False])
gl.to_csv(os.path.join(D, "P5T1_gene_level_detection.tsv"), sep="\t", index=False)

# ---------- 4. 轴级汇总
def seg(df, name):
    return dict(segment=name, n=len(df),
                ccle_med=df.detect_ccle.median(), ccle_mean=df.detect_ccle.mean(),
                tcga_med=df.detect_tcga.median(),
                frac_expressed_50=float((df.detect_ccle >= .5).mean()),
                delta_med=df.delta_detect.median())

rows = [seg(gl[gl.axis == "tissue-compositional"], "tissue-compositional (16)"),
        seg(gl[gl.axis == "tumour-intrinsic"], "tumour-intrinsic (292)"),
        seg(gl[gl.segment_primary == "SEG1 adipokine core"], "SEG1 adipokine core (22)"),
        seg(gl[gl.segment_primary == "SEG2 lipid pathway"], "SEG2 lipid pathway (247)"),
        seg(gl[gl.segment_primary == "SEG3 executor"], "SEG3 executor (39)")]
sm = pd.DataFrame(rows)
sm.to_csv(os.path.join(D, "P5T1_axis_summary.tsv"), sep="\t", index=False)
print("\n=== 轴级检出率 ===")
print(sm.to_string(index=False))

# ---------- 5. 置换检验：核心集 22 基因是否异常低
core = gl[gl.segment_primary == "SEG1 adipokine core"].detect_ccle.values
rest = gl[gl.segment_primary != "SEG1 adipokine core"].detect_ccle.values
rng = np.random.default_rng(20261023)
null = np.array([rng.choice(rest, size=len(core), replace=False).mean() for _ in range(10000)])
obs = core.mean()
p_emp = (null <= obs).mean()
print("\n核心集 CCLE 平均检出率 %.3f vs 其余基因 %.3f | 置换 p=%.4f (单侧, 10000 次)"
      % (obs, rest.mean(), p_emp))
intr = gl[gl.axis == "tumour-intrinsic"].detect_ccle.values
print("内在轴平均检出率 %.3f" % intr.mean())
u = stats.mannwhitneyu(core, intr, alternative="less")
print("核心 vs 内在 逐基因 Mann-Whitney U p=%.3g" % u.pvalue)

# ---------- 6. 跨样本相关：核心 vs 内在（细胞系 vs 实体瘤）
def axis_score(mat, genes):
    g = [x for x in genes if x in mat.columns]
    return np.log1p(mat[g].clip(lower=0)).mean(axis=1), len(g)

genes_core = gl[gl.segment_primary == "SEG1 adipokine core"].gene_symbol.tolist()
genes_intr = gl[gl.axis == "tumour-intrinsic"].gene_symbol.tolist()
sc_cc, n1 = axis_score(M, genes_core); si_cc, n2 = axis_score(M, genes_intr)
sc_tc, _ = axis_score(T, genes_core); si_tc, _ = axis_score(T, genes_intr)
r_cc = stats.spearmanr(sc_cc, si_cc)
r_tc = stats.spearmanr(sc_tc, si_tc)
print("\n核心×内在（细胞系 n=%d 系，用 %d/%d 基因）: rho=%.3f p=%.3g"
      % (len(M), n1, n2, r_cc.correlation, r_cc.pvalue))
print("核心×内在（TCGA n=%d 样本）: rho=%.3f p=%.3g"
      % (len(T), r_tc.correlation, r_tc.pvalue))

# ---------- 7. 细胞系“纯瘤”前提验证（对照基因）
fp = os.path.join(D, "ccle_melanoma_controlgens_rpkm.tsv")
if os.path.exists(fp):
    ctl = pd.read_csv(fp, sep="\t")
    ctl["val"] = np.where(LOG_CC, 2 ** ctl.rpkm - 1, ctl.rpkm)
    print("\n=== 对照基因（细胞系中位表达 RPKM）===")
    s = ctl.groupby(["group", "symbol"]).val.median().reset_index()
    for grp in ["melanocyte", "stromal", "immune", "housekeeping"]:
        sub = s[s.group == grp].sort_values("val", ascending=False)
        print(f"  {grp:13s}", ", ".join(f"{r.symbol}={r.val:.0f}" for r in sub.itertuples()))
    s.to_csv(os.path.join(D, "P5T1_control_genes_ccle.tsv"), sep="\t", index=False)

# ---------- 8. 输出名单
ex = gl[(gl.segment_primary == "SEG1 adipokine core") & (gl.detect_ccle >= .5)]
no = gl[(gl.segment_primary == "SEG1 adipokine core") & (gl.detect_ccle < .2)]
print("\n核心集在细胞系中表达(≥50%%系) %d 个: %s" % (len(ex), ", ".join(ex.gene_symbol)))
print("核心集在细胞系中不表达(<20%%) %d 个: %s" % (len(no), ", ".join(no.gene_symbol)))
ca = gl[(gl.axis == "tumour-intrinsic") & (gl.detect_ccle >= .8)]
print("内在轴中细胞系普遍表达(≥80%%) %d 个" % len(ca))
gl[["gene_symbol", "detect_ccle", "detect_tcga"]].to_csv(
    os.path.join(D, "P5T1_cell_autonomous_intrinsic_list.tsv"), sep="\t", index=False)
json.dump({"p_perm": float(p_emp), "core_mean": float(obs), "rest_mean": float(rest.mean()),
           "intr_mean": float(intr.mean()), "rho_ccle": float(r_cc.correlation),
           "p_ccle": float(r_cc.pvalue), "rho_tcga": float(r_tc.correlation),
           "p_tcga": float(r_tc.pvalue), "n_lines": int(len(M)),
           "n_tcga": int(len(T)), "log_ccle": bool(LOG_CC)},
          io.open(os.path.join(D, "P5T1_stats.json"), "w"), indent=1)

# ---------- 9. 阈值敏感性
print("\n=== 检出阈值敏感性（核心 vs 内在，细胞系检出率中位）===")
for th in [0.5, 1.0, 2.0, 5.0]:
    d_cc = (M >= th).mean(axis=0)
    c = gl.set_index("gene_symbol").loc[d_cc.index]
    core_t = d_cc[c.segment_primary == "SEG1 adipokine core"]
    intr_t = d_cc[c.axis == "tumour-intrinsic"]
    print(f"  阈值 RPKM≥{th}: 核心 {core_t.median():.3f} | 内在 {intr_t.median():.3f} "
          f"| Δ {core_t.median()-intr_t.median():+.3f}")

# ---------- 10. 内在轴中「细胞系不表达」的成员（提示基质来源误分类）
low = gl[(gl.axis == "tumour-intrinsic") & (gl.detect_ccle < 0.2)].sort_values("detect_ccle")
low.to_csv(os.path.join(D, "P5T1_intrinsic_low_detection.tsv"), sep="\t", index=False)
print("\n内在轴中细胞系检出<20%% 的成员 %d 个:" % len(low))
print("  " + ", ".join(low.gene_symbol.tolist()))
hi = gl[(gl.axis == "tumour-intrinsic") & (gl.detect_ccle >= 0.9)]
print("内在轴中细胞系检出≥90%% 的成员 %d 个（细胞自主子集）" % len(hi))
hi[["gene_symbol", "detect_ccle", "detect_tcga", "segment_primary"]].to_csv(
    os.path.join(D, "P5T1_cell_autonomous_intrinsic_subset.tsv"), sep="\t", index=False)

print("\n保存完毕:", D)
