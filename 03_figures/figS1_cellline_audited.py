# -*- coding: utf-8 -*-
"""Figure S1 v2 —— 按谱系审计后的细胞系面板重绘

相对 v1 的变更：
  a 三组对照（50 黑素瘤 / 7 成纤维 / 472 TCGA），不再把污染当背景
  b 轴中位数同三组
  c 核心—内在耦合：肿瘤 +0.60、黑素瘤系 +0.06（旧图画的是含污染的 −0.19）
  d 166 细胞自主子集：肿瘤 −0.10、黑素瘤系 −0.29（旧图 −0.46 为污染值）
  e 对照转录本改为按谱系并列（黑素细胞标记只在 50 株，基质标记只在 7 株）
  f 新增：核心 22 基因成纤维/黑素瘤倍数（log2），显示核心按细胞来源可拆
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 6.4,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
T2 = os.path.join(ROOT, "P5-T2_跨平台秩校准v1", "data")
FIG = os.path.join(ROOT, "figures")
C = dict(blue="#2166a5", red="#b8232f", orange="#e08214", teal="#2a9d8f",
         grey="#7f7f7f", purple="#6b4c9a", light="#d9d9d9", ink="#40566b")

# ---------- 数据
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
hi = pd.read_csv(os.path.join(D, "P5T1_cell_autonomous_intrinsic_subset_v2.tsv"), sep="\t")
ca = [g for g in hi.gene_symbol if g in T.columns and g in P.columns]
print("lines", len(mel), len(fib), "| tumours", T.shape)


def det(sub, genes):
    return (sub[[g for g in genes if g in sub.columns]] >= 1).mean(axis=0)


def score(sub, genes):
    g = [g for g in genes if g in sub.columns]
    return np.log1p(sub[g].clip(lower=0)).mean(axis=1)


def rho(a, b):
    r = stats.spearmanr(a, b)
    return float(r.correlation), float(r.pvalue)


# ---------- 绘图
IN = 1 / 25.4
W, H = 180, 196
fig = plt.figure(figsize=(W * IN, H * IN))


def txt(x, y, s, **kw):
    fig.text(x / W, 1 - y / H, s, **kw)


def block(letter, title, x, y, subs=()):
    txt(x, y, letter, fontsize=8.4, fontweight="bold", color="black")
    txt(x + 5.0, y, title, fontsize=6.8, fontweight="bold", color="black")
    for k, s in enumerate(subs):
        txt(x, y + 3.4 + 2.6 * k, s, fontsize=4.8, color=C["ink"], style="italic")


def A(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def thin(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


# ================= a 逐基因检出
d_mel = det(P.loc[mel], core)
d_fib = det(P.loc[fib], core)
d_tc = det(T, core)
o = np.argsort(d_mel.values)
genes = [core[i] for i in o]
axa = A(28, 32, 56, 52)
y = np.arange(len(genes))
axa.hlines(y, d_mel.values[o], d_tc.values[o], color=C["light"], lw=1.3, zorder=1)
axa.scatter(d_tc.values[o], y, s=8, facecolor="white", edgecolor=C["blue"], lw=0.6, zorder=2)
axa.scatter(d_fib.values[o], y, s=12, color=C["orange"], marker="D", zorder=3)
axa.scatter(d_mel.values[o], y, s=8, color=C["red"], zorder=4)
axa.set_yticks(y)
axa.set_yticklabels(genes, fontsize=4.7, fontstyle="italic")
axa.set_xlim(-0.04, 1.05)
axa.set_ylim(-0.9, len(genes) - 0.1)
axa.set_xlabel("fraction of samples with RPKM / RSEM \u2265 1", fontsize=5.1, labelpad=2)
axa.tick_params(axis="x", labelsize=4.7, length=2, pad=1)
axa.tick_params(axis="y", length=0, pad=2)
thin(axa)
txt(4, 18.4, "\u25cf", fontsize=5.2, color=C["blue"]); txt(7, 18.4, "tumours (472)", fontsize=4.7)
txt(34, 18.4, "\u25cf", fontsize=5.2, color=C["red"]); txt(37, 18.4, "melanoma lines (50)", fontsize=4.7)
txt(74, 18.4, "\u25c6", fontsize=4.6, color=C["orange"]); txt(77, 18.4, "fibroblast lines (7)", fontsize=4.7)
block("a", "The core is a fibroblast product, not a melanoma one", 3, 10,
      ["detection per gene and per dataset; no gene universe enters"])

# ================= b 轴中位数（三组）
sm = pd.read_csv(os.path.join(D, "P5T1_axis_summary.tsv"), sep="\t")
seg = {"tissue-compositional (16)": "comp.\n16", "tumour-intrinsic (292)": "intr.\n288",
       "SEG1 adipokine core (22)": "core\n22", "SEG2 lipid pathway (243)": "lipid\n243",
       "SEG3 executor (39)": "exec.\n39"}
s1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v2.csv"))
axs = dict(zip(s1.gene_symbol, s1.axis)); sg = dict(zip(s1.gene_symbol, s1.segment_primary))
groups = [("comp.", [g for g in s1.gene_symbol if axs.get(g) == "tissue-compositional"]),
          ("intr.", [g for g in s1.gene_symbol if axs.get(g) == "tumour-intrinsic"]),
          ("core", [g for g in s1.gene_symbol if sg.get(g) == "SEG1 adipokine core"]),
          ("lipid", [g for g in s1.gene_symbol if sg.get(g) == "SEG2 lipid pathway"]),
          ("exec.", [g for g in s1.gene_symbol if sg.get(g) == "SEG3 executor"])]
axb = A(106, 32, 30, 52)
xx = np.arange(len(groups))
vm = [float(det(P.loc[mel], g).median()) for _, g in groups]
vf = [float(det(P.loc[fib], g).median()) for _, g in groups]
vt = [float(det(T, g).median()) for _, g in groups]
axb.vlines(xx, vm, vt, color=C["light"], lw=1.8, zorder=1)
axb.scatter(xx - 0.13, vt, s=15, facecolor="white", edgecolor=C["blue"], lw=0.7, zorder=2)
axb.scatter(xx + 0.06, vf, s=13, color=C["orange"], marker="D", zorder=3)
axb.scatter(xx + 0.17, vm, s=14, color=C["red"], zorder=4)
axb.set_xticks(xx)
axb.set_xticklabels([n for n, _ in groups], fontsize=4.6)
axb.set_ylim(-0.03, 1.06)
axb.set_xlim(-0.6, len(groups) - 0.4)
axb.set_ylabel("median fraction detected", fontsize=5.1, labelpad=2)
axb.tick_params(labelsize=4.6, length=2, pad=1)
thin(axb)
block("b", "Axis-wise medians", 106, 10)

# ================= c/d 散点四联
panels = [("tumours (472)\nall 288 intrinsic", T, intr, C["blue"]),
          ("melanoma lines (50)\nall 288 intrinsic", P.loc[mel], intr, C["red"]),
          ("tumours (472)\n171 cell-autonomous", T, ca, C["teal"]),
          ("melanoma lines (50)\n171 cell-autonomous", P.loc[mel], ca, C["purple"])]
xs = [26, 49, 73, 96]
for k, (lab, sub, genes, col) in enumerate(panels):
    ax = A(xs[k], 108, 20, 24)
    a_, b_ = score(sub, core), score(sub, genes)
    r, p = rho(a_, b_)
    ax.scatter(a_, b_, s=1.5, color=col, alpha=.55, lw=0)
    ps = "p = %.3f" % p if p >= 1e-3 else "p = %.0e" % p
    ax.set_title(f"{lab}\n\u03c1 = {r:+.2f}, {ps}", fontsize=4.5, pad=2)
    ax.tick_params(labelsize=3.8, length=1.4, pad=1)
    ax.set_xlabel("core score", fontsize=4.4, labelpad=1)
    ax.set_ylabel("intrinsic score", fontsize=4.4, labelpad=1)
    thin(ax)
block("c", "Coupling with the intrinsic axis is a property of tissue, not of melanoma cells",
      3, 98, ["scores are means of log(1 + expression) over the genes of each set"])

# ================= e 对照转录本（按谱系）
ctl = pd.read_csv(os.path.join(D, "ccle_melanoma_controlgens_rpkm.tsv"), sep="\t")
CG = ctl.pivot_table(index="sampleId", columns="symbol", values="rpkm", aggfunc="mean")
ordg = ["S100B", "TYR", "SOX10", "MLANA", "PMEL", "DCT", "COL1A1", "COL1A2", "DCN",
        "LUM", "THY1", "PDGFRB", "PTPRC", "CD3D", "CD8A"]
ordg = [g for g in ordg if g in CG.columns]
axe = A(122, 108, 52, 24)
xm = CG.loc[[s for s in mel if s in CG.index], ordg].median()
xf = CG.loc[[s for s in fib if s in CG.index], ordg].median()
yy = np.arange(len(ordg))
axe.hlines(yy, np.log10(xm.clip(lower=.02)), np.log10(xf.clip(lower=.02)),
           color=C["light"], lw=1.3, zorder=1)
axe.scatter(np.log10(xm.clip(lower=.02)), yy, s=8, color=C["red"], zorder=3)
axe.scatter(np.log10(xf.clip(lower=.02)), yy, s=11, color=C["orange"], marker="D", zorder=4)
axe.set_yticks(yy)
axe.set_yticklabels(ordg, fontsize=4.0, fontstyle="italic")
axe.axvline(0, color=C["grey"], lw=0.4, ls=":")
axe.set_xlabel("log$_{10}$ RPKM (median)", fontsize=4.9, labelpad=1)
axe.tick_params(labelsize=4.0, length=1.6, pad=1)
thin(axe)
for y0, y1, lab, col in [(0.0, 5.0, "melanocytic", C["purple"]),
                         (6.0, 11.0, "collagen", C["orange"]),
                         (12.0, 14.0, "lymphoid", C["teal"])]:
    axe.annotate("", xy=(-3.55, y0), xytext=(-3.55, y1),
                 arrowprops=dict(arrowstyle="-", lw=1.4, color=col))
    axe.text(-3.62, (y0 + y1) / 2, lab, fontsize=4.0, color=col, style="italic",
             rotation=90, ha="right", va="center")
axe.set_xlim(-4.3, 3.95)
block("e", "Lineage audit of the same 57 lines", 119, 98,
      ["melanocytic markers survive only in the 50, collagen only in the 7"])

# ================= f 核心基因 成纤维/黑素瘤 倍数
G = pd.read_csv(os.path.join(T2, "P5T2_core_genes_fibroblast_vs_melanoma.tsv"), sep="\t")
G = G.sort_values("fib_over_mel", ascending=False)
axf = A(26, 156, 146, 26)
xxf = np.arange(len(G))
cols = [C["orange"] if v > 0 else C["red"] for v in G.fib_over_mel]
axf.bar(xxf, G.fib_over_mel, color=cols, width=0.72, lw=0)
axf.axhline(0, color="black", lw=0.5)
axf.set_xticks(xxf)
axf.set_xticklabels(G.gene, rotation=90, fontsize=4.2, fontstyle="italic")
axf.set_ylabel("log$_2$ (fibroblast / melanoma)", fontsize=4.9, labelpad=2)
axf.set_xlim(-0.7, len(G) - 0.3)
axf.set_ylim(-4.2, 9.0)
axf.tick_params(axis="y", labelsize=4.2, length=2, pad=1)
axf.tick_params(axis="x", length=0, pad=1)
thin(axf)
axf.annotate("", xy=(0.0, 8.3), xytext=(0.0, 9.4), arrowprops=dict(arrowstyle="-", lw=0.5))
axf.text(1.2, 8.3, "SERPINE1 978 vs 2.5 RPKM", fontsize=4.0, va="bottom")
axf.text(1.2, 7.2, "IL6 6.5 vs 0.3 RPKM", fontsize=4.0, va="bottom")
block("f", "The core splits by cell of origin: secreted members are fibroblast-derived, "
           "receptor members tumour-derived", 3, 150,
      ["orange, higher in the seven fibroblast lines; red, higher in the 50 melanoma lines"])

for ext in ["pdf", "png"]:
    fig.savefig(os.path.join(FIG, f"figS1_cellline_audited.{ext}"), dpi=300)
print("saved ->", FIG)
print("rho c:", rho(score(P.loc[mel], core), score(P.loc[mel], intr)))
