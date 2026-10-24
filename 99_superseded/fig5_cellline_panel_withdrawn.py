# -*- coding: utf-8 -*-
"""图：细胞系验证候选面板（v2，重排版式）"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 6.4,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
FIG = os.path.join(ROOT, "figures")
C = dict(blue="#2166a5", red="#b8232f", orange="#e08214", teal="#2a9d8f",
         grey="#7f7f7f", purple="#6b4c9a", light="#d9d9d9", ink="#40566b")

gl = pd.read_csv(os.path.join(D, "P5T1_gene_level_detection.tsv"), sep="\t")
core = gl[gl.segment_primary == "SEG1 adipokine core"].sort_values("detect_ccle")
sm = pd.read_csv(os.path.join(D, "P5T1_axis_summary.tsv"), sep="\t")
ctl = pd.read_csv(os.path.join(D, "P5T1_control_genes_ccle.tsv"), sep="\t")
sc = {"t_all": ("P5T1_scatter_all.tsv", "bulk tumours (n = 472)", "+0.60", C["blue"]),
      "c_all": ("P5T1_scatter_ccle_all.tsv", "cell lines (n = 57)", "\u22120.19", C["red"]),
      "t_sub": ("P5T1_scatter_cell_autonomous.tsv", "bulk tumours (n = 472)", "\u22120.10", C["teal"])}
DS = {k: pd.read_csv(os.path.join(D, v[0]), sep="\t") for k, v in sc.items()}

IN = 1 / 25.4
W, H = 180, 150
fig = plt.figure(figsize=(W * IN, H * IN))


def txt(x, y, s, **kw):
    fig.text(x / W, 1 - y / H, s, **kw)


def block(letter, title, x, y, subs=()):
    txt(x, y, letter, fontsize=8.4, fontweight="bold", color="black")
    txt(x + 5.0, y, title, fontsize=6.8, fontweight="bold", color="black")
    for k, s in enumerate(subs):
        txt(x, y + 3.5 + 2.7 * k, s, fontsize=4.8, color=C["ink"], style="italic")


def A(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def thin(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


# ---------------- a
axa = A(26, 30, 58, 56)
y = np.arange(len(core))
axa.hlines(y, core.detect_ccle, core.detect_tcga, color=C["light"], lw=1.5, zorder=1)
axa.scatter(core.detect_tcga, y, s=9, facecolor="white", edgecolor=C["blue"], lw=0.7, zorder=2)
axa.scatter(core.detect_ccle, y, s=9, color=C["red"], zorder=3)
axa.set_yticks(y)
axa.set_yticklabels(core.gene_symbol, fontsize=4.9, fontstyle="italic")
axa.set_xlim(-0.04, 1.04)
axa.set_ylim(-0.9, len(core) - 0.1)
axa.set_xlabel("fraction detected (RPKM / RSEM \u2265 1)", fontsize=5.2, labelpad=2)
axa.tick_params(axis="x", labelsize=4.8, length=2, pad=1)
axa.tick_params(axis="y", length=0, pad=2)
thin(axa)
txt(60, 15.6, "\u25cf", fontsize=5.4, color=C["blue"])
txt(63, 15.6, "bulk tumours", fontsize=4.8, color="black")
txt(60, 18.4, "\u25cf", fontsize=5.4, color=C["red"])
txt(63, 18.4, "melanoma cell lines", fontsize=4.8, color="black")
block("a", "The adipokine core is largely not transcribed by melanoma cells", 3, 12,
      ["detection is computed within each dataset, so the",
       "comparison is per gene and not a score"])

# ---------------- b
axb = A(103, 30, 34, 56)
xx = np.arange(len(sm))
axb.vlines(xx, sm.ccle_med, sm.tcga_med, color=C["light"], lw=2.0, zorder=1)
axb.scatter(xx - 0.10, sm.tcga_med, s=16, facecolor="white", edgecolor=C["blue"], lw=0.8, zorder=2)
axb.scatter(xx + 0.10, sm.ccle_med, s=16, color=C["red"], zorder=3)
axb.set_xticks(xx)
axb.set_xticklabels(["comp.\n16", "intr.\n288", "core\n22", "lipid\n243", "exec.\n39"], fontsize=4.6)
axb.set_ylim(-0.03, 1.06)
axb.set_xlim(-0.6, len(sm) - 0.4)
axb.set_ylabel("median fraction detected", fontsize=5.2, labelpad=2)
axb.tick_params(labelsize=4.8, length=2, pad=1)
thin(axb)
block("b", "Axis-wise medians", 103, 12, ["same data, grouped by axis"])

# ---------------- c（两张散点）
for k, (key, x0) in enumerate([("t_all", 26), ("c_all", 50)]):
    ax = A(x0, 106, 21, 32)
    d = DS[key]
    ax.scatter(d.core, d.intrinsic, s=1.6, color=sc[key][3], alpha=.55, lw=0)
    ax.set_title(f"{sc[key][1]}\n\u03c1 = {sc[key][2]}", fontsize=4.8, pad=2)
    ax.tick_params(labelsize=4.0, length=1.5, pad=1)
    ax.set_xlabel("core score", fontsize=4.6, labelpad=1)
    if k == 0:
        ax.set_ylabel("intrinsic score", fontsize=4.6, labelpad=1)
    thin(ax)
block("c", "Core\u2013intrinsic coupling exists only in tissue", 3, 94,
      ["all 288 intrinsic genes: strongly positive in tumours,",
       "absent in pure malignant cells"])

# ---------------- d（细胞自主子集）
axd = A(80, 106, 21, 32)
d = DS["t_sub"]
axd.scatter(d.core, d.intrinsic, s=1.6, color=C["teal"], alpha=.55, lw=0)
axd.set_title("cell-autonomous subset (166)\nbulk tumours, \u03c1 = \u22120.10", fontsize=4.8, pad=2)
axd.set_xlabel("core score", fontsize=4.6, labelpad=1)
axd.set_ylabel("intrinsic score", fontsize=4.6, labelpad=1)
axd.tick_params(labelsize=4.0, length=1.5, pad=1)
thin(axd)
block("d", "Cell-autonomous subset", 74, 94,
      ["the tumour-level positive correlation disappears"])

# ---------------- e（对照基因）
axe = A(112, 106, 56, 32)
order = ["housekeeping", "melanocyte", "stromal", "immune"]
cols = {"housekeeping": C["grey"], "melanocyte": C["purple"], "stromal": C["orange"],
        "immune": C["teal"]}
lab = {"housekeeping": "housekeeping", "melanocyte": "melanocytic", "stromal": "stromal",
       "immune": "immune"}
rng = np.random.default_rng(7)
for i, g in enumerate(order):
    sub = ctl[ctl.group == g].val.clip(lower=.01)
    axe.scatter(np.full(len(sub), i) + rng.normal(0, .045, len(sub)),
                np.log10(sub), s=10, color=cols[g], alpha=.8, lw=0)
    axe.hlines(np.log10(sub).median(), i - .26, i + .26, color="black", lw=1.0, zorder=4)
axe.set_xticks(range(4))
axe.set_xticklabels([lab[o] for o in order], fontsize=4.8)
axe.set_ylabel("log$_{10}$ RPKM", fontsize=5.2, labelpad=2)
axe.set_ylim(-2.4, 3.7)
axe.tick_params(labelsize=4.6, length=2, pad=1)
axe.axhline(0, color=C["grey"], lw=0.4, ls=":")
axe.text(-0.42, 0.12, "1 RPKM", fontsize=4.2, color=C["grey"], ha="left")
thin(axe)
block("e", "No stroma, no immune signal", 112, 94,
      ["black bars = medians; the premise that makes panel a interpretable"])

for ext in ["pdf", "png"]:
    fig.savefig(os.path.join(FIG, f"fig5_panel_cellline_candidate.{ext}"), dpi=300)
print("saved")
