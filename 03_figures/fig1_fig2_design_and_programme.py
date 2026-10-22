# -*- coding: utf-8 -*-
"""图 1（研究设计与预登记时间轴）＋ 图 2（ALMP 构建与双轴拆分）"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
OUT = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output/figures"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.6,
    "axes.linewidth": 0.7, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300,
})

C_BLUE = "#2B5C8A"; C_TEAL = "#2A7F7F"; C_ORANGE = "#C1602A"
C_GREY = "#8A8F98"; C_RED = "#A8322A"; C_LGREY = "#EDEFF2"

# ------------------------------------------------------------------ 图 1
fig = plt.figure(figsize=(180 / 25.4, 150 / 25.4))
gs = fig.add_gridspec(2, 1, height_ratios=[1.55, 1.0], hspace=0.20,
                      left=0.02, right=0.985, top=0.945, bottom=0.045)

# ---- 1a 流程
ax = fig.add_subplot(gs[0]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
ax.text(0, 100, "a", fontsize=9, fontweight="bold", va="top")
ax.text(4.2, 99.4, "Evidence hierarchy: four layers, one pre-specified audit", fontsize=7.2,
        fontweight="bold", va="top")

def box(x, y, w, h, title, lines, fc, ec, fs=6.2, tfs=6.5, dy=3.4):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.30,rounding_size=1.3",
                                facecolor=fc, edgecolor=ec, linewidth=0.9))
    ax.text(x + w / 2, y + h - 2.4, title, ha="center", va="top", fontsize=tfs,
            fontweight="bold", color=ec)
    ax.text(x + 1.7, y + h - 2.4 - dy, "\n".join(lines), ha="left", va="top", fontsize=fs,
            color="#1A1A1A", linespacing=1.45)

def arrow(x1, y1, x2, y2, c=C_GREY):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7,
                                 linewidth=0.9, color=c, shrinkA=0, shrinkB=0))

box(1, 68, 30.5, 25, "I.  Bulk transcriptomes",
    ["TCGA-SKCM  n = 472", "368 metastases / 103 primary", "primary = exploratory layer",
     "GSE65904  n = 214"], C_LGREY, C_BLUE)
box(34.7, 68, 30.5, 25, "II.  ICB cohorts",
    ["Gide 91 (73 with response)", "Liu 122  ·  Riaz 107 (49)", "UCLA 27  ·  Van Allen 40",
     "Van Allen = post hoc only"], C_LGREY, C_TEAL)
box(68.4, 68, 30.5, 25, "III.  Single cell / spatial",
    ["GSE72056  4,645 cells / 19 pts", "neutrophils not captured",
     "public spatial: supplement", "no peritumoural adipose"], C_LGREY, C_ORANGE)

box(1, 49, 45, 13.5, "Scoring",
    ["ALMP 308 genes  →  ssGSEA / GSVA on full transcriptome",
     "six segments  ·  axes: compositional 16 / intrinsic 292"], "#F4F7FA", C_BLUE)
box(49, 49, 50, 13.5, "Deconvolution",
    ["quanTIseq (absolute mode)  ·  MCP-counter (no-QN, QN)",
     "concordance reported at the lower of the two variants"], "#F2F8F8", C_TEAL)

box(1, 23, 30.5, 23, "Aim 1  —  decomposition",
    ["adipokine core = stromal readout", "4/4 cohorts, ρ +0.51 to +0.69",
     "intrinsic axis: stroma-independent", "orthogonal to TMB (36 tests)"], "#EFF6F1", "#2E6B3E")
box(34.7, 23, 30.5, 23, "Aim 2  —  causal audit",
    ["MR: leptin / adiponectin / resistin", "cis-SMR + HEIDI, 32 combinations",
     "mediation via leukocyte traits", "no supporting genetic evidence"], "#FAF0EE", C_ORANGE)
box(68.4, 23, 30.5, 23, "Aim 3/4  —  exploratory",
    ["ICB response: 48 tests, q > 0.05", "survival: one external cohort",
     "FADS2: four tissues, n.s.", "exploratory only"], "#F5F2F8", "#5B3F7A")

arrow(16, 67.6, 14, 62.9)
arrow(50, 67.6, 38, 62.9)
arrow(84, 67.6, 62, 62.9)
arrow(84, 67.6, 80, 62.9)
arrow(14, 48.6, 12, 46.4)
arrow(28, 48.6, 44, 46.4)
arrow(66, 48.6, 56, 46.4)
arrow(86, 48.6, 84, 46.4)

ax.add_patch(Rectangle((1, 3.4), 97.9, 17.0, facecolor="#FBF7EC", edgecolor="#B08A3E",
                       linewidth=0.9, linestyle=(0, (3, 2))))
ax.text(2.6, 18.5, "Pre-registration envelope  —  ten dated clauses, four raising a threshold and none lowering one",
        fontsize=6.4, fontweight="bold", color="#7A5A1E", va="top")
ax.text(2.6, 15.3,
        "multiplicative random-effects IVW fixed as the primary analysis before any MR was run  ·  ICB downgrade criteria registered before cohort expansion\n"
        "SMR + HEIDI made decisive in place of colocalisation (v1.6, correction)  ·  outcome hierarchy sealed and the unavailable meta-analysis disclosed (v1.8)\n"
        "seeded random gene-set null (v1.9)  ·  effect-scale correction with MDE obligation (v1.10)  ·  gene-universe consistency rule (v1.11)",
        fontsize=5.2, color="#3A3A3A", va="top", linespacing=1.55)
for x in (14, 50, 84):
    arrow(x, 20.8, x, 22.6, c="#B08A3E")

# ---- 1b 时间轴
ax2 = fig.add_subplot(gs[1])
ax2.set_xlim(0, 10.4); ax2.set_ylim(-1.05, 1.80); ax2.axis("off")
ax2.text(0, 1.78, "b", fontsize=9, fontweight="bold", va="top")
ax2.text(0.42, 1.74, "Pre-registration sequence and the points at which claims were withdrawn",
         fontsize=7.2, fontweight="bold", va="top")

CLAUSES = [
    ("v1.3", "frozen", "Sealing"),
    ("v1.4", "MR specification", "Raising"),
    ("v1.5", "ICB downgrade gate", "Raising"),
    ("v1.5.1", "single external cohort", "Raising"),
    ("v1.6", "SMR+HEIDI decisive", "Correction"),
    ("v1.7", "CDS2 demoted", "Lowering"),
    ("v1.8", "outcome hierarchy", "Sealing"),
    ("v1.9", "random null, seeded", "Raising"),
    ("v1.10", "effect-scale + MDE", "Raising"),
    ("v1.11", "gene universe", "Raising"),
]
COL = {"Sealing": C_BLUE, "Raising": "#2E6B3E", "Correction": C_ORANGE, "Lowering": C_RED}
ax2.plot([0.35, 10.15], [0.55, 0.55], color=C_GREY, linewidth=1.0, zorder=1)
for i, (v, lab, typ) in enumerate(CLAUSES):
    x = 0.55 + i * 0.97
    up = (i % 2 == 0)
    ax2.plot([x, x], [0.55, 0.92 if up else 0.18], color=COL[typ], linewidth=0.9, zorder=2)
    ax2.scatter([x], [0.55], s=15, color=COL[typ], zorder=3, edgecolor="white", linewidth=0.5)
    ax2.text(x, 0.96 if up else 0.13, v, fontsize=6.0, fontweight="bold", ha="center",
             va="bottom" if up else "top", color=COL[typ])
    ax2.text(x, 1.11 if up else -0.02, lab, fontsize=5.1, ha="center",
             va="bottom" if up else "top", color="#333333")

ax2.text(0.35, -0.52,
         "green  raising a threshold      blue  sealing      orange  correction (PP.H4 > 0.8 restored to the frozen value)      red  claim withdrawn",
         fontsize=5.3, color="#555555", va="top")
ax2.text(0.35, -0.78,
         "Dates from the SAP ledger: v1.5 2026-10-14 · v1.5.1 and v1.6 2026-10-19 · v1.8 2026-10-20 · v1.9–v1.11 2026-10-21.\n"
         "The two retractions reported in Discussion 4.4 (413-gene archived scores; the 89% positive-sign statistic) occurred after these clauses were sealed and are declared post hoc.",
         fontsize=5.1, color="#555555", va="top", linespacing=1.6)

fig.savefig(os.path.join(OUT, "fig1_design_hierarchy.pdf"))
fig.savefig(os.path.join(OUT, "fig1_design_hierarchy.png"), dpi=300)
plt.close(fig)
print("saved fig1")

# ------------------------------------------------------------------ 图 2
rows = list(csv.DictReader(open(os.path.join(ROOT, "output/P2-T1_基因集v1/ALMP_gene_axis_assignment_v1.csv"),
                               encoding="utf-8-sig")))
fig = plt.figure(figsize=(180 / 25.4, 82 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.20, 1.02], wspace=0.40,
                      left=0.075, right=0.982, top=0.885, bottom=0.175)

# 2a 组成
axa = fig.add_subplot(gs[0])
segs = [("Adipokine core", 22, C_ORANGE), ("Lipid pathway", 262, C_BLUE), ("Executor", 39, C_TEAL)]
left = 0
for nm, n, c in segs:
    axa.barh(1.35, n, left=left, height=0.40, color=c, edgecolor="white", linewidth=0.8)
    if n > 40:
        axa.text(left + n / 2, 1.35, str(n), ha="center", va="center", color="white",
                 fontsize=7.0, fontweight="bold")
    left += n
axa.annotate("22", (11, 1.55), fontsize=6.4, fontweight="bold", color=C_ORANGE,
             ha="center", va="bottom")
axa.set_xlim(0, 330); axa.set_ylim(0.15, 2.35); axa.set_yticks([])
axa.set_xlabel("Genes", fontsize=6.4)
axa.set_title("308 genes, six segments", fontsize=7.0, fontweight="bold", pad=5)
for yy, nm, n, c in zip((1.02, 0.78, 0.54), [s_[0] for s_ in segs], [22, 262, 39],
                        [C_ORANGE, C_BLUE, C_TEAL]):
    axa.add_patch(Rectangle((2, yy - 0.045), 9, 0.09, facecolor=c, edgecolor="none"))
    axa.text(13, yy, f"{nm} ({n})", fontsize=5.8, va="center")
axa.text(0, 0.30, "two axes: compositional 16 / intrinsic 292", fontsize=5.3,
         color="#555555", va="top")
for sp in ("top", "right"):
    axa.spines[sp].set_visible(False)

# 2b HPA 散点
axb = fig.add_subplot(gs[1])
xs, ys, cs, names = [], [], [], []
for r in rows:
    try:
        a = float(r["hpa_adipose"]); q = float(r["hpa_ratio"])
    except (TypeError, ValueError):
        continue
    if a <= 0 or q <= 0:
        continue
    xs.append(np.log10(a)); ys.append(np.log10(q))
    cs.append(C_ORANGE if r["axis"] == "组织成分轴" else C_BLUE)
    names.append(r["gene"])
xs = np.array(xs); ys = np.array(ys)
axb.scatter(xs, ys, s=5.0, c=cs, alpha=0.60, linewidth=0)
axb.axhline(np.log10(3), color=C_GREY, linewidth=0.7, linestyle=(0, (2, 2)))
axb.axvline(np.log10(10), color=C_GREY, linewidth=0.7, linestyle=(0, (2, 2)))
axb.text(-0.92, np.log10(3) + 0.12, "ratio = 3", fontsize=5.2, color=C_GREY)
axb.text(np.log10(10) + 0.06, -1.90, "nTPM = 10", fontsize=4.8, color=C_GREY,
         ha="left", va="bottom")
axb.set_xlabel("log₁₀ adipose nTPM (HPA)", fontsize=6.4)
axb.set_ylabel("log₁₀ adipose / max(non-adipose)", fontsize=6.3)
axb.set_title("Adipose-specificity criteria", fontsize=7.0, fontweight="bold", pad=5)
axb.set_xlim(-1.05, 3.55); axb.set_ylim(-2.05, 3.95)
axb.legend([Line2D([], [], marker="o", linestyle="", color=C_ORANGE, markersize=4),
            Line2D([], [], marker="o", linestyle="", color=C_BLUE, markersize=4)],
           ["Compositional axis (16)", "Intrinsic axis (292)"], fontsize=5.6, frameon=False,
           loc="upper left", handletextpad=0.2, borderpad=0.1)
OFF = {"LEP": (6, 2), "ADIPOQ": (5, -11), "FABP4": (-34, 4), "CD36": (-30, -10)}
for g, off in OFF.items():
    if g in names:
        i = names.index(g)
        axb.annotate(g, (xs[i], ys[i]), fontsize=5.0, xytext=off, textcoords="offset points",
                     color="#222222",
                     arrowprops=dict(arrowstyle="-", color="#999999", linewidth=0.5,
                                     shrinkA=0, shrinkB=2))
for sp in ("top", "right"):
    axb.spines[sp].set_visible(False)

# 2c 跨算法一致性（保守端）
axc = fig.add_subplot(gs[2])
pairs = [("CD8 T", 0.825, 0.174), ("M2 / mono", 0.396, 0.313), ("NK", 0.347, 0.312),
         ("Tregs / T", 0.306, 0.240), ("B cells", 0.234, 0.217),
         ("Neutrophils", 0.011, 0.106), ("DC", -0.154, 0.079)]
yy = np.arange(len(pairs))[::-1]
for y, (nm, a, b) in zip(yy, pairs):
    axc.plot([a, b], [y, y], color="#C9CDD4", linewidth=1.4, zorder=1, solid_capstyle="round")
    axc.scatter([a], [y], s=16, color=C_BLUE, zorder=3, linewidth=0)
    axc.scatter([b], [y], s=16, color=C_ORANGE, zorder=3, linewidth=0)
axc.axvline(0, color=C_GREY, linewidth=0.6)
axc.set_yticks(yy); axc.set_yticklabels([p_[0] for p_ in pairs], fontsize=5.7)
axc.set_xlim(-0.25, 0.95); axc.set_ylim(-0.9, len(pairs) - 0.35)
axc.set_xlabel("quanTIseq vs MCP-counter  ρ", fontsize=6.4)
axc.set_title("Concordance vs normalisation", fontsize=7.0, fontweight="bold", pad=5)
axc.legend([Line2D([], [], marker="o", linestyle="", color=C_BLUE, markersize=4),
            Line2D([], [], marker="o", linestyle="", color=C_ORANGE, markersize=4)],
           ["MCP-counter, no QN", "MCP-counter, QN"], fontsize=5.5, frameon=False,
           loc="lower right", handletextpad=0.2, borderpad=0.1)
axc.annotate("0.825 → 0.174", (0.50, yy[0]), xytext=(0, 12), textcoords="offset points",
             fontsize=5.4, color=C_RED, ha="center",
             arrowprops=dict(arrowstyle="-", color=C_RED, linewidth=0.6))
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)

fig.savefig(os.path.join(OUT, "fig2_almp_construction.pdf"))
fig.savefig(os.path.join(OUT, "fig2_almp_construction.png"), dpi=300)
plt.close(fig)
print("saved fig2")
