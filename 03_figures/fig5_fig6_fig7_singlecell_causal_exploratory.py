# -*- coding: utf-8 -*-
"""图 4a（跨方法一致性热图）／图 5（单细胞）／图 7b-d（因果审计）／图 8（探索层）"""
import os, sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

sys.stdout.reconfigure(encoding="utf-8")
ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
OUT = os.path.join(ROOT, "output", "figures")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.4, "axes.linewidth": 0.7,
                     "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300})
C_BLUE = "#2B5C8A"; C_TEAL = "#2A7F7F"; C_ORANGE = "#C1602A"
C_GREY = "#8A8F98"; C_RED = "#A8322A"; C_GREEN = "#2E6B3E"; C_PURPLE = "#5B3F7A"

TISS = {"皮肤日晒": "sun-exposed skin", "皮肤非日晒": "suprapubic skin",
        "皮下脂肪": "subcutaneous adipose", "全血": "whole blood"}

def tsv(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig"), delimiter="\t"))

def csvr(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))

# ============================================================ 图 4a
d = tsv(os.path.join(ROOT, "output/P2-T5b_多去卷积一致性v1/data/P2T5b_crossmethod_corr.tsv"))
labels = [r["feature"] for r in d]
M = np.array([[float(r[c]) for c in labels] for r in d])
methods = [l.split(".")[0] for l in labels]
blocks, cur, start = [], methods[0], 0
for i in range(1, len(methods) + 1):
    if i == len(methods) or methods[i] != cur:
        blocks.append((cur, start, i)); start = i
        if i < len(methods):
            cur = methods[i]

fig = plt.figure(figsize=(118 / 25.4, 112 / 25.4))
ax = fig.add_axes([0.315, 0.325, 0.485, 0.505])
im = ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=1.0),
               interpolation="nearest")
ax.set_xticks(range(len(labels)))
ax.set_yticks(range(len(labels)))
ax.set_xticklabels([l.split(".", 1)[1].replace(".", " ") for l in labels], rotation=90, fontsize=4.2)
ax.set_yticklabels([l.split(".", 1)[1].replace(".", " ") for l in labels], fontsize=4.2)
for (nm, a, b) in blocks:
    ax.add_patch(Rectangle((a - 0.5, a - 0.5), b - a, b - a, fill=False,
                           edgecolor="#111111", linewidth=0.9))
ax.set_xticks(np.arange(len(labels) + 1) - 0.5, minor=True)
ax.set_yticks(np.arange(len(labels) + 1) - 0.5, minor=True)
ax.grid(which="minor", color="white", linewidth=0.30)
ax.tick_params(which="minor", length=0)
ax.tick_params(axis="both", length=1.6, width=0.5)
for b in ax.spines.values():
    b.set_linewidth(0.7)
cax = fig.add_axes([0.815, 0.325, 0.016, 0.505])
cb = fig.colorbar(im, cax=cax)
cb.set_label("Spearman ρ", fontsize=5.4, labelpad=3)
cb.ax.tick_params(labelsize=4.8, width=0.5, length=1.5)
fig.text(0.012, 0.975, "a", fontsize=9, fontweight="bold", va="top")
fig.text(0.055, 0.968,
         "Cross-framework deconvolution concordance\n(quanTIseq · MCP-counter · MCP-counter QN)",
         fontsize=6.6, fontweight="bold", va="top", linespacing=1.5)
fig.savefig(os.path.join(OUT, "fig4a_crossmethod.pdf"))
fig.savefig(os.path.join(OUT, "fig4a_crossmethod.png"), dpi=300)
plt.close(fig)
print("fig4a")


# ============================================================ 图 5
cs = tsv(os.path.join(ROOT, "output/P4-T1_scRNA_ALMP定位v1/results/celltype_summary.tsv"))
pp = tsv(os.path.join(ROOT, "output/P4-T1_scRNA_ALMP定位v1/results/patient_paired.tsv"))
seg1 = [r for r in cs if r["geneset"] == "SEG1_ADIPOKINE_CORE"]
celltypes = []
for r in seg1:
    if r["celltype"] not in celltypes:
        celltypes.append(r["celltype"])
counts = {r["celltype"]: int(r["n"]) for r in seg1}

fig = plt.figure(figsize=(180 / 25.4, 62 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[0.85, 1.35, 1.35], wspace=0.52,
                      left=0.105, right=0.985, top=0.845, bottom=0.185)

axa = fig.add_subplot(gs[0])
order = sorted(celltypes, key=lambda c: counts[c])
axa.barh(range(len(order)), [counts[c] for c in order], color=C_BLUE, height=0.62,
         edgecolor="white", linewidth=0.5)
axa.set_yticks(range(len(order))); axa.set_yticklabels(order, fontsize=5.6)
axa.set_xscale("log"); axa.set_xlabel("Cells (log scale)", fontsize=6.2)
axa.set_title("4,645 cells / 19 patients", fontsize=6.8, fontweight="bold", pad=5)
for i, c in enumerate(order):
    axa.text(counts[c] * 1.15, i, str(counts[c]), fontsize=5.0, va="center")
axa.set_xlim(10, 4000)
for s in ("top", "right"):
    axa.spines[s].set_visible(False)

axb = fig.add_subplot(gs[1])
x = np.arange(len(celltypes))
for meth, col, off in (("UCell", C_ORANGE, -0.17), ("AUCell", C_TEAL, 0.17)):
    vals = []
    for c in celltypes:
        hit = [r for r in seg1 if r["celltype"] == c and r["method"] == meth]
        vals.append(float(hit[0]["median"]) if hit else np.nan)
    axb.bar(x + off, vals, width=0.32, color=col, label=meth, edgecolor="white", linewidth=0.4)
    for xi, v in zip(x + off, vals):
        if not np.isnan(v) and v > 0.0015:
            axb.text(xi, v + 0.002, "%.3f" % v, fontsize=3.9, rotation=90,
                     ha="center", va="bottom", color="#333333")
axb.set_xticks(x); axb.set_xticklabels(celltypes, fontsize=5.4, rotation=30, ha="right")
axb.set_ylabel("Median score (adipokine core)", fontsize=6.2)
axb.set_title("Absent from lymphoid compartments", fontsize=6.8, fontweight="bold", pad=5)
axb.legend(fontsize=5.2, frameon=False, loc="upper left", handlelength=1.0, borderpad=0.1)
axb.set_ylim(0, max(max(float(r["median"]) for r in seg1) * 1.42, 0.02))
for s in ("top", "right"):
    axb.spines[s].set_visible(False)

axc = fig.add_subplot(gs[2])
p1 = [r for r in pp if r["geneset"] == "SEG1_ADIPOKINE_CORE" and r["method"] == "UCell"]
p1 = sorted(p1, key=lambda r: float(r["median_delta_malignant_minus_other"]))
yv = np.arange(len(p1))
for y, r in zip(yv, p1):
    dv = float(r["median_delta_malignant_minus_other"]); q = float(r["q"])
    col = C_RED if q < 0.05 else C_GREY
    axc.plot([0, dv], [y, y], color="#CFD3D9", linewidth=1.0, zorder=1)
    axc.scatter([dv], [y], s=18, color=col, zorder=3, linewidth=0)
    axc.text(dv + (0.004 if dv >= 0 else -0.004), y, "q=%.3f" % q, fontsize=4.4,
             va="center", ha="left" if dv >= 0 else "right",
             color=col if q < 0.05 else "#777777")
axc.axvline(0, color="#444444", linewidth=0.7)
axc.set_yticks(yv)
axc.set_yticklabels(["%s (n=%s)" % (r["vs"], r["n_patients"]) for r in p1], fontsize=5.4)
axc.set_xlabel("Δ malignant − other (within patient)", fontsize=6.2)
axc.set_title("Paired comparison, 19 patients", fontsize=6.8, fontweight="bold", pad=5)
axc.set_xlim(min(float(r["median_delta_malignant_minus_other"]) for r in p1) - 0.045, 0.115)
for s in ("top", "right"):
    axc.spines[s].set_visible(False)
axc.text(0.99, 0.03, "red = q < 0.05", fontsize=5.0, color=C_RED, ha="right",
         transform=axc.transAxes)

fig.savefig(os.path.join(OUT, "fig5_singlecell.pdf"))
fig.savefig(os.path.join(OUT, "fig5_singlecell.png"), dpi=300)
plt.close(fig)
print("fig5")

# ============================================================ 图 7b-d
smr = tsv(os.path.join(ROOT, "output/P2-T13_GTEx定量SMR与共定位v1/results/smr_heidi_coloc_matrix.tsv"))
genes = sorted({r["gene"] for r in smr})
tissues = []
for r in smr:
    k = TISS.get(r["tissue_cn"], r["tissue_cn"])
    if k not in tissues:
        tissues.append(k)

def cat(r):
    if r["gate1"].strip() == "未评":
        return "not evaluable"
    if r["gate2"].strip().startswith("过"):
        return "SMR gate not met, HEIDI compatible"
    return "HEIDI rejected or too few loci"

CAT = {"both gates passed": C_GREEN,
       "SMR gate not met, HEIDI compatible": C_ORANGE,
       "HEIDI rejected or too few loci": "#7A4A2E",
       "not evaluable": "#DDE0E5"}

fig = plt.figure(figsize=(180 / 25.4, 82 / 25.4))
gs = fig.add_gridspec(1, 3, width_ratios=[1.18, 0.95, 0.90], wspace=0.46,
                      left=0.085, right=0.982, top=0.825, bottom=0.245)

axb = fig.add_subplot(gs[0])
for i, g in enumerate(genes):
    for j, t in enumerate(tissues):
        hit = [r for r in smr if r["gene"] == g and TISS.get(r["tissue_cn"], r["tissue_cn"]) == t]
        if not hit:
            continue
        k = cat(hit[0])
        axb.add_patch(Rectangle((j - 0.45, i - 0.45), 0.90, 0.90, facecolor=CAT[k],
                               edgecolor="white", linewidth=0.9))
        pv = hit[0]["smr_p"]
        if pv and pv != "NA":
            axb.text(j, i, "%.2f" % float(pv), ha="center", va="center", fontsize=4.5,
                     color="white" if k != "not evaluable" else "#8A8F98")
axb.set_xlim(-0.5, len(tissues) - 0.5); axb.set_ylim(len(genes) - 0.5, -0.5)
axb.set_xticks(range(len(tissues)))
axb.set_xticklabels([t.replace("subcutaneous", "subcut.").replace("sun-exposed", "sun-exp.")
                     .replace(" ", "\n", 1) for t in tissues], fontsize=5.0)
axb.set_yticks(range(len(genes))); axb.set_yticklabels(genes, fontsize=5.8)
axb.set_title("32 tissue × locus combinations (cell = SMR p)", fontsize=6.6,
              fontweight="bold", pad=5)
for sp in axb.spines.values():
    sp.set_visible(False)
axb.tick_params(length=0)
axb.legend([Rectangle((0, 0), 1, 1, fc=CAT[k], ec="white") for k in
            ("SMR gate not met, HEIDI compatible", "HEIDI rejected or too few loci", "not evaluable")],
           ["evaluable, HEIDI compatible (12)", "HEIDI rejected / < 3 loci (3)", "not evaluable (17)"],
           fontsize=5.0, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.235),
           ncol=1, handlelength=0.9, handleheight=0.9, borderpad=0.1, labelspacing=0.35)

ci = tsv(os.path.join(ROOT, "output/P2-T13_GTEx定量SMR与共定位v1/results/coloc_informativeness.tsv"))
axc = fig.add_subplot(gs[1])
xs, ys = [], []
for r in ci:
    try:
        p_ = float(r["smr_p"]); h = float(r["pp_h4"])
    except (TypeError, ValueError):
        continue
    xs.append(-np.log10(p_)); ys.append(h)
axc.scatter(xs, ys, s=15, color=C_BLUE, alpha=0.80, linewidth=0)
axc.axhline(0.8, color=C_RED, linewidth=0.7, linestyle=(0, (3, 2)))
axc.text(0.06, 0.825, "PP.H4 = 0.8", fontsize=5.0, color=C_RED)
cds2 = [(x_, y_, r_) for x_, y_, r_ in zip(xs, ys, ci) if r_["gene"] == "CDS2"]
if cds2:
    x_, y_, r_ = max(cds2, key=lambda t: t[1])
    axc.annotate("CDS2", (x_, y_), xytext=(8, -11), textcoords="offset points",
                 fontsize=5.4, color=C_ORANGE, fontweight="bold",
                 arrowprops=dict(arrowstyle="-", color=C_ORANGE, linewidth=0.6))
axc.set_xlabel("−log₁₀ SMR p", fontsize=6.2)
axc.set_ylabel("PP.H4 (colocalisation)", fontsize=6.2)
axc.set_title("Colocalisation adds no discrimination", fontsize=6.6, fontweight="bold", pad=5)
axc.set_ylim(-0.04, 1.04); axc.set_xlim(-0.05, 1.65)
axc.text(0.03, 0.28, "ρ = 0.122\n(p = 0.665)", fontsize=5.2, ha="left",
         transform=axc.transAxes, color="#333333", linespacing=1.5)
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)

axd = fig.add_subplot(gs[2]); axd.axis("off")
axd.set_xlim(0, 100); axd.set_ylim(0, 100)
axd.text(0, 100, "d", fontsize=9, fontweight="bold", va="top")
axd.text(0, 94, "Mediation via leukocyte traits", fontsize=6.6, fontweight="bold", va="top")
def nd(x, y, w, h, txt, ec):
    axd.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5,rounding_size=2.0",
                                 facecolor="white", edgecolor=ec, linewidth=0.9))
    axd.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=5.0, color="#222222")
nd(2, 58, 34, 16, "Adipokine", C_BLUE)
nd(64, 58, 34, 16, "Melanoma", C_RED)
nd(33, 20, 34, 16, "Leukocyte\ncounts", C_TEAL)
axd.add_patch(FancyArrowPatch((36, 66), (64, 66), arrowstyle="-|>", mutation_scale=8,
                              linewidth=0.9, color=C_GREY))
axd.add_patch(FancyArrowPatch((26, 57), (42, 37), arrowstyle="-|>", mutation_scale=8,
                              linewidth=0.9, color=C_GREY))
axd.add_patch(FancyArrowPatch((58, 37), (74, 57), arrowstyle="-|>", mutation_scale=8,
                              linewidth=0.9, color=C_GREY))
axd.text(50, 71.5, "OR 0.67\n(0.42–1.07)", fontsize=4.7, ha="center", va="bottom",
         color="#333333", linespacing=1.4)
axd.text(17, 45, "a: null", fontsize=4.7, color="#777777", rotation=-48)
axd.text(83, 45, "b: null", fontsize=4.7, color="#777777", rotation=48)
axd.text(50, 8, "Indirect effect 95% CI includes 0\nBarton 2021, UK Biobank only\n(no overlap with MVP)",
         fontsize=4.7, ha="center", va="top", color="#444444", linespacing=1.55)

fig.savefig(os.path.join(OUT, "fig7bcd_causal_audit.pdf"))
fig.savefig(os.path.join(OUT, "fig7bcd_causal_audit.png"), dpi=300)
plt.close(fig)
print("fig7bcd")

# ============================================================ 图 8
surv = tsv(os.path.join(ROOT, "output/P2-T11_外部验证v1/data/P2T11_axis_survival_assoc.tsv"))
resp = csvr(os.path.join(ROOT, "output/P2-T6b-R_五队列基因宇宙复查v1/data/P2T6bR_response_FULL.csv"))

fig = plt.figure(figsize=(180 / 25.4, 92 / 25.4))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.58, left=0.165, right=0.965,
                      top=0.895, bottom=0.10)

axa = fig.add_subplot(gs[0, 0])
rows = [r for r in surv if r["feature"] == "SEG1_ADIPOKINE_CORE"]
labs = []
for i, r in enumerate(rows):
    hr = float(r["HR"]); lo = float(r["ci_lo"]); hi = float(r["ci_hi"]); p_ = float(r["p"])
    axa.plot([lo, hi], [i, i], color=C_GREY, linewidth=1.0)
    axa.scatter([hr], [i], s=16, color=C_BLUE, zorder=3, linewidth=0)
    axa.text(1.06, i, "%.2f p=%.3f" % (hr, p_), fontsize=4.2, va="center",
             transform=axa.get_yaxis_transform())
    short = "split" if "median" in r["model"] else "cont."
    labs.append("%s %s %s" % (r["readout"], r["endpoint"], short))
axa.axvline(1, color="#444444", linewidth=0.7, linestyle=(0, (3, 2)))
axa.set_yticks(range(len(rows))); axa.set_yticklabels(labs, fontsize=4.9)
axa.set_xscale("log"); axa.set_xlim(0.30, 1.95)
axa.minorticks_off()
from matplotlib.ticker import NullFormatter, NullLocator
axa.xaxis.set_minor_locator(NullLocator())
axa.xaxis.set_minor_formatter(NullFormatter())
axa.set_xticks([0.3, 0.5, 0.7, 1.0, 1.5, 2.0])
axa.set_xticklabels(["0.3", "0.5", "0.7", "1.0", "1.5", "2.0"], fontsize=5.2)
axa.set_xlabel("HR (95% CI) per axis score", fontsize=5.9)
axa.set_title("Survival: direction not reproducible", fontsize=6.6, fontweight="bold", pad=5)
axa.set_ylim(-0.7, len(rows) - 0.3)
for sp in ("top", "right"):
    axa.spines[sp].set_visible(False)

axb = fig.add_subplot(gs[0, 1])
qs = np.array([float(r["q"]) for r in resp])
axb.scatter(np.random.RandomState(1).uniform(-0.30, 0.30, len(qs)), qs, s=8,
            color=C_PURPLE, alpha=0.75, linewidth=0)
axb.axhline(0.05, color=C_RED, linewidth=0.8, linestyle=(0, (3, 2)))
axb.text(0.30, 0.075, "q = 0.05", fontsize=5.0, color=C_RED, ha="right")
axb.set_xlim(-0.5, 0.5); axb.set_ylim(0, 1.05); axb.set_xticks([])
axb.set_ylabel("BH q", fontsize=5.9)
axb.set_title("ICB response: %d tests, none q < 0.05" % len(qs), fontsize=6.6,
              fontweight="bold", pad=5)
axb.text(0.03, 0.22, "minimum q = %.2f" % qs.min(), fontsize=5.2, va="top",
         transform=axb.transAxes)
for sp in ("top", "right"):
    axb.spines[sp].set_visible(False)

axc = fig.add_subplot(gs[1, 0])
fa = [r for r in smr if r["gene"] == "FADS2"]
order_t = [TISS.get(r["tissue_cn"], r["tissue_cn"]) for r in fa]
SHORT = {"whole blood": "whole blood", "subcutaneous adipose": "subcut. adipose",
         "suprapubic skin": "suprapubic skin", "sun-exposed skin": "sun-exp. skin"}
for i, r in enumerate(fa):
    orv = float(r["smr_or"]); se = float(r["smr_se"]); p_ = float(r["smr_p"])
    lo, hi = np.exp(np.log(orv) - 1.96 * se), np.exp(np.log(orv) + 1.96 * se)
    axc.plot([lo, hi], [i, i], color=C_ORANGE, linewidth=1.1)
    axc.scatter([orv], [i], s=16, color=C_ORANGE, zorder=3, linewidth=0)
    axc.text(1.335, i, "p=%.3f" % p_, fontsize=4.6, va="center")
axc.axvline(1, color="#444444", linewidth=0.7, linestyle=(0, (3, 2)))
axc.set_yticks(range(len(fa)))
axc.set_yticklabels([SHORT.get(t, t) for t in order_t], fontsize=5.0)
axc.set_xlim(0.90, 1.44)
axc.set_xticks([0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
axc.set_xlabel("OR per predicted expression unit", fontsize=5.9)
axc.set_title("FADS2: same direction, none significant", fontsize=6.6, fontweight="bold", pad=5)

for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)

axd = fig.add_subplot(gs[1, 1])
cd = tsv(os.path.join(ROOT, "output/P2-T11_外部验证v1/data/P2T11_CDS2_survival.tsv"))
for i, r in enumerate(cd):
    hr = float(r["HR"]); lo = float(r["ci_lo"]); hi = float(r["ci_hi"])
    axd.plot([lo, hi], [i, i], color=C_GREY, linewidth=1.0)
    axd.scatter([hr], [i], s=15, color=C_TEAL, zorder=3, linewidth=0)
axd.axvline(1, color="#444444", linewidth=0.7, linestyle=(0, (3, 2)))
axd.set_yticks(range(len(cd)))
axd.set_yticklabels(["%s %s" % (r["endpoint"], "split" if "median" in r["model"] else "continuous")
                     for r in cd], fontsize=5.0)
axd.set_xscale("log"); axd.set_xlim(0.55, 2.2)
axd.set_xticks([0.6, 0.8, 1.0, 1.5, 2.0])
axd.set_xticklabels(["0.6", "0.8", "1.0", "1.5", "2.0"], fontsize=5.2)
axd.set_xlabel("HR (95% CI)", fontsize=5.9)
axd.set_title("CDS2: demoted to an ordinary member", fontsize=6.6, fontweight="bold", pad=14)
axd.text(0.5, 1.03, "SMR p = 0.033 · HEIDI 0.998 · PP.H4 0.816 · but q = 0.110",
         fontsize=4.6, color="#444444", transform=axd.transAxes, ha="center", va="bottom")
for sp in ("top", "right"):
    axd.spines[sp].set_visible(False)

fig.savefig(os.path.join(OUT, "fig8_exploratory.pdf"))
fig.savefig(os.path.join(OUT, "fig8_exploratory.png"), dpi=300)
plt.close(fig)
print("fig8")
