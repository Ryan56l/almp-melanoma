# -*- coding: utf-8 -*-
"""Figure S6：GSE115978 第二队列的轴格局与患者内配对差（干净版）"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "font.size": 6.4,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T3_GSE115978复现v1", "data")
FIG = os.path.join(ROOT, "figures")
C = dict(blue="#2166a5", red="#b8232f", orange="#e08214", teal="#2a9d8f",
         grey="#7f7f7f", purple="#6b4c9a", light="#d9d9d9", ink="#40566b")

sc = pd.read_csv(os.path.join(D, "P5T3_cell_scores.tsv"), sep="\t", index_col=0)
ORD = ["Mal", "Macrophage", "CAF", "Endo.", "T.CD8", "T.CD4", "T.cell", "B.cell", "NK"]
sc = sc[sc.cell_type.isin(ORD)]
CTS = ["T.CD8", "T.CD4", "B.cell", "NK", "Macrophage", "CAF", "Endo."]

IN = 1 / 25.4
W, H = 180, 150
fig = plt.figure(figsize=(W * IN, H * IN))


def txt(x, y, s, **kw):
    fig.text(x / W, 1 - y / H, s, **kw)


def block(letter, title, x, y, subs=()):
    txt(x, y, letter, fontsize=8.4, fontweight="bold", color="black")
    if title:
        txt(x + 5.0, y, title, fontsize=6.8, fontweight="bold", color="black")
    for k, s in enumerate(subs):
        txt(x + 5.0, y + 3.4 + 2.6 * k, s, fontsize=4.8, color=C["ink"], style="italic")


def A(x, y, w, h):
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def thin(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)


# ---------- a
axa = A(24, 32, 60, 42)
xx = np.arange(len(ORD))
SER = [("core", "core (22)", C["red"], "o"), ("intrinsic", "intrinsic (273)", C["blue"], "s"),
       ("intrinsic166", "intrinsic, 166", C["teal"], "^")]
axa.axvspan(0.6, 2.4, color=C["light"], alpha=.6, lw=0, zorder=0)
for k, (col, lab, cc, mk) in enumerate(SER):
    v = [float(sc.loc[sc.cell_type == g, col].median()) for g in ORD]
    axa.plot(xx + (k - 1) * 0.18, v, mk, ms=3.0, color=cc, label=lab, ls="none")
axa.set_xticks(xx)
axa.set_xticklabels(ORD, rotation=40, ha="right", fontsize=4.6, fontstyle="italic")
axa.set_ylabel("median UCell-style score", fontsize=5.1, labelpad=2)
axa.tick_params(axis="y", labelsize=4.6, length=2, pad=1)
axa.tick_params(axis="x", length=2, pad=1)
axa.legend(frameon=False, fontsize=4.4, loc="upper center", ncol=3, handletextpad=0.3,
           borderpad=0.1, columnspacing=0.8, bbox_to_anchor=(0.5, 1.16))
thin(axa)
axa.text(1.5, 0.4735, "non-lymphoid", fontsize=4.4, color=C["grey"], style="italic", ha="center", va="bottom")
block("a", "Cell-type medians in GSE115978 (7,186 cells, 32 patients)", 3, 12,
      ["each score is a per-cell mean of full-transcriptome percentile rank"])

# ---------- 配对差
rows = []
for coh in ["ALL", "New"]:
    for axis in ["core", "intrinsic"]:
        sub = sc if coh == "ALL" else sc[sc.cohort == coh]
        pv = sub.pivot_table(index="patient", columns="cell_type", values=axis, aggfunc="median")
        for ct in CTS:
            if "Mal" not in pv.columns or ct not in pv.columns:
                continue
            d = (pv["Mal"] - pv[ct]).dropna()
            if len(d) < 5:
                continue
            rows.append(dict(cohort=coh, axis=axis, ct=ct, n=len(d), med=float(d.median()),
                             p=float(stats.wilcoxon(d).pvalue)))
DD = pd.DataFrame(rows)
DD.to_csv(os.path.join(D, "P5T3_paired_deltas.tsv"), sep="\t", index=False)

pairs = [("b", "ALL", "core", C["red"], "core score: malignant \u2212 comparator"),
         ("c", "ALL", "intrinsic", C["blue"], "intrinsic score: malignant \u2212 comparator")]
for k, (letter, coh, axis, col, title) in enumerate(pairs):
    ax = A(112, 28 + k * 56, 58, 34)
    d = DD[(DD.cohort == coh) & (DD.axis == axis)].set_index("ct").reindex(
        [c for c in CTS if c in set(DD.ct)])
    yy = np.arange(len(d))
    ax.barh(yy, d.med, color=[col if m > 0 else C["grey"] for m in d.med], height=0.64, lw=0)
    lim = max(abs(d.med.min()), abs(d.med.max())) * 1.75
    for i, (m, p, n) in enumerate(zip(d.med, d.p, d.n)):
        s = "q<0.001" if p < 1e-3 else f"p={p:.3f}"
        ax.text(m + (0.001 if m > 0 else -0.001), i, f"{s} (n={int(n)})", fontsize=3.6,
                va="center", ha="left" if m > 0 else "right")
    ax.axvline(0, color="black", lw=0.5)
    ax.set_yticks(yy)
    ax.set_yticklabels(d.index, fontsize=4.4, fontstyle="italic")
    ax.tick_params(axis="x", labelsize=4.2, length=2, pad=1)
    ax.tick_params(axis="y", length=0, pad=1)
    ax.set_xlim(-lim, lim)
    ax.set_xlabel("paired difference in score", fontsize=4.8, labelpad=1)
    thin(ax)
    block(letter, title, 108, 20 + k * 56,
          ["within-patient medians; Wilcoxon signed-rank"])

# ---------- d 两队列方向不一致
txt(108, 125, "d", fontsize=8.4, fontweight="bold", color="black")
txt(113, 125, "Two cohorts, two directions",
    fontsize=6.4, fontweight="bold", color="black")
nm = DD[(DD.cohort == "New") & (DD.axis == "intrinsic")].set_index("ct")
am = DD[(DD.cohort == "ALL") & (DD.axis == "intrinsic")].set_index("ct")
m1 = nm.loc["Macrophage"]
txt(111, 129.6, f"malignant \u2212 macrophage: {m1.med:+.4f} (2018, n={int(m1.n)}, p={m1.p:.3f}) "
    f"vs {am.loc['Macrophage','med']:+.4f}", fontsize=4.8, color=C["ink"], style="italic")
txt(111, 132.6, f"(all {int(am.loc['Macrophage','n'])} patients, "
    f"p={am.loc['Macrophage','p']:.2f}): the sign is not stable between cohorts",
    fontsize=4.8, color=C["ink"], style="italic")
txt(111, 135.6, "attribution rests on the lymphoid contrast and the cell-line data,",
    fontsize=4.8, color=C["ink"], style="italic")
txt(111, 138.6, "not on a comparison with macrophages", fontsize=4.8, color=C["ink"], style="italic")

for ext in ["pdf", "png"]:
    fig.savefig(os.path.join(FIG, f"figS6_two_cohort_scRNA.{ext}"), dpi=300)
print("saved")
print(DD.round(5).to_string(index=False))
