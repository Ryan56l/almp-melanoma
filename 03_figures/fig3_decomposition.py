# -*- coding: utf-8 -*-
"""
y52 v2: 图 3 出图（指令五十一）—— 版式修正版
a 五队列 SEG1 x Fibro/Endo/Neutro 森林（全转录组口径）+ CD8 对照
b 主轴线合并 vs ICB 应答（4 队列 + 合并菱形，d=-0.095, I2=7.2%）—— 取代原 40 格全景
c ALMP x TMB 正交性（36 条）
d 各集合自有零分布对照（段内横向比较，规避"规模不同不可比百分位"陷阱）
  + 两套零模型并列（均为 GSVA，与观测 0.410 同算法）

输出：output/figures/fig3_niche_coupling.pdf / .png (300 dpi)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
D6R = os.path.join(ROOT, "output/P2-T6b-R_五队列基因宇宙复查v1/data")
D12 = os.path.join(ROOT, "output/P2-T12_随机基因集阴性对照v1/results")
FIG = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7.0,
    "axes.linewidth": 0.6, "axes.titlesize": 7.8, "axes.labelsize": 7.2,
    "xtick.labelsize": 6.6, "ytick.labelsize": 6.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.6, "ytick.major.size": 2.6,
    "legend.fontsize": 6.2, "legend.frameon": False,
    "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42,
})

C = dict(blue="#2166a5", red="#b8232f", orange="#e08214", teal="#2a9d8f",
         purple="#6b4c9a", grey="#7f7f7f", light="#d9d9d9")
MM = 1 / 25.4
W, H = 180 * MM, 232 * MM


def fisher_ci(r, n, z=1.959964):
    r = np.clip(r, -0.999, 0.999)
    zr = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return np.tanh(zr - z * se), np.tanh(zr + z * se)


coup = pd.read_csv(os.path.join(D6R, "P2T6bR_coupling_FULL.csv"))
cmp_ = pd.read_csv(os.path.join(D6R, "P2T6bR_coupling_comparison.csv"))
tmb = pd.read_csv(os.path.join(D6R, "P2T6bR_tmb_FULL.csv"))
seg = pd.read_csv(os.path.join(D6R, "P2T6bR_segment_specificity_FULL.csv"))
tcga = pd.read_csv(os.path.join(D6R, "P2T6bR_tcga_full_coupling.csv"))
null = pd.read_csv(os.path.join(D12, "y46_random_control_main.tsv"), sep="\t")
nulls = pd.read_csv(os.path.join(D12, "y46_random_control_sensitivity.tsv"), sep="\t")

COH = ["TCGA-SKCM", "Gide2019", "Liu2019", "Riaz2017", "UCLA2016"]
NN = {"TCGA-SKCM": 472, "Gide2019": 91, "Liu2019": 122, "Riaz2017": 107,
      "UCLA2016": 27, "VanAllen2015": 40}
POP = ["Fibroblasts", "Endothelial cells", "Neutrophils", "CD8 T cells"]
PCOL = {"Fibroblasts": C["blue"], "Endothelial cells": C["teal"],
        "Neutrophils": C["orange"], "CD8 T cells": C["red"]}
PLAB = {"Fibroblasts": "Fibroblasts", "Endothelial cells": "Endothelial cells",
        "Neutrophils": "Neutrophils", "CD8 T cells": "CD8$^+$ T cells"}

recs = []
for pop in POP:
    for c in COH[1:]:
        r = coup[(coup.cohort == c) & (coup.feature == pop)].iloc[0]
        recs.append(dict(cohort=c, pop=pop, n=int(r.n), rho=r.rho, q=r.q))
    t = tcga[(tcga.algorithm == "ssGSEA") & (tcga.population == pop)].iloc[0]
    recs.append(dict(cohort="TCGA-SKCM", pop=pop, n=int(t.n), rho=t.rho, q=t.q))
    v = seg[(seg.cohort == "VanAllen2015") & (seg.feature == pop) &
            (seg.set == "SEG1_ADIPOKINE_CORE")]
    if len(v):
        recs.append(dict(cohort="VanAllen2015", pop=pop, n=40,
                         rho=v.iloc[0].rho, q=v.iloc[0].q))
FA = pd.DataFrame(recs)

fig = plt.figure(figsize=(W, H))
gs = GridSpec(3, 2, figure=fig, height_ratios=[1.22, 0.92, 1.00],
              width_ratios=[0.88, 1.12], hspace=1.30, wspace=0.46,
              left=0.155, right=0.972, top=0.935, bottom=0.055)

# ================================================================ a
axa = fig.add_subplot(gs[0, :])
labs = ["TCGA-SKCM", "Gide2019", "Liu2019", "Riaz2017", "UCLA2016", "VanAllen2015"]
nlabs = [f"{l} (n={NN[l]})" if l != "VanAllen2015" else f"{l} \u2020 (n=40)"
         for l in labs]
for k, pop in enumerate(POP):
    off = (k - 1.5) * 0.20
    col = PCOL[pop]
    for i, c in enumerate(labs):
        row = FA[(FA.cohort == c) & (FA["pop"] == pop)]
        if not len(row):
            continue
        row = row.iloc[0]
        lo, hi = fisher_ci(row.rho, row.n)
        sig = (row.q == row.q) and (row.q < 0.05)
        yy = i + off
        axa.plot([lo, hi], [yy, yy], "-", color=col, lw=0.85,
                 solid_capstyle="butt", zorder=2)
        if c == "VanAllen2015":
            axa.plot(row.rho, yy, "o", mfc="white", mec=col, mew=0.9, ms=3.5, zorder=3)
        else:
            axa.plot(row.rho, yy, "o", mfc=col if sig else "white", mec=col,
                     mew=0.9, ms=3.6, zorder=3)
axa.axhspan(-0.62, 5.28, color="#f4f4f4", zorder=0)
axa.axhspan(5.28, 5.88, color="#fbf3e4", zorder=0)
axa.set_yticks(range(len(labs)))
axa.set_yticklabels(nlabs, fontsize=6.3)
axa.invert_yaxis()
axa.set_ylim(5.95, -0.70)
axa.axvline(0, color="black", lw=0.6, ls=(0, (3, 2)))
axa.set_xlim(-0.60, 1.00)
axa.set_xlabel("Spearman's $\\rho$  (ALMP adipokine-core score vs MCP-counter population)",
               labelpad=2)
axa.set_title("a   Adipokine-core coupling with stromal and immune populations, five cohorts",
              loc="left", fontweight="bold", pad=6)
axa.spines[["top", "right"]].set_visible(False)
axa.text(0.0, -0.250,
         "Fibroblast, endothelial and neutrophil coefficients are of similar magnitude\n"
         "in every cohort: this is the fingerprint of total specimen cellularity, and is\n"
         "therefore not evidence for any compartment-specific association. CD8$^+$ T cells\n"
         "are shown as a boundary case; their inverse association is confined to TCGA-SKCM\n"
         "(4/4 ICB cohorts non-significant). \u2020 VanAllen2015: post hoc sensitivity\n"
         "cohort without response annotation.",
         transform=axa.transAxes, ha="left", va="top", fontsize=5.2,
         color="#555555", linespacing=1.30)
handles = [Line2D([], [], marker="o", ls="-", color=c, mfc=c, mec=c, lw=0.85, ms=3.4,
                  label=PLAB[p]) for p, c in PCOL.items()]
handles += [Line2D([], [], marker="o", ls="none", mfc="black", mec="black", ms=3.6,
                   label="q $<$ 0.05"),
            Line2D([], [], marker="o", ls="none", mfc="white", mec="black", ms=3.6,
                   label="q $\\geq$ 0.05")]
axa.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.155),
           ncol=6, handletextpad=0.35, columnspacing=1.0, borderpad=0.1)

# ================================================================ b
DB = os.path.join(ROOT, "output/P2-T6b_ICB队列扩面v1/data")
bys = pd.read_csv(os.path.join(DB, "P2T6b_axis_by_cohort.tsv"), sep="\t")
pool = pd.read_csv(os.path.join(DB, "P2T6b_pooled_meta.tsv"), sep="\t")

axb = fig.add_subplot(gs[1, 0])
sub = bys[(bys.readout == "ssGSEA") & (bys.feature == "SEG1_ADIPOKINE_CORE")].copy()
order_b = [c for c in ["Gide2019", "Liu2019", "Riaz2017", "UCLA2016"]
           if c in set(sub.cohort)]
sub = sub.set_index("cohort").loc[order_b].reset_index()
n1, n2 = sub.n_R.astype(float), sub.n_NR.astype(float)
sub["se"] = np.sqrt((n1 + n2) / (n1 * n2) + sub.cohen_d ** 2 / (2 * (n1 + n2)))
sub["lo"] = sub.cohen_d - 1.959964 * sub.se
sub["hi"] = sub.cohen_d + 1.959964 * sub.se
pm = pool[(pool.readout == "ssGSEA") &
          (pool.feature == "SEG1_ADIPOKINE_CORE")].iloc[0]

for i, r in enumerate(sub.itertuples()):
    axb.plot([r.lo, r.hi], [i, i], "-", color=C["grey"], lw=0.9, zorder=2)
    axb.plot(r.cohen_d, i, "o", mfc="white", mec=C["blue"], mew=1.0, ms=3.8, zorder=3)
    axb.text(0.985, i, f"{r.cohen_d:+.2f}", transform=axb.get_yaxis_transform(),
             ha="right", va="center", fontsize=5.8, color="#444444")

yp = len(sub) + 0.35
axb.fill([pm.ci_lo, pm.d_pooled, pm.ci_hi, pm.d_pooled],
         [yp, yp - 0.30, yp, yp + 0.30], color=C["red"], alpha=0.85,
         lw=0.7, edgecolor=C["red"], zorder=3)
axb.plot([pm.ci_lo, pm.ci_hi], [yp, yp], "-", color=C["red"], lw=0.8, zorder=2)
axb.text(0.985, yp, f"{pm.d_pooled:+.3f}", transform=axb.get_yaxis_transform(),
         ha="right", va="center", fontsize=5.8, color=C["red"], fontweight="bold")

axb.axvline(0, color="black", lw=0.6, ls=(0, (3, 2)))
axb.set_yticks(list(range(len(sub))) + [yp])
axb.set_yticklabels([f"{r.cohort} ({int(r.n_R)}/{int(r.n_NR)})" for r in sub.itertuples()]
                    + [f"pooled (k = {int(pm.k)})"], fontsize=6.0)
axb.invert_yaxis()
axb.set_ylim(yp + 0.75, -0.75)
axb.set_xlim(-1.35, 1.35)
axb.set_xlabel("Cohen's $d$, responders vs non-responders\n"
               f"pooled $d$ = {pm.d_pooled:.3f} ({pm.ci_lo:.3f} to {pm.ci_hi:.3f}); "
               f"$I^2$ = {pm.I2:.1f}%; p = {pm.p:.2f}", fontsize=6.2, labelpad=1)
axb.set_title("b   Axis and response, pooled",
              loc="left", fontweight="bold", pad=5)
axb.spines[["top", "right"]].set_visible(False)
axb.text(0.0, -0.34, "All four cohort estimates non-significant, with inconsistent direction.\n"
         "Cohorts are single-arm: predictive cannot be separated from prognostic.",
         transform=axb.transAxes, fontsize=5.3, color="#555555", va="top")


# ================================================================ c
axc = fig.add_subplot(gs[1, 1])
SETS = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
        "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]
SHORT = {"SEG1_ADIPOKINE_CORE": "SEG1 core", "SEG2_LIPID_PATHWAY": "SEG2 lipid",
         "SEG3_EXECUTOR": "SEG3 executor", "ALMP_COMPOSITIONAL_AXIS": "Compositional",
         "ALMP_TUMOR_INTRINSIC_AXIS": "Intrinsic", "ALMP_ALL_308": "ALL_308"}
mk = {"GSVA": "o", "ssGSEA": "^"}
coc = {"Liu2019": C["blue"], "Riaz2017": C["teal"], "UCLA2016": C["orange"]}
for i, s in enumerate(SETS):
    for _, r in tmb[tmb.feature == s].iterrows():
        lag = {"GSVA": -0.14, "ssGSEA": 0.14}[r.readout]
        axc.plot(r.rho, i + lag, mk[r.readout], ms=3.3,
                 mfc=coc.get(r.cohort, C["grey"]), mec="white", mew=0.4,
                 alpha=0.95, zorder=3)
axc.axvline(0, color="black", lw=0.6, ls=(0, (3, 2)))
for xv in (-0.364, 0.364):
    axc.axvline(xv, color=C["red"], lw=0.7, ls=":")
axc.set_yticks(range(len(SETS)))
axc.set_yticklabels([SHORT[s] for s in SETS], fontsize=6.2)
axc.invert_yaxis()
axc.set_xlim(-0.45, 0.45)
axc.set_xlabel("$\\rho$ (ALMP score vs TMB)\n"
               "min q = 0.770 ; |$\\rho$|max = 0.364 (dotted) ; 0/36 significant",
               fontsize=6.4, labelpad=1)
axc.set_title("c   Orthogonality to TMB, 36 tests", loc="left", fontweight="bold", pad=16)
axc.spines[["top", "right"]].set_visible(False)
h = [Line2D([], [], marker=mk[a], ls="none", mfc="none", mec="black", ms=3.6, label=a)
     for a in ("GSVA", "ssGSEA")]
h += [Line2D([], [], marker="o", ls="none", mfc=c, mec="white", ms=3.6,
            label=k.replace("2019", "").replace("2017", "").replace("2016", ""))
      for k, c in coc.items()]
axc.legend(handles=h, loc="lower right", bbox_to_anchor=(1.02, 0.995), ncol=5,
           handletextpad=0.22, labelspacing=0.2, borderpad=0.1, fontsize=5.2,
           columnspacing=0.50)

# ================================================================ d
gsd = gs[2, :].subgridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.30)
axd1 = fig.add_subplot(gsd[0])
neu = null[null.target == "PanglaoDB::Neutrophils"].copy()
order = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
         "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]
neu["o"] = neu.set.map({s: i for i, s in enumerate(order)})
neu = neu.sort_values("o")
for i, (_, r) in enumerate(neu.iterrows()):
    axd1.plot([r.null_p95, r.null_p95], [i - 0.17, i + 0.17], color=C["grey"],
              lw=1.7, alpha=0.55, solid_capstyle="butt", zorder=2)
    axd1.plot(r.null_mean, i, "|", color=C["grey"], ms=5, mew=1.0, zorder=3)
    axd1.plot(r.rho_real, i, "o", mfc=C["red"] if r.pctile > 90 else "white",
              mec=C["red"], mew=1.1, ms=4.2, zorder=4)
    axd1.text(0.70, i, f"{r.pctile:4.1f}%", va="center", ha="right", fontsize=6.0,
              color=C["red"] if r.pctile > 90 else "#555555")
axd1.axvline(0, color="black", lw=0.6, ls=(0, (3, 2)))
axd1.set_yticks(range(len(neu)))
axd1.set_yticklabels([SHORT[s] for s in neu.set], fontsize=6.2)
axd1.invert_yaxis()
axd1.set_ylim(5.9, -0.9)
axd1.set_xlim(-0.32, 0.76)
axd1.set_xlabel("$\\rho$ vs PanglaoDB neutrophils (TCGA-SKCM)\n"
                "each set compared with its own size- and expression-matched null",
                fontsize=6.4, labelpad=1)
axd1.set_title("", pad=2)
axd1.set_title("d   Each segment vs its own 1000 matched random gene sets (TCGA-SKCM)",
               loc="left", fontweight="bold", pad=16, fontsize=7.6)
axd1.spines[["top", "right"]].set_visible(False)
h = [Line2D([], [], marker="o", ls="none", mfc="black", mec="black", ms=3.8,
            label="observed $\\rho$"),
     Line2D([], [], color=C["grey"], lw=1.7, label="null P95"),
     Line2D([], [], marker="|", ls="none", color=C["grey"], ms=5, label="null mean")]
axd1.legend(handles=h, loc="upper left", fontsize=5.8, ncol=1,
            labelspacing=0.22, borderpad=0.15, handletextpad=0.5,
            bbox_to_anchor=(-0.02, 1.14))
axd1.text(0.70, -0.62, "percentile of\nobserved $\\rho$", ha="right", va="center",
          fontsize=5.4, color="#777777")

axd2 = fig.add_subplot(gsd[1])
xx = np.linspace(-0.55, 0.95, 600)
a_ = null[(null.target == "PanglaoDB::Neutrophils") &
          (null.set == "SEG1_ADIPOKINE_CORE")].iloc[0]
b_ = nulls[(nulls.scenario == "B_mean_x_var_30_GSVA") &
           (nulls.target == "PanglaoDB::Neutrophils")].iloc[0]
for mu, sd, p95, col, lab in (
        (a_.null_mean, a_.null_sd, a_.null_p95, C["blue"], "A  mean-expression decile"),
        (b_.null_mean, b_.null_sd, b_.null_p95, C["teal"], "B  mean $\\times$ variance (30 strata)")):
    yy = np.exp(-0.5 * ((xx - mu) / sd) ** 2)
    axd2.plot(xx, yy, color=col, lw=1.1, label=lab)
    axd2.fill_between(xx, 0, yy, color=col, alpha=0.10)
    axd2.plot(p95, np.exp(-0.5 * ((p95 - mu) / sd) ** 2), "v", color=col, ms=3.2,
              mec="white", mew=0.3)
axd2.axvline(a_.rho_real, color=C["red"], lw=1.1, zorder=1)
axd2.annotate("observed\n$\\rho$ = 0.410", xy=(a_.rho_real, 1.02), xytext=(0.58, 1.24),
              color=C["red"], fontsize=5.8, va="top", ha="left",
              arrowprops=dict(arrowstyle="-", color=C["red"], lw=0.6))
axd2.text(-0.53, 0.06, f"A: {a_.pctile:.1f}th pct (p$_{{emp}}$ = {a_.p_emp:.3f})\n"
          f"B: {b_.pctile:.1f}th pct (p$_{{emp}}$ = {b_.p_emp:.3f})",
          fontsize=5.6, color="#444444", va="bottom")
axd2.set_xlim(-0.55, 0.95)
axd2.set_ylim(0, 1.34)
axd2.set_yticks([])
axd2.set_xlabel("$\\rho$", labelpad=1)
axd2.set_title("")
axd2.spines[["top", "right", "left"]].set_visible(False)
axd2.legend(loc="upper left", fontsize=5.6, labelspacing=0.22, borderpad=0.15,
            bbox_to_anchor=(-0.02, 1.14))
axd2.text(0.0, -0.40, "Variance matching is not a neutral null in bulk matrices "
          "(matched random sets reach 1.9$\\times$ the variance of the real set);\n"
          "we therefore report both and make no claim of gene-set specificity.",
          transform=axd2.transAxes, fontsize=5.3, color="#555555", va="top")

for p in (axa, axb, axc, axd1, axd2):
    for s in p.spines.values():
        s.set_linewidth(0.6)

fig.savefig(os.path.join(FIG, "fig3_axis_decomposition.pdf"))
fig.savefig(os.path.join(FIG, "fig3_axis_decomposition.png"), dpi=300)
print("[ok] fig3 v5 ->", os.path.join(FIG, "fig3_axis_decomposition.pdf"))
print("null A pctile", a_.pctile, "| null B pctile", b_.pctile)
