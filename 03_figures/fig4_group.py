# -*- coding: utf-8 -*-
"""
y68 v2: 图 4 组版（a–e 五面板合一）—— 指令五十七 A 案重构后定稿
版式修正（相对 v1）：
  1. 面板抬头改为「以坐标区顶边为基准、va=bottom 向上堆叠」，副标题不再压进数据区；
  2. b 面板两个热图共享底部 x 轴标签（上图不画），方法名贴左侧竖排，消除 ssGSEA 标签与 GSVA 标题相撞；
  3. c / d 的图例移到坐标区上方留白带，不再覆盖数据点；
  4. 全图 180 × 224 mm，纵向预算逐带核算。
输出：manuscript/output/figures/fig4_group.pdf / .png
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
R50 = r"C:/Rtmp/p2t12/results50"
SRC = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output/P2-T20_纯度去混杂v1"
FIG = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output/figures"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 6.4,
    "axes.linewidth": 0.6, "axes.titlesize": 7.2, "axes.labelsize": 6.2,
    "xtick.labelsize": 5.6, "ytick.labelsize": 5.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.2, "ytick.major.size": 2.2,
    "legend.fontsize": 5.4, "legend.frameon": False,
    "savefig.dpi": 300, "figure.dpi": 150, "pdf.fonttype": 42, "ps.fonttype": 42,
})
C = dict(blue="#2166a5", red="#b8232f", orange="#e08214", teal="#2a9d8f",
         grey="#7f7f7f", purple="#6b4c9a", light="#d9d9d9", ink="#40566b")
MM = 1 / 25.4

SETS = ["SEG1_ADIPOKINE_CORE", "SEG2_LIPID_PATHWAY", "SEG3_EXECUTOR",
        "ALMP_COMPOSITIONAL_AXIS", "ALMP_TUMOR_INTRINSIC_AXIS", "ALMP_ALL_308"]
SHORT = {"SEG1_ADIPOKINE_CORE": "SEG1 adipokine core (22)",
         "SEG2_LIPID_PATHWAY": "SEG2 lipid pathway (262)",
         "SEG3_EXECUTOR": "SEG3 executor (39)",
         "ALMP_COMPOSITIONAL_AXIS": "Compositional axis (16)",
         "ALMP_TUMOR_INTRINSIC_AXIS": "Intrinsic axis (292)",
         "ALMP_ALL_308": "ALL_308 (308)"}

BSHORT = ["SEG1 core (22)", "SEG2 lipid (262)", "SEG3 exec (39)",
          "Composition (16)", "Intrinsic (292)", "ALL_308"]

W, H = 180.0, 252.0
fig = plt.figure(figsize=(W * MM, H * MM))


def rect(x, y, w, h):
    return [x / W, y / H, w / W, h / H]


def txt(x, y, s, **kw):
    kw.setdefault("va", "bottom")
    fig.text(x / W, y / H, s, **kw)


def head(letter, title, top, x, subs=()):
    """抬头：以坐标区顶边 top(mm) 为基准向上堆叠。subs 依次向下排。"""
    txt(x, top + 10.4, letter, fontsize=9.0, fontweight="bold")
    txt(x + 5.6, top + 10.4, title, fontsize=7.0, fontweight="bold")
    for k, s in enumerate(subs):
        txt(x, top + 6.5 - 2.9 * k, s, fontsize=5.0, color=C["ink"], style="italic")


# ================================================================ a
d = pd.read_csv(os.path.join(ROOT, "output/P2-T5b_多去卷积一致性v1/data/"
                                   "P2T5b_crossmethod_corr.tsv"), sep="\t")
labels = [r for r in d["feature"]]
M = d[labels].values.astype(float)
methods = [l.split(".")[0] for l in labels]
blocks, cur, start = [], methods[0], 0
for i in range(1, len(methods) + 1):
    if i == len(methods) or methods[i] != cur:
        blocks.append((cur, start, i)); start = i
        if i < len(methods):
            cur = methods[i]

A_TOP, A_H = 238, 58
axa = fig.add_axes(rect(21, A_TOP - A_H, A_H, A_H))
im = axa.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-0.8, vcenter=0, vmax=1.0),
                interpolation="nearest")
axa.set_xticks(range(len(labels)))
axa.set_yticks(range(len(labels)))
axa.set_xticklabels([l.split(".", 1)[1].replace(".", " ") for l in labels],
                    rotation=90, fontsize=3.9)
axa.set_yticklabels([l.split(".", 1)[1].replace(".", " ") for l in labels], fontsize=3.9)
for (nm, a_, b_) in blocks:
    axa.add_patch(plt.Rectangle((a_ - 0.5, a_ - 0.5), b_ - a_, b_ - a_, fill=False,
                                edgecolor="#111111", linewidth=0.75))
axa.set_xticks(np.arange(len(labels) + 1) - 0.5, minor=True)
axa.set_yticks(np.arange(len(labels) + 1) - 0.5, minor=True)
axa.grid(which="minor", color="white", linewidth=0.22)
axa.tick_params(which="minor", length=0)
axa.tick_params(axis="both", length=1.2, width=0.45, pad=1.0)
for b_ in axa.spines.values():
    b_.set_linewidth(0.55)
caxa = fig.add_axes(rect(81.2, A_TOP - A_H + 14, 3.0, A_H - 28))
cba = fig.colorbar(im, cax=caxa)
cba.set_label("Spearman's $\\rho$", fontsize=5.2, labelpad=2)
cba.ax.tick_params(labelsize=4.6, width=0.45, length=1.3)
cba.outline.set_linewidth(0.5)
head("a", "Cross-framework deconvolution concordance", A_TOP, 12.0,
     ["31 cell-type features; three scoring variants (quanTIseq,", "MCP-counter, MCP-counter QN)."])

# ================================================================ b（两热图共享 x 标签）
coup = pd.read_csv(os.path.join(R50, "P2T5_coupling_fulltranscriptome.tsv"), sep="\t")
imm = list(dict.fromkeys(coup.immune))
lab = ["::".join(s.split("::")[1:]) for s in imm]
grp = {"PanglaoDB": (0, 10), "CellMarker": (11, 18), "Hallmark": (19, 26)}

B1_TOP, B1_H = 238, 19
B2_TOP, B2_H = 213, 19
for meth, top, hh, show_x in [("ssGSEA", B1_TOP, B1_H, False),
                              ("GSVA", B2_TOP, B2_H, True)]:
    axb = fig.add_axes(rect(122, top - hh, 46, hh))
    sub = coup[coup.method == meth].pivot_table(index="almp_set", columns="immune",
                                                values="rho")
    sub = sub.loc[SETS, imm]
    imb = axb.imshow(sub.values, cmap="RdBu_r", vmin=-0.45, vmax=0.55, aspect="auto")
    axb.set_xticks(range(len(imm)))
    if show_x:
        axb.set_xticklabels(lab, rotation=90, fontsize=3.2)
        axb.tick_params(axis="x", length=1.2, pad=1.0)
    else:
        axb.set_xticklabels([])
        axb.tick_params(axis="x", length=0)
    axb.set_yticks(range(len(SETS)))
    axb.set_yticklabels(BSHORT if meth == "ssGSEA" else [], fontsize=4.3)
    for xx in (10.5, 20.5):
        axb.axvline(xx, color="black", lw=0.65)
    axb.set_xlim(-0.5, len(imm) + 6.4)
    axb.set_ylim(len(SETS) - 0.5, -0.5)
    axb.tick_params(axis="y", length=1.2, pad=1.0)
    med = sub.median(axis=1)
    for i in range(len(SETS)):
        axb.text(len(imm) + 1.1, i, f"{med.iloc[i]:+.3f}", fontsize=4.3, va="center",
                 ha="left", color=C["red"] if i == 0 else "#444444",
                 fontweight="bold" if i == 0 else "normal")
    txt(104.0, top - hh / 2.0, meth, fontsize=5.2, fontweight="bold", va="center",
        rotation=90)

txt(168.0, B1_TOP + 1.0, "median $\\rho$", fontsize=4.0, va="bottom", ha="right",
    color=C["red"], fontweight="bold")
for gname, (a_, b_) in grp.items():
    txt(122 + (a_ + b_) / 2.0 * 46.0 / 31.0, B1_TOP + 1.0, gname, fontsize=4.0,
        va="bottom", ha="center", color=C["blue"])
caxb = fig.add_axes(rect(169.9, B2_TOP + 2, 2.4, B1_TOP - B2_TOP - 4))
cbb = fig.colorbar(imb, cax=caxb)
cbb.set_label("$\\rho$", fontsize=4.8, labelpad=1.0)
cbb.ax.tick_params(labelsize=4.2, width=0.45, length=1.2)
cbb.outline.set_linewidth(0.5)
head("b", "Within-set comparison across six segments", B1_TOP, 96.0,
     ["TCGA-SKCM, n = 472; full-transcriptome input, same algorithm both sides.",
      "Descriptive only; segment ordering is not supported (see Results §3.3)."])

# ================================================================ c
prune = pd.read_csv(os.path.join(R50, "shared_gene_pruning_SEG1.tsv"), sep="\t")
pr = prune[prune.n_shared > 0].copy().sort_values("delta")
pr["lab"] = ["::".join(s.split("::")[1:]) for s in pr.immune]

C_TOP, C_H = 148, 22
axc = fig.add_axes(rect(30, C_TOP - C_H, 128, C_H))
for i, (_, r) in enumerate(pr.iterrows()):
    axc.plot([r.rho_orig, r.rho_pruned], [i, i], "-", color=C["light"], lw=1.8, zorder=1)
    axc.plot(r.rho_orig, i, "o", mfc="white", mec=C["grey"], mew=0.85, ms=3.4, zorder=2)
    axc.plot(r.rho_pruned, i, "o", mfc=C["red"], mec=C["red"], ms=3.4, zorder=3)
    axc.text(0.612, i, f"{int(r.n_shared)}", fontsize=4.4, va="center", ha="left",
             color="#777777")
axc.set_yticks(range(len(pr)))
axc.set_yticklabels(pr.lab, fontsize=5.0)
axc.set_xlim(0.10, 0.66)
axc.set_ylim(-0.7, len(pr) - 0.3)
axc.set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6])
axc.tick_params(axis="both", length=1.6, pad=1.2)
axc.set_xlabel("Spearman's $\\rho$ (SEG1 core vs immune signature)", labelpad=1.5)
axc.text(0.612, len(pr) - 0.55, "shared", fontsize=4.0, va="bottom", ha="left",
         color="#777777")
for sp in ("top", "right"):
    axc.spines[sp].set_visible(False)
axc.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=C["grey"],
                           mew=0.85, ms=3.2, label="before pruning"),
                    Line2D([], [], marker="o", ls="", mfc=C["red"], mec=C["red"],
                           ms=3.2, label="after pruning")],
           loc="lower left", ncol=2, handletextpad=0.4, columnspacing=1.4,
           borderpad=0.0, bbox_to_anchor=(0.0, 1.01))
head("c", "Shared-gene pruning", C_TOP, 12.0,
     ["Coupling persists or strengthens in every set after shared genes are removed from both sides."])

# ================================================================ d
F = pd.read_csv(os.path.join(SRC, "P2T20_F_fibroblast_mediation.tsv"), sep="\t")
D = pd.read_csv(os.path.join(SRC, "P2T20_D_covariate_ladder.tsv"), sep="\t")
s = F[(F.algorithm == "ssGSEA") & (F.mcp == "noQN") &
      (F.almp == "SEG1_ADIPOKINE_CORE") & (F.target != "Fibroblasts")].copy()
s = s.sort_values("rho", ascending=True).reset_index(drop=True)

D_TOP, D_H = 102, 22
axd = fig.add_axes(rect(36, D_TOP - D_H, 124, D_H))
yd = np.arange(len(s))
for i, r in s.iterrows():
    axd.plot([r.rho_c_fib, r.rho], [i, i], color=C["light"], lw=1.3, zorder=1,
             solid_capstyle="round")
axd.scatter(s.rho_c_fib, yd, s=15, facecolor="white", edgecolor=C["grey"],
            linewidth=0.65, zorder=3)
axd.scatter(s.rho, yd, s=15, facecolor=C["red"], edgecolor="none", zorder=4)
axd.axvline(0, color="#4d4d4d", lw=0.5, ls=(0, (2, 1.6)), zorder=0)
axd.set_yticks(yd)
axd.set_yticklabels(s.target, fontsize=5.0)
axd.set_ylim(-0.65, len(s) - 0.35)
axd.set_xlim(-0.32, 0.50)
axd.set_xticks([-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5])
axd.tick_params(axis="both", length=1.6, pad=1.2)
axd.set_xlabel("Spearman's $\\rho$ (SEG1 core vs inferred population)", labelpad=1.5)
for sp in ("top", "right"):
    axd.spines[sp].set_visible(False)
for i, r in s.iterrows():
    keep = r.rho_c_fib / r.rho * 100 if abs(r.rho) > 1e-9 else np.nan
    if r.rho >= 0:
        axd.annotate("%+.0f%%" % keep, xy=(r.rho + 0.011, i), va="center", ha="left",
                     fontsize=4.7, color="#666666")
    else:
        axd.annotate("%+.0f%%" % keep, xy=(min(r.rho, r.rho_c_fib) - 0.011, i),
                     va="center", ha="right", fontsize=4.7, color="#666666")
axd.legend(handles=[Line2D([], [], marker="o", ls="", mfc="white", mec=C["grey"],
                           mew=0.65, ms=3.2, label="adjusted for fibroblast content"),
                    Line2D([], [], marker="o", ls="", mfc=C["red"], mec="none",
                           ms=3.2, label="unadjusted")],
           loc="lower left", ncol=2, handletextpad=0.4, columnspacing=1.4,
           borderpad=0.0, bbox_to_anchor=(0.0, 1.01))
head("d", "Composition adjustment", D_TOP, 12.0,
     ["TCGA-SKCM, n = 470; ssGSEA, MCP-counter without quantile normalisation. Single covariates only:",
      "stromal and immune content correlate at r = 0.81, so no joint model is reported."])

# ================================================================ e
fib = F[(F.algorithm == "ssGSEA") & (F.mcp == "noQN") & (F.target == "Fibroblasts")]
fib = fib.set_index("almp").loc[SETS].reset_index()
seg1 = D[(D.algorithm == "ssGSEA") & (D.mcp == "noQN") &
         (D.target == "Fibroblasts")].iloc[0]


def star(p):
    return "***" if p < 1e-3 else ("**" if p < 1e-2 else ("*" if p < 0.05 else "n.s."))


E_TOP, E_H = 56, 30
axe1 = fig.add_axes(rect(28, E_TOP - E_H, 66, E_H))
cols = [C["red"] if a == "SEG1_ADIPOKINE_CORE" else
        (C["blue"] if a == "ALMP_TUMOR_INTRINSIC_AXIS" else C["grey"]) for a in fib.almp]
xe = np.arange(len(fib))
axe1.bar(xe, fib.rho.values, width=0.66, color=cols, edgecolor="none", zorder=2)
axe1.axhline(0, color="#4d4d4d", lw=0.5)
axe1.set_xticks(xe)
axe1.set_xticklabels(["SEG1 core (22)", "SEG2 lipid (262)", "SEG3 exec (39)",
                      "Composition (16)", "Intrinsic (292)", "ALL_308"],
                     rotation=32, ha="right", fontsize=4.4)
axe1.set_ylabel("$\\rho$ vs fibroblast content", fontsize=5.6, labelpad=1.5)
axe1.set_ylim(-0.18, 0.62)
for sp in ("top", "right"):
    axe1.spines[sp].set_visible(False)
for xi, r in zip(xe, fib.itertuples()):
    axe1.annotate(star(r.p), xy=(xi, r.rho + (0.022 if r.rho >= 0 else -0.022)),
                  ha="center", va="bottom" if r.rho >= 0 else "top", fontsize=4.6,
                  color="#333333")
axe1.tick_params(axis="both", length=1.6, pad=1.2)

axe2 = fig.add_axes(rect(112, E_TOP - E_H, 58, E_H))
lab2 = ["unadjusted", "+ endothelial", "+ malignant", "+ immune", "+ stromal\n& malignant"]
val2 = [seg1.rho_raw, seg1.rho_c_str, seg1.rho_c_pur, seg1.rho_c_imm, seg1.rho_c_str_pur]
x2 = np.arange(len(val2))
greys = [C["grey"], C["orange"], C["orange"], C["orange"], C["purple"]]
axe2.bar(x2, val2, width=0.66, color=greys, edgecolor="none", zorder=2)
axe2.set_xticks(x2)
axe2.set_xticklabels(lab2, rotation=32, ha="right", fontsize=4.4)
axe2.set_ylabel("SEG1 $\\rho$ vs fibroblast content", fontsize=5.6, labelpad=1.5)
axe2.set_ylim(0, 0.63)
for sp in ("top", "right"):
    axe2.spines[sp].set_visible(False)
base = val2[0]
for xi, v in zip(x2, val2):
    axe2.annotate("%.0f%%" % (v / base * 100), xy=(xi, v + 0.016), ha="center",
                  va="bottom", fontsize=4.6, color="#333333")
axe2.annotate("retention 82–93%", xy=(0.5, 0.97), xycoords="axes fraction",
              ha="center", va="top", fontsize=4.8, color="#8c8c8c", style="italic")
axe2.tick_params(axis="both", length=1.6, pad=1.2)
head("e", "Axis decomposition", E_TOP, 12.0,
     ["Left: each score against fibroblast content. Right: SEG1 under each single covariate."])

fig.savefig(os.path.join(FIG, "fig4_group.pdf"))
fig.savefig(os.path.join(FIG, "fig4_group.png"), dpi=300)
plt.close(fig)
print("OK -> fig4_group.pdf / .png  (%.0f x %.0f mm)" % (W, H))
