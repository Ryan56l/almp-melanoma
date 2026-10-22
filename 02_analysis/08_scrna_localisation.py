# -*- coding: utf-8 -*-
"""
P4-T1~T3：单细胞分辨率下 ALMP 的细胞归属（指令四十三）
=========================================================
输入
  data/GSE72056.txt.gz   Tirosh 2016（GSE72056）黑色素瘤 scRNA，TPM 矩阵
                         文件自带原作者细胞标签：malignant(0/1) 与 non-malignant type
                         1=T 2=B 3=Macrophage 4=Endothelial 5=CAF 6=NK
  output/P2-T1_基因集v1/ALMP_gene_axis_assignment_v1.csv   308 基因＋轴归属

方法（含三处须写入 Methods 的诚实声明）
  ① 不做 de novo 聚类与 UMAP/Harmony 整合 —— 原始研究已提供人工注释的细胞类型；
     Harmony/scVI 的目的是"为聚类消除批次"，本问不需要，故不引入。
  ② 打分用 UCell（rank-based，maxRank=1500）与 AUCell（top 5% 恢复曲线）双法 ——
     两者都只依赖"每细胞内的基因排序"，对文库大小与测序深度天然稳健，适合跨患者合并。
  ③ 不做配受体推断（CellChat/NicheNet）—— 脂肪细胞在 scRNA 中捕获率近零，
     配受体推断无可信度（用户裁定②）。

统计
  - Kruskal-Wallis：全体细胞类型间分数差异
  - pairwise Mann-Whitney vs malignant（BH 校正）
  - 患者内配对：每患者每细胞类型取中位分数，Wilcoxon signed-rank（malignant vs 各型）
    —— 用于排除"患者构成差异"造成的假象（防 Simpson 悖论）

输出（results/）
  almp_scores_per_cell.tsv  每细胞三套分数
  celltype_summary.tsv      细胞类型 × 基因集的描述统计
  kruskal.tsv / pairwise_vs_malignant.tsv / patient_paired.tsv
  fig_score_by_celltype.png
"""
import os
import sys
import io
import gzip
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
DATA = os.path.join(HERE, "..", "data", "GSE72056.txt.gz")
RES = os.path.join(HERE, "..", "results")
os.makedirs(RES, exist_ok=True)
GS = os.path.join(ROOT, "output", "P2-T1_基因集v1", "ALMP_gene_axis_assignment_v1.csv")

MAXRANK = 1500
TOPFRAC = 0.05
TYPE = {1: "T", 2: "B", 3: "Macrophage", 4: "Endothelial", 5: "CAF", 6: "NK"}


# ---------------- 读取 GSE72056 ----------------
def load_gse72056(path):
    # ★ 实测布局（2026-10-21 核对）：共 4 行注释，其后为基因
    #   row0 = "Cell"          + 细胞名
    #   row1 = "tumor"         + 患者编号
    #   row2 = "malignant(1=no,2=yes,0=unresolved)" + 代码
    #   row3 = "non-malignant cell type(1=T,2=B,3=Macro,4=Endo,5=CAF,6=NK)" + 代码
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        r0 = f.readline().rstrip("\n").split("\t")
        r1 = f.readline().rstrip("\n").split("\t")
        r2 = f.readline().rstrip("\n").split("\t")
        r3 = f.readline().rstrip("\n").split("\t")
        print(f"[0] 行标签：{r0[0]!r} / {r1[0]!r} / {r2[0]!r} / {r3[0]!r}")
        cells = [c.replace('"', "") for c in r0[1:]]
        patient = [x.replace('"', "") for x in r1[1:]]
        mal = np.array([int(x) for x in r2[1:]], dtype=int)   # 2 = malignant
        typ = np.array([int(x) for x in r3[1:]], dtype=int)
        n_cells = len(cells)
        genes, rows = [], []
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) != n_cells + 1:
                continue
            genes.append(p[0].replace('"', ""))
            rows.append(p[1:])
    print(f"[0] 基因 {len(genes)} × 细胞 {n_cells}；malignant 代码分布 {dict(zip(*np.unique(mal, return_counts=True)))}")
    X = np.array(rows, dtype=np.float32)
    return X, genes, cells, mal, typ, patient


# ---------------- 打分 ----------------
def rank_matrix(X, i0, i1):
    """返回 cells × genes 的秩（1 = 表达最高）"""
    sub = X[:, i0:i1]                      # genes × cells
    order = np.argsort(-sub, axis=0)       # 每列降序的基因下标
    r = np.empty_like(order, dtype=np.uint16)
    for j in range(order.shape[1]):
        r[order[:, j], j] = np.arange(1, order.shape[0] + 1, dtype=np.uint16)
    return r.T                             # cells × genes


def ucell(ranks, idx):
    n = min(len(idx), MAXRANK)
    sub = np.sort(ranks[:, idx], axis=1)[:, :n].astype(np.float64)
    s = np.where(sub <= MAXRANK, sub - 1, MAXRANK).sum(axis=1)
    return np.clip(1.0 - s / (MAXRANK * n), 0.0, 1.0)


def aucell(ranks, idx, n_genes):
    n_top = int(round(TOPFRAC * n_genes))
    k = len(idx)
    pos = np.sort(ranks[:, idx], axis=1)
    pos = np.where(pos <= n_top, pos, n_top + 1)
    pos = np.sort(pos, axis=1)
    # AUC：Σ (n_top - pos_i + 1)，除以理想情形 Σ (n_top - i + 1)
    num = np.where(pos <= n_top, n_top - pos + 1, 0).sum(axis=1, dtype=np.float64)
    den = (n_top * k) - (k * (k + 1) / 2.0) + k * 1.0
    return np.clip(num / den, 0.0, 1.0)


def main():
    X, genes, cells, mal, typ, patient = load_gse72056(DATA)
    gs = pd.read_csv(GS, encoding="utf-8-sig")
    gs["gene"] = gs["gene"].str.upper().str.strip()
    up = [g.upper().strip() for g in genes]

    sets = {
        "SEG1_ADIPOKINE_CORE": gs[gs["segments"].str.contains("脂肪因子核心", na=False)]["gene"].tolist(),
        "AXIS_INTRINSIC": gs[gs["axis"] == "肿瘤-免疫内在轴"]["gene"].tolist(),
        "AXIS_TISSUE": gs[gs["axis"] == "组织成分轴"]["gene"].tolist(),
    }
    idxs, cover = {}, {}
    for k, gl in sets.items():
        found = [i for i, g in enumerate(up) if g in set(gl)]
        idxs[k] = found
        cover[k] = (len(found), len(gl))
        print(f"[1] {k}: 命中 {len(found)}/{len(gl)}")

    # QC
    detected = (X > 0).sum(axis=0)
    mito = np.array([up[i].startswith("MT-") for i in range(len(up))])
    mito_frac = X[mito, :].sum(axis=0) / np.maximum(X.sum(axis=0), 1e-9)
    keep = (detected >= 1000) & (mito_frac < 0.20)
    print(f"[2] QC：细胞 {len(cells)} → 保留 {keep.sum()}（检出基因≥1000 且 线粒体<20%）")
    Xk = X[:, keep]
    cells_k = [c for c, kp in zip(cells, keep) if kp]
    mal_k, typ_k = mal[keep], typ[keep]

    # 细胞类型标签
    labels = []
    for m, t in zip(mal_k, typ_k):
        if m == 2:                      # 2 = malignant
            labels.append("Malignant")
        elif m == 0:
            labels.append("Unresolved")
        else:
            labels.append(TYPE.get(int(t), f"other_{int(t)}"))
    labels = np.array(labels)
    print("[3] 细胞类型分布：", dict(pd.Series(labels).value_counts()))

    # 打分（分块，控制内存）
    scores = {k: np.zeros(len(cells_k), dtype=np.float64) for k in sets}
    aucs = {k: np.zeros(len(cells_k), dtype=np.float64) for k in sets}
    step = 500
    for s in range(0, len(cells_k), step):
        e = min(s + step, len(cells_k))
        rk = rank_matrix(Xk, s, e)
        for k in sets:
            scores[k][s:e] = ucell(rk, idxs[k])
            aucs[k][s:e] = aucell(rk, idxs[k], Xk.shape[0])
        if s % 2000 == 0:
            print(f"    打分进度 {e}/{len(cells_k)}")

    pat_k = [pt for pt, kp in zip(patient, keep) if kp]
    df = pd.DataFrame({"cell": cells_k, "celltype": labels, "patient": pat_k,
                       "malignant_code": mal_k, "nonmal_type_code": typ_k})
    for k in sets:
        df[f"UCell_{k}"] = scores[k]
        df[f"AUCell_{k}"] = aucs[k]
    df.to_csv(os.path.join(RES, "almp_scores_per_cell.tsv"), sep="\t", index=False)

    # 描述统计
    rows = []
    for k in sets:
        for method in ("UCell", "AUCell"):
            col = f"{method}_{k}"
            g = df.groupby("celltype")[col]
            for ct, v in g:
                rows.append({"geneset": k, "method": method, "celltype": ct,
                             "n": len(v), "median": float(np.median(v)),
                             "mean": float(np.mean(v)), "q25": float(np.percentile(v, 25)),
                             "q75": float(np.percentile(v, 75))})
    summ = pd.DataFrame(rows)
    summ.to_csv(os.path.join(RES, "celltype_summary.tsv"), sep="\t", index=False)

    # Kruskal + pairwise
    kr, pw = [], []
    cts = [c for c in df["celltype"].unique() if (df["celltype"] == c).sum() >= 20]
    for k in sets:
        for method in ("UCell", "AUCell"):
            col = f"{method}_{k}"
            groups = [df.loc[df["celltype"] == c, col].values for c in cts]
            if len(groups) >= 2:
                h, p = stats.kruskal(*groups)
                kr.append({"geneset": k, "method": method, "k_groups": len(groups),
                           "H": float(h), "p": float(p)})
            if "Malignant" in cts:
                base = df.loc[df["celltype"] == "Malignant", col].values
                for c in cts:
                    if c == "Malignant":
                        continue
                    v = df.loc[df["celltype"] == c, col].values
                    u, p = stats.mannwhitneyu(base, v, alternative="two-sided")
                    auc = u / (len(base) * len(v))
                    pw.append({"geneset": k, "method": method, "vs": c,
                               "n_malignant": len(base), "n_other": len(v),
                               "median_malignant": float(np.median(base)),
                               "median_other": float(np.median(v)),
                               "cliff_delta": float(2 * auc - 1), "p": float(p)})
    pd.DataFrame(kr).to_csv(os.path.join(RES, "kruskal.tsv"), sep="\t", index=False)

    pwd = pd.DataFrame(pw)
    if len(pwd):
        from scipy.stats import false_discovery_control
        pwd["q"] = false_discovery_control(pwd["p"].values)
    pwd.to_csv(os.path.join(RES, "pairwise_vs_malignant.tsv"), sep="\t", index=False)

    # 患者内配对（防患者构成混杂）
    # GSE72056 的细胞名形如 "Cy72_1"；患者号取首段
    prs = []
    for k in sets:
        for method in ("UCell", "AUCell"):
            col = f"{method}_{k}"
            piv = df.groupby(["patient", "celltype"])[col].median().unstack()
            if "Malignant" not in piv.columns:
                continue
            for c in piv.columns:
                if c == "Malignant":
                    continue
                sub = piv[["Malignant", c]].dropna()
                if len(sub) < 5:
                    continue
                w, p = stats.wilcoxon(sub["Malignant"], sub[c])
                prs.append({"geneset": k, "method": method, "vs": c,
                            "n_patients": len(sub),
                            "median_delta_malignant_minus_other": float((sub["Malignant"] - sub[c]).median()),
                            "p": float(p)})
    prd = pd.DataFrame(prs)
    if len(prd):
        from scipy.stats import false_discovery_control
        prd["q"] = false_discovery_control(prd["p"].values)
    prd.to_csv(os.path.join(RES, "patient_paired.tsv"), sep="\t", index=False)

    # 图
    order = [c for c in ["Malignant", "T", "NK", "B", "Macrophage", "Endothelial", "CAF"] if c in cts]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2), sharey=False)
    for ax, k in zip(axes, sets):
        data = [df.loc[df["celltype"] == c, f"UCell_{k}"].values for c in order]
        ax.boxplot(data, tick_labels=order, showfliers=False)
        ax.set_title(f"UCell · {k}", fontsize=10)
        ax.tick_params(axis="x", rotation=45, labelsize=8)
        ax.set_ylabel("UCell score")
    plt.tight_layout()
    plt.savefig(os.path.join(RES, "fig_score_by_celltype.png"), dpi=200)
    print("[4] 图已保存 results/fig_score_by_celltype.png")

    # 摘要
    print("\n=== Kruskal-Wallis ===")
    print(pd.DataFrame(kr).to_string(index=False))
    print("\n=== 患者内配对（前 12 行） ===")
    print(prd.sort_values("p").head(12).to_string(index=False) if len(prd) else "（样本不足）")


if __name__ == "__main__":
    main()
