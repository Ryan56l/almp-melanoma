# -*- coding: utf-8 -*-
"""S5 v2：按谱系审计口径重建逐基因检出表

变更（相对 v1）：
  1. 57 株按标记审计拆为 50 黑素瘤 + 7 成纤维（HS600T/HS688AT/HS834T/HS839T/HS895T/HS934T/HS940T）
  2. 新增 detection_fraction_50_melanoma / detection_fraction_7_fibroblast / fold_fibroblast_over_melanoma
  3. 保留 v1 的 57 株合计列与 TCGA 列，便于逐列对账
"""
import os, json
import numpy as np
import pandas as pd

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
S5 = os.path.join(ROOT, "supp", "S5_ALMP_gene_detection_celllines_v1.csv")

gm = json.load(open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gm["entrez2symbol"].items()}
cc = pd.read_csv(os.path.join(D, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t")
cc["symbol"] = cc.entrezGeneId.astype(int).map(e2s)
P = cc.dropna(subset=["symbol"]).pivot_table(index="sampleId", columns="symbol",
                                             values="rpkm", aggfunc="mean")
audit = pd.read_csv(os.path.join(ROOT, "P5-T2_跨平台秩校准v1", "data",
                                 "P5T2_lineage_audit.tsv"), sep="\t")
call = dict(zip(audit.sampleId, audit["call"]))
mel = [s for s in P.index if call.get(s) == "melanoma"]
fib = [s for s in P.index if call.get(s) == "fibroblast"]
print("melanoma", len(mel), "fibroblast", len(fib), "total", len(P))

old = pd.read_csv(S5)
rows = []
for _, r in old.iterrows():
    g = r.gene_symbol
    if g in P.columns:
        v = P[g]
        vm, vf = v.loc[mel], v.loc[fib]
        dm = float((vm >= 1).mean()); df = float((vf >= 1).mean()); da = float((v >= 1).mean())
        mm = float(vm.median()); mf = float(vf.median())
        fold = (mf / mm) if mm > 0 else np.nan
    else:
        dm = df = da = mm = mf = fold = np.nan
    rows.append(dict(
        gene_symbol=g, axis=r.axis, segment=r.segment,
        meets_adipose_criterion=r.meets_adipose_criterion,
        n_lines_deposited_as_melanoma=int(len(P)),
        lineage_call_applied="50 melanoma + 7 fibroblast (marker audit)",
        detection_fraction_50_melanoma=dm,
        detection_fraction_7_fibroblast=df,
        detection_fraction_all_57=da,
        median_RPKM_50_melanoma=mm,
        median_RPKM_7_fibroblast=mf,
        fold_fibroblast_over_melanoma=fold,
        detection_fraction_472_TCGA_tumours=r.detection_fraction_472_TCGA_tumours,
        median_RSEM_tumours=r.median_RSEM_tumours,
    ))
new = pd.DataFrame(rows)
outp = os.path.join(ROOT, "supp", "S5_ALMP_gene_detection_celllines_v2.csv")
new.to_csv(outp, index=False, encoding="utf-8-sig")
print("wrote", outp, new.shape)

# 对账：v1 的 57 列应等于新 all_57 列
m = new.dropna(subset=["detection_fraction_all_57"])
d = (m.detection_fraction_all_57.astype(float) -
     old.set_index("gene_symbol").loc[m.gene_symbol,
        "detection_fraction_57_melanoma_cell_lines"].values.astype(float))
print("max|Δ| vs v1 all-57 column =", float(np.abs(d).max()))

# 关键块核对
for g in ["SERPINE1", "IL6", "ANGPTL4", "VEGFA", "ADIPOR1", "ADIPOR2", "NAMPT", "GPNMB",
          "ADIPOQ", "LEP", "RETN", "SERPINA12", "ITLN1"]:
    r = new[new.gene_symbol == g].iloc[0]
    print(f"{g:10s} mel={r.detection_fraction_50_melanoma:.3f} "
          f"fib={r.detection_fraction_7_fibroblast:.3f} "
          f"medR_mel={r.median_RPKM_50_melanoma:8.2f} medR_fib={r.median_RPKM_7_fibroblast:8.2f} "
          f"fold={r.fold_fibroblast_over_melanoma}")

# 轴层面汇总
core = new[new.segment == "SEG1 adipokine core"]
intr = new[new.axis == "tumour-intrinsic"]
print("core median det (mel) =", round(core.detection_fraction_50_melanoma.median(), 4),
      "(fib) =", round(core.detection_fraction_7_fibroblast.median(), 4))
print("intrinsic median det (mel) =", round(intr.detection_fraction_50_melanoma.median(), 4),
      "(fib) =", round(intr.detection_fraction_7_fibroblast.median(), 4))
