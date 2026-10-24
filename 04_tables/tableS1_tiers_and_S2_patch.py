# -*- coding: utf-8 -*-
"""S1 v2：加细胞自主分档列；S2 v2：更正 S2f 为 171 基因口径"""
import os, pandas as pd, numpy as np
D = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output/supp"

# ---------- S1 v2 ----------
s1 = pd.read_csv(os.path.join(D, "S1_ALMP_308gene_list_and_axis_v1.csv"))
s5 = pd.read_csv(os.path.join(D, "S5_ALMP_gene_detection_celllines_v2.csv"))
d = s5.set_index("gene_symbol").detection_fraction_50_melanoma
tier, sub = [], []
for _, r in s1.iterrows():
    if r.axis != "tumour-intrinsic":
        tier.append("not_applicable"); sub.append("no")
        continue
    v = d.get(r.gene_symbol, np.nan)
    if pd.isna(v):
        tier.append("not_assessed"); sub.append("no")
    elif v >= 0.90:
        tier.append("supported"); sub.append("yes")
    elif v < 0.20:
        tier.append("non_attributable"); sub.append("no")
    else:
        tier.append("intermediate"); sub.append("no")
s1["cell_autonomous_support"] = tier
s1["cell_autonomous_subset"] = sub
s1.to_csv(os.path.join(D, "S1_ALMP_308gene_list_and_axis_v2.csv"), index=False, encoding="utf-8-sig")
print("S1 v2", s1.shape)
print(pd.crosstab(s1.axis, s1.cell_autonomous_support))
print("subset=yes:", int((s1.cell_autonomous_subset == "yes").sum()))

# ---------- S2 v2 ----------
s2 = pd.read_csv(os.path.join(D, "S2_sensitivity_analyses_v1.csv"))
NEW = {  # dataset -> (rho, p, n_genes)
    "472 TCGA-SKCM tumours": ("-0.097", "0.0349"),
    "50 melanoma cell lines (marker audit)": ("-0.288", "0.0426"),
    "7 fibroblast lines (positive control)": ("+0.214", "0.6445"),
}
m = s2.comparison.str.contains("166 cell-autonomous", na=False)
print("\nS2 rows to patch:", int(m.sum()))
for i in s2.index[m]:
    lvl = s2.at[i, "level"]
    rho, p = NEW[lvl]
    old = (s2.at[i, "comparison"], s2.at[i, "value"], s2.at[i, "note"])
    s2.at[i, "comparison"] = "core vs intrinsic, 171 cell-autonomous genes"
    s2.at[i, "value"] = rho
    s2.at[i, "note"] = f"p = {p}; 171 genes"
    print("  ", lvl, old, "->", (s2.at[i, "comparison"], rho, s2.at[i, "note"]))
s2.to_csv(os.path.join(D, "S2_sensitivity_analyses_v2.csv"), index=False, encoding="utf-8-sig")
print("S2 v2", s2.shape)
