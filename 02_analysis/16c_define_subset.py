# -*- coding: utf-8 -*-
"""生成 171 基因细胞自主子集（分母更正：50 黑素瘤系，而非 57 株含 7 成纤维系）"""
import os, pandas as pd, numpy as np
ROOT = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output"
D = os.path.join(ROOT,"P5-T1_细胞系内在轴验证v1","data")
S5 = pd.read_csv(os.path.join(ROOT,"supp","S5_ALMP_gene_detection_celllines_v2.csv"))
intr = S5[(S5.axis=="tumour-intrinsic")&(S5.detection_fraction_50_melanoma>=0.90)].copy()
out = pd.DataFrame({
  "gene_symbol": intr.gene_symbol.values,
  "detect_ccle": intr.detection_fraction_50_melanoma.values,
  "detect_tcga": intr.detection_fraction_472_TCGA_tumours.values,
  "segment_primary": intr.segment.values,
  "detection_fraction_all_57": intr.detection_fraction_all_57.values,
}).sort_values("gene_symbol").reset_index(drop=True)
p = os.path.join(D,"P5T1_cell_autonomous_intrinsic_subset_v2.tsv")
out.to_csv(p,sep="\t",index=False)
print("wrote",p,out.shape)
v1 = pd.read_csv(os.path.join(D,"P5T1_cell_autonomous_intrinsic_subset.tsv"),sep="\t")
print("v1",v1.shape,"v2",out.shape)
print("added:",sorted(set(out.gene_symbol)-set(v1.gene_symbol)))
print("removed:",sorted(set(v1.gene_symbol)-set(out.gene_symbol)))
