# -*- coding: utf-8 -*-
"""166 vs 171：细胞自主子集的分母更正影响评估"""
import os, json
import numpy as np, pandas as pd
from scipy import stats

ROOT = r"C:/Users/刘宸瑞/Desktop/脂肪因子生信/manuscript/output"
D  = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
S5 = pd.read_csv(os.path.join(ROOT,"supp","S5_ALMP_gene_detection_celllines_v2.csv"))

gm = json.load(open(os.path.join(D,"gene_map.json"))); e2s={int(k):v for k,v in gm["entrez2symbol"].items()}
cc = pd.read_csv(os.path.join(D,"ccle_melanoma_308genes_rpkm.tsv"),sep="\t")
cc["symbol"]=cc.entrezGeneId.astype(int).map(e2s)
P = cc.dropna(subset=["symbol"]).pivot_table(index="sampleId",columns="symbol",values="rpkm",aggfunc="mean")
gl = pd.read_csv(os.path.join(D,"P5T1_gene_level_detection.tsv"),sep="\t")
core=[g for g in gl[gl.segment_primary=="SEG1 adipokine core"].gene_symbol if g in P.columns]
intr=[g for g in gl[gl.axis=="tumour-intrinsic"].gene_symbol if g in P.columns]

ctl = pd.read_csv(os.path.join(D,"ccle_melanoma_controlgens_rpkm.tsv"),sep="\t")
C = ctl.pivot_table(index="sampleId",columns="symbol",values="rpkm",aggfunc="mean")
mel_mk=["MLANA","TYR","SOX10","S100B","PMEL","DCT"]; fib_mk=["COL1A1","COL1A2","DCN","LUM","THY1","PDGFRB"]
aud=pd.DataFrame({"melanocytic_max":C[mel_mk].max(axis=1),"fibroblast_median":C[fib_mk].median(axis=1)}).loc[P.index]
aud["call"]=np.where((aud.melanocytic_max<1)&(aud.fibroblast_median>=10),"fibroblast","melanoma")
fib=aud.index[aud.call=="fibroblast"].tolist(); mel=aud.index[aud.call=="melanoma"].tolist()
print("melanoma",len(mel),"fibroblast",len(fib))

tc=pd.read_csv(os.path.join(os.path.expanduser("~"),"WorkBuddy","2026-10-07-22-15-36",".workbuddy","data","cbio_skcm","mrna_data.tsv"),sep="\t")
tc["symbol"]=tc.entrezGeneId.astype(int).map(e2s)
T=tc.dropna(subset=["symbol"]).pivot_table(index="sampleId",columns="symbol",values="value",aggfunc="mean")
print("T",T.shape,"P",P.shape)

hi=pd.read_csv(os.path.join(D,"P5T1_cell_autonomous_intrinsic_subset.tsv"),sep="\t")
set166=[g for g in hi.gene_symbol if g in T.columns and g in P.columns]
intr57=S5[(S5.axis=="tumour-intrinsic")&(S5.detection_fraction_all_57>=0.90)].gene_symbol.tolist()
intr50=S5[(S5.axis=="tumour-intrinsic")&(S5.detection_fraction_50_melanoma>=0.90)].gene_symbol.tolist()
set171=[g for g in intr50 if g in T.columns and g in P.columns]
print("set166 n=%d  set171(from S5 50-line) n=%d"%(len(set166),len(set171)))
print("171 minus 166:",sorted(set(set171)-set(set166)))
print("166 minus 171:",sorted(set(set166)-set(set171)))

def rr(sub,genes):
    genes=[g for g in genes if g in sub.columns]
    co=[g for g in core if g in sub.columns]
    a=np.log1p(sub[co].clip(lower=0)).mean(axis=1); b=np.log1p(sub[genes].clip(lower=0)).mean(axis=1)
    r=stats.spearmanr(a,b); return float(r.correlation),float(r.pvalue),len(genes)

rows=[]
for tag,sub in [("tumours(TCGA)",T),("melanoma50",P.loc[mel]),("fibroblast7",P.loc[fib])]:
    for gname,genes in [("all intrinsic 288",intr),("166 frozen",set166),("171 corrected",set171)]:
        r,p,k=rr(sub,genes); rows.append(dict(dataset=tag,set=gname,n=k,rho=r,p=p))
Q=pd.DataFrame(rows); print("\n=== 2x3 ==="); print(Q.round(4).to_string(index=False))
Q.to_csv(os.path.join(D,"P5T1_rho_166_vs_171.tsv"),sep="\t",index=False)
print("\nsaved ->",D)
