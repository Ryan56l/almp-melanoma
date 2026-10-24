# -*- coding: utf-8 -*-
"""②b 全转录组分位验证：检验「308 基因宇宙」结论是否依赖受限宇宙。

取 12 株黑色素瘤细胞系 + 12 例 TCGA-SKCM，抓全基因向量，
在每个样本内按**全部表达基因**（value>0）计算分位，比较核心段与内在轴。
"""
import sys, os, json, time
import numpy as np
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\tmp")
from netlib import post, API

ROOT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"
D = os.path.join(ROOT, "P5-T1_细胞系内在轴验证v1", "data")
OD = os.path.join(ROOT, "P5-T2_跨平台秩校准v1", "data")
rng = np.random.default_rng(20261023)

s1 = pd.read_csv(os.path.join(ROOT, "supp", "S1_ALMP_308gene_list_and_axis_v1.csv"))
gm = json.load(open(os.path.join(D, "gene_map.json")))
e2s = {int(k): v for k, v in gm["entrez2symbol"].items()}
sym2e = {v: k for k, v in e2s.items()}
ax = dict(zip(s1.gene_symbol, s1.axis))
seg = dict(zip(s1.gene_symbol, s1.segment_primary))
core = {g for g in s1.gene_symbol if seg.get(g) == "SEG1 adipokine core"}
intr = {g for g in s1.gene_symbol if ax.get(g) == "tumour-intrinsic"}

mel = sorted(json.load(open(os.path.join(D, "ccle_melanoma_lines.json"))))
tc = pd.read_csv(r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36\.workbuddy\data\cbio_skcm\mrna_data.tsv",
                 sep="\t")
tc_samples = sorted(tc.sampleId.unique())

PICK = {"CCLE": list(rng.choice(mel, 12, replace=False)),
        "TCGA": list(rng.choice(tc_samples, 12, replace=False))}
PROF = {"CCLE": "ccle_broad_2019_rna_seq_mrna",
        "TCGA": "skcm_tcga_pan_can_atlas_2018_rna_seq_v2_mrna"}
print(PICK, flush=True)


def job(args):
    ds, sid = args
    body = {"sampleIds": [sid]}
    try:
        d = post(f"{API}/molecular-profiles/{PROF[ds]}/molecular-data/fetch?projection=SUMMARY", body,
                 tries=5, base=4)
    except Exception as e:
        return ds, sid, None
    if not d:
        return ds, sid, None
    v = pd.Series({e2s.get(int(x["entrezGeneId"]), str(x["entrezGeneId"])): float(x["value"])
                   for x in d})
    v = v[~v.index.duplicated()]
    return ds, sid, v


tasks = [("CCLE", s) for s in PICK["CCLE"]] + [("TCGA", s) for s in PICK["TCGA"]]
res = []
with ThreadPoolExecutor(5) as ex:
    for k, r in enumerate(ex.map(job, tasks)):
        print(f"  {k+1}/{len(tasks)} {r[0]} {r[1]} n={0 if r[2] is None else len(r[2])}", flush=True)
        if r[2] is not None:
            res.append(r)

rows, detail = [], []
for ds, sid, v in res:
    expr = v[v > 0]
    n_all, n_expr = len(v), len(expr)
    pct = expr.rank(pct=True)          # 在「该样本全部表达基因」中的分位
    for gene, grp in [(g, "core") for g in core] + [(g, "intrinsic") for g in intr]:
        if gene in pct.index:
            detail.append(dict(dataset=ds, sample=sid, gene=gene, group=grp,
                               pct=float(pct[gene]), value=float(v[gene])))
    rows.append(dict(dataset=ds, sample=sid, n_genes_all=n_all, n_genes_expressed=n_expr))
    print(f"{ds:5s} {sid:26s} genes={n_all} expressed={n_expr}", flush=True)

det = pd.DataFrame(detail)
det.to_csv(os.path.join(OD, "P5T2_fulltranscriptome_gene_pct.tsv"), sep="\t", index=False)
pd.DataFrame(rows).to_csv(os.path.join(OD, "P5T2_fulltranscriptome_samples.tsv"), sep="\t", index=False)

summ = (det.groupby(["dataset", "group"]).pct
        .agg(["median", "count",
              lambda x: float((x >= .5).mean())]).reset_index())
summ.columns = ["dataset", "group", "median_pct", "n", "frac_top_half"]
by_sample = det.groupby(["dataset", "sample", "group"]).pct.median().reset_index()
contrast = by_sample.pivot_table(index=["dataset", "sample"], columns="group",
                                 values="pct").reset_index()
contrast["delta_core_minus_intrinsic"] = contrast["core"] - contrast["intrinsic"]
print("\n=== 全转录组分位（基因级） ===")
print(summ.to_string(index=False))
print("\n=== 每样本核心−内在 差值 ===")
print(contrast.groupby("dataset").delta_core_minus_intrinsic.describe().to_string())
summ.to_csv(os.path.join(OD, "P5T2_fulltranscriptome_summary.tsv"), sep="\t", index=False)
contrast.to_csv(os.path.join(OD, "P5T2_fulltranscriptome_contrast.tsv"), sep="\t", index=False)
print("saved")
