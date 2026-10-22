#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P2-T13 Step A：从 GTEx v8 eQTL 归档抽取 8 基因座 × 4 组织的 cis-eQTL 明细
================================================================================
输入（本地，已解包）：
  C:/Rtmp/smr/gtex/<Tissue>.v8.egenes.txt.gz                    ← 每基因 top 变异 + 阈值
  C:/Rtmp/smr/gtex/<Tissue>.v8.signif_variant_gene_pairs.txt.gz ← 显著 SNP-基因对明细

输出：
  data/gtex_egenes_top.tsv   每 (tissue,gene) 的 top 变异、基因坐标、阈值
  data/gtex_snp_pairs.tsv    每 (tissue,gene,SNP) 的 slope / slope_se / p / maf
  data/gtex_loci.tsv         每基因的 cis 窗口（±1 Mb，用于抓结局）

★ 关键口径
  - GTEx `slope` 为 **alt 等位基因**（= egenes 的 alt 列）对表达的效果；结局文件的
    effect_allele 须与之对齐，本脚本只在 Step C 做对齐，此处原样保留 ref/alt。
  - `signif_variant_gene_pairs` 只含 p 低于该基因阈值的 SNP（不含量化层全部位点），
    因此 HEIDI 可用位点数为「显著对数」，不是 cis 区域全位点数 —— 须在报告中声明。
  - 归档不含 allpairs → 无从补齐无效位点 → coloc 只能作支持性证据（SAP v1.6 勘误 1）。
"""
import gzip
import os
import sys
import csv

GTEX = r"C:/Rtmp/smr/gtex"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
os.makedirs(OUT, exist_ok=True)

TISSUES = [
    ("Adipose_Subcutaneous", "皮下脂肪"),
    ("Skin_Not_Sun_Exposed_Suprapubic", "皮肤非日晒"),
    ("Skin_Sun_Exposed_Lower_leg", "皮肤日晒"),
    ("Whole_Blood", "全血"),
]
GENES = ["LEP", "ADIPOQ", "NAMPT", "FABP4", "CD36", "FADS2", "CDS2", "SERPINE1"]
CIS_HALF = 1_000_000  # GTEx cis 窗口 ±1 Mb（与本项目 coloc/HEIDI 窗口一致）


def log(*a):
    print(*a, flush=True)


# ---------- 1. egenes：取基因坐标 + top 变异 + 阈值 ----------
eg_top, loci = [], {}
for tis, cn in TISSUES:
    path = os.path.join(GTEX, f"{tis}.v8.egenes.txt.gz")
    if not os.path.exists(path):
        log(f"[FATAL] 缺文件 {path}")
        sys.exit(1)
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        need = ["gene_id", "gene_name", "gene_chr", "gene_start", "gene_end", "num_var",
                "variant_id", "tss_distance", "chr", "variant_pos", "ref", "alt",
                "rs_id_dbSNP151_GRCh38p7", "maf", "pval_nominal", "slope", "slope_se",
                "qval", "pval_nominal_threshold"]
        for c in need:
            if c not in ix:
                log(f"[FATAL] {tis} egenes 缺列 {c}；实际列={hdr}")
                sys.exit(2)
        n = 0
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            g = p[ix["gene_name"]]
            if g not in GENES:
                continue
            chrom = p[ix["chr"]].lstrip("chr")
            eg_top.append({
                "tissue": tis, "tissue_cn": cn, "gene": g,
                "gene_id": p[ix["gene_id"]], "gene_chr": p[ix["gene_chr"]].lstrip("chr"),
                "gene_start": p[ix["gene_start"]], "gene_end": p[ix["gene_end"]],
                "num_var": p[ix["num_var"]],
                "top_variant_id": p[ix["variant_id"]],
                "top_rsid": p[ix["rs_id_dbSNP151_GRCh38p7"]],
                "top_chr": chrom, "top_pos": p[ix["variant_pos"]],
                "top_ref": p[ix["ref"]], "top_alt": p[ix["alt"]],
                "top_maf": p[ix["maf"]],
                "top_pval": p[ix["pval_nominal"]],
                "top_slope": p[ix["slope"]], "top_slope_se": p[ix["slope_se"]],
                "tss_distance": p[ix["tss_distance"]],
                "qval": p[ix["qval"]],
                "pval_threshold": p[ix["pval_nominal_threshold"]],
            })
            loci.setdefault(g, (chrom, int(p[ix["gene_start"]]), int(p[ix["gene_end"]])))
            n += 1
        log(f"[A1] {tis:34s} 命中 eGene {n}/8")

with open(os.path.join(OUT, "gtex_egenes_top.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(eg_top[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(eg_top)
log(f"[A1] → gtex_egenes_top.tsv  {len(eg_top)} 行")

# ---------- 2. signif_variant_gene_pairs：全部显著 SNP ----------
pairs = []
for tis, cn in TISSUES:
    path = os.path.join(GTEX, f"{tis}.v8.signif_variant_gene_pairs.txt.gz")
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(hdr)}
        # 该文件无 gene_name，只有 gene_id → 需用 egenes 的 gene_id→symbol 映射
        gid2sym = {r["gene_id"]: r["gene"] for r in eg_top if r["tissue"] == tis}
        if not gid2sym:
            log(f"[A2] {tis}: 无目标基因，跳过")
            continue
        # ★ 实测列名（2026-10-19 核验）：无 samples_num / af，有 ma_samples / ma_count /
        #   pval_nominal_threshold / min_pval_nominal —— 按实返列名取，勿凭记忆
        cols = ["variant_id", "gene_id", "tss_distance", "maf", "ma_samples", "ma_count",
                "slope", "slope_se", "pval_nominal", "pval_nominal_threshold"]
        for c in cols:
            if c not in ix:
                log(f"[FATAL] {tis} signif 缺列 {c}；实际列={hdr}")
                sys.exit(3)
        n = 0
        for ln in fh:
            p = ln.rstrip("\n").split("\t")
            if len(p) < len(hdr):
                continue
            gid = p[ix["gene_id"]]
            g = gid2sym.get(gid)
            if g is None:
                continue
            vid = p[ix["variant_id"]]
            # variant_id = chr_pos_ref_alt_b38
            parts = vid.split("_")
            if len(parts) < 4:
                continue
            chrom, pos, ref, alt = parts[0].lstrip("chr"), parts[1], parts[2], parts[3]
            pairs.append({
                "tissue": tis, "gene": g, "variant_id": vid,
                "snp_chr": chrom, "snp_pos": pos, "snp_ref": ref, "snp_alt": alt,
                "maf": p[ix["maf"]],
                "pval_threshold": p[ix["pval_nominal_threshold"]],
                "slope": p[ix["slope"]], "slope_se": p[ix["slope_se"]],
                "pval_nominal": p[ix["pval_nominal"]],
                "tss_distance": p[ix["tss_distance"]],
            })
            n += 1
        log(f"[A2] {tis:34s} 显著 SNP-基因对 {n}")

with open(os.path.join(OUT, "gtex_snp_pairs.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(pairs[0].keys()), delimiter="\t")
    w.writeheader()
    w.writerows(pairs)
log(f"[A2] → gtex_snp_pairs.tsv  {len(pairs)} 行")

# ---------- 3. 基因座窗口（用于抓结局） ----------
with open(os.path.join(OUT, "gtex_loci.tsv"), "w", encoding="utf-8", newline="") as f:
    f.write("gene\tchrom\tgene_start\tgene_end\twin_beg\twin_end\n")
    for g in GENES:
        c, s, e = loci[g]
        f.write(f"{g}\t{c}\t{s}\t{e}\t{max(1, s - CIS_HALF)}\t{e + CIS_HALF}\n")
log("[A3] → gtex_loci.tsv")
for g in GENES:
    c, s, e = loci[g]
    log(f"      {g:9s} chr{c}:{s}-{e}  窗口 {max(1,s-CIS_HALF)}-{e+CIS_HALF}")

# ---------- 4. 自检：每基因×组织的显著 SNP 数与 top SNP 是否在显著对内 ----------
log("\n[自检] 每 (tissue,gene) 显著 SNP 数 / top SNP 是否在其中")
import collections
cnt = collections.Counter((r["tissue"], r["gene"]) for r in pairs)
bad = 0
for r in eg_top:
    k = (r["tissue"], r["gene"])
    c = cnt.get(k, 0)
    inper = any(q["tissue"] == r["tissue"] and q["gene"] == r["gene"]
                and q["variant_id"] == r["top_variant_id"] for q in pairs)
    flag = "" if inper else "  ← top SNP 不在显著对内(异常)"
    if not inper:
        bad += 1
    log(f"   {r['tissue']:34s} {r['gene']:9s} n_sig={c:5d} top_in_pairs={inper}{flag}")
log(f"\n[自检结论] top SNP 未命中显著对的组合数 = {bad}（应为 0）")
