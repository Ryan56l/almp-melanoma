#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P2-T13 Step C：SMR + HEIDI + coloc（SAP v1.5 条款 F，含 v1.6 勘误 1）
================================================================================
判据（v1.6 生效）：
  闸 0（前置）  该 (组织, 基因) 必须是显著 eGene（GTEx qval < 0.05），否则登记「该组织不可评估」
  闸 1（决定性）SMR p < 0.05 且 组织内 BH q < 0.05
  闸 2（决定性）HEIDI p > 0.05 且 可评位点数 ≥ 3；不足 3 个 → 判未通过
  闸 3（支持性）coloc PP.H4 > 0.8 且 PP.H3 < 0.5（显著位点子集 → 上偏，仅作正面支持）
  补偿性强化：禁止用低 PP.H4 作「无共定位」的反面证据

实现要点（均须在 Methods 中如实写明）：
  * 结局 `standard_error` 全为 NA → 由 OR 与 p 反推 SE：Z=sign(lnOR)·Φ⁻¹(1−p/2)，SE=|lnOR|/|Z|
  * GTEx `slope` 为 alt 等位基因对表达的效果；结局 effect_allele 与之对齐，必要时翻转符号
  * 回文位点（A/T、C/G）在 HEIDI/coloc 中剔除（无 LD 参考、无法可靠定相），单独计数
  * SMR 基于 top eQTL 位点（单 SNP）；HEIDI 用显著 SNP 对；coloc 用显著 SNP 对交集
  * coloc：V1=sdY²=0.0225、W1/W2 双设定（0.15/0.2 与 0.1/0.1）作敏感性；先验 p1=p2=1e-4、p12=1e-5
"""
import os
import sys
import csv
import math
import collections
import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

TISSUE_CN = {
    "Adipose_Subcutaneous": "皮下脂肪",
    "Skin_Not_Sun_Exposed_Suprapubic": "皮肤非日晒",
    "Skin_Sun_Exposed_Lower_leg": "皮肤日晒",
    "Whole_Blood": "全血",
}
# 组织优先序（条款 F §2B.3 第 1 条）：皮肤（日晒/非日晒）> 皮下脂肪 > 全血
TISSUE_ORDER = {"Skin_Sun_Exposed_Lower_leg": 1, "Skin_Not_Sun_Exposed_Suprapubic": 2,
                "Adipose_Subcutaneous": 3, "Whole_Blood": 4}
GENES = ["LEP", "ADIPOQ", "NAMPT", "FABP4", "CD36", "FADS2", "CDS2", "SERPINE1"]

PALIN = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def log(*a):
    print(*a, flush=True)


# ---------------- 读入 ----------------
top = list(csv.DictReader(open(os.path.join(DATA, "gtex_egenes_top.tsv"), encoding="utf-8"), delimiter="\t"))
pairs = list(csv.DictReader(open(os.path.join(DATA, "gtex_snp_pairs.tsv"), encoding="utf-8"), delimiter="\t"))

outcome = {}
n_out_tot = 0
n_p0 = 0
OUTFILE = sys.argv[1] if len(sys.argv) > 1 else "outcome_loci_v2.tsv"
log(f"[C0] 结局输入 = {OUTFILE}")
with open(os.path.join(DATA, OUTFILE), encoding="utf-8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for ln in f:
        p = ln.rstrip("\n").split("\t")
        if len(p) < len(ix):
            continue
        n_out_tot += 1
        try:
            chrom = p[ix["chromosome"]].lstrip("chr")
            pos = int(p[ix["base_pair_location"]])
            ea = p[ix["effect_allele"]].upper()
            oa = p[ix["other_allele"]].upper()
            orv = float(p[ix["odds_ratio"]])
            pv = float(p[ix["p_value"]])
            eaf = float(p[ix["effect_allele_frequency"]]) if p[ix["effect_allele_frequency"]] not in ("NA", "", "nan") else float("nan")
            nc = float(p[ix["num_cases"]]); nk = float(p[ix["num_controls"]])
        except Exception:
            continue
        if not (orv > 0) or not (0 < pv <= 1):
            n_p0 += 1
            continue
        bz = math.log(orv)
        if pv >= 1 or bz == 0:
            se = float("nan")
        else:
            z = stats.norm.isf(pv / 2.0)
            se = abs(bz) / z if z > 0 else float("nan")
        outcome[(chrom, pos)] = dict(ea=ea, oa=oa, beta=bz, se=se, p=pv, eaf=eaf,
                                     ncase=nc, nctrl=nk, rsid=p[ix["rsid"]])
log(f"[C0] 结局行 {n_out_tot}；不可用（OR/p 异常）{n_p0}；入索引 {len(outcome)}")


def align_pairs(g_alt, g_ref, o):
    """返回 (beta_out_aligned, 说明)。GTEx 效果等位基因为 alt。"""
    ea, oa = o["ea"], o["oa"]
    if g_alt == ea and g_ref == oa:
        return o["beta"], "same"
    if g_alt == oa and g_ref == ea:
        return -o["beta"], "flip"
    # 反向互补（理论上一律正向链，仍做兜底）
    ca, cr = COMP.get(g_alt), COMP.get(g_ref)
    if ca == ea and cr == oa:
        return o["beta"], "same_complement"
    if ca == oa and cr == ea:
        return -o["beta"], "flip_complement"
    return None, "allele_mismatch"


# ---------------- 分组合并 ----------------
bykey = collections.defaultdict(list)
for r in pairs:
    bykey[(r["tissue"], r["gene"])].append(r)

rows_out = []
detail = {}

for r in top:
    tis, gene = r["tissue"], r["gene"]
    qval = float(r["qval"])
    key = (tis, gene)
    pset = bykey.get(key, [])
    rec = dict(tissue=tis, tissue_cn=TISSUE_CN[tis], gene=gene, qval=qval,
               n_sig=len(pset), n_matched=0, n_palin=0, n_mismatch=0,
               top_rsid=r["top_rsid"], top_variant=r["top_variant_id"])
    # --- 闸 0 ---
    if not (qval < 0.05):
        rec.update(status="闸0未过：非显著 eGene → 该组织不可评估",
                   smr_b=float("nan"), smr_se=float("nan"), smr_p=float("nan"),
                   heidi_p=float("nan"), heidi_n=0,
                   pp_h4=float("nan"), pp_h3=float("nan"),
                   gate1="未评", gate2="未评", gate3="未评", verdict=rec.get("verdict", ""))
        rec["verdict"] = "不可评估（非显著 eGene）"
        rows_out.append(rec)
        continue
    # --- 取位点 ---
    snps = []
    for q in pset:
        k = (q["snp_chr"], int(q["snp_pos"]))
        o = outcome.get(k)
        if o is None:
            continue
        if (q["snp_alt"].upper(), q["snp_ref"].upper()) in PALIN:
            rec["n_palin"] += 1
            continue
        b, how = align_pairs(q["snp_alt"].upper(), q["snp_ref"].upper(), o)
        if b is None:
            rec["n_mismatch"] += 1
            continue
        if not np.isfinite(o["se"]) or o["se"] <= 0:
            continue
        snps.append(dict(vid=q["variant_id"], pos=k[1], bzx=float(q["slope"]), sezx=float(q["slope_se"]),
                         pzx=float(q["pval_nominal"]), bzy=b, sezy=o["se"], pzy=o["p"],
                         maf=float(q["maf"]), eaf=o["eaf"], ncase=o["ncase"], nctrl=o["nctrl"],
                         rsid=o["rsid"], how=how))
    rec["n_matched"] = len(snps)
    # --- top SNP ---
    tkey = (r["top_chr"], int(r["top_pos"]))
    o = outcome.get(tkey)
    top_ok = False
    b_smr = se_smr = p_smr = float("nan")
    fail = ""
    if o is None:
        fail = "top SNP 不在结局数据中"
    elif (r["top_alt"].upper(), r["top_ref"].upper()) in PALIN:
        fail = "top SNP 为回文位点（A/T、C/G），无 LD 参考不可定相"
    else:
        bz, how = align_pairs(r["top_alt"].upper(), r["top_ref"].upper(), o)
        bzx, sezx = float(r["top_slope"]), float(r["top_slope_se"])
        if bz is None:
            fail = "等位基因不符（GTEx alt/ref 与结局 effect/other 无法对齐）"
        elif not (np.isfinite(o["se"]) and o["se"] > 0):
            fail = "结局 SE 不可反推（p 或 OR 异常）"
        elif bzx == 0:
            fail = "eQTL 斜率恰为 0"
        else:
            b_smr = bz / bzx
            se_smr = o["se"] / abs(bzx)
            p_smr = 2 * stats.norm.sf(abs(b_smr / se_smr))
            top_ok = True
            rec["top_harmonise"] = how
    rec["top_fail"] = fail
    if top_ok:
        rec["smr_or"] = math.exp(b_smr)
    rec.update(top_ok=top_ok, top_in_outcome=(o is not None),
               smr_b=b_smr, smr_se=se_smr, smr_p=p_smr)
    if not top_ok:
        rec["status"] = "SMR 不可评（top SNP 缺失/回文/等位不符）"
        rec["gate1"] = "未评"; rec["gate2"] = "未评"; rec["gate3"] = "未评"
        rec["heidi_p"] = float("nan"); rec["heidi_n"] = 0
        rec["pp_h4"] = float("nan"); rec["pp_h3"] = float("nan")
        rec["verdict"] = "不可评估（top SNP 不可用）"
        rows_out.append(rec)
        continue
    rec["status"] = "已计算"

    # --- HEIDI ---
    # 阈值：SMR 默认 p_eQTL < 1.57e-3；本项目显著对全部低于该值，仍显式过滤
    heidi_snps = [s for s in snps if s["pzx"] < 1.57e-3 and s["bzx"] != 0]
    # SMR 官方过滤：剔除 b_smr²·se_zx² ≥ se_zy² 的位点（eQTL 估计过粗）
    keep = [s for s in heidi_snps if (b_smr ** 2) * (s["sezx"] ** 2) < (s["sezy"] ** 2)]
    rec["heidi_n_raw"] = len(heidi_snps)
    rec["heidi_n"] = len(keep)
    if len(keep) >= 3:
        w = np.array([1.0 / (s["sezy"] ** 2 + b_smr ** 2 * s["sezx"] ** 2) for s in keep])
        d = np.array([s["bzy"] - b_smr * s["bzx"] for s in keep])
        dbar = float(np.sum(w * d) / np.sum(w))
        T = float(np.sum(w * (d - dbar) ** 2))
        df = len(keep) - 1
        p_h = float(stats.chi2.sf(T, df))
        rec.update(heidi_T=T, heidi_df=df, heidi_p=p_h)
    else:
        rec.update(heidi_T=float("nan"), heidi_df=max(0, len(keep) - 1), heidi_p=float("nan"))

    # --- coloc（显著位点交集；ABF） ---
    W_SET = [(0.15, 0.2), (0.10, 0.10)]
    pp = {}
    for wi, (W1, W2) in enumerate(W_SET):
        V1 = 0.15 ** 2
        # V2：由 se2² = V2/(2f(1-f)Neff) 反推，取中位数
        neff = [4.0 / (1.0 / s["ncase"] + 1.0 / s["nctrl"]) for s in snps if s["ncase"] > 0 and s["nctrl"] > 0]
        Ne = float(np.median(neff)) if neff else np.nan
        v2s = [s["sezy"] ** 2 * 2 * s["maf"] * (1 - s["maf"]) * Ne for s in snps if 0 < s["maf"] < 1]
        V2 = float(np.median(v2s)) if v2s else np.nan
        if not (np.isfinite(V2) and V2 > 0):
            pp[wi] = (float("nan"),) * 4
            continue
        l1, l2 = [], []
        for s in snps:
            z1 = s["bzx"] / s["sezx"]; z2 = s["bzy"] / s["sezy"]
            l1.append(0.5 * math.log(V1 / (V1 + W1)) + 0.5 * z1 * z1 * W1 / (V1 + W1))
            l2.append(0.5 * math.log(V2 / (V2 + W2)) + 0.5 * z2 * z2 * W2 / (V2 + W2))
        l1 = np.array(l1); l2 = np.array(l2)
        s1 = np.logaddexp.reduce(l1)
        s2 = np.logaddexp.reduce(l2)
        # ΣΣ over i≠j
        tot = np.logaddexp.reduce(l1[:, None] + l2[None, :], axis=None)
        ii = np.logaddexp.reduce(l1 + l2)
        s12 = tot + math.log(max(1e-300, 1 - math.exp(min(0.0, ii - tot))))
        p1 = p2 = 1e-4; p12 = 1e-5
        lh0 = math.log(1 - p1) + math.log(1 - p2)
        lh1 = math.log(p1) + math.log(1 - p2) + s1
        lh2 = math.log(1 - p1) + math.log(p2) + s2
        lh3 = math.log(p1) + math.log(p2) + s12
        lh4 = math.log(p12) + tot
        den = np.logaddexp.reduce(np.array([lh0, lh1, lh2, lh3, lh4]))
        pp[wi] = tuple(math.exp(x - den) for x in (lh0, lh1, lh2, lh3, lh4))
    h0, h1, h2, h3, h4 = pp[0]
    rec.update(pp_h1=h1, pp_h2=h2, pp_h3=h3, pp_h4=h4,
               pp_h4_sens=pp[1][4] if len(pp[1]) == 5 else float("nan"),
               coloc_n=len(snps))
    detail[(tis, gene)] = dict(snps=snps, keep=keep)
    rows_out.append(rec)

# ---------------- BH 校正（组织内跨 8 基因） ----------------
by_tis = collections.defaultdict(list)
for rec in rows_out:
    if np.isfinite(rec.get("smr_p", float("nan"))):
        by_tis[rec["tissue"]].append(rec)
for tis, lst in by_tis.items():
    ps = np.array([r["smr_p"] for r in lst])
    order = np.argsort(ps)
    m = len(ps)
    q = np.empty(m)
    prev = 1.0
    for rank in range(m - 1, -1, -1):
        i = order[rank]
        val = ps[i] * m / (rank + 1)
        prev = min(prev, val)
        q[i] = prev
    for r, qq in zip(lst, q):
        r["smr_q"] = float(qq)

# ---------------- 判定 ----------------
for r in rows_out:
    r.setdefault("smr_q", float("nan"))
    if r["status"].startswith("已计算"):
        g1 = np.isfinite(r["smr_p"]) and r["smr_p"] < 0.05 and np.isfinite(r["smr_q"]) and r["smr_q"] < 0.05
        g2 = np.isfinite(r["heidi_p"]) and r["heidi_p"] > 0.05 and r["heidi_n"] >= 3
        g3 = np.isfinite(r["pp_h4"]) and r["pp_h4"] > 0.8 and r["pp_h3"] < 0.5
        r["gate1"] = "过" if g1 else ("未过" if np.isfinite(r["smr_p"]) else "未评")
        if not np.isfinite(r["heidi_p"]):
            r["gate2"] = "未过（位点<3，不可评即未通过）"
        else:
            r["gate2"] = "过" if g2 else "未过"
        r["gate3"] = ("支持" if g3 else "未支持") + "（支持性）"
        if g1 and g2 and g3:
            r["verdict"] = "两闸过＋共定位支持"
        elif g1 and g2:
            r["verdict"] = "两闸过（共定位未支持）"
        elif g1:
            r["verdict"] = "仅 SMR 过（HEIDI 否决/不可评）"
        else:
            r["verdict"] = "两闸均未过"
    r.setdefault("gate1", "未评")
    r.setdefault("gate2", "未评")
    r.setdefault("gate3", "未评")

# ---------------- 输出 ----------------
cols = ["tissue", "tissue_cn", "gene", "qval", "n_sig", "n_matched", "n_palin", "n_mismatch",
        "top_rsid", "top_variant", "top_ok", "top_in_outcome", "top_fail",
        "smr_b", "smr_se", "smr_or", "smr_p", "smr_q",
        "heidi_n_raw", "heidi_n", "heidi_T", "heidi_df", "heidi_p",
        "pp_h1", "pp_h2", "pp_h3", "pp_h4", "pp_h4_sens", "coloc_n",
        "gate1", "gate2", "gate3", "verdict", "status"]
rows_out.sort(key=lambda r: (TISSUE_ORDER[r["tissue"]], GENES.index(r["gene"])))
with open(os.path.join(RES, "smr_heidi_coloc_matrix.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for r in rows_out:
        rr = {}
        for c in cols:
            v = r.get(c, "")
            if isinstance(v, float):
                rr[c] = "NA" if not np.isfinite(v) else f"{v:.6g}"
            else:
                rr[c] = v
        w.writerow(rr)
log(f"[C1] → results/smr_heidi_coloc_matrix.tsv  {len(rows_out)} 行")

# ---------------- 控制台摘要 ----------------
log("\n=== 判定摘要（32 组合）===")
log(f"{'组织':10s}{'基因':10s}{'qval':>10s}{'n':>5s}{'SMRp':>10s}{'q':>9s}{'HEIDIp':>10s}{'nH':>5s}{'PPH4':>8s}{'判定':>26s}")
for r in rows_out:
    f_ = lambda x, fm="%.3g": ("NA" if not np.isfinite(x) else fm % x)
    log(f"{r['tissue_cn']:10s}{r['gene']:10s}{r['qval']:10.2e}{r['n_matched']:5d}"
        f"{f_(r['smr_p'],'%.4f'):>10s}{f_(r.get('smr_q',float('nan')),'%.3f'):>9s}"
        f"{f_(r['heidi_p'],'%.3f'):>10s}{r['heidi_n']:5d}{f_(r['pp_h4'],'%.2f'):>8s}{r['verdict']:>26s}")
# ★ 共定位区分度诊断：PP.H4 是否随 SMR 强度变化
diag = [r for r in rows_out if np.isfinite(r.get("pp_h4", float("nan"))) and np.isfinite(r.get("smr_p", float("nan")))]
if len(diag) >= 5:
    x = np.array([-math.log10(max(r["smr_p"], 1e-300)) for r in diag])
    y = np.array([r["pp_h4"] for r in diag])
    rho = stats.spearmanr(x, y)
    log("\n=== 共定位区分度诊断（PP.H4 与 SMR 显著性是否同向）===")
    log(f"  n={len(diag)}  Spearman rho={rho.statistic:.3f}  p={rho.pvalue:.3f}")
    hi = [r for r in diag if r["pp_h4"] > 0.8]
    log(f"  PP.H4>0.8 的组合数={len(hi)}（其中 SMR p<0.05 者={sum(1 for r in hi if r['smr_p']<0.05)}）")
    with open(os.path.join(RES, "coloc_informativeness.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["tissue_cn", "gene", "smr_p", "pp_h4", "pp_h3", "pp_h4_gt08"])
        for r in sorted(diag, key=lambda z: z["smr_p"]):
            w.writerow([r["tissue_cn"], r["gene"], f"{r['smr_p']:.4g}", f"{r['pp_h4']:.4f}",
                        f"{r['pp_h3']:.4g}", r["pp_h4"] > 0.8])
    log("  → results/coloc_informativeness.tsv")

n_pass = sum(1 for r in rows_out if r["verdict"].startswith("两闸过"))
log(f"\n两闸通过者：{n_pass} / 32")
with open(os.path.join(RES, "run_log.txt"), "w", encoding="utf-8") as f:
    f.write(f"n_outcome_variants={len(outcome)}\n")
    f.write(f"two_gate_pass={n_pass}\n")
