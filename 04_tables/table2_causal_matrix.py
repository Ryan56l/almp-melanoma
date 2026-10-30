# -*- coding: utf-8 -*-
"""表 2：32 组组织×基因 cis-SMR / HEIDI / coloc 判定矩阵（离线，不改动任何分析结果）"""
import csv, io, os, collections

ROOT = "C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36"
SRC = os.path.join(ROOT, "output/P2-T13_GTEx定量SMR与共定位v1/results/smr_heidi_coloc_matrix.tsv")
OUT_MD = os.path.join(ROOT, "output/表2_因果层32组合判定矩阵_v1.md")
OUT_CSV = os.path.join(ROOT, "output/表2_因果层32组合判定矩阵_v1.csv")


def cn(s):
    """表格文本字段：半角逗号/分号 -> 中文标点（用户约定）"""
    if s is None:
        return ""
    return s.replace(",", "，").replace(";", "；")


def num(v, nd=4):
    if v in ("", "NA", None):
        return "—"
    try:
        return ("%%.%df" % nd) % float(v)
    except ValueError:
        return cn(v)


rows = list(csv.DictReader(io.open(SRC, encoding="utf-8"), delimiter="\t"))
assert len(rows) == 32, len(rows)

TISSUE_ORDER = ["皮肤日晒", "皮肤非日晒", "皮下脂肪", "全血"]
GENE_ORDER = []
for r in rows:
    if r["gene"] not in GENE_ORDER:
        GENE_ORDER.append(r["gene"])

# ---------- 派生统计 ----------
VC = collections.Counter(r["verdict"] for r in rows)
G1 = collections.Counter(r["gate1"] for r in rows)
G2 = collections.Counter(r["gate2"] for r in rows)
G3 = collections.Counter(r["gate3"] for r in rows)

evaluable = [r for r in rows if r["verdict"] == "两闸均未过"]
# 逐组织计数
per_tissue = collections.OrderedDict()
for t in TISSUE_ORDER:
    sub = [r for r in rows if r["tissue_cn"] == t]
    per_tissue[t] = {
        "n": len(sub),
        "noneg": sum(1 for r in sub if r["verdict"] == "不可评估（非显著 eGene）"),
        "notop": sum(1 for r in sub if r["verdict"] == "不可评估（top SNP 不可用）"),
        "fail": sum(1 for r in sub if r["verdict"] == "两闸均未过"),
    }

# ---------- 写 markdown ----------
L = []
L.append("# 表 2　因果层 32 组「组织 × 基因」cis-SMR／HEIDI／共定位判定矩阵（v1）")
L.append("")
L.append("> **日期**：2026-10-22　|　**数据源**：`output/P2-T13_GTEx定量SMR与共定位v1/results/smr_heidi_coloc_matrix.tsv`（32 行，逐行直读，未重算）")
L.append("> **判据**：SAP v1.6 条款 F——**决定性判据＝SMR＋HEIDI 两闸**；coloc（PP.H4）为**支持性**，不作判定依据。")
L.append("> **重要**：`smr_p` 未做多重校正前的名义 p；闸 1 用 **BH q<0.05**（`smr_q`）。两闸通过数 **0/32**。")
L.append("")

L.append("## A　判定分布")
L.append("")
L.append("| 判定值 | 组合数 | 含义 |")
L.append("|---|---|---|")
L.append("| 两闸均未过 | **%d** | eGene 显著、可算 SMR，但闸 1（q<0.05）未过 → **不成立** |" % VC["两闸均未过"])
L.append("| 不可评估（非显著 eGene） | **%d** | 该组织该基因在 GTEx 中 eGene q≥0.05，无可用 cis 工具 → **不判为阴性** |" % VC["不可评估（非显著 eGene）"])
L.append("| 不可评估（top SNP 不可用） | **%d** | top eSNP 缺失或无法与结局侧匹配 → **不判为阴性** |" % VC["不可评估（top SNP 不可用）"])
L.append("| 通过两闸 | **0** | — |")
L.append("")
L.append("**闸门分解**：")
L.append("")
L.append("- 闸 1（SMR q<0.05）：未评 %d｜未过 %d｜过 **0**" % (G1["未评"], G1["未过"]))
L.append("- 闸 2（HEIDI p>0.05）：未评 %d｜过 %d｜未过（位点<3，不可评即未通过）**%d**" % (G2["未评"], G2["过"], G2["未过（位点<3，不可评即未通过）"]))
L.append("- 闸 3（coloc PP.H4>0.8，**支持性**）：%s" % "；".join("%s %d" % (k, v) for k, v in G3.items()))
L.append("")
L.append("**逐组织计数**：")
L.append("")
L.append("| 组织 | 组合数 | 非显著 eGene | top SNP 不可用 | 两闸均未过 |")
L.append("|---|---|---|---|---|")
for t, d in per_tissue.items():
    L.append("| %s | %d | %d | %d | %d |" % (t, d["n"], d["noneg"], d["notop"], d["fail"]))
L.append("")
L.append("> 读法：**「不可评估」不是阴性结果**。论文正文与图 7b 一律不得把 17 组不可评估计为「未发现关联」；")
L.append("> 该 17 组只能表述为「本设计下无判定力」。这是 Methods 与 Discussion §4.3 的口径底线。")
L.append("")

L.append("## B　逐组合明细（32 行）")
L.append("")
L.append("| # | 组织 | 基因 | eGene q | 显著位点 | 已配对 | 等位基因同向 | SMR b | OR | SMR p | SMR q | HEIDI p | PP.H4 | 闸1 | 闸2 | 判定 |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
i = 0
for t in TISSUE_ORDER:
    for g in GENE_ORDER:
        for r in rows:
            if r["tissue_cn"] != t or r["gene"] != g:
                continue
            i += 1
            L.append("| %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                i, r["tissue_cn"], r["gene"],
                num(r["qval"], 4), r["n_sig"] or "0", r["n_matched"] or "0", r["n_palin"] or "0",
                num(r["smr_b"], 3), num(r["smr_or"], 3), num(r["smr_p"]), num(r["smr_q"]),
                num(r["heidi_p"], 4), num(r["pp_h4"], 3),
                cn(r["gate1"]), cn(r["gate2"]), cn(r["verdict"]),
            ))
assert i == 32, i
L.append("")

L.append("## C　两个必须在正文中出现的说明")
L.append("")
L.append("1. **SMR 的 z 恒等于结局侧 z**。eQTL 斜率 b 在 SMR 统计量中相消（`z_SMR = b_xy / se ≈ z_outcome`），")
L.append("   因此本表的 SMR p 与**结局 GWAS 自身的 p** 同源，**不构成**「基因表达 → 黑色素瘤」的独立证据增量。")
L.append("   这一点已在 Discussion §4.4 明确披露，表注不得省略。")
L.append("2. **coloc 无区分力**。ρ(PP.H4, −log₁₀ SMR p) = **0.122（p=0.665）**，即 PP.H4 与「SMR 有多显著」近乎正交，")
L.append("   所以 PP.H4 不能用来补救闸 1 的失败。CDS2 即典型案例：PP.H4 0.816 看似支持，q=0.110 仍未过闸 1 → **降为普通成员**。")
L.append("")
L.append("## D　版本")
L.append("")
L.append("| 版本 | 日期 | 变更 |")
L.append("|---|---|---|")
L.append("| v1 | 2026-10-22 | 首版。自 `smr_heidi_coloc_matrix.tsv` 直读生成，未重跑分析、未改动任何数值。表内文本字段半角逗号/分号按工作区约定替换为中文标点。 |")
L.append("")

io.open(OUT_MD, "w", encoding="utf-8").write("\n".join(L))

# ---------- CSV ----------
with io.open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["序号", "组织", "基因", "eGene q", "显著位点数", "已配对位点", "回文位点（已剔除）",
                "SMR b", "SMR OR", "SMR p", "SMR q(BH)", "HEIDI p", "PP.H4",
                "闸1(SMR q<0.05)", "闸2(HEIDI p>0.05)", "判定"])
    i = 0
    for t in TISSUE_ORDER:
        for g in GENE_ORDER:
            for r in rows:
                if r["tissue_cn"] != t or r["gene"] != g:
                    continue
                i += 1
                w.writerow([i, cn(r["tissue_cn"]), cn(r["gene"]), num(r["qval"], 4), r["n_sig"] or "0",
                            r["n_matched"] or "0", r["n_palin"] or "0", num(r["smr_b"], 3),
                            num(r["smr_or"], 3), num(r["smr_p"]), num(r["smr_q"]),
                            num(r["heidi_p"], 4), num(r["pp_h4"], 3),
                            cn(r["gate1"]), cn(r["gate2"]), cn(r["verdict"])])

print("rows:", i)
print("verdict:", dict(VC))
print("md bytes:", os.path.getsize(OUT_MD), "| csv bytes:", os.path.getsize(OUT_CSV))
