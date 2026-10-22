#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
P2-T13 Step B：抓取结局 GCST90475578 在 8 基因座窗口内的全部变异
================================================================================
窗口来源：data/grab_windows.tsv（= 该基因全部 GTEx 显著位点 + top SNP 的最小跨度 ±50 kb）
★ 为何不用整个 ±1 Mb：GTEx v8 归档无 allpairs，coloc 本来就只能在显著位点上跑，
  故显著位点跨度之外的数据对本层分析无用；且远端抓取每 64 KB 块约 4 s，窗口越大越慢。
可续跑：data/outcome_loci.tsv.progress 记录已完成窗口（gene:chrom:beg-end）。
"""
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36/.workbuddy/tmp")
from grab_90475578 import _stats  # noqa: E402
sys.path.insert(0, HERE)
from bgzf_reader_v2 import grab_region_v2 as grab_region  # noqa: E402  ★ carry-over 修正版

HDR = ("chromosome\tbase_pair_location\teffect_allele\tother_allele\todds_ratio\t"
       "standard_error\teffect_allele_frequency\tp_value\trsid\tci_upper\tci_lower\talt\t"
       "n\tcase_af\tnum_cases\tcontrol_af\tnum_controls\tr2\tq_pval\ti2\tdirection\t"
       "hm_coordinate_conversion\thm_code\tvariant_id")

win_file = os.path.join(ROOT, "data", "grab_windows.tsv")
out_file = os.path.join(ROOT, "data", "outcome_loci_v2.tsv")
prog_file = out_file + ".progress"

wins = []
with open(win_file, encoding="utf-8") as f:
    next(f)
    for ln in f:
        g, c, b, e, _sp = ln.rstrip("\n").split("\t")
        wins.append((g, c, int(b), int(e)))

done = set()
if os.path.exists(prog_file):
    with open(prog_file, encoding="utf-8") as f:
        done = {ln.strip() for ln in f if ln.strip()}

new = not os.path.exists(out_file)
fo = open(out_file, "a", encoding="utf-8")
if new:
    fo.write(HDR + "\n")
    fo.flush()

for g, c, b, e in wins:
    key = f"{g}:{c}:{b}-{e}"
    if key in done:
        print(f"[skip] {key}", flush=True)
        continue
    t0 = time.time()
    try:
        rows = grab_region(int(c), b, e, verbose=True)
    except Exception as ex:
        print(f"[ERR] {key}: {type(ex).__name__} {ex}", flush=True)
        continue
    n_bad = 0
    for r in rows:
        if len(r.split("\t")) != len(HDR.split("\t")):
            n_bad += 1
            continue
        fo.write(r + "\n")
    fo.flush()
    with open(prog_file, "a", encoding="utf-8") as fp:
        fp.write(key + "\n")
    print(f"[OK] {key}  rows={len(rows)} bad={n_bad}  {time.time()-t0:.0f}s", flush=True)

fo.close()
print(f"\n>>> done. total requests={_stats['requests']}  bytes={_stats['bytes']/1e6:.1f} MB", flush=True)
