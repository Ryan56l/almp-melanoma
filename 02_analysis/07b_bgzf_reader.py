#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
BGZF 区域读取器 v2（跨块 carry-over 修正版）
================================================================================
★ 缺陷背景（2026-10-19 自检发现）
  旧实现（grab_90475578.grab_region）逐 BGZF 块解压后各自 split('\\n')：
  块尾的半行被丢弃，下一块块首的续行也被丢弃 → **每个块边界净丢 1 行真实变异**。
  实测：LEP 窗口 19 块丢 14 行、chr20 窗口 19 块丢 3 行（`bad` 计数即此）。
  宽窗与紧窗抓取结果因此不一致（紧窗 1310 vs 宽窗 1308），违反永久自检铁律。

★ 修正
  1. 首块：丢弃至第一个 '\\n'（chunk 起点可能落在行中间）。
  2. 其后各块：把「最后一个未终止换行的片段」保留为 buf，与下一块解压结果拼接。
  3. 按行解析，仅在读到 p > end 时停止（不因块边界提前截断）。

一致性要求：同一区间「宽窗」与「紧窗」两次抓取的坐标集合必须**完全相同**。
"""
import os
import sys

sys.path.insert(0, r"C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36/.workbuddy/tmp")
from tbi import read_tbi                      # noqa: E402
from grab_90475578 import fetch_range, bgzf_bsize, multi_member_decompress, TBI  # noqa: E402


def grab_region_v2(chrom, beg, end, max_blocks=20000, verbose=False):
    """抓取 1-based 闭区间 [beg,end] 的全部行；跨块行不丢失。"""
    tbi = read_tbi(TBI)
    cs = str(chrom)
    if cs not in tbi["names"]:
        return []
    ref = tbi["refs"][tbi["names"].index(cs)]
    win = (beg - 1) >> 14
    if win >= len(ref["intvs"]):
        return []
    v = ref["intvs"][win]
    if v == 0:
        k = win
        while k > 0 and ref["intvs"][k] == 0:
            k -= 1
        v = ref["intvs"][k]
    if v == 0:
        return []
    pos = v >> 16

    out = {}
    buf = b""
    first = True
    blocks = 0
    stop = False
    while blocks < max_blocks and not stop:
        head = fetch_range(pos, pos + 65535)
        bs = bgzf_bsize(head, 0)
        if bs is None:
            break
        bsize = bs + 1
        if bsize > len(head):
            head = fetch_range(pos, pos + bsize - 1)
        txt = multi_member_decompress(head[:bsize])
        if not txt:
            break
        data = buf + txt
        if first:
            i = data.find(b"\n")
            data = data[i + 1:] if i >= 0 else b""
            first = False
        last_nl = data.rfind(b"\n")
        if last_nl < 0:
            buf = data
            pos += bsize
            blocks += 1
            continue
        buf = data[last_nl + 1:]
        for line in data[:last_nl].split(b"\n"):
            if not line:
                continue
            f = line.decode("utf-8", "replace").rstrip("\r").split("\t")
            if len(f) < 3 or f[0].lstrip("chr") != cs:
                continue
            try:
                p = int(f[1])
            except ValueError:
                continue
            if p < beg:
                continue
            if p > end:
                stop = True
                break
            out.setdefault(p, []).append(line.decode("utf-8", "replace").rstrip("\r"))
        pos += bsize
        blocks += 1
        if verbose and blocks % 10 == 0:
            print(f"    ... {blocks} blocks, {sum(len(v) for v in out.values())} rows", flush=True)

    rows = []
    for p in sorted(out):
        rows.extend(out[p])
    if verbose:
        print(f"    done v2: {blocks} blocks, {len(rows)} rows", flush=True)
    return rows
