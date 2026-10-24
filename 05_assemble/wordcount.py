# -*- coding: utf-8 -*-
"""权威词数机扫（Word 口径近似）v2
令牌化规则：按空白切分；剔除纯破折号/标点令牌（— – -）；连字符复合词算 1 词。
只计正文散文行：剔除 # 标题行、> 引用块、| 表格行、--- 分隔线、** 小标题行、``` 代码块。
"""
import io, os, re, sys

D = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output"

DASHES = set("—–-−")


def is_prose(line):
    s = line.strip()
    if not s:
        return False
    if s.startswith("#") or s.startswith(">") or s.startswith("|"):
        return False
    if set(s) <= set("-—=*_ "):
        return False
    if s.startswith("**") and s.endswith("**"):
        return False
    if s.startswith("```") or s.startswith("    "):
        return False
    return True


def count(text):
    """strict：全部空白令牌，含独立破折号 —— 与系统 `wc -w` 等价，**权威口径**"""
    return len(text.split())


def count_adj(text):
    """参考口径：剔除纯破折号/标点令牌（不用于对外报数）"""
    return len([t for t in text.split()
                if not all((ch in DASHES) or ch in ".,;:()[]" for ch in t)])


def sections(path, cut_markers=()):
    """返回 [(heading, [段落文本, ...]), ...]"""
    raw = io.open(path, encoding="utf-8").read()
    # 截断到第一个 cut marker（自查表/变更日志等不计入）
    for m in cut_markers:
        i = raw.find(m)
        if i > 0:
            raw = raw[:i]
    out, cur, buf = [], "_preamble", []
    for line in raw.splitlines():
        if line.startswith("#"):
            if buf:
                out.append((cur, buf))
            cur, buf = line.lstrip("# ").strip(), []
            continue
        if is_prose(line):
            buf.append(line.strip())
        else:
            if buf:
                out.append((cur, buf))
                buf = []
                cur = cur  # keep
    if buf:
        out.append((cur, buf))
    # 合并同名 heading 的连续块
    merged = []
    for h, b in out:
        if merged and merged[-1][0] == h:
            merged[-1][1].extend(b)
        else:
            merged.append((h, list(b)))
    return merged


def report(tag, path, keep_prefix, cut_markers=()):
    print("=" * 70)
    print(tag, os.path.basename(path))
    total = 0
    for h, blocks in sections(path, cut_markers):
        if not h.startswith(keep_prefix):
            continue
        n = sum(count(b) for b in blocks)
        if n == 0:
            continue
        per = " / ".join(str(count(b)) for b in blocks)
        print(f"  {h[:70]:<72} {n:>5}   [{per}]")
        total += n
    print(f"  {'-- 小计':<72} {total:>5}")
    return total


A = os.path.join(D, "论文初稿_Abstract与Introduction_v1.md")
R = os.path.join(D, "论文初稿_Results_v1.md")
DISC = os.path.join(D, "论文初稿_Discussion_v1.md")
M = os.path.join(D, "论文初稿_Methods_v1.md")

abs_n = report("ABSTRACT", A, "二、Abstract")
intro_n = report("INTRODUCTION", A, "三、Introduction")
res_n = report("RESULTS", R, "3.1") + 0  # placeholder overridden below
# Results: 3.1–3.5 五节，逐节取
print("=" * 70)
print("RESULTS")
res_n = 0
for h, blocks in sections(R, ("## 附",)):
    if not re.match(r"^3\.[1-6] ", h):
        continue
    n = sum(count(b) for b in blocks)
    print(f"  {h[:70]:<72} {n:>5}   [{' / '.join(str(count(b)) for b in blocks)}]")
    res_n += n
print(f"  {'-- 小计':<72} {res_n:>5}")

print("=" * 70)
print("DISCUSSION")
dis_n = 0
for h, blocks in sections(DISC, ("## 4.8",)):
    if not re.match(r"^4\.[1-7] ", h):
        continue
    n = sum(count(b) for b in blocks)
    print(f"  {h[:70]:<72} {n:>5}   [{' / '.join(str(count(b)) for b in blocks)}]")
    dis_n += n
print(f"  {'-- 小计':<72} {dis_n:>5}")

print("=" * 70)
print("METHODS")
met_n = 0
for h, blocks in sections(M, ("## 9.",)):
    if not re.match(r"^2\.[1-8] ", h):
        continue
    n = sum(count(b) for b in blocks)
    print(f"  {h[:70]:<72} {n:>5}   [{' / '.join(str(count(b)) for b in blocks)}]")
    met_n += n
print(f"  {'-- 小计':<72} {met_n:>5}")

imrad = intro_n + res_n + dis_n
print("=" * 70)
print(f"Abstract            {abs_n:>5}   (上限 300)")
print(f"IMRaD 三节合计      {imrad:>5}   (预算 5,000 → 差 {imrad - 5000:+d})")
print(f"Methods（单列）     {met_n:>5}")
print(f"全部含摘要与方法    {abs_n + imrad + met_n:>5}")
print("=" * 70)
print(f"段落数：Intro {len([b for h, bs in sections(A, ('## 四、',)) if h.startswith('三、Introduction') for b in bs])} 段")
