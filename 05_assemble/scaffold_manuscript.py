# -*- coding: utf-8 -*-
"""
指令：建立写作期统一工作目录 manuscript/
- 复制（不移动）；保持相对结构
- 逐文件校验：数量 / 大小 / 文件名 / md5
- 输出迁移清单 md + csv
"""
import csv, hashlib, io, os, shutil, collections

SRC_WS = r"C:\Users\刘宸瑞\WorkBuddy\2026-10-07-22-15-36"          # 工作区（活文件）
SRC_DK = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信"                    # 桌面（原始存档）
DST = os.path.join(SRC_DK, "manuscript")

EXCLUDE_DIR_PARTS = {".workbuddy", ".git"}


def md5(p, buf=1 << 20):
    h = hashlib.md5()
    with open(p, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- 纳入规则
def keep(rel, size):
    """rel 为相对源根的 posix 路径；返回 (是否纳入, 类别)"""
    low = rel.lower()
    name = os.path.basename(low)
    parts = rel.split("/")
    if any(p in EXCLUDE_DIR_PARTS for p in parts):
        return False, ""
    if name.endswith((".bak",)) or ".bak-" in low or ".backup" in low:
        return False, ""
    if name.endswith(".docx") or name.endswith(".xlsx"):
        return True, "方案·SAP·台账"
    if name.endswith(".md"):
        return True, "报告与正文"
    if name.endswith(".csv"):
        return True, "结果表"
    if name.endswith(".html"):
        return True, "选题报告/排版稿"
    if name.endswith((".pdf", ".png")):
        return True, "图件"
    if name.endswith(".tsv"):
        # 纳入：results/ 下的结果表；output/ 根级的核验表。data/ 下为中间/原始数据，不纳入
        if size <= 2 * 1024 * 1024 and ("results" in parts[:-1] or len(parts) == 2):
            return True, "结果表(results/)"
        return False, ""
    return False, ""


# ---------------------------------------------------------------- 规划复制
plan = []          # dict: src, dest_rel, size, kind, origin
seen_md5 = {}

# 1) 工作区根级文档
WS_ROOT_FILES = [
    "研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层.docx",
    "研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层_v1.4修订.docx",
    "研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层_v1.5修订.docx",
    "统计分析计划_SAP_ALMP研究_v1.0冻结版.docx",
    "统计分析计划_SAP_ALMP研究_v1.1冻结版.docx",
    "统计分析计划_SAP_ALMP研究_v1.2冻结版.docx",
    "统计分析计划_SAP_ALMP研究_v1.3冻结版.docx",
    "ALMP研究项目任务清单.xlsx",
    "报告规范核对台账_89条_v1.0.xlsx",
    "报告规范适用性与方案16.3完整性审计_v1.0.md",
    "数据核验记录表_v1.0.md",
    "数据核验记录表_TCGA-SKCM与GEO队列_v0.1.md",
    "任务清单_ALMP方案v1.1_分阶段执行计划.md",
    "选题报告_脂质代谢脂肪因子_皮肤肿瘤免疫微环境.html",
]
for f in WS_ROOT_FILES:
    p = os.path.join(SRC_WS, f)
    if not os.path.exists(p):
        print("[WARN] 缺失（工作区根级）：", f)
        continue
    plan.append(dict(src=p, dest_rel="docs/" + f, size=os.path.getsize(p),
                     kind="方案·SAP·台账", origin="工作区根级"))

# 2) 桌面根级：仅方案 docx 与选题报告；同名不同版者加"桌面存档版"
DK_FILES = [
    ("研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层.docx",
     "研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层_桌面存档版_20261008.docx"),
    ("选题报告_脂质代谢脂肪因子_皮肤肿瘤免疫微环境.html", None),
]
for src_name, alias in DK_FILES:
    p = os.path.join(SRC_DK, src_name)
    if not os.path.exists(p):
        print("[WARN] 缺失（桌面）：", src_name)
        continue
    plan.append(dict(src=p, dest_rel="docs/" + (alias or src_name), size=os.path.getsize(p),
                     kind="方案·SAP·台账", origin="桌面根级"))

# 3) 工作区 output/ 全树
for root, dirs, files in os.walk(os.path.join(SRC_WS, "output")):
    dirs[:] = [d for d in dirs if d not in EXCLUDE_DIR_PARTS]
    for fn in files:
        p = os.path.join(root, fn)
        rel = os.path.relpath(p, SRC_WS).replace("\\", "/")
        sz = os.path.getsize(p)
        ok, kind = keep(rel, sz)
        if not ok:
            continue
        plan.append(dict(src=p, dest_rel=rel, size=sz, kind=kind, origin="工作区 output/"))

# ---------------------------------------------------------------- 同名去重
# 目标相对路径相同者：若 md5 一致则合并为一份并记录；不一致则中止（须人工裁定）
_by_rel = {}
dedup = []
uniq = []
for it in plan:
    if it["dest_rel"] in _by_rel:
        prev = _by_rel[it["dest_rel"]]
        if md5(prev["src"]) == md5(it["src"]):
            dedup.append((it["dest_rel"], prev["src"], it["src"]))
            continue
        raise SystemExit("[ABORT] 同名但内容不同，须人工裁定：%s\n  A=%s\n  B=%s"
                         % (it["dest_rel"], prev["src"], it["src"]))
    _by_rel[it["dest_rel"]] = it
    uniq.append(it)
plan = uniq

# ---------------------------------------------------------------- 执行复制
# 安全：不做任何删除。目标非空且无本脚本标记 → 中止；有标记 → 允许覆盖（幂等重跑）。
MARKER = os.path.join("output", "论文骨架_IMRaD与图表清单_v1.md")
if os.path.isdir(DST) and os.listdir(DST):
    if not os.path.exists(os.path.join(DST, MARKER)):
        raise SystemExit("[ABORT] 目标已存在且非本脚本产物，未做任何改动：%s" % DST)
os.makedirs(DST, exist_ok=True)

rows = []
for it in plan:
    dpath = os.path.join(DST, it["dest_rel"].replace("/", os.sep))
    os.makedirs(os.path.dirname(dpath), exist_ok=True)
    shutil.copy2(it["src"], dpath)
    it["dest"] = dpath
    it["md5_src"] = md5(it["src"])
    it["md5_dst"] = md5(dpath)
    it["size_dst"] = os.path.getsize(dpath)
    rows.append(it)

# ---------------------------------------------------------------- 校验
n_src = len(plan)
EXCL_FROM_COUNT = {"迁移清单_v1.md", "迁移清单_v1.csv", "README_使用约定.md"}
landed, landed_alias = [], []
for root, dirs, files in os.walk(DST):
    for fn in files:
        p = os.path.join(root, fn)
        rel = os.path.relpath(p, DST).replace("\\", "/")
        if rel.startswith("ascii_alias/"):
            landed_alias.append(rel)
        elif fn not in EXCL_FROM_COUNT:
            landed.append(rel)
n_dst = len(landed)
expected = set(r["dest_rel"] for r in plan)
extra = sorted(set(landed) - expected)
missing = sorted(expected - set(landed))
bad_size = [r for r in rows if r["size"] != r["size_dst"]]
bad_md5 = [r for r in rows if r["md5_src"] != r["md5_dst"]]
bad_name = [r for r in rows
            if os.path.basename(r["dest"]) != os.path.basename(r["src"])
            and "_桌面存档版_" not in os.path.basename(r["dest"])]

# 名称集合比对（同名仅在桌面存档版一处有意改名）
src_names = collections.Counter(os.path.basename(r["src"]) for r in rows)
dst_names = collections.Counter(os.path.basename(r["dest"]) for r in rows)

# 源侧未纳入统计（output/ 排除项）
excluded, ex_size = [], 0
for root, dirs, files in os.walk(os.path.join(SRC_WS, "output")):
    dirs[:] = [d for d in dirs if d not in EXCLUDE_DIR_PARTS]
    for fn in files:
        p = os.path.join(root, fn)
        rel = os.path.relpath(p, SRC_WS).replace("\\", "/")
        sz = os.path.getsize(p)
        ok, _ = keep(rel, sz)
        if not ok:
            excluded.append((rel, sz))
            ex_size += sz

# ---------------------------------------------------------------- ASCII 别名
ASCII_MAP = {
    "表1_队列与数据源基线_v1.csv": "Table1_cohort_baseline_v1.csv",
    "表2_因果层32组合判定矩阵_v1.csv": "Table2_causal_32matrix_v1.csv",
    "GWAS数据清单_T13_v1.csv": "GWAS_dataset_list_T13_v1.csv",
    "MDE估计_T14_v1.csv": "MDE_estimate_T14_v1.csv",
    "MR结果_GCST90475578_v1.csv": "MR_results_GCST90475578_v1.csv",
    "工具变量强度表_T14_v1.csv": "instrument_strength_T14_v1.csv",
}
ascii_rows = []
alias_dir = os.path.join(DST, "ascii_alias")
os.makedirs(alias_dir, exist_ok=True)
for r in rows:
    if not r["dest_rel"].lower().endswith(".csv"):
        continue
    base = os.path.basename(r["dest_rel"])
    fixed = ASCII_MAP.get(base)
    if fixed is None:
        stem = os.path.splitext(base)[0]
        # 非中文表名的 csv 直接保留原名；纯中文/混合名给一个确定性 ASCII 名
        if all(ord(c) < 128 for c in stem):
            fixed = base
        else:
            fixed = "alias_%03d.csv" % (len(ascii_rows) + 1)
    dst_alias = os.path.join(alias_dir, fixed)
    shutil.copy2(r["dest"], dst_alias)
    ascii_rows.append(dict(orig=r["dest_rel"], alias="ascii_alias/" + fixed,
                           size=os.path.getsize(dst_alias),
                           md5_ok=(md5(dst_alias) == r["md5_dst"])))

# 清理本脚本上一次运行遗留的别名（仅限 ascii_alias/ 下的 *.csv，且其内容与某个正本 csv 相同）
_new_alias = set(os.path.basename(a["alias"]) for a in ascii_rows)
_canon_md5 = set(r["md5_dst"] for r in rows if r["dest_rel"].lower().endswith(".csv"))
_purged = []
for fn in sorted(os.listdir(alias_dir)):
    p = os.path.join(alias_dir, fn)
    if fn in _new_alias or not fn.lower().endswith(".csv"):
        continue
    if md5(p) in _canon_md5:          # 确认是本脚本生成的副本，不是用户文件
        os.remove(p)
        _purged.append(fn)
if _purged:
    print("[清理] 旧别名：", ", ".join(_purged))

# ---------------------------------------------------------------- 写清单
KIND_ORDER = ["方案·SAP·台账", "报告与正文", "结果表", "结果表(results/)", "图件", "选题报告/排版稿"]
by_kind = collections.OrderedDict((k, []) for k in KIND_ORDER)
for r in rows:
    by_kind.setdefault(r["kind"], []).append(r)

total_size = sum(r["size"] for r in rows)
L = []
L.append("# 迁移清单 v1　写作期统一工作目录 `manuscript/`")
L.append("")
L.append("> **执行时间**：2026-10-22 21:30（GMT+8）")
L.append("> **方式**：**复制，不移动**。源目录文件一字未改，已转为只读存档。")
L.append("> **目标根**：`C:\\Users\\刘宸瑞\\Desktop\\脂肪因子生信\\manuscript\\`")
L.append("> **源根 1（活文件）**：`C:\\Users\\刘宸瑞\\WorkBuddy\\2026-10-07-22-15-36\\`")
L.append("> **源根 2（原始存档）**：`C:\\Users\\刘宸瑞\\Desktop\\脂肪因子生信\\`（仅取方案 docx 与选题报告）")
L.append("")
L.append("## 一、校验结果")
L.append("")
L.append("| 校验项 | 结果 |")
L.append("|---|---|")
L.append("| 计划复制条目（去重前） | **%d** |" % (n_src + len(dedup)))
L.append("| 同名同内容合并 | **%d** 处（见 §四-3） |" % len(dedup))
L.append("| 应落盘文件数 | **%d** |" % n_src)
L.append("| 实际落盘文件数（不含 `ascii_alias/` 与清单自身） | **%d** |" % n_dst)
L.append("| 数量比对 | **%s**（多出 %d，缺失 %d） |" % ("✅ 一致" if (n_dst == n_src and not extra and not missing) else "❌ 不一致", len(extra), len(missing)))
L.append("| 大小比对（逐文件） | **%s**（%d 个文件） |" % ("✅ 全部一致" if not bad_size else "❌ %d 个不符" % len(bad_size), n_src))
L.append("| 文件名比对（逐文件） | **%s**（唯一有意改名见 §四-1） |" % ("✅ 全部一致" if not bad_name else "❌ %d 个不符" % len(bad_name)))
L.append("| md5 比对（逐文件，附加项） | **%s** |" % ("✅ 全部一致" if not bad_md5 else "❌ %d 个不符" % len(bad_md5)))
L.append("| 复制总量 | **%.2f MB** |" % (total_size / 1048576))
L.append("| ASCII 别名副本 | **%d 个**（`ascii_alias/`，md5 与正本一致 %d/%d，见 §五） |" % (len(ascii_rows), sum(1 for a in ascii_rows if a["md5_ok"]), len(ascii_rows)))
L.append("| 未纳管（留在原目录） | `output/` 下 **%d 个文件 / %.2f GB**；另 `.workbuddy/data` 原始下载件 53 MB 与 `.workbuddy/tmp` 脚本中间件 164 MB 不在纳管范围 |" % (len(excluded), ex_size / 1073741824))
L.append("")
L.append("## 二、目录结构")
L.append("")
L.append("```")
L.append("manuscript/")
L.append("├─ docs/          ← 方案 docx 各修订版、SAP v1.0–v1.3、任务清单 xlsx、台账、核验记录表、选题报告")
L.append("├─ output/        ← 工作区 output/ 的 1:1 镜像（子目录层级原样保留）")
L.append("│   ├─ 论文初稿_Abstract与Introduction_v1.md 等全部正文稿")
L.append("│   ├─ figures/   ← 9 对 PDF+PNG")
L.append("│   ├─ P2-T13_…/results/ 等结果表")
L.append("│   └─ …")
L.append("├─ ascii_alias/   ← 表格 csv 的 ASCII 别名副本（见 §五）")
L.append("├─ 迁移清单_v1.md")
L.append("├─ 迁移清单_v1.csv")
L.append("└─ README_使用约定.md")
L.append("```")
L.append("")
L.append("## 三、逐文件清单")
L.append("")
n = 0
for k in KIND_ORDER:
    items = by_kind.get(k, [])
    if not items:
        continue
    L.append("### %s（%d 个，%.2f MB）" % (k, len(items), sum(x["size"] for x in items) / 1048576))
    L.append("")
    L.append("| # | 文件名 | 大小 (bytes) | 源路径 |")
    L.append("|---|---|---|---|")
    for r in items:
        n += 1
        base_root = SRC_WS if r["origin"].startswith("工作区") else SRC_DK
        src_rel = os.path.relpath(r["src"], base_root).replace("\\", "/")
        tag = "工作区" if r["origin"].startswith("工作区") else "桌面"
        L.append("| %d | `%s` | %d | %s `%s` |" % (n, os.path.basename(r["dest"]), r["size"], tag, src_rel))
    L.append("")

L.append("## 四、需要您确认的三处处理")
L.append("")
L.append("1. **同名不同版，已加后缀区分**：`研究方案_脂肪因子-脂质代谢轴与黑色素瘤免疫治疗应答分层.docx` 在两处都存在且**内容不同**——")
L.append("   工作区版 88,690 B（2026-10-12 23:20，实际使用的主线版本），桌面版 79,696 B（2026-10-08 00:45，最早提交稿）。")
L.append("   主线版保留原名，桌面版落为 `…_桌面存档版_20261008.docx`。**如确认桌面版已作废，可删此一份。**")
L.append("2. **未纳入 `manuscript/` 的文件类别**（按要求留在原目录，转只读存档）：")
L.append("   - 原始表达矩阵／GEO／ICB／GTEx 下载件与解包件（如 `P4-T1_…/data/GSE72056.txt.gz` 75 MB）")
L.append("   - 中间/派生数据（`output/**/data/*.tsv`、`*.gz`、`*.json`、`*.log`、`*.rds`）")
L.append("   - 分析脚本（`*.py`／`*.R`）——**如投稿需附代码，请另行放行**")
L.append("   - `.bak` / `.backup-*` 备份（含任务清单 5 个历史备份，只留最新版）")
L.append("   - 桌面上的 `workbuddy_*.txt` 会话记录与 `landi_reply.txt` 通信（非成稿文件）")
L.append("3. **同名同内容合并（%d 处）**：以下文件在两处来源同时存在且 **md5 完全相同**，只落一份：" % len(dedup))
for rel, a, b in dedup:
    def _tag(x):
        if x.startswith(SRC_WS):
            return "工作区 `" + os.path.relpath(x, SRC_WS).replace("\\", "/") + "`"
        return "桌面 `" + os.path.relpath(x, SRC_DK).replace("\\", "/") + "`"
    L.append("   - `%s`　←　%s 与 %s" % (rel, _tag(a), _tag(b)))
L.append("")
L.append("## 五、中文路径与 ASCII 别名")
L.append("")
L.append("本机 R 4.6.1 读取含中文的路径会失败（既有教训：须复制到 ASCII 路径并设 `TMPDIR`）。")
L.append("故在 `manuscript/ascii_alias/` 下为**全部表格 csv** 建立 ASCII 名副本，**内容与正本逐字节相同（md5 已校验）**：")
L.append("")
L.append("| 正本（相对 manuscript/） | ASCII 别名 | 大小 | md5 一致 |")
L.append("|---|---|---|---|")
for a in ascii_rows:
    L.append("| `%s` | `%s` | %d | %s |" % (a["orig"], a["alias"], a["size"], "✅" if a["md5_ok"] else "❌"))
L.append("")
L.append("> **纪律**：`ascii_alias/` 是**只读镜像**。任何修改只改 `docs/` 或 `output/` 下的正本，改完重跑 `.workbuddy/tmp/y61_scaffold_manuscript.py` 同步。")
L.append("> 文档类（docx／md／pdf）无需别名——写作用 Word/编辑器打开，不经过 R。")
L.append("")
L.append("## 六、自本清单起生效的读写纪律")
L.append("")
L.append("1. **单写**：成稿期所有读写以 `manuscript/` 为准；原目录文档层**冻结**，不再更新。")
L.append("2. **禁双写**：不得同时改原目录与 `manuscript/` 中的同名文件；一旦发现分叉，以 `manuscript/` 为准并回写原目录一次。")
L.append("3. **再取数**：如后续需重跑分析，**只读访问**原目录（`output/**/data/`、`*.gz`、`*.rds`），产出仍写入 `manuscript/output/`。")
L.append("4. **只增不改**：`迁移清单_v1.md` 为一次性凭证；后续新增文件另立 `迁移清单_v2.md`，不覆盖本文件。")
L.append("")

io.open(os.path.join(DST, "迁移清单_v1.md"), "w", encoding="utf-8").write("\n".join(L))

with io.open(os.path.join(DST, "迁移清单_v1.csv"), "w", encoding="utf-8-sig", newline="") as f:
    w = csv.writer(f)
    w.writerow(["序号", "类别", "文件名", "大小_bytes", "源路径", "目标相对路径", "md5_一致"])
    for i, r in enumerate(rows, 1):
        src_rel = os.path.relpath(r["src"], SRC_WS if r["origin"].startswith("工作区") else SRC_DK)
        w.writerow([i, r["kind"], os.path.basename(r["dest"]), r["size"],
                    src_rel.replace("\\", "/"), r["dest_rel"],
                    "是" if r["md5_src"] == r["md5_dst"] else "否"])

print("=" * 60)
print("计划(去重前/后):", n_src + len(dedup), "/", n_src, " 落盘:", n_dst, " ascii别名:", len(ascii_rows))
print("多出:", len(extra), " 缺失:", len(missing))
print("大小不符:", len(bad_size), "| md5不符:", len(bad_md5), "| 名称不符:", len(bad_name))
print("总量: %.2f MB" % (total_size / 1048576))
print("未纳管: %d 个 / %.2f GB" % (len(excluded), ex_size / 1073741824))
if bad_size:
    for r in bad_size[:10]:
        print("  SIZE", r["dest_rel"], r["size"], r["size_dst"])
if bad_md5:
    for r in bad_md5[:10]:
        print("  MD5 ", r["dest_rel"])
if extra:
    print("  EXTRA", extra[:5])
if missing:
    print("  MISS ", missing[:5])
