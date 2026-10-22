#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
指令二十八：eQTLGen Z-score 方向一致性筛查（8 个 ALMP 关键基因座）
================================================================
★ 纪律（用户指令，强制）：
   ① 本输出一律称「Z-score 方向一致性筛查」，**不得**以 "SMR" 单独修饰或命名。
      理由：无 beta/SE/MAF → 无法计算基于比值的 SMR 点估计与标准误，更无法做 HEIDI。
   ② 每基因座登记：cis 位点数、Z 方向一致性比例、最强位点（|Z_zx| 最大）的 Z_zy 与 Z_zx。
      HEIDI 与共定位一律标注「待 GTEx」。
   ③ CDS2 单列汇总。

输入：
  - eQTLGen cis-eQTL（FDR<0.05，322.8 MB；仅 Zscore，无 beta/SE/MAF）
  - 结局 GCST90475578（harmonised hg38；rsid 列可用；standard_error 全 NA → Z 由 OR+p 反推）
输出：
  - P2T12_eQTL_cis_rows.tsv      8 基因座全部 cis 位点（eQTLGen 原始）
  - P2T12_outcome_rows.tsv       对应区间结局行
  - P2T12_locus_summary.tsv      每基因座汇总（3 项必填 + 描述性符号检验）
  - P2T12_merged_snps.tsv        逐位点合并表（含等位基因对齐标记）
"""
import os, sys, io, ssl, gzip, time, math, json, csv, collections
import urllib.request, urllib.parse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE

ROOT = 'C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36'
TMPDIR = os.path.join(ROOT, '.workbuddy', 'tmp')
OUT = 'C:/Rtmp/smr/out'; os.makedirs(OUT, exist_ok=True)
EQTL = 'C:/Rtmp/smr/eQTLGen_cis_FDR005.txt.gz'

GENES = ['LEP', 'ADIPOQ', 'NAMPT', 'FABP4', 'CD36', 'FADS2', 'CDS2', 'SERPINE1']
WIN = 500_000   # cis 窗：基因 ± 500 kb

sys.path.insert(0, TMPDIR)
import grab_90475578 as G   # noqa: E402


# ---------------------------------------------------------------- 工具
def norm(p):
    """正态分布上尾 p → Z（单侧绝对值）"""
    if p is None or not (0 < p <= 1):
        return None
    # 二分求 Φ⁻¹(1 - p/2)
    lo, hi = 0.0, 40.0
    target = 1 - p / 2
    for _ in range(200):
        mid = (lo + hi) / 2
        if 0.5 * (1 + math.erf(mid / math.sqrt(2))) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


COMP = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def align(z_eqtl, a1, a2, eff, oth):
    """把 eQTL 的 Z 对齐到结局的效应等位基因方向。返回 (z_aligned, status)。"""
    a1, a2, eff, oth = (x.upper() if x else '' for x in (a1, a2, eff, oth))
    if not all([a1, a2, eff, oth]):
        return None, 'missing_allele'
    palin = (a1, a2) in [('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C')]
    if palin:
        return None, 'palindromic_excluded'
    if (a1, a2) == (eff, oth):
        return z_eqtl, 'same'
    if (a1, a2) == (oth, eff):
        return -z_eqtl, 'swapped'
    if (COMP.get(a1), COMP.get(a2)) == (eff, oth):
        return z_eqtl, 'complement'
    if (COMP.get(a1), COMP.get(a2)) == (oth, eff):
        return -z_eqtl, 'complement_swapped'
    return None, 'allele_mismatch'


# ---------------------------------------------------------------- 1. 基因坐标（hg38 用于抓结局；hg19 备用）
def gene_pos(genes):
    q = ' OR '.join(f'symbol:{g}' for g in genes)
    u = ('https://mygene.info/v3/query?q=' + urllib.parse.quote(q) +
         '&fields=symbol,genomic_pos,genomic_pos_hg19&species=human&size=100')
    d = json.load(urllib.request.urlopen(u, context=CTX, timeout=60))
    out = {}
    for h in d['hits']:
        s = h.get('symbol')
        if s not in genes or s in out:
            continue
        g38 = h.get('genomic_pos'); g19 = h.get('genomic_pos_hg19')
        if isinstance(g38, list): g38 = g38[0]
        if isinstance(g19, list): g19 = g19[0]
        if g38:
            out[s] = dict(chr=str(g38['chr']), start=g38['start'], end=g38['end'],
                          hg19_start=(g19['start'] if g19 else None))
    return out


GP = gene_pos(GENES)
print('[1] hg38 坐标：')
for g in GENES:
    if g in GP:
        p = GP[g]
        print(f'    {g:10s} chr{p["chr"]}:{p["start"]:,}-{p["end"]:,}')
    else:
        print(f'    {g:10s} ✗ 未取到')
assert len(GP) == len(GENES), '有基因未取到坐标'

# ---------------------------------------------------------------- 2. 流式筛 eQTLGen
def scan_eqtl():
    cache = os.path.join(OUT, 'P2T12_eQTL_cis_rows.tsv')
    if os.path.exists(cache) and os.path.getsize(cache) > 1000:
        print(f'[2] 复用缓存 {cache}')
        return cache
    want = set(GENES)
    t0 = time.time(); n = 0; kept = 0
    with gzip.open(EQTL, 'rt', encoding='utf-8', errors='replace') as f, \
         open(cache, 'w', encoding='utf-8', newline='') as o:
        hdr = f.readline().rstrip('\n').split('\t')
        idx = {k: i for i, k in enumerate(hdr)}
        cols = ['Pvalue', 'SNP', 'SNPChr', 'SNPPos', 'AssessedAllele', 'OtherAllele',
                'Zscore', 'GeneSymbol', 'NrSamples']
        o.write('\t'.join(cols) + '\n')
        ci = [idx[c] for c in cols]
        for ln in f:
            n += 1
            c = ln.rstrip('\n').split('\t')
            if len(c) <= max(ci):
                continue
            if c[ci[7]] not in want:
                continue
            o.write('\t'.join(c[i] for i in ci) + '\n'); kept += 1
            if n % 2_000_000 < 1:
                print(f'    扫描 {n/1e6:.0f}M 行，命中 {kept}（{time.time()-t0:.0f}s）', flush=True)
    print(f'[2] eQTLGen 扫描完成：{n:,} 行，命中 {kept:,} 行（{time.time()-t0:.0f}s）')
    return cache


EQ = scan_eqtl()

# ---------------------------------------------------------------- 3. 抓结局区间
def grab_outcome():
    cache = os.path.join(OUT, 'P2T12_outcome_rows.tsv')
    if os.path.exists(cache) and os.path.getsize(cache) > 1000:
        print(f'[3] 复用缓存 {cache}')
        return cache
    rows = set()
    with open(cache, 'w', encoding='utf-8') as o:
        for g in GENES:
            p = GP[g]
            c = int(p['chr'])
            a = max(1, p['start'] - WIN); b = p['end'] + WIN
            print(f'[3] {g}: chr{c}:{a:,}-{b:,}', flush=True)
            rs = G.grab_region(c, a, b)
            print(f'    -> {len(rs)} 行', flush=True)
            for r in rs:
                if r not in rows:
                    rows.add(r); o.write(r + '\n')
    print(f'[3] 结局行合计 {len(rows):,}')
    return cache


OC = grab_outcome()

# ---------------------------------------------------------------- 4. 合并 + 对齐 + 汇总
eq = {}
with open(EQ, encoding='utf-8') as f:
    rd = csv.DictReader(f, delimiter='\t')
    for r in rd:
        eq.setdefault(r['SNP'], []).append(r)
print(f'[4] eQTL 唯一 rsID {len(eq):,}')

HDR = ['chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'odds_ratio',
       'standard_error', 'effect_allele_frequency', 'p_value', 'rsid', 'ci_upper', 'ci_lower',
       'alt', 'n', 'case_af', 'num_cases', 'control_af', 'num_controls', 'r2', 'q_pval',
       'i2', 'direction', 'hm_coordinate_conversion', 'hm_code', 'variant_id']
oc = {}
with open(OC, encoding='utf-8') as f:
    rd = csv.DictReader(f, delimiter='\t', fieldnames=HDR)
    for r in rd:
        rs_ = (r.get('rsid') or '').strip()
        if not rs_ or rs_ == 'rsid':      # 跳过表头行/空行
            continue
        oc[rs_] = r
print(f'[4] 结局唯一 rsID {len(oc):,}')

merged = []
for rs, lst in eq.items():
    o = oc.get(rs)
    if not o:
        continue
    try:
        orv = float(o['odds_ratio']); pv = float(o['p_value'])
    except Exception:
        continue
    zy = norm(pv)
    if zy is None:
        continue
    zy = zy if orv > 1 else -zy          # OR>1 → 风险升高
    for e in lst:
        try:
            zx = float(e['Zscore'])
        except Exception:
            continue
        za, st = align(zx, e['AssessedAllele'], e['OtherAllele'],
                       o['effect_allele'], o['other_allele'])
        if za is None:
            continue
        merged.append(dict(gene=e['GeneSymbol'], rsid=rs, z_eqtl=za, z_out=zy,
                           z_eqtl_raw=zx, align=st,
                           chr=e['SNPChr'], pos_hg19=e['SNPPos'],
                           a1=e['AssessedAllele'], a2=e['OtherAllele'],
                           eff=o['effect_allele'], oth=o['other_allele'],
                           orv=orv, p=pv, n_eqtl=e['NrSamples']))
print(f'[4] 合并可用位点 {len(merged):,}')

with open(os.path.join(OUT, 'P2T12_merged_snps.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(merged[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(merged)

# 汇总
summ = []
for g in GENES:
    ms = [m for m in merged if m['gene'] == g]
    if not ms:
        n_eqtl = sum(1 for r in csv.DictReader(open(EQ, encoding='utf-8'), delimiter='\t')
                     if r['GeneSymbol'] == g)
        note = ('eQTLGen 无该基因 cis-eQTL（全血检测不到）' if n_eqtl == 0
                else f'eQTLGen 有 {n_eqtl} 条 cis 位点，但全部在等位基因对齐/结局匹配阶段被剔除')
        summ.append(dict(gene=g, n_cis=0, prop_concordant=None,
                         top_rsid=None, top_z_eqtl=None, top_z_out=None,
                         top_ratio=None, direction=None, sign_test_p=None, note=note))
        continue
    conc = sum(1 for m in ms if m['z_eqtl'] * m['z_out'] > 0)
    top = max(ms, key=lambda m: abs(m['z_eqtl']))
    ratio = top['z_out'] / top['z_eqtl']
    # 描述性符号检验（双侧；对数空间计算，避免大 n 组合数溢出）
    k, n = conc, len(ms)
    def tail_le(kk):
        if kk < 0: return 0.0
        if kk >= n: return 1.0
        lg = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + n * math.log(0.5)
              for i in range(0, kk + 1)]
        m = max(lg)
        return math.exp(m) * sum(math.exp(x - m) for x in lg)
    p_lo = tail_le(k)            # P(X <= k)
    p_hi = 1.0 - tail_le(k - 1)  # P(X >= k)
    p = min(1.0, 2 * min(p_lo, p_hi))
    if abs(top['z_out']) < 1.96:
        direction = '不可判（top 位点与结局无关联 |Z_out|<1.96）'
    else:
        direction = ('表达升高→风险升高' if ratio > 0 else '表达升高→风险降低')
    summ.append(dict(gene=g, n_cis=n, prop_concordant=round(conc / n, 3),
                     top_rsid=top['rsid'], top_z_eqtl=round(top['z_eqtl'], 3),
                     top_z_out=round(top['z_out'], 3), top_ratio=round(ratio, 4),
                     direction=direction, sign_test_p=round(p, 4), note=''))

with open(os.path.join(OUT, 'P2T12_locus_summary.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(summ[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(summ)

print('\n===== Z-score 方向一致性筛查 · 每基因座汇总 =====')
print(f'{"gene":10s} {"n_cis":>6s} {"concord":>8s} {"top_rsid":>14s} {"Z_eQTL":>9s} {"Z_out":>8s} {"ratio":>8s} {"sign_p":>8s}  direction')
for s in summ:
    if s['n_cis'] == 0:
        print(f'{s["gene"]:10s} {0:>6d}     —            —         —        —        —        —    无可用位点')
        continue
    print(f'{s["gene"]:10s} {s["n_cis"]:>6d} {s["prop_concordant"]:>8.3f} {s["top_rsid"]:>14s} '
          f'{s["top_z_eqtl"]:>9.3f} {s["top_z_out"]:>8.3f} {s["top_ratio"]:>8.4f} {s["sign_test_p"]:>8.4f}  {s["direction"]}')

print('\n★ CDS2 单列：')
for s in summ:
    if s['gene'] == 'CDS2':
        print('   ', s)
print('\n===== 位点剔除归因（透明登记）=====')
attr = []
for g in GENES:
    row = dict(gene=g)
    rows = [e for r in eq.values() for e in r if e['GeneSymbol'] == g]
    row['n_eqtlgen_rows'] = len(rows)
    row['no_outcome_match'] = sum(1 for e in rows if e['SNP'] not in oc)
    st = collections.Counter()
    for e in rows:
        o = oc.get(e['SNP'])
        if not o:
            continue
        try:
            orv = float(o['odds_ratio']); pv = float(o['p_value'])
        except Exception:
            st['bad_or_p'] += 1; continue
        if norm(pv) is None:
            st['bad_or_p'] += 1; continue
        za, s2 = align(float(e['Zscore']), e['AssessedAllele'], e['OtherAllele'],
                       o['effect_allele'], o['other_allele'])
        st[s2] += 1
    row['palindromic'] = st.get('palindromic_excluded', 0)
    row['allele_mismatch'] = st.get('allele_mismatch', 0)
    row['missing_allele'] = st.get('missing_allele', 0)
    row['bad_or_p'] = st.get('bad_or_p', 0)
    row['used'] = st.get('same', 0) + st.get('swapped', 0) + st.get('complement', 0) + st.get('complement_swapped', 0)
    attr.append(row)
for a in attr:
    print('   ', a)
with open(os.path.join(OUT, 'P2T12_snp_attrition.tsv'), 'w', encoding='utf-8', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(attr[0].keys()), delimiter='\t')
    w.writeheader(); w.writerows(attr)

print('\nHEIDI：待 GTEx　｜　共定位（coloc PP4）：待 GTEx')
print('DONE')
