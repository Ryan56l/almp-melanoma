#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
指令二十一 · v3：补抓 Riaz 2017（MSK, anti-PD-1, n=107）与 UCLA/Hugo 2016（iAtlas 协调版, n=38）。
基因集沿用 v2 生成的 symbol_entrez_map_union.tsv（ALMP ∪ MCP-counter 111 标记）。

Hugo/UCLA 与 GSE78220 同源（Hugo et al. Cell 2016）→ ★ 两者不得作为独立队列合并计入。
"""
import os, sys, io, csv, json, ssl, time
import urllib.request, urllib.parse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
ROOT = 'C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36'
OUT = os.path.join(ROOT, '.workbuddy', 'data', 'icb_v2')
API = 'https://www.cbioportal.org/api'

STUDIES = [
    ('mel_iatlas_riaz_nivolumab_2017', 'mel_iatlas_riaz_nivolumab_2017_rna_seq_mrna', 'Riaz2017'),
    ('mel_ucla_2016',                  'mel_ucla_2016_rna_seq_mrna',                  'UCLA2016'),
]


def call(url, body=None, tries=4):
    data = json.dumps(body).encode() if body is not None else None
    req = urllib.request.Request(url, data=data, headers={
        'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json',
        'Content-Type': 'application/json'})
    last = None
    for _ in range(tries):
        try:
            return json.load(urllib.request.urlopen(req, context=CTX, timeout=240))
        except Exception as e:
            last = e; time.sleep(3)
    raise last


# ---------- 基因清单：沿用并集映射 ----------
e2s = {}
with open(os.path.join(OUT, 'symbol_entrez_map_union.tsv'), encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        if r['entrez']:
            e2s[r['symbol']] = r['entrez']
entrez_all = sorted({int(v) for v in e2s.values()})
print(f'[0] 并集 entrez = {len(entrez_all)}')

report = {}
for sid, pid, tag in STUDIES:
    print('=' * 88); print(f'### {tag} ({sid})')
    # 临床
    clin = call(f'{API}/studies/{sid}/clinical-data?clinicalDataType=SAMPLE&projection=SUMMARY')
    attrs = sorted({c['clinicalAttributeId'] for c in clin})
    bys = {}
    for c in clin:
        bys.setdefault(c['sampleId'], {})[c['clinicalAttributeId']] = c['value']
    rows = [dict(sampleId=s, **v) for s, v in bys.items()]
    hdr = ['sampleId'] + [a for a in attrs if a != 'sampleId']
    cp = os.path.join(OUT, f'{tag}_clinical.tsv')
    with open(cp, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f'[1] clinical: {len(rows)} samples x {len(hdr)-1} attrs -> {cp}')
    for k in ['RESPONSE', 'TREATMENT_RESPONSE', 'SAMPLE_TREATMENT', 'TMB_NONSYNONYMOUS']:
        if k in attrs:
            import collections
            cc = collections.Counter(str(r.get(k, '')) for r in rows)
            print(f'      {k}: {dict(sorted(cc.items(), key=lambda t: -t[1])[:8])}')

    # 表达
    mat, samples = {}, []
    B = 90
    for i in range(0, len(entrez_all), B):
        batch = entrez_all[i:i+B]
        d = call(f'{API}/molecular-profiles/{pid}/molecular-data/fetch',
                 {'sampleListId': pid, 'entrezGeneIds': batch})
        for rec in d:
            s = rec['sampleId']; e = str(rec['entrezGeneId'])
            mat.setdefault(e, {})[s] = rec['value']
            if s not in samples:
                samples.append(s)
        print(f'[2] batch {i//B+1}/{-(-len(entrez_all)//B)}: 累积基因 {len(mat)}')
        time.sleep(0.4)
    samples = sorted(samples)
    ep = os.path.join(OUT, f'{tag}_expr_union.tsv')
    genes = sorted(mat.keys(), key=lambda x: int(x))
    with open(ep, 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['entrez'] + samples)
        for g in genes:
            w.writerow([g] + [mat[g].get(s, '') for s in samples])
    print(f'[2] expr_union: {len(genes)} genes x {len(samples)} samples -> {ep}')
    report[tag] = dict(genes=len(genes), samples=len(samples))

print('\nSUMMARY', json.dumps(report, ensure_ascii=False))
print('DONE ->', OUT)
