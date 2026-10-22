#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""指令十八① —— 免疫签名基因表达获取（MyGene 映射 + cBioPortal 批量拉取）"""
import json, os, ssl, sys, io, csv, time, urllib.request
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
TMP = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(TMP, 'immune_sets')

syms = [s.strip() for s in open(os.path.join(TMP, 'immune_union_symbols.txt'), encoding='utf-8') if s.strip()]
print('symbols:', len(syms), flush=True)

# ---- 1. MyGene 批量映射 ----
def mygene(batch):
    q = ' OR '.join(f'symbol:{s}' for s in batch)
    url = ('https://mygene.info/v3/query?q=' + urllib.parse.quote(f'symbol:({q})') +
           '&fields=symbol,entrezgene&size=1000&species=human')
    for i in range(3):
        try:
            return json.loads(urllib.request.urlopen(url, timeout=90, context=ctx).read().decode())
        except Exception:
            time.sleep(2)
    return None


import urllib.parse
sym2ent = {}
B = 40
for i in range(0, len(syms), B):
    b = syms[i:i + B]
    r = mygene(b)
    if not r:
        print(f'  batch {i} FAILED', flush=True); continue
    for hit in r.get('hits', []):
        s = (hit.get('symbol') or '').upper()
        e = hit.get('entrezgene')
        if s and e and s in b:
            sym2ent.setdefault(s, str(e))
    print(f'  map {i+len(b)}/{len(syms)} → {len(sym2ent)}', flush=True)
json.dump(sym2ent, open(os.path.join(OUT, 'immune_sym2entrez.json'), 'w', encoding='utf-8'), ensure_ascii=False)
print('mapped symbols:', len(sym2ent), flush=True)

# ---- 2. cBioPortal 拉表达 ----
B2 = 'https://www.cbioportal.org/api'
def req(url, data=None, timeout=180):
    h = {'User-Agent': 'Mozilla/5.0', 'Accept': 'application/json', 'Content-Type': 'application/json'}
    r = urllib.request.Request(url, data=(json.dumps(data).encode() if data is not None else None),
                               headers=h, method=('POST' if data is not None else 'GET'))
    return json.loads(urllib.request.urlopen(r, timeout=timeout, context=ctx).read().decode())

ents = sorted(set(int(v) for v in sym2ent.values()))
print('entrez ids:', len(ents), flush=True)
# 读已有 305 基因（避免重复拉取）
have = set()
with open('C:/Users/刘宸瑞/WorkBuddy/2026-10-07-22-15-36/.workbuddy/data/cbio_skcm/mrna_data.tsv', encoding='utf-8') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        have.add(int(r['entrezGeneId']))
todo = [e for e in ents if e not in have]
print(f'already have {len(have)}; need {len(todo)} new', flush=True)

allrows = []
BS = 40
for i in range(0, len(todo), BS):
    chunk = todo[i:i + BS]
    ok = False
    for att in range(3):
        try:
            d = req(f'{B2}/molecular-profiles/skcm_tcga_rna_seq_v2_mrna/molecular-data/fetch?projection=SUMMARY',
                    {'entrezGeneIds': chunk, 'sampleListId': 'skcm_tcga_rna_seq_v2_mrna'})
            ok = True; break
        except Exception as e:
            print(f'  retry {i} {type(e).__name__}', flush=True); time.sleep(3)
    if not ok:
        print(f'  batch {i} FAILED', flush=True); continue
    allrows.extend(d)
    if (i // BS) % 5 == 0:
        print(f'  fetch {i+len(chunk)}/{len(todo)} → {len(allrows)} rows', flush=True)

with open(os.path.join(OUT, 'immune_mrna_data.tsv'), 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['entrezGeneId', 'sampleId', 'patientId', 'value'],
                       delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(allrows)
print('TOTAL new rows', len(allrows), flush=True)
print('→', os.path.join(OUT, 'immune_mrna_data.tsv'), flush=True)
