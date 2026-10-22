#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GPL8432 平台表提取：SOFT 文件 896 MB，但平台表位于文件最前 → gzip 部分解压即可（不必全量下载）。"""
import ssl, urllib.request, zlib, sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
U = 'https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL8nnn/GPL8432/soft/GPL8432_family.soft.gz'

for span in (13631487, 20971519, 29360127):
    raw = urllib.request.urlopen(
        urllib.request.Request(U, headers={'Range': f'bytes=0-{span}', 'User-Agent': 'Mozilla/5.0'}),
        context=CTX, timeout=1200).read()
    txt = zlib.decompressobj(31).decompress(raw).decode('utf-8', 'replace')
    lines = txt.split('\n')
    b = [i for i, l in enumerate(lines) if l.startswith('!platform_table_begin')]
    e = [i for i, l in enumerate(lines) if l.startswith('!platform_table_end')]
    print(f'range {span/1e6:.1f} MB -> 解出 {len(txt):,} 字符; end={e[:1]}', flush=True)
    if not b:
        continue
    hdr = lines[b[0] + 1].rstrip().split('\t')
    stop = e[0] if e else len(lines)
    rows = [r for r in lines[b[0] + 2:stop] if r.strip() and not r.startswith('!')]
    out = [f'probe\tsymbol']
    for r in rows:
        c = r.rstrip().split('\t')
        if len(c) >= 6 and c[0].startswith('ILMN'):
            out.append(f'{c[0]}\t{c[5]}')
    open('GPL8432_probe2symbol.tsv', 'w', encoding='utf-8').write('\n'.join(out))
    print(f'  行数 {len(rows)}，写出 {len(out)-1} 映射 -> GPL8432_probe2symbol.tsv', flush=True)
    if e:
        print('  ✓ 平台表完整取到', flush=True)
        break
print('DONE', flush=True)
