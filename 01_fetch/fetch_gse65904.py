#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""并行分片下载（GEO FTP 单连接仅 ~11 KB/s，但支持 Range → 分片并行提速）。"""
import os, sys, io, ssl, time, json
import urllib.request
from concurrent.futures import ThreadPoolExecutor

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
CTX = ssl.create_default_context(); CTX.check_hostname = False; CTX.verify_mode = ssl.CERT_NONE
DST = os.path.dirname(os.path.abspath(__file__))
if not os.path.isdir(DST) or DST == '':
    DST = os.getcwd()

TASKS = [
    ('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65904/matrix/GSE65904_series_matrix.txt.gz',
     'GSE65904_series_matrix.txt.gz', 6),
    ('https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz',
     'GPL10558.annot.gz', 4),
]


def size_of(u):
    r = urllib.request.Request(u, method='HEAD')
    return int(urllib.request.urlopen(r, context=CTX, timeout=90).headers['Content-Length'])


def fetch(url, out, nchunk):
    total = size_of(url)
    part = out + '.parts'
    os.makedirs(part, exist_ok=True)
    bounds = []
    step = total // nchunk
    for i in range(nchunk):
        a = i * step
        b = total - 1 if i == nchunk - 1 else (a + step - 1)
        bounds.append((i, a, b))
    done = [0]

    def one(t):
        i, a, b = t
        pf = os.path.join(part, f'{os.path.basename(out)}.{i:03d}')
        if os.path.exists(pf) and os.path.getsize(pf) == b - a + 1:
            done[0] += b - a + 1
            return pf
        for attempt in range(6):
            try:
                r = urllib.request.Request(url, headers={'Range': f'bytes={a}-{b}'})
                o = urllib.request.urlopen(r, context=CTX, timeout=600)
                buf = b''.join(iter(lambda: o.read(262144), b''))
                if len(buf) == b - a + 1:
                    with open(pf, 'wb') as f:
                        f.write(buf)
                    done[0] += len(buf)
                    return pf
            except Exception as e:
                if attempt == 3:
                    print(f'   chunk {i} FAIL {e}', flush=True)
                time.sleep(4 * (attempt + 1))
        return None

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=nchunk) as ex:
        res = list(ex.map(one, bounds))
    ok = [p for p in res if p]
    print(f'   {len(ok)}/{nchunk} 分片完成，{done[0]/1e6:.1f} MB，{time.time()-t0:.0f}s', flush=True)
    if len(ok) != nchunk:
        raise RuntimeError(f'{out}: 仅 {len(ok)}/{nchunk} 分片成功')
    with open(out, 'wb') as f:
        for i in range(nchunk):
            with open(os.path.join(part, f'{os.path.basename(out)}.{i:03d}'), 'rb') as g:
                f.write(g.read())
    n = os.path.getsize(out)
    assert n == total, f'{out} 大小不符 {n} != {total}'
    print(f'   ✓ {out}: {n:,} bytes（应为 {total:,}）', flush=True)


if __name__ == '__main__':
    for url, out, nc in TASKS:
        p = os.path.join(DST, out)
        print(f'=== {out} ===', flush=True)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            try:
                if os.path.getsize(p) == size_of(url):
                    print('   已完整，跳过', flush=True); continue
            except Exception:
                pass
        fetch(url, p, nc)
    print('ALL DONE', flush=True)
