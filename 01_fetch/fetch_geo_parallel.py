# -*- coding: utf-8 -*-
"""
通用 GEO 分片并发下载器（断点续传）——P4-T1 用
=================================================
用法： python geo_parallel_dl.py <url> <out_path> [线程数]
机制：先 Range 0-0 探总长 → 按 1 MiB 切片 → 并发取 → 合并 → 校验大小
背景：GEO 按 IP 限流（单线程 ~10 KB/s）；HTTP 支持 Range → 多线程可提速数倍
"""
import os
import sys
import ssl
import time
import io
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")
ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE
UA = {"User-Agent": "Mozilla/5.0 (academic-research)"}

CHUNK = 1 << 20  # 1 MiB


def probe(url):
    req = urllib.request.Request(url, headers={**UA, "Range": "bytes=0-0"})
    with urllib.request.urlopen(req, timeout=60, context=ctx) as r:
        cr = r.headers.get("Content-Range", "")
        if "/" in cr:
            return int(cr.split("/")[-1])
        cl = r.headers.get("Content-Length")
        return int(cl) if cl else None


def fetch(url, st, en, path, tries=4):
    want = en - st + 1
    if os.path.exists(path) and os.path.getsize(path) == want:
        return True, want, 0.0
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={**UA, "Range": f"bytes={st}-{en}"})
            t0 = time.time()
            with urllib.request.urlopen(req, timeout=90, context=ctx) as r:
                buf = bytearray()
                while len(buf) < want:
                    b = r.read(65536)
                    if not b:
                        break
                    buf += b
            if len(buf) != want:
                raise IOError(f"short {len(buf)}/{want}")
            with open(path, "wb") as f:
                f.write(buf)
            return True, want, time.time() - t0
        except Exception as e:
            if a == tries - 1:
                return False, str(e)[:70], 0.0
            time.sleep(1.5 * (a + 1))


def main():
    url, out = sys.argv[1], sys.argv[2]
    nthread = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    total = probe(url)
    if total is None:
        raise SystemExit("[FATAL] 无法取得总长度（服务器未返回 Content-Range）")
    print(f"[1] 总长 {total:,} B ({total/1048576:.1f} MiB)；线程 {nthread}")

    parts = os.path.join(os.path.dirname(os.path.abspath(out)), "_parts")
    os.makedirs(parts, exist_ok=True)
    chunks, s, i = [], 0, 0
    while s < total:
        e = min(s + CHUNK - 1, total - 1)
        chunks.append((i, s, e))
        s = e + 1
        i += 1
    print(f"[2] 切片 {len(chunks)} 个")

    ok = fail = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=nthread) as ex:
        futs = {ex.submit(fetch, url, st, en, os.path.join(parts, "c%05d.part" % idx)): idx
                for idx, st, en in chunks}
        for n, f in enumerate(as_completed(futs), 1):
            good, info, dt = f.result()
            if good:
                ok += 1
            else:
                fail += 1
                print(f"    ✗ 片 {futs[f]} 失败：{info}")
            if n % 20 == 0:
                print(f"    进度 {n}/{len(chunks)}  用时 {time.time()-t0:.0f}s")

    print(f"[3] 成功 {ok} / 失败 {fail}；总耗时 {time.time()-t0:.0f}s ({total/1048576/max(1e-9,(time.time()-t0)):.2f} MiB/s)")
    if fail:
        raise SystemExit("[FATAL] 存在失败片，重跑本命令可续传")

    with open(out, "wb") as fo:
        for idx, _, _ in chunks:
            with open(os.path.join(parts, "c%05d.part" % idx), "rb") as fi:
                fo.write(fi.read())
    sz = os.path.getsize(out)
    print(f"[4] 合并 → {out} ({sz:,} B)；大小一致={sz == total}")
    if sz == total:
        for idx, _, _ in chunks:
            os.remove(os.path.join(parts, "c%05d.part" % idx))
        os.rmdir(parts)


if __name__ == "__main__":
    main()
