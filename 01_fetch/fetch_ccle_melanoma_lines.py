# -*- coding: utf-8 -*-
"""取 CCLE 黑色素瘤细胞系 308 基因表达（cBioPortal API）"""
import json, os, sys, time, io
import urllib.request as ur
import pandas as pd

API = "https://www.cbioportal.org/api"
HDR = {"Content-Type": "application/json", "Accept": "application/json"}
OUT = r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output\P5-T1_细胞系内在轴验证v1\data"
os.makedirs(OUT, exist_ok=True)


def get(url):
    req = ur.Request(url, headers={"Accept": "application/json"})
    with ur.urlopen(req, timeout=90) as r:
        return json.loads(r.read().decode())


def post(url, body):
    req = ur.Request(url, data=json.dumps(body).encode(), headers=HDR)
    with ur.urlopen(req, timeout=300) as r:
        return json.loads(r.read().decode())


# ---- 1. 黑色素瘤细胞系
rows = get(f"{API}/studies/ccle_broad_2019/clinical-data"
           "?clinicalDataType=SAMPLE&attributeId=CANCER_TYPE_DETAILED")
mel = sorted({r["sampleId"] for r in rows if (r["value"] or "").strip().lower() == "cutaneous melanoma"})
print("melanoma lines:", len(mel))
json.dump(mel, io.open(os.path.join(OUT, "ccle_melanoma_lines.json"), "w"), indent=1)

hist = get(f"{API}/studies/ccle_broad_2019/clinical-data"
           "?clinicalDataType=SAMPLE&attributeId=HISTOLOGY")
h = {r["sampleId"]: r["value"] for r in hist if r["sampleId"] in set(mel)}
print("histology 非空:", sum(1 for v in h.values() if v))

# ---- 2. 308 基因 → entrez
s1 = pd.read_csv(r"C:\Users\刘宸瑞\Desktop\脂肪因子生信\manuscript\output\supp"
                 r"\S1_ALMP_308gene_list_and_axis_v1.csv")
sym = [s.strip() for s in s1["gene_symbol"].tolist()]
print("genes:", len(sym))
m = {}
for k, sg in enumerate(sym):
    try:
        x = get(f"{API}/genes/{sg}")
        m[sg] = int(x["entrezGeneId"])
    except Exception:
        pass
    if k % 50 == 0:
        print(f"  gene {k}/{len(sym)}", flush=True)
miss = [s for s in sym if s not in m]
print("mapped:", len(m), "missing:", miss)
entrez = sorted(set(m.values()))

# ---- 3. 取表达（按基因分批，避免单次返回上限）
recs = []
B = 40
el = sorted(entrez)
for i in range(0, len(el), B):
    chunk = el[i:i + B]
    d = post(f"{API}/molecular-profiles/ccle_broad_2019_rna_seq_mrna/molecular-data/fetch",
             {"sampleIds": mel, "entrezGeneIds": chunk})
    recs.extend(d)
    print(f"  genes {i + len(chunk)}/{len(el)} -> {len(d)} 条", flush=True)
    time.sleep(0.3)

df = pd.DataFrame(recs)
df = df.rename(columns={"value": "rpkm"})
df = df[["sampleId", "entrezGeneId", "rpkm"]]
print("rows", len(df), "pairs", len(mel) * len(el),
      "覆盖 %.1f%%" % (100 * len(df) / (len(mel) * len(el))))
df.to_csv(os.path.join(OUT, "ccle_melanoma_308genes_rpkm.tsv"), sep="\t", index=False)
print("saved")
