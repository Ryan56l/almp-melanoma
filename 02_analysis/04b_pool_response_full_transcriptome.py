#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""04b_pool_response_full_transcriptome.py

Pool the four-cohort ICB response contrasts on **full-transcriptome** ALMP scores
and write the pooled block of Table S3 (``S3c_pooled``).

Why this script exists
----------------------
The archived pooling (``04_icb_cohort_expansion.R``) was computed on the 413-gene
subset universe. After the gene-universe rule (SAP v1.11) the primary scores are
full-transcriptome, so the pooled estimates had to be recomputed on that universe.
This script is the recomputation and it is the source of the numbers printed in
Figure 3b and in ``S3c_pooled`` of Additional file 3.

Input
-----
``P2T6bR_response_FULL.csv`` — 48 rows, one per cohort x collection x algorithm,
with columns: cohort, readout, feature, n_R, n_NR, med_R, med_NR, delta,
cohen_d, auc, p, q, universe, q_pool.

Method
------
1. Variance of Cohen's d from the two group sizes:
       var(d) = (n_R + n_NR) / (n_R * n_NR) + d^2 / (2 * (n_R + n_NR))
2. DerSimonian–Laird random effects:
       tau^2 = max(0, (Q - (k-1)) / C),  C = sum(w) - sum(w^2)/sum(w)
       pooled d = sum(w* d) / sum(w*),  w* = 1 / (var + tau^2)
       95% CI  = pooled d +/- 1.96 * sqrt(1 / sum(w*))
3. I^2 = max(0, (Q - (k-1)) / Q) * 100; Q p-value from chi-square with k-1 df
4. Benjamini–Hochberg across the 12 pooled tests (6 collections x 2 algorithms)

Usage
-----
    python 04b_pool_response_full_transcriptome.py \
        --input data/P2T6bR_response_FULL.csv \
        --out   S3c_pooled_response.csv

Output columns match the ``S3c_pooled`` block of Table S3.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys

try:
    from scipy.stats import chi2, norm
except ImportError:  # pragma: no cover
    sys.exit("scipy is required: pip install scipy")


COLLECTIONS = [
    "SEG1_ADIPOKINE_CORE",
    "SEG2_LIPID_PATHWAY",
    "SEG3_EXECUTOR",
    "ALMP_COMPOSITIONAL_AXIS",
    "ALMP_TUMOR_INTRINSIC_AXIS",
    "ALMP_ALL_308",
]
COHORT_ORDER = ["Gide2019", "Liu2019", "Riaz2017", "UCLA2016"]


def var_d(d: float, n1: float, n2: float) -> float:
    """Variance of Cohen's d from two independent group sizes."""
    return (n1 + n2) / (n1 * n2) + d * d / (2.0 * (n1 + n2))


def pool(rows: list[dict]) -> dict:
    """DerSimonian–Laird random-effects pooling of one collection x algorithm."""
    d = [float(r["cohen_d"]) for r in rows]
    v = [var_d(x, float(r["n_R"]), float(r["n_NR"])) for x, r in zip(d, rows)]
    k = len(d)

    w = [1.0 / x for x in v]
    sw = sum(w)
    fixed = sum(wi * di for wi, di in zip(w, d)) / sw
    q = sum(wi * (di - fixed) ** 2 for wi, di in zip(w, d))
    c = sw - sum(wi * wi for wi in w) / sw
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0

    ws = [1.0 / (x + tau2) for x in v]
    sws = sum(ws)
    est = sum(wi * di for wi, di in zip(ws, d)) / sws
    se = (1.0 / sws) ** 0.5
    lo, hi = est - 1.959964 * se, est + 1.959964 * se
    p = 2.0 * norm.sf(abs(est / se))
    i2 = max(0.0, (q - (k - 1)) / q) * 100.0 if q > 0 else 0.0
    qp = chi2.sf(q, k - 1) if k > 1 else float("nan")

    n_pos = sum(1 for x in d if x > 0)
    return {
        "k": k,
        "d_pooled": est,
        "ci_lo": lo,
        "ci_hi": hi,
        "p": p,
        "I2": i2,
        "Qp": qp,
        "n_pos": n_pos,
        "n_neg": k - n_pos,
        "tau2": tau2,
        "cohorts": "+".join(r["cohort"] for r in rows),
    }


def bh(pvals: list[float]) -> list[float]:
    """Benjamini–Hochberg adjusted p-values."""
    m = len(pvals)
    order = sorted(range(m), key=lambda i: pvals[i])
    out = [1.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        val = min(prev, pvals[i] * m / rank)
        out[i] = val
        prev = val
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="data/P2T6bR_response_FULL.csv")
    ap.add_argument("--out", default="S3c_pooled_response.csv")
    args = ap.parse_args()

    rows = list(csv.DictReader(io.open(args.input, encoding="utf-8-sig")))
    if not rows:
        sys.exit("empty input: %s" % args.input)

    universes = {r.get("universe", "") for r in rows}
    if universes - {"FULL"}:
        print("[warn] input contains non-FULL rows: %s" % sorted(universes - {"FULL"}))

    results = []
    for coll in COLLECTIONS:
        for alg in ("ssGSEA", "GSVA"):
            sub = [r for r in rows if r["feature"] == coll and r["readout"] == alg]
            if not sub:
                continue
            order = {c: i for i, c in enumerate(COHORT_ORDER)}
            sub.sort(key=lambda r: order.get(r["cohort"], 99))
            res = pool(sub)
            res.update({"readout": alg, "feature": coll})
            results.append(res)

    qs = bh([r["p"] for r in results])
    for r, qv in zip(results, qs):
        r["q"] = qv

    cols = ["readout", "feature", "k", "d_pooled", "ci_lo", "ci_hi", "p",
            "I2", "Qp", "n_neg", "n_pos", "cohorts", "q", "tau2"]
    with io.open(args.out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        for r in results:
            w.writerow([
                r["readout"], r["feature"], r["k"],
                "%.4f" % r["d_pooled"], "%.4f" % r["ci_lo"], "%.4f" % r["ci_hi"],
                "%.4f" % r["p"], "%.1f" % r["I2"], "%.4f" % r["Qp"],
                r["n_neg"], r["n_pos"], r["cohorts"], "%.4f" % r["q"], "%.4f" % r["tau2"],
            ])

    print("wrote %s (%d pooled estimates)" % (args.out, len(results)))
    for r in results:
        if r["feature"] == "SEG1_ADIPOKINE_CORE":
            print("  %-6s d=%+.3f (%+.3f to %+.3f)  I2=%.1f%%  p=%.4f  q=%.4f"
                  % (r["readout"], r["d_pooled"], r["ci_lo"], r["ci_hi"],
                     r["I2"], r["p"], r["q"]))


if __name__ == "__main__":
    main()
