# Adipokine–lipid metabolic programme (ALMP) in cutaneous melanoma — analysis code

Archival analysis code for the manuscript *"A decomposable adipokine–lipid programme in
cutaneous melanoma: a stromal adipokine core and a stroma-independent tumour-intrinsic
axis"*.

This repository contains **the scripts that produced every number, figure and table** in the
manuscript. It is a **transparency archive, not a turnkey pipeline**: scripts are kept in the
form in which they were run, and file paths are hard-coded to the author's local directory
layout (see *Running the code* below).

---

## What the study does

A 308-gene adipokine–lipid metabolic programme (ALMP) is defined by an **external** criterion
(HPA adipose nTPM ≥ 10 AND adipose / max(non-adipose) ≥ 3), not by differential expression in
the data at hand. The programme is then decomposed into a 22-gene adipokine core, a 247-gene
lipid-pathway segment and a 39-gene executor segment, and along a second axis into
**16 tissue-compositional** and **292 tumour-intrinsic** members.

The paper reports (i) that the core behaves as a stromal readout rather than a tumour-cell
signal, (ii) that the intrinsic axis is the only component attributable to malignant cells,
and (iii) a pre-registered causal audit (two-sample MR, cis-SMR with HEIDI, colocalisation)
that is **negative** across all layers. Four claims were withdrawn during the analysis and the
withdrawals are reported in the manuscript.

## Repository layout

| Directory | Contents |
|---|---|
| `01_fetch/` | data acquisition from cBioPortal, GEO, GTEx and CCLE |
| `02_analysis/` | deconvolution, scoring, coupling, composition adjustment, MR/SMR, single-cell, cell-line validation |
| `03_figures/` | figure-generating scripts (Figures 1–7) |
| `04_tables/` | Tables 1–2 |
| `05_assemble/` | word counting and manuscript assembly |
| `MANIFEST.tsv` | file-by-file description and MD5 |
| `DATA_SOURCES.md` | every accession, version and retrieval date |

## Analysis order

```
01_fetch/            →  raw expression, genotypes, cohort tables
02_analysis/01–05    →  scores, deconvolution, coupling
02_analysis/06       →  random gene-set negative control
02_analysis/07       →  cis-SMR / HEIDI / colocalisation
02_analysis/08       →  single-cell localisation
02_analysis/09–11    →  external validation, direction screen, LASSO-Cox (exploratory)
02_analysis/12       →  full-transcriptome recount (gene-universe discipline)
02_analysis/13*      →  composition adjustment and identifiability  ← the P0 analysis
02_analysis/14*      →  melanoma cell-line validation
03–05                →  figures, tables, assembly
```

## Two rules that the code enforces

1. **Gene-universe rule.** Enrichment scores are computed on the **full transcriptome**, never
   on a subset matrix, and both sides of any comparison use the same universe
   (`02_analysis/12_full_transcriptome_recount.py`, `13_purity_deconfound.py`).
2. **Pre-specified model.** The primary MR estimate is a multiplicative random-effects IVW
   (`02_analysis/07_smr_heidi_coloc.py` and the SAP, Supplementary Table S4b). Alternative
   estimators are reported as sensitivity analyses only.

## Running the code

* Python 3.13 with `numpy`, `scipy`, `pandas`, `matplotlib` (see `requirements.txt`);
  R 4.6.1 with `GSVA`, `quantiseqr`, `survival`, `glmnet` (see `environment.txt`).
* Each script begins with a path block. Replace the two roots
  (`ROOT`, and where present `DATA` / `OUT`) with your local checkout before running.
* Scripts are ordered rather than idempotent: later scripts read tables written by earlier
  ones. `MANIFEST.tsv` gives the intended order.

## Data

No primary data are redistributed here. Every input is public and listed with its accession,
version and retrieval date in `DATA_SOURCES.md`. TCGA-SKCM and CCLE expression were retrieved
through the cBioPortal API; GTEx and the four ICB cohorts from their original repositories.

## Licence

Code: MIT (`LICENSE`). Data remain under the terms of their original repositories.

## Contact

Liu Chenrui — Kunming University of Science and Technology — chenruiliu@stu.kust.edu.cn
