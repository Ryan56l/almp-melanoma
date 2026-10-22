# Data sources, versions and retrieval dates

All data are public. Nothing in this repository redistributes primary data.

## Bulk tumour transcriptomes

| Resource | Accession | Used for | Retrieved |
|---|---|---|---|
| TCGA-SKCM (n = 472; 470 with matched composition) | cBioPortal `skcm_tcga` mRNA (RSEM) | main bulk layer; coupling, composition adjustment, cell-line comparison | 2026-10-14 / 2026-10-23 |
| GSE65904 (n = 214) | GEO | external prognostic validation (the only external cohort used, per SAP v1.5.1) | 2026-10-20 |
| GSE59455 (141 cores, 23 patients) | GEO | array–sample integrity check | 2026-10-21 |
| GSE72056 (23,686 cells, 19 patients) | GEO | single-cell localisation (4,645 malignant cells analysed) | 2026-10-19 |
| GSE98394 | GEO | exploratory only | 2026-10-19 |

## Immune checkpoint blockade cohorts

| Cohort | Source | Note |
|---|---|---|
| Gide 2019 | cBioPortal / GEO | pre-treatment samples, CR/PR vs SD/PD |
| Riaz 2017 | cBioPortal / GEO | pre-treatment samples |
| Liu 2019 | cBioPortal / GEO | pre-treatment samples |
| Hugo 2016 (UCLA) | cBioPortal / GEO | pre-treatment samples |

Response was harmonised as CR/PR versus SD/PD on pre-treatment samples. All four are
single-arm, so no treatment-by-score interaction is estimable.

## Malignant-cell-only layer

| Resource | Accession | Used for | Retrieved |
|---|---|---|---|
| CCLE cutaneous melanoma cell lines (n = 57 with mRNA) | cBioPortal `ccle_broad_2019_rna_seq_mrna` (RPKM) | direct test of whether the programme is transcribed by malignant cells | 2026-10-23 |
| TCGA-SKCM, same visit | cBioPortal `skcm_tcga` | comparator for detection fractions | 2026-10-23 |

## Causal layer

| Resource | Accession | Used for |
|---|---|---|
| Melanoma outcome (primary) | GCST90475578 (MVP, 10,468 cases) | two-sample MR outcome |
| Replication outcome | ieu-b-4969 | MR replication |
| Sensitivity outcomes | finn-b-C3_SKIN, ieu-b-4959 | MR sensitivity |
| Leptin instrument | GCST90007319 (5 SNPs) | exposure |
| Adiponectin instrument | GCST90007313 (11 SNPs reported; see SAP) | exposure |
| Resistin instrument | GCST90007319-series | exposure |
| BMI instrument | GCST90000501 | negative-control exposure |
| Lymphocyte / WBC counts | GCST90000507 / GCST90000504 | negative-control exposures |
| GTEx v8 | four tissues (skin, subcutaneous adipose, visceral adipose, whole blood) | cis-SMR and colocalisation |

Note: the Landi 2020 melanoma GWAS was **not** obtained; dbGaP access was not initiated. This
is disclosed in the manuscript and in SAP v1.8.

## Software

| Tool | Version | Role |
|---|---|---|
| GSVA | R/Bioconductor | ssGSEA and GSVA scoring on the full transcriptome |
| quantiseqr / quanTIseq | R/Bioconductor | absolute-mode deconvolution |
| MCP-counter | official R implementation | marker-mean deconvolution |
| SMR + HEIDI | 1.3.x | cis-SMR and heterogeneity test |
| coloc | R | colocalisation (supportive only, per SAP v1.6) |
| UCell | R | single-cell signature scoring |
| TwoSampleMR | R | multiplicative random-effects IVW (primary) |

OpenGWAS API tokens used for instrument retrieval expired on 2026-10-25; the retrieved
instrument tables are archived with the project and are reproduced from the accessions above.
