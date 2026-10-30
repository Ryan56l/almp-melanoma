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
| Riaz 2017 | GEO (GSE91061) | pre-treatment samples |
| Liu 2019 | cBioPortal / GEO | pre-treatment samples |
| Hugo 2016 (UCLA) | GEO (GSE78220) | pre-treatment samples |

Response was harmonised as CR/PR versus SD/PD on pre-treatment samples. All four are
single-arm, so no treatment-by-score interaction is estimable.

Expression matrices were taken from the source giving the fuller gene coverage for that cohort
and merged on sample identifiers; the manuscript therefore cites the repositories collectively
rather than a single interface.

## Malignant-cell-only layer

| Resource | Accession | Used for | Retrieved |
|---|---|---|---|
| CCLE cutaneous melanoma cell lines (n = 57 with mRNA) | cBioPortal `ccle_broad_2019_rna_seq_mrna` (RPKM) | direct test of whether the programme is transcribed by malignant cells | 2026-10-23 |
| TCGA-SKCM, same visit | cBioPortal `skcm_tcga` | comparator for detection fractions | 2026-10-23 |

## Causal layer

Accessions below are exactly those stated in the manuscript's Data Availability section and were
verified against source records between 8 and 22 October 2026.

| Resource | Accession | Role |
|---|---|---|
| Leptin (BMI-adjusted) instrument | `ebi-a-GCST90007319` | exposure, 5 SNPs |
| Adiponectin instrument | `ieu-a-1` | exposure, 11 SNPs |
| Resistin instrument | `ebi-a-GCST90012034` | exposure, 13 SNPs |
| Body mass index | `ieu-b-40` (primary), `ebi-a-GCST90029007` (sensitivity) | negative-control exposure |
| Lymphocyte count | `ebi-a-GCST90025984` | mediation trait (Barton 2021, UK Biobank) |
| White blood cell count | `ebi-a-GCST90025985` | mediation trait (Barton 2021, UK Biobank) |
| Cutaneous melanoma (primary outcome) | `GCST90475578` (Million Veteran Program; 10,468 cases, 436,283 controls) | outcome |
| Replication outcome | `ieu-b-4969` (3,751 cases, 372,016 controls) | outcome |
| Sensitivity outcomes | `finn-b-C3_SKIN` (10,384 cases, 208,408 controls); `ieu-b-4959` (23,694 cases, 372,016 controls) | outcome |
| GTEx v8 | skin, subcutaneous adipose, visceral adipose, whole blood | cis-SMR and colocalisation |

**Sourcing note.** The primary outcome was obtained from the EBI GWAS Catalog FTP archive rather
than from OpenGWAS. The Landi 2020 melanoma meta-analysis was not obtained and dbGaP access was
not initiated; this is disclosed in the manuscript and in SAP v1.8. OpenGWAS tokens used during
retrieval expired on 2026-10-25.

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
