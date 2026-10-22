#!/usr/bin/env Rscript
# P2-T4: ALMP 基因集 GSVA / ssGSEA 评分（TCGA-SKCM）
# 输入: mrna_data.tsv (entrez × sample × value), ALMP_gene_set_v1.gmt, ALMP_ID_mapping_v1.csv
# 输出: P2T4_scores_ssgsea.tsv, P2T4_scores_gsva.tsv, P2T4_cutpoints.tsv, P2T4_summary.txt
suppressMessages({library(GSVA)})
setwd('C:/Rtmp/p2t4')
w <- function(...) cat(sprintf(...), '\n', sep='')

## ---- 1. 表达矩阵 ----
m <- read.delim('mrna_data.tsv', stringsAsFactors=FALSE, check.names=FALSE)
ent <- as.character(m$entrezGeneId); smp <- m$sampleId
genes <- unique(ent); samples <- unique(smp)
mat <- matrix(NA_real_, length(genes), length(samples), dimnames=list(genes, samples))
mat[cbind(match(ent, genes), match(smp, samples))] <- m$value
w('expression: %d genes x %d samples', nrow(mat), ncol(mat))
w('  NA cells: %d (%.2f%%)', sum(is.na(mat)), 100*mean(is.na(mat)))
mat <- log2(mat + 1)          # RSEM-like → log2(x+1)
w('  after log2(x+1): min=%.3f median=%.3f max=%.3f', min(mat,na.rm=TRUE),
  median(mat,na.rm=TRUE), max(mat,na.rm=TRUE))

## ---- 2. 基因集（symbol → entrez）----
map <- read.delim('map_symbol2entrez.tsv', stringsAsFactors=FALSE, check.names=FALSE)
sym2ent <- setNames(as.character(map$entrez), as.character(map$gene_symbol))
w('mapping: %d symbols', length(sym2ent))
gs <- list(); gs_meta <- data.frame()
for (ln in readLines('ALMP_gene_set_v1.gmt')) {
  p <- strsplit(ln, '\t')[[1]]
  name <- p[1]; syms <- unique(p[-c(1,2)])
  e <- unique(sym2ent[syms]); e <- e[!is.na(e)]; e <- e[e %in% rownames(mat)]
  gs[[name]] <- e
  gs_meta <- rbind(gs_meta, data.frame(set=name, n_gmt=length(syms), n_used=length(e)))
}
w('\n--- gene sets ---'); print(gs_meta, row.names=FALSE)

## ---- 3. GSVA / ssGSEA ----
run_safe <- function(expr) tryCatch(expr, error=function(e){ message('ERR: ', conditionMessage(e)); NULL })
w('\n--- ssGSEA ---')
res_ss <- run_safe(gsva(ssgseaParam(mat, gs, normalize=TRUE), verbose=FALSE))
if (is.null(res_ss)) { w('  retry ssGSEA normalize=FALSE'); res_ss <- run_safe(gsva(ssgseaParam(mat, gs, normalize=FALSE), verbose=FALSE)) }
if (is.null(res_ss)) { w('  retry ssGSEA alpha=0.25 normalize=FALSE'); res_ss <- run_safe(gsva(ssgseaParam(mat, gs, alpha=0.25, normalize=FALSE), verbose=FALSE)) }
if (!is.null(res_ss)) w('  ssGSEA dim: %s  NA=%d', paste(dim(res_ss), collapse=' x '), sum(is.na(res_ss)))

w('--- GSVA (Gaussian kernel) ---')
res_gv <- run_safe(gsva(gsvaParam(mat, gs, kcdf='Gaussian'), verbose=FALSE))
if (!is.null(res_gv)) w('  GSVA dim: %s  NA=%d', paste(dim(res_gv), collapse=' x '), sum(is.na(res_gv)))

## ---- 4. 切点：中位数 + 三分位 ----
cutpoints <- function(res, tagname) {
  out <- data.frame()
  for (s in rownames(res)) {
    x <- res[s, ]; x <- x[!is.na(x)]
    if (length(x) < 10) { w('  [skip] %s/%s (n=%d)', tagname, s, length(x)); next }
    med <- median(x)
    q <- quantile(x, c(1/3, 2/3), names=FALSE)
    out <- rbind(out, data.frame(method=tagname, set=s,
      min=min(x), q1=unname(q[1]), median=med, q3=unname(q[2]), max=max(x),
      n_high_median=sum(x > med), n_low_median=sum(x <= med),
      n_low=sum(x <= q[1]), n_mid=sum(x > q[1] & x <= q[2]), n_high=sum(x > q[2])))
  }
  out
}
cp <- rbind(cutpoints(res_gv, 'GSVA'), if (!is.null(res_ss)) cutpoints(res_ss, 'ssGSEA') else NULL)
w('\n--- cut points ---'); print(cp[,c('method','set','min','median','max','n_low','n_mid','n_high')], row.names=FALSE)

## ---- 5. 输出 ----
fmt <- function(x) formatC(x, format='g', digits=8)
write_tsv <- function(res, file) {
  d <- data.frame(set=rownames(res), res, check.names=FALSE)
  write.table(d, file, sep='\t', quote=FALSE, row.names=FALSE)
}
write_tsv(res_gv, 'P2T4_scores_gsva.tsv')
if (!is.null(res_ss)) write_tsv(res_ss, 'P2T4_scores_ssgsea.tsv')
write.table(cp, 'P2T4_cutpoints.tsv', sep='\t', quote=FALSE, row.names=FALSE)
write.table(gs_meta, 'P2T4_geneset_meta.tsv', sep='\t', quote=FALSE, row.names=FALSE)

sink('P2T4_summary.txt')
w('P2-T4 GSVA/ssGSEA 评分汇总')
w('样本: %d  (TCGA-SKCM rna_seq_v2_mrna)', ncol(mat))
w('基因: %d  (ALMP 308 中的表达可用子集)', nrow(mat))
w('表达变换: log2(RSEM+1)')
w('')
w('基因集（使用到的基因数 / GMT 符号数）：')
print(gs_meta, row.names=FALSE)
w('')
w('GSVA  得分范围：'); print(summary(as.vector(res_gv)))
if (!is.null(res_ss)) { w('ssGSEA 得分范围：'); print(summary(as.vector(res_ss))) }
w('')
w('切点统计：'); print(cp, row.names=FALSE)
sink()
w('\nDONE → P2T4_scores_gsva.tsv%s / P2T4_cutpoints.tsv',
  if (!is.null(res_ss)) ' + P2T4_scores_ssgsea.tsv' else '')
