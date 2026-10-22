#!/usr/bin/env Rscript
# 指令十八① —— ALMP 双轴评分 × 免疫浸润耦合分析（三套独立签名框架一致性）
suppressMessages({library(GSVA)})
setwd('C:/Rtmp/p2t4')
w <- function(...) cat(sprintf(...), '\n', sep='')

## ---- 1. 合并表达（ALMP 305 + 免疫签名基因）----
read_expr <- function(f) {
  m <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE)
  list(ent=as.character(m$entrezGeneId), smp=m$sampleId, val=m$value)
}
a <- read_expr('mrna_data.tsv')
e <- read_expr('immune_mrna_data.tsv')
allent <- unique(c(a$ent, e$ent)); allsmp <- unique(c(a$smp, e$smp))
mat <- matrix(NA_real_, length(allent), length(allsmp), dimnames=list(allent, allsmp))
mat[cbind(match(c(a$ent, e$ent), allent), match(c(a$smp, e$smp), allsmp))] <- c(a$val, e$val)
w('combined matrix: %d genes x %d samples', nrow(mat), ncol(mat))
w('  immune-only genes added: %d', length(setdiff(e$ent, a$ent)))
keep <- rowSums(is.na(mat)) < 0.2 * ncol(mat)   # 基因级：缺失 <20% 保留
w('  genes kept (NA<20%%): %d / %d', sum(keep), nrow(mat))
mat <- mat[keep, , drop=FALSE]
for (i in seq_len(nrow(mat))) { v <- mat[i, ]; if (any(is.na(v))) mat[i, is.na(v)] <- median(v, na.rm=TRUE) }
mat <- log2(mat + 1)

## ---- 2. 基因集 ----
read_gmt <- function(f) { out <- list(); for (ln in readLines(f)) { p <- strsplit(ln,'\t')[[1]]; out[[p[1]]] <- unique(p[-c(1,2)]) }; out }
sym2ent <- read.delim('map_symbol2entrez_all.tsv', stringsAsFactors=FALSE)
s2e <- setNames(as.character(sym2ent$entrez), sym2ent$gene_symbol)
w('symbol->entrez map: %d symbols', nrow(sym2ent))
imm <- read_gmt('immune_signatures.gmt')
gs <- list()
for (nm in names(imm)) {
  e2 <- unique(s2e[imm[[nm]]]); e2 <- e2[!is.na(e2) & e2 %in% rownames(mat)]
  if (length(e2) >= 5) gs[[nm]] <- e2
}
w('\nimmune signatures usable: %d / %d', length(gs), length(imm))
for (nm in names(gs)) w('  %-42s %3d genes', nm, length(gs[[nm]]))

## ---- 3. 免疫评分（ssGSEA + GSVA）----
res_ss <- tryCatch(gsva(ssgseaParam(mat, gs, normalize=FALSE), verbose=FALSE), error=function(e){ message('ssGSEA ERR: ', conditionMessage(e)); NULL })
res_gv <- tryCatch(gsva(gsvaParam(mat, gs, kcdf='Gaussian'), verbose=FALSE), error=function(e){ message('GSVA ERR: ', conditionMessage(e)); NULL })
if (!is.null(res_ss)) w('ssGSEA: %s', paste(dim(res_ss), collapse=' x '))
if (!is.null(res_gv)) w('GSVA  : %s', paste(dim(res_gv), collapse=' x '))

## ---- 4. ALMP 双轴评分（复用 P2-T4 结果）----
almp <- list()
for (fx in c('gsva','ssgsea')) {
  f <- sprintf('P2T4_scores_%s.tsv', fx)
  if (!file.exists(f)) next
  s <- read.delim(f, stringsAsFactors=FALSE, check.names=FALSE); rownames(s) <- s$set; s$set <- NULL
  for (st in rownames(s)) almp[[paste0(fx, '::', st)]] <- setNames(as.numeric(s[st, ]), colnames(s))
}
w('\nALMP score sets: %d', length(almp))

## ---- 5. Spearman 耦合（ssGSEA 与 GSVA 两种算法 × 三套签名框架）----
comsmp <- Reduce(intersect, lapply(almp, names))
w('common samples: %d', length(comsmp))
res <- data.frame()
mk <- function(imm_mat, method) {
  for (im in rownames(imm_mat)) {
    for (am in names(almp)) {
      x <- imm_mat[im, comsmp]; y <- almp[[am]][comsmp]
      if (sum(!is.na(x)) < 3 || sum(!is.na(y)) < 3) next
      if (isTRUE(sd(x, na.rm=TRUE) == 0) || isTRUE(sd(y, na.rm=TRUE) == 0)) next
      ct <- suppressWarnings(cor.test(x, y, method='spearman'))
      res <<- rbind(res, data.frame(immune=im, almp=am, method=method,
        framework=sub('::.*','',im), almp_fx=sub('::.*','',am), almp_set=sub('.*::','',am),
        rho=unname(ct$estimate), p=ct$p.value, n=length(comsmp)))
    }
  }
}
if (!is.null(res_ss)) mk(res_ss, 'ssGSEA')
if (!is.null(res_gv)) mk(res_gv, 'GSVA')
res$q <- p.adjust(res$p, method='BH')
write.table(res, 'P2T4_coupling_spearman.tsv', sep='\t', quote=FALSE, row.names=FALSE)
w('→ P2T4_coupling_spearman.tsv  (rows=%d)', nrow(res))

## ---- 5b. 算法/框架稳健性：按 (ALMP评分 × 免疫签名) 汇总方向一致率 ----
res$key <- paste(res$almp, res$immune, sep=' | ')
dirsum <- aggregate(list(n_pos=res$rho), by=list(key=res$key, almp_set=res$almp_set),
                    FUN=function(v) sum(v>0))
dirsum$n_total <- as.numeric(table(res$key)[dirsum$key])
dirsum$frac_pos <- dirsum$n_pos / dirsum$n_total

## ---- 6. 框架/算法一致性汇总 ----
w('\n--- 按 ALMP 评分看各框架 × 算法 的方向一致性（rho 中位数） ---')
for (a1 in unique(res$almp_set)) {
  sub <- res[res$almp_set==a1, ]
  fw <- sapply(split(sub$rho, paste(sub$framework, sub$method, sep='/')),
               function(v) sprintf('%+.2f', median(v)))
  w('  %-30s %s', a1, paste(names(fw), fw, sep='=', collapse='  '))
}
w('\n--- 强耦合（|rho|>=0.5 且 q<0.05）条目数：%d / %d ---', sum(abs(res$rho)>=0.5 & res$q<0.05, na.rm=TRUE), nrow(res))
w('--- 名义显著（q<0.05）条目数：%d / %d ---', sum(res$q<0.05, na.rm=TRUE), nrow(res))

sink('P2T4_coupling_summary.txt')
w('ALMP × 免疫浸润耦合（Spearman），算法 ssGSEA + GSVA；框架 PanglaoDB / CellMarker / Hallmark')
w('样本 %d；免疫签名 %d；ALMP 评分 %d；条目 %d', length(comsmp), length(unique(res$immune)), length(unique(res$almp)), nrow(res))
w('')
w('显著（q<0.05）条目数：%d / %d', sum(res$q<0.05, na.rm=TRUE), nrow(res))
w('强耦合（|rho|>=0.5 且 q<0.05）：%d / %d', sum(abs(res$rho)>=0.5 & res$q<0.05, na.rm=TRUE), nrow(res))
w('')
w('--- 按 |rho| 降序前 30 ---')
print(head(res[order(-abs(res$rho)), c('immune','almp','method','rho','p','q')], 30), row.names=FALSE, digits=3)
w('')
w('--- 各 ALMP 评分 × 免疫签名：跨 6 格方向一致率（frac_pos，1=全正） ---')
ds <- dirsum[order(-dirsum$frac_pos), ]
print(head(ds, 30), row.names=FALSE, digits=3)
sink()
w('→ P2T4_coupling_summary.txt')
