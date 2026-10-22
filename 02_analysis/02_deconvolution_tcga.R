#!/usr/bin/env Rscript
# =============================================================================
# P2-T5b v2 : 多去卷积算法一致性（TCGA-SKCM 472 样本） 【尺度/矩阵范围修正版】
#
# v1 缺陷（已由对照实验证实）：
#   把「仅 250 个签名基因」的子集喂给 run_quantiseq -> fixMixture() 的
#   t(t(x)*1e6/colSums(x)) 以 250 基因的子集列和归一化 -> b 被放大上百倍
#   -> 约束最小二乘被顶到单纯形顶点 -> 每样本只有一个细胞类型=1，其余=0（角点解）
#   证据：对 quantiseqr 自带 demo 做对照——全矩阵输入 median=0.011/max=0.71，
#         仅签名子集输入 median=0/max=1.000、11 类中 9~10 类恒为 0。
#
# v2 修正（全转录组归一化铁律）：
#   b = 签名基因的线性 TPM ÷ 【全转录组列和】× 1e6，再直接调用内部核心
#   quantiseqr:::quanTIseq()（绕开 fixMixture 的子集重归一化）。
#   已在 demo 上与 run_quantiseq(全矩阵) 逐样本比对，max|Δ|=2.5e-3（等价）。
#   MCP-counter 的分位归一化同样改为在全转录组上做。
# =============================================================================
suppressMessages({library(quantiseqr); library(GSVA); library(preprocessCore); library(data.table)})
setwd('C:/Rtmp/deconv')
w <- function(...) cat(sprintf(...), '\n', sep='')
OUT <- 'C:/Rtmp/deconv/out_v2'; dir.create(OUT, showWarnings=FALSE)

## ================= 1. 签名基因 =================
til <- read.delim('C:/Rtmp/quantiseqr/inst/extdata/TIL10_signature.txt', check.names=FALSE)
til_genes <- til$ID
mcp <- fread('MCPcounter_genes.txt', header=TRUE); setnames(mcp, 1:3, c('symbol','population','entrez'))
mcp <- mcp[symbol != '' & population != '']
mcp_list <- split(mcp$symbol, mcp$population)
need <- unique(c(til_genes, mcp$symbol))
w('[1] TIL10 %d genes; MCP-counter %d pops (%d markers); union need %d',
  length(til_genes), length(mcp_list), nrow(mcp), length(need))

## ================= 2. TPM 矩阵（注意：gdc-hub star_tpm 实为 log2(TPM+1)）=================
tpm <- fread('TCGA-SKCM.star_tpm.tsv.gz'); setnames(tpm, 1, 'Ensembl_ID')
tpm[, ens := sub('\\..*$', '', Ensembl_ID)]
sc <- read.delim('C:/Rtmp/p2t4/P2T4_scores_gsva.tsv', check.names=FALSE)
cbio_ids <- colnames(sc)[-1]
tpm_smps <- setdiff(names(tpm), c('Ensembl_ID','ens'))
sel <- tpm_smps[substr(tpm_smps,1,15) %in% substr(cbio_ids,1,15)]
w('[2] TPM %d genes x %d samples; matched %d', nrow(tpm), length(tpm_smps), length(sel))

## ★ 全转录组列和（线性 TPM 尺度），分块计算以控内存
cs_full <- numeric(length(sel)); names(cs_full) <- sel
ch <- 5000
for (s in seq(1, nrow(tpm), by = ch)) {
  e <- min(s + ch - 1, nrow(tpm))
  blk <- as.matrix(tpm[s:e, ..sel])
  cs_full <- cs_full + colSums(2^blk - 1)
}
w('[2] full-transcriptome colSums (linear TPM): median=%.0f  min=%.0f  max=%.0f',
  median(cs_full), min(cs_full), max(cs_full))

## ================= 3. Ensembl -> symbol =================
smap <- fread('symbol_ensembl_map.tsv'); smap <- smap[!is.na(ensembl) & ensembl != '']
fix <- fread('symbol_alias_fix.tsv')
smap <- rbind(smap[, .(symbol, ensembl)], fix[, .(symbol, ensembl)])
m <- rbindlist(lapply(seq_len(nrow(smap)), function(i)
  data.table(ens = strsplit(smap$ensembl[i], '|', fixed=TRUE)[[1]], symbol = smap$symbol[i])))
m <- unique(m[ens %in% tpm$ens & symbol %in% need])
expr <- as.matrix(tpm[match(m$ens, tpm$ens), ..sel]); colnames(expr) <- sel
m[, mu := rowMeans(expr)]; setorder(m, symbol, -mu); m1 <- m[, .SD[1], by = symbol]
w('[3] resolved symbols: %d / %d needed', nrow(m1), length(need))
w('[3] missing: %s', paste(setdiff(need, m1$symbol), collapse=', '))
mat <- expr[match(m1$ens, m$ens), , drop=FALSE]; rownames(mat) <- m1$symbol
mat <- mat[!duplicated(rownames(mat)), , drop=FALSE]
w('[3] signature-gene matrix: %d genes x %d samples (log2 scale)', nrow(mat), ncol(mat))

## ================= 4. quanTIseq（★ 全转录组列和归一化后调用内部核心）=================
sigraw <- as.matrix(read.delim('C:/Rtmp/quantiseqr/inst/extdata/TIL10_signature.txt',
                               row.names=1, check.names=FALSE))
rmg <- as.vector(read.table('C:/Rtmp/quantiseqr/inst/extdata/TIL10_rmgenes.txt', header=FALSE, sep='\t')[,1])
abg <- as.vector(read.table('C:/Rtmp/quantiseqr/inst/extdata/TIL10_TCGA_aberrant_immune_genes.txt',
                            header=FALSE, sep='\t')[,1])
sig <- sigraw[!rownames(sigraw) %in% rmg, , drop=FALSE]
sig <- sig[!rownames(sig) %in% abg, , drop=FALSE]
msc <- read.table('C:/Rtmp/quantiseqr/inst/extdata/TIL10_mRNA_scaling.txt',
                  header=FALSE, sep='\t', stringsAsFactors=FALSE)
colnames(msc) <- c('celltype','scaling')
mRNA <- as.vector(as.matrix(msc$scaling[match(colnames(sig), msc$celltype)]))
w('[4] signature after removals: %d genes x %d cell types', nrow(sig), ncol(sig))

b <- mat                                                  # log2 尺度
b <- 2^b - 1                                              # -> 线性 TPM
b <- t(t(b) * 1e6 / cs_full)                              # ★ 全转录组列和归一化
res1 <- quantiseqr:::quanTIseq(sig, b, scaling = mRNA, method = 'lsei')
i <- which(colnames(sig) == 'T.cells.CD4')
res2 <- quantiseqr:::quanTIseq(sig[, -i], b, scaling = mRNA[-i], method = 'lsei')
ind <- which(res1[, 'Tregs'] < 0.02)
if (length(ind) > 0) {
  res1[ind, 'Tregs'] <- (res2[ind, 'Tregs'] + res1[ind, 'Tregs']) / 2
  res1[ind, 'T.cells.CD4'] <- pmax(0, res1[ind, 'T.cells.CD4'] -
                                     (res2[ind, 'Tregs'] + res1[ind, 'Tregs']) / 2)
}
qt <- as.data.frame(t(res1 / rowSums(res1)))              # cell types x samples
w('[4] quanTIseq: %d cell types x %d samples', nrow(qt), ncol(qt))
for (ct in rownames(qt)) { v <- as.numeric(qt[ct,])
  w('      %-18s med=%.4f mean=%.4f max=%.4f  ~1:%d  ==0:%d', ct, median(v), mean(v), max(v), sum(v>0.99), sum(v==0)) }
write.table(data.frame(cell_type = rownames(qt), qt, check.names = FALSE),
            file.path(OUT, 'P2T5b_quantiseq_tcga.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

## ================= 5. MCP-counter（★ 分位归一化在全转录组上做）=================
lg_full <- as.matrix(tpm[, ..sel]); rownames(lg_full) <- tpm$ens
lg_qn <- preprocessCore::normalize.quantiles(lg_full); dimnames(lg_qn) <- dimnames(lg_full)
rm(lg_full); gc(verbose=FALSE)
appendSignatures <- function(xp, markers) as.data.frame(do.call(cbind, lapply(markers, function(x)
  apply(xp[intersect(rownames(xp), x), , drop=FALSE], 2, mean, na.rm=TRUE))))
mcp_calc <- function(lgmat, markers, ensmap) {
  # markers 为 symbol 列表；映射到 ENSG 后取行
  ml <- lapply(markers, function(s) unique(ensmap$ens[ensmap$symbol %in% s]))
  names(ml) <- names(markers)
  as.data.frame(t(appendSignatures(lgmat, ml)))
}
mc    <- mcp_calc(lg_qn, mcp_list, m1)                    # 变体 A：全转录组分位归一化（QN）
mcraw <- as.data.frame(t(appendSignatures(mat, mcp_list))) # 变体 B：log2(TPM+1) 不做 QN（MCP-counter 常规 RNA-seq 输入）
g_ok <- sum(sapply(mcp_list, function(s) sum(m1$symbol %in% s)))
w('[5] MCP-counter: %d populations x %d samples (markers %d/%d)', nrow(mc), ncol(mc), g_ok, nrow(mcp))
write.table(data.frame(population = rownames(mc), mc, check.names = FALSE),
            file.path(OUT, 'P2T5b_mcpcounter_tcga_QN.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
write.table(data.frame(population = rownames(mcraw), mcraw, check.names = FALSE),
            file.path(OUT, 'P2T5b_mcpcounter_tcga_noQN.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
qn_chk <- data.frame(population = rownames(mc),
                     rho_QN_vs_noQN = round(sapply(rownames(mc), function(p)
                       suppressWarnings(cor(as.numeric(mc[p,]), as.numeric(mcraw[p,]), method='spearman'))), 4))
write.table(qn_chk, file.path(OUT, 'P2T5b_mcpcounter_QN_sensitivity.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
w('[5] QN vs no-QN rho: %s', paste(sprintf('%s=%.3f', qn_chk$population, qn_chk$rho_QN_vs_noQN), collapse=' '))

## ================= 6. P2-T4 签名评分对齐 =================
ss <- read.delim('C:/Rtmp/p2t4/P2T4_scores_ssgsea.tsv', check.names=FALSE)
gv <- read.delim('C:/Rtmp/p2t4/P2T4_scores_gsva.tsv', check.names=FALSE)
rownames(ss) <- ss$set; ss$set <- NULL
rownames(gv) <- gv$set; gv$set <- NULL
mapid <- setNames(substr(colnames(ss),1,15), colnames(ss))
key <- substr(colnames(mat),1,15)
ss2 <- ss[, match(key, mapid), drop=FALSE]; colnames(ss2) <- colnames(mat)
gv2 <- gv[, match(key, mapid), drop=FALSE]; colnames(gv2) <- colnames(mat)
w('[6] aligned P2-T4 scores: ssGSEA %d sets, GSVA %d sets', nrow(ss2), nrow(gv2))

## ================= 7. 跨方法一致性 =================
common <- data.frame()
pairs <- list(c('B.cells','B lineage'), c('Macrophages.M1','Monocytic lineage'),
              c('Macrophages.M2','Monocytic lineage'), c('Monocytes','Monocytic lineage'),
              c('Neutrophils','Neutrophils'), c('NK.cells','NK cells'),
              c('T.cells.CD8','CD8 T cells'), c('T.cells.CD4','T cells'),
              c('Tregs','T cells'), c('Dendritic.cells','Myeloid dendritic cells'))
mcp_vars <- list(noQN = mcraw, QN = mc)   # 主 = noQN（MCP-counter 常规 RNA-seq 输入）
for (vn in names(mcp_vars)) {
  mv <- mcp_vars[[vn]]
  for (p in pairs) {
    if (!(p[1] %in% rownames(qt)) || !(p[2] %in% rownames(mv))) next
    x <- as.numeric(qt[p[1],]); y <- as.numeric(mv[p[2],])
    if (sum(is.finite(x)) < 5 || sum(is.finite(y)) < 5 || sd(x)==0 || sd(y)==0) {
      common <- rbind(common, data.frame(mcp_variant=vn, method_a='quanTIseq', type_a=p[1],
                                         method_b='MCPcounter', type_b=p[2], rho=NA, p=NA)); next }
    ct <- suppressWarnings(cor.test(x, y, method='spearman'))
    common <- rbind(common, data.frame(mcp_variant=vn, method_a='quanTIseq', type_a=p[1],
                                       method_b='MCPcounter', type_b=p[2],
                                       rho=round(unname(ct$estimate),3), p=signif(ct$p.value,4)))
  }
}
common$q <- signif(p.adjust(common$p, 'BH'), 4)
w('\n[7a] quanTIseq vs MCP-counter（主 noQN / 敏感性 QN）:')
print(common[common$mcp_variant=='noQN', c('type_a','type_b','rho','p','q')], row.names=FALSE)
write.table(common, file.path(OUT, 'P2T5b_crossmethod_celltype.tsv'), sep='\t', quote=FALSE, row.names=FALSE)

keys <- list()
for (t in rownames(qt)) keys[[paste0('quanTIseq|', t)]] <- as.numeric(qt[t,])
for (t in rownames(mcraw)) keys[[paste0('MCPcounter|', t)]]    <- as.numeric(mcraw[t,])   # 主
for (t in rownames(mc))    keys[[paste0('MCPcounter_QN|', t)]] <- as.numeric(mc[t,])      # 敏感性
sig_pat <- c('Neutrophil','Macrophage','Monocyt','Myeloid','Dendritic','Treg','CD8','B cell','NK','Cytotoxic','T cell','Inflammatory')
for (s in rownames(ss2)) if (any(sapply(sig_pat, function(p) grepl(p, s, ignore.case=TRUE)))) keys[[paste0('ssGSEA|', s)]] <- as.numeric(ss2[s,])
for (s in rownames(gv2)) if (any(sapply(sig_pat, function(p) grepl(p, s, ignore.case=TRUE)))) keys[[paste0('GSVA|', s)]] <- as.numeric(gv2[s,])
K <- as.data.frame(keys); K <- K[, sapply(K, function(x) sd(x, na.rm=TRUE) > 0), drop=FALSE]
cm <- cor(K, method='spearman', use='pairwise.complete.obs')
write.table(data.frame(feature=rownames(cm), cm, check.names=FALSE),
            file.path(OUT, 'P2T5b_crossmethod_corr.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
w('[7b] cross-method correlation matrix: %d features', ncol(K))

coup <- data.frame()
for (nm in names(keys)) {
  x <- as.numeric(keys[[nm]])
  for (ax in rownames(ss2)) for (alg in c('ssGSEA','GSVA')) {
    y <- as.numeric(if (alg=='ssGSEA') ss2[ax,] else gv2[ax,])
    if (sum(is.finite(x)) < 5 || sum(is.finite(y)) < 5) next
    if (sd(x)==0 || sd(y, na.rm=TRUE)==0) next
    ct <- suppressWarnings(cor.test(x, y, method='spearman'))
    coup <- rbind(coup, data.frame(feature=nm, axis=ax, almp_method=alg,
                                   rho=round(unname(ct$estimate),3), p=signif(ct$p.value,4)))
  }
}
coup$q <- signif(p.adjust(coup$p, 'BH'), 4)
write.table(coup, file.path(OUT, 'P2T5b_seg_coupling_by_method.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
w('\n[7c] 免疫估计 x SEG1_ADIPOKINE_CORE（|rho|>=0.3 且 q<0.05）')
pr <- coup[coup$axis=='SEG1_ADIPOKINE_CORE' & abs(coup$rho)>=0.3 & coup$q<0.05, ]
pr <- pr[order(-abs(pr$rho)), ]; print(pr, row.names=FALSE)

saveRDS(list(qt=qt, mc=mc, mcraw=mcraw, ss=ss2, gv=gv2, cross=common, coup=coup, cs_full=cs_full),
        file.path(OUT, 'P2T5b_all.rds'))
writeLines(capture.output(sessionInfo()), file.path(OUT, 'P2T5b_sessionInfo.txt'))
w('\nDONE')
