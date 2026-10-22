#!/usr/bin/env Rscript
# 指令四十六 / SAP 条款 H 敏感性分析：
#   (S1) 更严的表达分布匹配 —— 均值十分位 × 方差三分位 交叉 30 层
#   (S2) 换算法 —— ssGSEA（与主分析口径一致）
# 仅针对主判据对（SEG1_ADIPOKINE_CORE 22 基因 × 3 个预指定签名）。
# 种子 20261021（与主分析同）。

suppressMessages({library(GSVA)})
setwd('C:/Rtmp/p2t12')
dir.create('results', showWarnings = FALSE)
SEED <- 20261021
N_RAND <- 1000
sink('results/y46_sensitivity_run.log', split = TRUE)
w <- function(...) cat(sprintf(...), '\n', sep = '')

m <- read.delim('full_mrna_data.tsv', stringsAsFactors = FALSE, check.names = FALSE)
genes <- unique(as.character(m$entrezGeneId)); samples <- unique(m$sampleId)
mat <- matrix(NA_real_, length(genes), length(samples), dimnames = list(genes, samples))
mat[cbind(match(m$entrezGeneId, genes), match(m$sampleId, samples))] <- m$value
keep <- rowSums(is.na(mat)) < 0.2 * ncol(mat)
mat <- mat[keep, , drop = FALSE]
for (i in seq_len(nrow(mat))) { v <- mat[i, ]; if (any(is.na(v))) mat[i, is.na(v)] <- median(v, na.rm = TRUE) }
mat <- log2(mat + 1)
w('矩阵: %d 基因 x %d 样本', nrow(mat), ncol(mat))

read_gmt <- function(f) { out <- list(); for (ln in readLines(f)) { p <- strsplit(ln, '\t')[[1]]; out[[p[1]]] <- unique(p[-c(1, 2)]) }; out }
s2e_df <- read.delim('map_symbol2entrez_all.tsv', stringsAsFactors = FALSE)
s2e <- setNames(as.character(s2e_df$entrez), s2e_df$gene_symbol)
to_ent <- function(x) { e <- unique(s2e[x]); e[!is.na(e) & e %in% rownames(mat)] }
ALMP <- lapply(read_gmt('ALMP_gene_set_v1.gmt'), to_ent)
IMM  <- lapply(read_gmt('immune_signatures.gmt'), to_ent); IMM <- IMM[sapply(IMM, length) >= 5]
TARGETS <- c('PanglaoDB::Neutrophils', 'PanglaoDB::Macrophages', 'Hallmark::Inflammatory Response')

gv <- gsva(gsvaParam(mat, c(ALMP, IMM), kcdf = 'Gaussian'), verbose = FALSE)
ss <- gsva(ssgseaParam(mat, c(ALMP, IMM), normalize = TRUE), verbose = FALSE)

gm <- rowMeans(mat); gv_ <- apply(mat, 1, var)
dec <- cut(gm, breaks = quantile(gm, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
ter <- cut(gv_, breaks = quantile(gv_, seq(0, 1, 1/3)), include.lowest = TRUE, labels = FALSE)
lvl30 <- paste(dec, ter, sep = '_'); names(lvl30) <- rownames(mat)
names(dec) <- rownames(mat)
w('30 层最小层基因数: %d', min(table(lvl30)))

draw <- function(real_genes, levels, n_draw) {
  tab <- table(levels[real_genes]); pool <- split(rownames(mat), levels)
  if (any(sapply(names(tab), function(k) length(pool[[k]]) < as.integer(tab[[k]])))) {
    w('  !! 有层基因数不足，将放回抽样')
  }
  lapply(seq_len(n_draw), function(b) {
    p <- character(0)
    for (lv in names(tab)) {
      n <- as.integer(tab[[lv]]); cand <- pool[[lv]]
      p <- c(p, if (length(cand) >= n) sample(cand, n) else sample(cand, n, replace = TRUE))
    }
    unique(p)
  })
}

SEG <- ALMP[['SEG1_ADIPOKINE_CORE']]
out <- data.frame()
for (scen in c('mean_decile', 'mean_x_var_30')) {
  set.seed(SEED)
  lv <- if (scen == 'mean_decile') dec else lvl30
  rs <- draw(SEG, lv, N_RAND)
  g1 <- gsva(gsvaParam(mat, rs, kcdf = 'Gaussian'), verbose = FALSE)
  s1 <- gsva(ssgseaParam(mat, rs, normalize = TRUE), verbose = FALSE)
  for (tgt in TARGETS) {
    for (mt in c('GSVA', 'ssGSEA')) {
      M <- if (mt == 'GSVA') g1 else s1; REF <- if (mt == 'GSVA') gv else ss
      rr <- apply(M, 1, function(z) cor(z, REF[tgt, ], method = 'spearman'))
      r0 <- cor(REF[['SEG1_ADIPOKINE_CORE'], ], REF[tgt, ], method = 'spearman')
      out <- rbind(out, data.frame(
        scenario = scen, method = mt, target = tgt,
        rho_real = r0, null_mean = mean(rr), null_p95 = quantile(rr, 0.95),
        pctile = 100 * mean(rr < r0), p_emp = (1 + sum(rr >= r0)) / (1 + N_RAND), N = N_RAND))
      w('%-14s %-7s %-40s rho=%+.4f null=%+.4f P95=%+.4f 分位=%5.1f%% p_emp=%.4f',
        scen, mt, tgt, r0, mean(rr), quantile(rr, 0.95), 100 * mean(rr < r0), (1 + sum(rr >= r0)) / (1 + N_RAND))
    }
  }
}
write.table(out, 'results/y46_random_control_sensitivity.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
sink(); w('DONE')
