#!/usr/bin/env Rscript
# 指令四十六 / SAP 条款 H：随机基因集阴性对照（主分析，合并一次 gsva 调用以节省预处理）
# 预登记参数见 SAP §2D（2026-10-21 18:13 CST）。
suppressMessages({library(GSVA)})
setwd('C:/Rtmp/p2t12')
dir.create('results', showWarnings = FALSE)
SEED <- 20261021; N_RAND <- 1000
sink('results/y46_run.log', split = TRUE)
w <- function(...) cat(sprintf(...), '\n', sep = '')
w('=== 指令四十六 随机基因集阴性对照（主分析）===')
w('R %s | GSVA %s | seed %d | N = %d', getRversion(), as.character(packageVersion('GSVA')), SEED, N_RAND)

## ---- 1. 全转录组矩阵 ----
m <- read.delim('full_mrna_data.tsv', stringsAsFactors = FALSE, check.names = FALSE)
genes <- unique(as.character(m$entrezGeneId)); samples <- unique(m$sampleId)
mat <- matrix(NA_real_, length(genes), length(samples), dimnames = list(genes, samples))
mat[cbind(match(m$entrezGeneId, genes), match(m$sampleId, samples))] <- m$value
w('矩阵: %d 基因 x %d 样本 | 缺失 %.2f%%', nrow(mat), ncol(mat), 100 * mean(is.na(mat)))
keep <- rowSums(is.na(mat)) < 0.2 * ncol(mat)
mat <- mat[keep, , drop = FALSE]
for (i in seq_len(nrow(mat))) { v <- mat[i, ]; if (any(is.na(v))) mat[i, is.na(v)] <- median(v, na.rm = TRUE) }
mat <- log2(mat + 1)
w('保留基因 %d | log2(x+1): median=%.2f', nrow(mat), median(mat))

## ---- 2. 基因集 ----
read_gmt <- function(f) { out <- list(); for (ln in readLines(f)) { p <- strsplit(ln, '\t')[[1]]; out[[p[1]]] <- unique(p[-c(1, 2)]) }; out }
s2e_df <- read.delim('map_symbol2entrez_all.tsv', stringsAsFactors = FALSE)
s2e <- setNames(as.character(s2e_df$entrez), s2e_df$gene_symbol)
to_ent <- function(x) { e <- unique(s2e[x]); e[!is.na(e) & e %in% rownames(mat)] }
ALMP <- lapply(read_gmt('ALMP_gene_set_v1.gmt'), to_ent)
IMM  <- lapply(read_gmt('immune_signatures.gmt'), to_ent); IMM <- IMM[sapply(IMM, length) >= 5]
TARGETS <- c('PanglaoDB::Neutrophils', 'PanglaoDB::Macrophages', 'Hallmark::Inflammatory Response')
w('ALMP 在场基因: %s', paste(names(ALMP), sapply(ALMP, length), sep = '=', collapse = ', '))

## ---- 3. 真实集 + 免疫签名（全转录组，一次调用）----
w('[3] 全转录组 GSVA / ssGSEA（真实集 + 27 免疫签名）...')
t0 <- Sys.time()
gv_real <- gsva(gsvaParam(mat, c(ALMP, IMM), kcdf = 'Gaussian'), verbose = FALSE)
ss_real <- gsva(ssgseaParam(mat, c(ALMP, IMM), normalize = TRUE), verbose = FALSE)
w('  完成 %d x %d，用时 %.1f s', nrow(gv_real), ncol(gv_real), as.numeric(difftime(Sys.time(), t0, units = 'secs')))

## ---- 4. 与归档（305 基因子集版）逐集对比 ----
old <- read.delim('P2T4_scores_gsva.tsv', stringsAsFactors = FALSE, check.names = FALSE)
rownames(old) <- old$set; old$set <- NULL
cmp <- data.frame()
for (st in intersect(rownames(old), names(ALMP))) {
  cs <- intersect(colnames(old), colnames(gv_real))
  ov <- suppressWarnings(as.numeric(old[st, cs]))
  cmp <- rbind(cmp, data.frame(
    set = st, n_used = length(ALMP[[st]]),
    n_NA_in_archived = sum(is.na(ov)),
    rho_archived_vs_full = if (all(is.na(ov))) NA else cor(ov, as.numeric(gv_real[st, cs]), method = 'spearman')
  ))
}
w('\n[4] 归档子集评分 vs 全转录组评分:')
print(cmp, row.names = FALSE, digits = 4)
write.table(cmp, 'results/y46_archived_vs_full_scores.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

## 真实耦合 rho：全转录组 vs 归档
oc <- read.delim('P2T4_coupling_spearman.tsv', stringsAsFactors = FALSE)
rt <- data.frame()
for (tgt in TARGETS) {
  for (st in names(ALMP)) {
    r_gv <- cor(gv_real[st, ], gv_real[tgt, ], method = 'spearman')
    r_ss <- cor(ss_real[st, ], ss_real[tgt, ], method = 'spearman')
    o_gv <- oc$rho[oc$immune == tgt & oc$almp == paste0('gsva::', st) & oc$method == 'GSVA']
    o_ss <- oc$rho[oc$immune == tgt & oc$almp == paste0('ssgsea::', st) & oc$method == 'ssGSEA']
    rt <- rbind(rt, data.frame(target = tgt, set = st, k = length(ALMP[[st]]),
      rho_full_GSVA = r_gv, rho_archived_GSVA = if (length(o_gv)) o_gv[1] else NA,
      rho_full_ssGSEA = r_ss, rho_archived_ssGSEA = if (length(o_ss)) o_ss[1] else NA))
  }
}
w('\n[4b] 真实 Spearman rho（全转录组 vs 归档）:')
print(rt, row.names = FALSE, digits = 4)
write.table(rt, 'results/y46_real_rho_full_vs_archived.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

## ---- 5. 表达分层 ----
gm <- rowMeans(mat); gvv <- apply(mat, 1, var)
dec <- cut(gm, quantile(gm, seq(0, 1, 0.1)), include.lowest = TRUE, labels = FALSE)
ter <- cut(gvv, quantile(gvv, seq(0, 1, 1/3)), include.lowest = TRUE, labels = FALSE)
names(dec) <- rownames(mat); names(ter) <- rownames(mat)
lvl30 <- paste(dec, ter, sep = '_'); names(lvl30) <- rownames(mat)
w('\n[5] 分层: 均值十分位(最小层 %d 基因) | 均值x方差 30 层(最小层 %d)', min(table(dec)), min(table(lvl30)))

draw <- function(real_genes, levels, n_draw) {
  tab <- table(levels[real_genes]); pool <- split(rownames(mat), levels)
  lapply(seq_len(n_draw), function(b) {
    p <- character(0)
    for (lv in names(tab)) {
      n <- as.integer(tab[[lv]]); cand <- pool[[lv]]
      if (is.null(cand) || !length(cand)) next
      p <- c(p, if (length(cand) >= n) sample(cand, n) else sample(cand, n, replace = TRUE))
    }
    unique(p)
  })
}

## ---- 6. 主分析：全部 4 个真实集 × 1000 随机集，一次 gsva 调用 ----
set.seed(SEED)
w('\n[6] 生成随机集 ...')
all_sets <- list(); meta <- data.frame()
for (st in names(ALMP)) {
  rs <- draw(ALMP[[st]], dec, N_RAND)
  nm <- sprintf('%s__rand%04d', st, seq_len(N_RAND))
  all_sets[nm] <- rs
  meta <- rbind(meta, data.frame(name = nm, real = st, k = length(ALMP[[st]]),
                                 rand_size = sapply(rs, length)))
}
w('  合计随机集 %d 个；规模校验: %s', length(all_sets),
  paste(tapply(meta$rand_size, meta$real, function(z) sprintf('%d-%d', min(z), max(z))), collapse = ' | '))
write.table(meta, 'results/y46_random_sets_meta.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

t0 <- Sys.time()
gv_rand <- gsva(gsvaParam(mat, all_sets, kcdf = 'Gaussian'), verbose = FALSE)
w('  1000x4 随机集 GSVA 用时 %.1f min', as.numeric(difftime(Sys.time(), t0, units = 'mins')))

## ---- 7. 统计 ----
res <- data.frame(); null_store <- list()
for (st in names(ALMP)) {
  nm <- meta$name[meta$real == st]
  M <- gv_rand[nm, , drop = FALSE]
  for (tgt in TARGETS) {
    rr <- apply(M, 1, function(z) cor(z, gv_real[tgt, ], method = 'spearman'))
    r0 <- cor(gv_real[st, ], gv_real[tgt, ], method = 'spearman')
    p_emp <- (1 + sum(rr >= r0)) / (1 + N_RAND)
    res <- rbind(res, data.frame(set = st, k = length(ALMP[[st]]), target = tgt, method = 'GSVA',
      match = 'mean_decile', rho_real = r0, null_mean = mean(rr), null_sd = sd(rr),
      null_p50 = median(rr), null_p95 = unname(quantile(rr, 0.95)), null_max = max(rr),
      pctile = 100 * mean(rr < r0), p_emp = p_emp, N = N_RAND))
    if (st == 'SEG1_ADIPOKINE_CORE') null_store[[tgt]] <- rr
    w('  %-26s x %-38s rho=%+.4f | null %+.4f±%.4f P95=%+.4f max=%+.4f | 分位 %5.1f%% | p_emp=%.4f',
      st, tgt, r0, mean(rr), sd(rr), unname(quantile(rr, 0.95)), max(rr), 100 * mean(rr < r0), p_emp)
  }
}
write.table(res, 'results/y46_random_control_main.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
if (length(null_store)) write.table(as.data.frame(null_store), 'results/y46_null_distribution_SEG1.tsv',
                                    sep = '\t', quote = FALSE, row.names = FALSE)
w('\n=> 主判据：SEG1_ADIPOKINE_CORE x PanglaoDB::Neutrophils (GSVA)  p_emp = %.4f',
  res$p_emp[res$set == 'SEG1_ADIPOKINE_CORE' & res$target == TARGETS[1]])
sink(); cat('DONE\n')
