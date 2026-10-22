#!/usr/bin/env Rscript
# =============================================================================
# P2-T6b：ICB 转化层队列扩面分析（指令二十一）
# -----------------------------------------------------------------------------
# 队列（cBioPortal 协调版，统一子集管线）：
#   Gide2019 / Liu2019 / Riaz2017 / UCLA2016   —— 主合并
#   VanAllen2015                                —— 无应答字段 → 仅 TMB/表达
#   （GSE78220 = UCLA/Hugo 同源；GSE91061 = BMS038）→ 独立一致性检验，不并入合并
#
# ★ 纪律 1：统一应答定义 = CR/PR → R（应答）；SD/PD → NR（无应答）
# ★ 纪律 2：quanTIseq 需全转录组 → 本层不跑（子集输入致角点解，见 P2-T5b 铁律）
# ★ 纪律 3：MCP-counter 仅标记均值；**不做子集 QN**（QN 必须在全转录组上做）
# ★ 纪律 4：尺度不得猜测 —— 先打印各队列数值范围，线性者显式 log2(x+1)
# =============================================================================
suppressMessages({library(GSVA); library(data.table); library(metafor)})
setwd('C:/Rtmp/icb_v2')
OUT <- 'C:/Rtmp/icb_v2/out'; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
DATA <- 'C:/Rtmp/icb_v2/data'
w <- function(...) cat(sprintf(...), '\n', sep = '')

## ================= 0. 基因集 / 标记 =================
gmt <- 'C:/Rtmp/p2t4/ALMP_gene_set_v1.gmt'
almp <- list()
for (ln in readLines(gmt)) {
  p <- strsplit(ln, '\t')[[1]]; almp[[p[1]]] <- unique(p[-c(1, 2)])
}
mp <- fread(file.path(DATA, 'symbol_entrez_map_union.tsv'))
mp <- mp[symbol != '' & entrez != '']
ez2sym <- setNames(mp$symbol, as.character(mp$entrez))
mcp_list <- split(mp$symbol[mp$mcp_population != ''], mp$mcp_population[mp$mcp_population != ''])
w('[0] ALMP 集 %d 个：%s', length(almp), paste(sprintf('%s(%d)', names(almp), sapply(almp, length)), collapse = ' '))
w('[0] MCP 群体 %d 类，标记 %d 个', length(mcp_list), length(unlist(mcp_list)))

## ================= 1. 读入表达式（entrez→symbol） =================
read_cohort <- function(tag) {
  d <- fread(file.path(DATA, paste0(tag, '_expr_union.tsv')), header = TRUE)
  setnames(d, 1, 'entrez')
  d[, entrez := as.character(entrez)]
  rn <- ez2sym[d$entrez]
  m <- as.matrix(d[, -1, drop = FALSE])
  ok <- !is.na(rn)
  m <- m[ok, , drop = FALSE]; rr <- rn[ok]
  mu <- rowMeans(m, na.rm = TRUE); o <- order(rr, -mu)
  m <- m[o, , drop = FALSE]; rr <- rr[o]
  k <- !duplicated(rr)
  m <- m[k, , drop = FALSE]; rownames(m) <- rr[k]
  # ★ 尺度自检
  rg <- range(m, na.rm = TRUE)
  logd <- if (max(m, na.rm = TRUE) > 30) { m <- log2(m + 1); TRUE } else FALSE
  w('[1] %-14s %4d 基因 x %4d 样本 | 原始范围 [%.2f, %.2f] | log2(+1) 变换=%s | 变换后 [%.2f, %.2f]',
    tag, nrow(m), ncol(m), rg[1], rg[2], logd, min(m, na.rm = TRUE), max(m, na.rm = TRUE))
  m
}
COH <- c('Gide2019', 'Liu2019', 'Riaz2017', 'UCLA2016', 'VanAllen2015')
MATS <- lapply(COH, read_cohort); names(MATS) <- COH

## ================= 2. 评分 =================
appendSignatures <- function(xp, markers) as.data.frame(do.call(cbind, lapply(markers, function(x)
  apply(xp[intersect(rownames(xp), x), , drop = FALSE], 2, mean, na.rm = TRUE))))

score_cohort <- function(mat, tag) {
  gs <- lapply(almp, function(s) intersect(s, rownames(mat)))
  w('[2] %-14s 覆盖: %s', tag, paste(sprintf('%s %d/%d', names(gs), sapply(gs, length), sapply(almp, length)), collapse = ' | '))
  gs <- gs[sapply(gs, length) >= 3]
  ss <- tryCatch(gsva(ssgseaParam(mat, gs, normalize = TRUE), verbose = FALSE), error = function(e) { w('   ssGSEA ERR %s', conditionMessage(e)); NULL })
  gv <- tryCatch(gsva(gsvaParam(mat, gs, kcdf = 'Gaussian'), verbose = FALSE), error = function(e) { w('   GSVA ERR %s', conditionMessage(e)); NULL })
  mc <- as.data.frame(t(appendSignatures(mat, mcp_list)))
  list(ss = ss, gv = gv, mc = mc)
}
SC <- lapply(COH, function(t) score_cohort(MATS[[t]], t)); names(SC) <- COH

## ================= 3. 应答统一（CR/PR vs SD/PD） =================
harm <- function(tag) {
  cl <- fread(file.path(DATA, paste0(tag, '_clinical.tsv')), header = TRUE)
  setnames(cl, 1, 'sampleId')
  bin_of <- function(v) {
    v <- trimws(as.character(v))
    ifelse(v %in% c('Complete Response', 'Partial Response', 'Responder', 'TRUE', 'True'), 'R',
    ifelse(v %in% c('Stable Disease', 'Progressive Disease', 'Non-responder', 'FALSE', 'False'), 'NR', NA))
  }
  if (tag == 'UCLA2016') {
    cl[, bin := bin_of(TREATMENT_RESPONSE)]; cl[, keep := !is.na(bin)]
  } else if (!'RESPONSE' %in% names(cl)) {
    w('[3] %-14s ★ 无 RESPONSE 字段 → 仅 TMB/表达，不入应答关联', tag)
    cl[, bin := NA_character_]; cl[, keep := FALSE]
  } else {
    cl[, bin := bin_of(RESPONSE)]
    if (!'SAMPLE_TREATMENT' %in% names(cl) || all(is.na(cl$SAMPLE_TREATMENT)) || all(cl$SAMPLE_TREATMENT == '')) {
      cl[, keep := !is.na(bin)]
    } else {
      cl[, keep := !is.na(bin) & (SAMPLE_TREATMENT %in% c('Pre', 'pre'))]
    }
  }
  cl[]
}
CL <- lapply(COH, harm); names(CL) <- COH
for (t in COH) {
  cl <- CL[[t]]
  if (all(is.na(cl$bin))) { w('[3] %-14s 无应答字段 → 排除出应答关联', t); next }
  w('[3] %-14s 可分析 %3d 例（R %d / NR %d）｜总样本 %d',
    t, sum(cl$keep, na.rm = TRUE), sum(cl$keep & cl$bin == 'R', na.rm = TRUE),
    sum(cl$keep & cl$bin == 'NR', na.rm = TRUE), nrow(cl))
}

## ================= 4. 检验 =================
wtest <- function(m, grp) {
  out <- data.frame()
  for (s in rownames(m)) {
    v <- as.numeric(m[s, ])
    x <- v[grp == 'R']; y <- v[grp == 'NR']
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x) < 3 || length(y) < 3 || isTRUE(sd(c(x, y)) == 0)) next
    wt <- suppressWarnings(wilcox.test(x, y))
    out <- rbind(out, data.frame(feature = s, n_R = length(x), n_NR = length(y),
      med_R = round(median(x), 3), med_NR = round(median(y), 3),
      delta = round(median(x) - median(y), 3), cohen_d = round((median(x) - median(y)) / sd(c(x, y)), 3),
      auc = round(unname(wt$statistic) / (length(x) * length(y)), 3),
      p = signif(wt$p.value, 4)))
  }
  if (nrow(out)) out$q <- signif(p.adjust(out$p, 'BH'), 4)
  out
}

main <- c('Gide2019', 'Liu2019', 'Riaz2017', 'UCLA2016')
RES <- list()
for (t in main) {
  cl <- CL[[t]][keep == TRUE]
  cl <- cl[sampleId %in% colnames(MATS[[t]])]
  grp <- setNames(cl$bin, cl$sampleId)
  for (rd in c('ssGSEA', 'GSVA', 'MCPcounter')) {
    m <- SC[[t]][[if (rd == 'ssGSEA') 'ss' else if (rd == 'GSVA') 'gv' else 'mc']]
    m <- m[, intersect(colnames(m), names(grp)), drop = FALSE]
    r <- wtest(m, grp[colnames(m)])
    if (nrow(r)) RES[[paste(t, rd)]] <- cbind(cohort = t, readout = rd, r)
  }
}
allres <- do.call(rbind, RES)
allres$q_pool <- signif(p.adjust(allres$p, 'BH'), 4)
fwrite(allres, file.path(OUT, 'P2T6b_response_assoc.tsv'), sep = '\t')

AXES <- c('SEG1_ADIPOKINE_CORE', 'SEG2_LIPID_PATHWAY', 'SEG3_EXECUTOR',
          'ALMP_COMPOSITIONAL_AXIS', 'ALMP_TUMOR_INTRINSIC_AXIS', 'ALMP_ALL_308')
axis_rows <- allres[allres$feature %in% AXES, ]
w('\n===== 主轴 × 应答（全部 %d 条主检验）=====', nrow(axis_rows))
print(axis_rows[, c('cohort','readout','feature','n_R','n_NR','med_R','med_NR','delta','cohen_d','p','q')], row.names = FALSE)
w('\n===== 任一 q<0.05 者 =====')
z <- allres[allres$q < 0.05, ]
if (nrow(z)) print(z, row.names = FALSE) else w('（无）')
w('\n===== 全表最小 p =====')
print(head(allres[order(allres$p), c('cohort','readout','feature','cohen_d','auc','p','q')], 10), row.names = FALSE)

## ================= 5. 跨队列合并（metafor，REML + HKSJ）=================
meta_in <- allres[!is.na(allres$cohen_d), ]
meta_out <- data.frame()
for (rd in unique(meta_in$readout)) for (ft in unique(meta_in$feature)) {
  s <- meta_in[meta_in$readout == rd & meta_in$feature == ft, ]
  if (nrow(s) < 2) next
  # 由中位数差与 Cohen's d 反推 SE(d)（近似）
  s$se_d <- abs(s$cohen_d) / pmax(abs(qnorm(s$p / 2)), 1e-6)
  s <- s[is.finite(s$se_d) & s$se_d > 0, ]
  if (nrow(s) < 2) next
  # ★ 必须显式写 s$cohen_d / s$se_d：裸列名在全局环境不存在 → tryCatch 会静默吞掉
  fit <- tryCatch(rma(yi = s$cohen_d, sei = s$se_d, method = 'REML', test = 'knha'),
                  error = function(e) { w('   meta ERR [%s %s] %s', rd, ft, conditionMessage(e)); NULL })
  if (is.null(fit)) next
  meta_out <- rbind(meta_out, data.frame(readout = rd, feature = ft, k = fit$k,
    d_pooled = round(fit$b, 3), ci_lo = round(fit$ci.lb, 3), ci_hi = round(fit$ci.ub, 3),
    p = signif(fit$pval, 4), I2 = round(fit$I2, 1), Qp = signif(fit$QEp, 4),
    n_neg = sum(s$cohen_d < 0), n_pos = sum(s$cohen_d > 0),
    cohorts = paste(s$cohort, collapse = '+')))
}
if (nrow(meta_out)) meta_out$q <- signif(p.adjust(meta_out$p, 'BH'), 4)
fwrite(meta_out, file.path(OUT, 'P2T6b_pooled_meta.tsv'), sep = '\t')
w('\n===== 跨队列合并（REML + HKSJ）—— 全部 %d 条 =====', nrow(meta_out))
print(meta_out[order(meta_out$p), ], row.names = FALSE)
w('\n===== 跨队列合并 —— 仅 6 条 ALMP 轴 =====')
print(meta_out[meta_out$feature %in% AXES, ], row.names = FALSE)

## ================= 5b. SEG1_ADIPOKINE_CORE × MCP 群体（生态位耦合，ICB 队列内）==========
coup2 <- data.frame()
for (t in main) {
  m <- SC[[t]]$mc
  ax <- SC[[t]]$ss['SEG1_ADIPOKINE_CORE', ]
  sh <- intersect(colnames(m), names(ax))
  for (ft in rownames(m)) {
    v <- as.numeric(m[ft, sh]); a <- as.numeric(ax[sh])
    ok <- is.finite(v) & is.finite(a)
    if (sum(ok) < 10 || isTRUE(sd(v[ok]) == 0)) next
    ct <- suppressWarnings(cor.test(a[ok], v[ok], method = 'spearman'))
    coup2 <- rbind(coup2, data.frame(cohort = t, feature = ft, n = sum(ok),
      rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 4)))
  }
}
if (nrow(coup2)) coup2$q <- signif(p.adjust(coup2$p, 'BH'), 4)
fwrite(coup2, file.path(OUT, 'P2T6b_seg1_mcp_coupling.tsv'), sep = '\t')
w('\n===== SEG1_ADIPOKINE_CORE × MCP 群体（Spearman，q<0.05）=====')
print(coup2[coup2$q < 0.05, ], row.names = FALSE)
w('（SEG1×MCP 检验 %d 条；q<0.05 者 %d）', nrow(coup2), sum(coup2$q < 0.05))

## ================= 6. TMB =================
tmb <- data.frame()
for (t in c('Liu2019', 'Riaz2017', 'UCLA2016', 'VanAllen2015')) {
  cl <- CL[[t]]
  if (!'TMB_NONSYNONYMOUS' %in% names(cl)) { w('[6] %s 无 TMB 字段', t); next }
  x <- suppressWarnings(as.numeric(cl$TMB_NONSYNONYMOUS))
  tmb <- rbind(tmb, data.frame(cohort = t, sampleId = cl$sampleId, tmb = x,
                               bin = cl$bin, keep = cl$keep))
}
fwrite(tmb, file.path(OUT, 'P2T6b_tmb.tsv'), sep = '\t')
w('\n===== TMB 描述（队列内）=====')
tt <- tmb[is.finite(tmb$tmb), ]
print(as.data.table(tt)[, .(n = .N, med = round(median(tmb), 2), p25 = round(quantile(tmb,.25), 2),
                            p75 = round(quantile(tmb,.75), 2), max = round(max(tmb), 1)), by = cohort], row.names = FALSE)

# 6b. ALMP 轴 × TMB 相关性（是否混杂）
tmbcor <- data.frame()
for (t in c('Liu2019', 'Riaz2017', 'UCLA2016')) {
  cl <- CL[[t]]
  x <- setNames(suppressWarnings(as.numeric(cl$TMB_NONSYNONYMOUS)), cl$sampleId)
  for (rd in c('ssGSEA', 'GSVA')) {
    m <- SC[[t]][[if (rd == 'ssGSEA') 'ss' else 'gv']]
    for (ft in intersect(c('SEG1_ADIPOKINE_CORE','SEG2_LIPID_PATHWAY','SEG3_EXECUTOR','ALMP_ALL_308'), rownames(m))) {
      v <- setNames(as.numeric(m[ft, ]), colnames(m))
      sh <- intersect(names(x), names(v)); sh <- sh[is.finite(x[sh]) & is.finite(v[sh])]
      if (length(sh) < 8) next
      ct <- suppressWarnings(cor.test(x[sh], v[sh], method = 'spearman'))
      tmbcor <- rbind(tmbcor, data.frame(cohort = t, readout = rd, feature = ft, n = length(sh),
        rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 4)))
    }
  }
}
if (nrow(tmbcor)) { tmbcor$q <- signif(p.adjust(tmbcor$p, 'BH'), 4) }
fwrite(tmbcor, file.path(OUT, 'P2T6b_axis_vs_tmb.tsv'), sep = '\t')
w('\n===== ALMP 轴 × TMB（Spearman）=====')
print(tmbcor, row.names = FALSE)

# 6c. 低 TMB 亚组：主切点 = 队列内中位数；敏感性 = 10 mut/Mb
lowres <- data.frame()
for (t in c('Liu2019', 'Riaz2017', 'UCLA2016')) {
  cl <- CL[[t]][keep == TRUE]
  cl <- cl[sampleId %in% colnames(MATS[[t]])]
  tv <- suppressWarnings(as.numeric(cl$TMB_NONSYNONYMOUS))
  cut_med <- median(tv, na.rm = TRUE)
  for (cutname in c('median', 'abs10')) {
    cutv <- if (cutname == 'median') cut_med else 10
    for (side in c('low', 'high')) {
      sel <- if (side == 'low') which(tv < cutv) else which(tv >= cutv)
      if (length(sel) < 8) next
      sub <- cl[sel]; grp <- setNames(sub$bin, sub$sampleId)
      if (sum(grp == 'R') < 3 || sum(grp == 'NR') < 3) next
      for (rd in c('ssGSEA', 'GSVA')) {
        m <- SC[[t]][[if (rd == 'ssGSEA') 'ss' else 'gv']]
        m <- m[, intersect(colnames(m), names(grp)), drop = FALSE]
        r <- wtest(m, grp[colnames(m)])
        if (nrow(r)) lowres <- rbind(lowres, data.frame(cohort = t, cut = cutname, side = side,
          cut_value = round(cutv, 2), readout = rd, r))
      }
    }
  }
}
if (nrow(lowres)) lowres$q <- signif(p.adjust(lowres$p, 'BH'), 4)
fwrite(lowres, file.path(OUT, 'P2T6b_tmb_subgroup.tsv'), sep = '\t')
w('\n===== 低 TMB 亚组 × 应答（仅列主轴，q<0.1）=====')
lr <- lowres[lowres$feature %in% c('SEG1_ADIPOKINE_CORE','SEG2_LIPID_PATHWAY','SEG3_EXECUTOR','ALMP_ALL_308'), ]
print(lr[lr$q < 0.1, c('cohort','cut','side','readout','feature','n_R','n_NR','med_R','med_NR','cohen_d','p','q')], row.names = FALSE)
w('（低 TMB 亚组检验总数 %d；q<0.05 者 %d）', nrow(lowres), sum(lowres$q < 0.05))

saveRDS(list(allres = allres, meta = meta_out, tmb = tmb, tmbcor = tmbcor, lowres = lowres,
             clin = CL, scores = SC), file.path(OUT, 'P2T6b_all.rds'))
writeLines(capture.output(sessionInfo()), file.path(OUT, 'P2T6b_sessionInfo.txt'))
w('\nDONE')

## ================= 7. 每队列 × 每轴 一览表（供报告直接用）=================
peraxis <- allres[allres$feature %in% AXES, c('cohort','readout','feature','n_R','n_NR',
                                             'med_R','med_NR','delta','cohen_d','auc','p','q')]
peraxis <- peraxis[order(peraxis$feature, peraxis$readout, peraxis$cohort), ]
fwrite(peraxis, file.path(OUT, 'P2T6b_axis_by_cohort.tsv'), sep = '\t')
w('\n===== 每队列 × 每轴 一览（%d 行）已写出 P2T6b_axis_by_cohort.tsv =====', nrow(peraxis))
w('★ 主检验族（ALMP 轴 × 应答）: 总 %d 条, p<0.05 者 %d, q<0.05 者 %d',
  nrow(peraxis), sum(peraxis$p < 0.05), sum(peraxis$q < 0.05))
