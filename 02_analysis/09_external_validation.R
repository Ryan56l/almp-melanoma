#!/usr/bin/env Rscript
# =============================================================================
# P2-T11：外部独立队列验证（指令二十二）
# 队列：GSE65904（Cirenajwis 2015，转移性黑色素瘤，Illumina HumanHT-12 v4 / GPL10558）
# 目标：验证 TCGA-SKCM 探索性信号 —— ssGSEA SEG1_ADIPOKINE_CORE 转移灶 HR 0.838（p=0.027）
#       方向能否在外部队列复现；CDS2 逐基因生存关联登记。
#
# ★ 纪律：① GSVA 过滤冻结规则（样本级检出率≥50%）；② 核心集仅在外部队列方向一致且名义显著
#        时方可冻结命名；③ 全轴全终点 BH-FDR，逐条登记，禁择优；④ 尺度直接看数值范围，不猜。
# =============================================================================
suppressMessages({library(GSVA); library(survival); library(data.table)})
setwd('C:/Rtmp/p2t11')
OUT <- 'C:/Rtmp/p2t11/out'; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(...) cat(sprintf(...), '\n', sep = '')

## ================= 0. 基因集 =================
almp <- list()
for (ln in readLines('C:/Rtmp/p2t4/ALMP_gene_set_v1.gmt')) {
  p <- strsplit(ln, '\t')[[1]]; almp[[p[1]]] <- unique(p[-c(1, 2)])
}
w('[0] ALMP 集 %d：%s', length(almp), paste(sprintf('%s(%d)', names(almp), sapply(almp, length)), collapse = ' '))

## ================= 1. 系列矩阵 =================
sm <- readLines(gzfile('GSE65904_series_matrix.txt.gz'))
gi <- grep('^!series_matrix_table_begin', sm); ge <- grep('^!series_matrix_table_end', sm)
w('[1] 系列矩阵 %d 行；表体 %d-%d', length(sm), gi, ge)
gv1 <- function(pat) { i <- grep(pat, sm); if (!length(i)) return(NULL)
  gsub('"', '', strsplit(gsub('^[^\t]*\t', '', sm[i[1]]), '\t')[[1]]) }
gvA <- function(pat) { i <- grep(pat, sm); if (!length(i)) return(NULL)
  lapply(i, function(k) gsub('"', '', strsplit(gsub('^[^\t]*\t', '', sm[k]), '\t')[[1]])) }
acc <- gv1('^!Sample_geo_accession'); ch <- gvA('^!Sample_characteristics_ch1')
w('[1] 样本 %d；characteristics %d 行', length(acc), length(ch))
for (k in seq_along(ch)) w('     char[%d] 例: %s', k, paste(head(unique(ch[[k]]), 4), collapse = ' | '))

tab <- sm[(gi + 1):(ge - 1)]
hdr <- gsub('"', '', strsplit(tab[1], '\t')[[1]])
rows <- strsplit(tab[-1], '\t')
m <- matrix(as.numeric(unlist(lapply(rows, function(r) r[-1]))), ncol = length(hdr) - 1, byrow = TRUE)
rownames(m) <- vapply(rows, function(r) gsub('"', '', r[1]), '')
colnames(m) <- hdr[-1]
qq <- quantile(m, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
w('[1] 原始矩阵 %d 探针 x %d 样本 | 分位 %.1f/%.1f/%.1f/%.1f/%.1f/%.1f', nrow(m), ncol(m), qq[1], qq[2], qq[3], qq[4], qq[5], qq[6])
## ★ 尺度自检：Illumina BeadArray 提交值多为线性荧光强度 → 需 log2；已 log 者值域 ~2-15
if (max(m, na.rm = TRUE) > 100 && median(m, na.rm = TRUE) > 20) {
  m <- log2(m); w('[1] ★ 判定为线性强度 → 已 log2(x)（变换后 %.2f/%.2f/%.2f）',
                  min(m, na.rm=TRUE), median(m, na.rm=TRUE), max(m, na.rm=TRUE))
} else w('[1] 判定已是 log 尺度 → 保持不变')

## ================= 2. 探针 → 基因 =================
an <- readLines(gzfile('GPL10558.annot.gz'))
b <- grep('^!platform_table_begin', an); e <- grep('^!platform_table_end', an)
ah <- strsplit(an[b + 1], '\t')[[1]]
w('[2] 平台注释列：%s', paste(ah, collapse = ' | '))
sp <- strsplit(an[(b + 2):(e - 1)], '\t')
gsc <- grep('Gene symbol', ah)[1]; w('[2] 探针列=%s；基因列=%s', ah[1], ah[gsc])
pr <- vapply(sp, function(x) x[1], '')
gs <- vapply(sp, function(x) if (length(x) >= gsc) x[gsc] else '', '')
map <- data.table(probe = pr, symbol = gs)[symbol != '' & !is.na(symbol)]
map[, symbol := tstrsplit(symbol, '///', fixed = TRUE)[[1]]][, symbol := trimws(symbol)]
w('[2] 注释 %d 条；唯一探针 %d；唯一符号 %d', nrow(map), uniqueN(map$probe), uniqueN(map$symbol))

m2 <- m[rownames(m) %in% map$probe, , drop = FALSE]
mm <- map[match(rownames(m2), map$probe)]
w('[2] 命中探针 %d / %d', nrow(m2), nrow(m))
mu <- rowMeans(m2, na.rm = TRUE); o <- order(mm$symbol, -mu)
m2 <- m2[o, , drop = FALSE]; sy <- mm$symbol[o]
keep <- !duplicated(sy)
expr <- m2[keep, , drop = FALSE]; rownames(expr) <- sy[keep]
w('[2] 基因矩阵 %d 基因 x %d 样本', nrow(expr), ncol(expr))

q10 <- quantile(expr, 0.10, na.rm = TRUE)
expr <- expr[rowMeans(expr > q10, na.rm = TRUE) >= 0.50, , drop = FALSE]
w('[2] 过滤后（检出率≥50%%）%d 基因；ALMP 覆盖：%s', nrow(expr),
  paste(sprintf('%s %d/%d', names(almp),
                sapply(almp, function(s) sum(s %in% rownames(expr))), sapply(almp, length)), collapse = ' | '))

## ================= 3. 临床终点 =================
## GSE65904 实际 characteristics（经 gzip 部分解压核验，214 样本）：
##   gender / age / tumor stage / tissue / distant metastasis free survival in days + (1=death, 0=alive)
##   / disease specific survival in days + (1=death, 0=alive)
clin <- data.frame(sample = acc, stringsAsFactors = FALSE)
pick <- function(pat) { for (k in seq_along(ch)) if (any(grepl(pat, ch[[k]], ignore.case = TRUE))) return(ch[[k]]); NULL }
cand <- list(
  DMFS_time = c('distant metastasis free survival in days'),
  DMFS_event= c('distant metastasis free survival \\(1='),
  DSS_time  = c('disease specific survival in days'),
  DSS_event = c('disease specific survival \\(1=')
)
for (nm in names(cand)) {
  hit <- NULL
  for (pat in cand[[nm]]) { v <- pick(pat)
    if (!is.null(v) && length(unique(v)) > 1) { clin[[nm]] <- v; hit <- pat; w('[3] %-10s <- "%s"', nm, pat); break } }
  if (is.null(hit)) w('[3] %-10s ✗ 未匹配', nm)
}
## 附加协变量（供敏感性）：性别 / 年龄 / 采样部位
for (nm in c('gender','age','tissue','tumor stage')) {
  v <- pick(paste0('^', nm, ':'))
  if (!is.null(v)) { clin[[gsub(' ', '_', nm)]] <- sub('^[^:]*:\\s*', '', v); w('[3] %-10s <- 协变量', nm) }
}
strip_lab <- function(x) sub('^[^:]*:[[:space:]]*', '', trimws(as.character(x)))
ev2num <- function(x) { x <- tolower(trimws(strip_lab(x)))
  ifelse(x %in% c('1','yes','dead','death','event','deceased','recurrence','metastasis','relapse','true','1:event'), 1L,
  ifelse(x %in% c('0','no','alive','censored','no event','disease free','false','0:censored'), 0L, NA_integer_)) }
t2num <- function(x) suppressWarnings(as.numeric(sub('^[^:]*:\\s*', '', trimws(as.character(x)))))
for (nm in c('DMFS','DSS')) {
  tc <- paste0(nm, '_time'); ec <- paste0(nm, '_event')
  if (tc %in% names(clin) && ec %in% names(clin)) {
    clin[[paste0(nm, '_t')]] <- t2num(clin[[tc]]); clin[[paste0(nm, '_e')]] <- ev2num(clin[[ec]])
    v <- is.finite(clin[[paste0(nm,'_t')]]) & !is.na(clin[[paste0(nm,'_e')]]) & clin[[paste0(nm,'_t')]] > 0
    w('[3] %s：可用 n=%d（事件 %d）｜时间中位 %.0f 天', nm, sum(v), sum(clin[[paste0(nm,'_e')]] == 1, na.rm = TRUE),
      median(clin[[paste0(nm,'_t')]][v]))
  }
}
fwrite(clin, file.path(OUT, 'P2T11_GSE65904_clinical_parsed.tsv'), sep = '\t')

## ================= 4. 评分 =================
gs <- lapply(almp, function(s) intersect(s, rownames(expr)))
gs <- gs[sapply(gs, length) >= 3]
ss  <- tryCatch(gsva(ssgseaParam(expr, gs, normalize = TRUE), verbose = FALSE),
                error = function(e2) { w('ssGSEA ERR %s', conditionMessage(e2)); NULL })
gv2 <- tryCatch(gsva(gsvaParam(expr, gs, kcdf = 'Gaussian'), verbose = FALSE),
                error = function(e2) { w('GSVA ERR %s', conditionMessage(e2)); NULL })
w('[4] ssGSEA %s | GSVA %s', paste(dim(ss), collapse = 'x'), paste(dim(gv2), collapse = 'x'))
fwrite(data.table(set = rownames(ss), as.data.frame(ss)),  file.path(OUT, 'P2T11_ssGSEA_scores.tsv'), sep = '\t')
fwrite(data.table(set = rownames(gv2), as.data.frame(gv2)), file.path(OUT, 'P2T11_GSVA_scores.tsv'), sep = '\t')

## ================= 5. Cox =================
cox_one <- function(score, clin, tcol, ecol, tag) {
  d <- data.frame(s = as.numeric(score), t = clin[[tcol]], e = clin[[ecol]])
  d <- d[is.finite(d$s) & is.finite(d$t) & !is.na(d$e), ]
  d <- d[d$t > 0, ]
  if (nrow(d) < 30 || sum(d$e == 1) < 10 || isTRUE(sd(d$s) == 0)) return(NULL)
  d$s <- as.numeric(scale(d$s))
  s1 <- summary(coxph(Surv(t, e) ~ s, data = d))
  a <- data.frame(endpoint = tag, model = 'continuous(per SD)', n = nrow(d), events = sum(d$e == 1),
                  HR = round(s1$conf.int[1, 'exp(coef)'], 4), ci_lo = round(s1$conf.int[1, 'lower .95'], 4),
                  ci_hi = round(s1$conf.int[1, 'upper .95'], 4), p = signif(s1$coefficients[1, 'Pr(>|z|)'], 4))
  dd <- d; dd$grp <- factor(ifelse(dd$s >= median(dd$s), 'high', 'low'), levels = c('low', 'high'))
  s2 <- summary(coxph(Surv(t, e) ~ grp, data = dd))
  b <- data.frame(endpoint = tag, model = 'median split(high vs low)', n = nrow(dd), events = sum(dd$e == 1),
                  HR = round(s2$conf.int[1, 'exp(coef)'], 4), ci_lo = round(s2$conf.int[1, 'lower .95'], 4),
                  ci_hi = round(s2$conf.int[1, 'upper .95'], 4), p = signif(s2$coefficients[1, 'Pr(>|z|)'], 4))
  rbind(a, b)
}
eps <- intersect(c('DMFS_t', 'DSS_t'), names(clin))
res <- data.frame()
for (ep in eps) { tag <- sub('_t$', '', ep)
  for (rd in c('ssGSEA', 'GSVA')) {
    sc <- if (rd == 'ssGSEA') ss else gv2; if (is.null(sc)) next
    for (ft in rownames(sc)) {
      if (!(ft %in% c(names(almp), 'ALMP_ALL_308'))) next
      r <- cox_one(sc[ft, ], clin, ep, paste0(tag, '_e'), tag)
      if (!is.null(r)) res <- rbind(res, cbind(readout = rd, feature = ft, r))
    }
  }
}
if (nrow(res)) res$q <- signif(p.adjust(res$p, 'BH'), 4)
fwrite(res, file.path(OUT, 'P2T11_axis_survival_assoc.tsv'), sep = '\t')
w('\n===== ALMP 轴 × 生存（全部 %d 条）；q<0.05 者 %d =====', nrow(res), sum(res$q < 0.05))
if (nrow(res)) print(res[order(res$p), ][1:min(20, nrow(res)), ], row.names = FALSE) else w('（无结果）')
w('\n===== ★ 主验证目标：ssGSEA SEG1_ADIPOKINE_CORE =====')
tv <- res[res$readout == 'ssGSEA' & res$feature == 'SEG1_ADIPOKINE_CORE', ]
if (nrow(tv)) print(tv, row.names = FALSE) else w('（未产出）')

## ================= 6. CDS2 逐基因生存 =================
w('\n===== CDS2 逐基因生存关联 =====')
cds <- data.frame()
if ('CDS2' %in% rownames(expr)) {
  for (ep in eps) { tag <- sub('_t$', '', ep)
    r <- cox_one(expr['CDS2', ], clin, ep, paste0(tag, '_e'), tag)
    if (!is.null(r)) cds <- rbind(cds, data.frame(cohort = 'GSE65904', gene = 'CDS2', r))
  }
  if (nrow(cds)) cds$q <- signif(p.adjust(cds$p, 'BH'), 4)
  fwrite(cds, file.path(OUT, 'P2T11_CDS2_survival.tsv'), sep = '\t')
  print(cds, row.names = FALSE)
} else w('CDS2 不在表达矩阵中')

## ================= 7. 主验证目标的协变量校正敏感性 =================
w('\n===== 协变量校正敏感性（主验证目标）=====')
adj <- data.frame()
for (ep in eps) { tag <- sub('_t$', '', ep)
  for (rd in c('ssGSEA', 'GSVA')) {
    sc <- if (rd == 'ssGSEA') ss else gv2
    for (ft in c('SEG1_ADIPOKINE_CORE', 'ALMP_TUMOR_INTRINSIC_AXIS')) {
      if (!ft %in% rownames(sc)) next
      d <- data.frame(s = as.numeric(scale(sc[ft, ])), t = clin[[ep]], e = clin[[paste0(tag, '_e')]],
                      age = suppressWarnings(as.numeric(clin$age)),
                      sex = factor(clin$gender), tis = factor(clin$tissue))
      d <- d[is.finite(d$s) & is.finite(d$t) & !is.na(d$e) & d$t > 0, ]
      if (nrow(d) < 30 || sum(d$e == 1) < 10) next
      f1 <- summary(coxph(Surv(t, e) ~ s, data = d))
      f2 <- summary(coxph(Surv(t, e) ~ s + age + sex, data = d))
      f3 <- summary(coxph(Surv(t, e) ~ s + age + sex + tis, data = d))
      g <- function(f) c(HR = round(f$conf.int['s','exp(coef)'],4), lo = round(f$conf.int['s','lower .95'],4),
                         hi = round(f$conf.int['s','upper .95'],4), p = signif(f$coefficients['s','Pr(>|z|)'],4))
      adj <- rbind(adj, data.frame(cohort='GSE65904', endpoint=tag, readout=rd, feature=ft, n=nrow(d),
        events=sum(d$e==1), model=c('crude','+age+sex','+age+sex+tissue'),
        rbind(g(f1), g(f2), g(f3))))
    }
  }
}
fwrite(adj, file.path(OUT, 'P2T11_covariate_adjusted.tsv'), sep = '\t')
print(adj, row.names = FALSE)

saveRDS(list(clin = clin, ss = ss, gv = gv2, expr = expr, res = res, cds2 = cds, adj = adj), file.path(OUT, 'P2T11_all.rds'))
writeLines(capture.output(sessionInfo()), file.path(OUT, 'P2T11_sessionInfo.txt'))
w('\nDONE')
