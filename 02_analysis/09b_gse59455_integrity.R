#!/usr/bin/env Rscript
# =============================================================================
# P2-T11b：第二外部验证队列 GSE59455（指令二十七）
# 队列：GSE59455（141 样本，FFPE，原发+转移混合，Illumina HumanRef-8 WG-DASL v3.0 / GPL8432）
# 终点：OS（周）
#
# ★★ 两个必须处理的数据完整性陷阱：
#   ① characteristics 各「行」并不按样本对齐 —— 提交者把不同字段写在不同行的不同列，
#      同一行在不同样本位置出现不同标签。→ 必须**按单元格标签驱动解析**，不能按行解析。
#   ② 同一样本标题带 _A/_B/_C 等后缀 —— 同一患者多块蜡块重复取样。
#      → 必须**折叠到患者级**，否则构成伪重复（pseudo-replication），CI 被系统性低估。
# =============================================================================
suppressMessages({library(GSVA); library(survival); library(data.table)})
setwd('C:/Rtmp/gse59455')
OUT <- 'C:/Rtmp/gse59455/out'; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
w <- function(...) cat(sprintf(...), '\n', sep = '')

## ================= 0. 基因集 =================
almp <- list()
for (ln in readLines('C:/Rtmp/p2t4/ALMP_gene_set_v1.gmt')) {
  p <- strsplit(ln, '\t')[[1]]; almp[[p[1]]] <- unique(p[-c(1, 2)])
}
w('[0] ALMP 集 %d：%s', length(almp), paste(sprintf('%s(%d)', names(almp), sapply(almp, length)), collapse = ' '))

## ================= 1. 系列矩阵 =================
sm <- readLines(gzfile('GSE59455_series_matrix.txt.gz'))
gi <- grep('^!series_matrix_table_begin', sm); ge <- grep('^!series_matrix_table_end', sm)
gv1 <- function(pat) { i <- grep(pat, sm); if (!length(i)) return(NULL)
  gsub('"', '', strsplit(gsub('^[^\t]*\t', '', sm[i[1]]), '\t')[[1]]) }
gvA <- function(pat) { i <- grep(pat, sm); if (!length(i)) return(NULL)
  lapply(i, function(k) gsub('"', '', strsplit(gsub('^[^\t]*\t', '', sm[k]), '\t')[[1]])) }
acc <- gv1('^!Sample_geo_accession'); ttl <- gv1('^!Sample_title'); ch <- gvA('^!Sample_characteristics_ch1')
w('[1] 样本 %d；characteristics 行 %d（★ 每行标签不一致，须按单元格解析）', length(acc), length(ch))

## ★ 标签驱动解析：逐单元格取 "label: value"
recs <- vector('list', length(acc))
for (j in seq_along(acc)) recs[[j]] <- list()
for (k in seq_along(ch)) {
  v <- ch[[k]]
  if (length(v) < length(acc)) next
  for (j in seq_along(acc)) {
    cell <- trimws(v[j])
    if (is.na(cell) || cell == '') next
    if (grepl(':', cell)) {
      lab <- trimws(sub(':.*$', '', cell))
      val <- trimws(sub('^[^:]*:', '', cell))
      if (nzchar(val)) recs[[j]][[lab]] <- val
    }
  }
}
labs <- sort(unique(unlist(lapply(recs, names))))
w('[1] 解析出 %d 个字段标签：', length(labs)); w('     %s', paste(labs, collapse = ' | '))

mk <- function(lab) vapply(recs, function(r) if (!is.null(r[[lab]])) r[[lab]] else NA_character_, '')
clin <- data.table(sample = acc, title = ttl)
for (lab in labs) {
  key <- gsub('[^A-Za-z0-9]', '_', tolower(lab))
  clin[[key]] <- mk(lab)
}
## 关键字段重命名（按实际标签模糊匹配）
ren <- function(pat, new) {
  hit <- names(clin)[grepl(pat, names(clin))][1]
  if (!is.na(hit) && !is.null(hit)) { setnames(clin, hit, new); w('[1] %-14s <- "%s"', new, hit) } else w('[1] %-14s ✗ 未匹配 /%s/', new, pat)
}
ren('overall_survival', 'os_w')
ren('age_at_diagosis|age_at_diagnosis', 'age_dx')
ren('age_at_death', 'age_death')
ren('gender|sex', 'sex')
ren('primary_metastatic', 'site_type')
ren('breslow', 'breslow')
ren('solar_elastosis', 'elastosis')
ren('sun_exposure', 'sun_score')
ren('braf_nras_wt', 'mutation')

## 患者级折叠键：标题形如 4531335087_A
clin[, patient := sub('_[A-Za-z0-9]+$', '', title)]
w('[1] 样本 %d → 患者 %d（每患者 %.2f 块）', nrow(clin), uniqueN(clin$patient), nrow(clin)/uniqueN(clin$patient))
print(table(table(clin$patient)))

## ================= 2. 表达矩阵 =================
tab <- sm[(gi + 1):(ge - 1)]
hdr <- gsub('"', '', strsplit(tab[1], '\t')[[1]])
rows <- strsplit(tab[-1], '\t')
m <- matrix(as.numeric(unlist(lapply(rows, function(r) r[-1]))), ncol = length(hdr) - 1, byrow = TRUE)
rownames(m) <- vapply(rows, function(r) gsub('"', '', r[1]), '')
colnames(m) <- hdr[-1]
qq <- quantile(m, c(0, .25, .5, .75, .99, 1), na.rm = TRUE)
w('[2] 原始 %d 探针 x %d 样本 | 分位 %.2f/%.2f/%.2f/%.2f/%.1f/%.1f',
  nrow(m), ncol(m), qq[1], qq[2], qq[3], qq[4], qq[5], qq[6])
if (max(m, na.rm = TRUE) > 100 && median(m, na.rm = TRUE) > 20) {
  m <- log2(m); w('[2] ★ 判定线性强度 → log2(x)（%.2f/%.2f/%.2f）', min(m,na.rm=T), median(m,na.rm=T), max(m,na.rm=T))
} else w('[2] 判定已是 log 尺度 → 不变')

pm <- fread('GPL8432_probe2symbol.tsv', header = TRUE)
pm <- pm[!is.na(symbol) & symbol != '' & symbol != '---']
w('[2] 平台映射 %d 条 / 唯一符号 %d', nrow(pm), uniqueN(pm$symbol))
m2 <- m[rownames(m) %in% pm$probe, , drop = FALSE]
mm <- pm[match(rownames(m2), probe)]
mu <- rowMeans(m2, na.rm = TRUE); o <- order(mm$symbol, -mu)
m2 <- m2[o, , drop = FALSE]; sy <- mm$symbol[o]
keep <- !duplicated(sy)
expr <- m2[keep, , drop = FALSE]; rownames(expr) <- sy[keep]
w('[2] 基因矩阵 %d 基因 x %d 样本', nrow(expr), ncol(expr))
q10 <- quantile(expr, 0.10, na.rm = TRUE)
expr <- expr[rowMeans(expr > q10, na.rm = TRUE) >= 0.50, , drop = FALSE]
w('[2] 过滤（检出率≥50%%）后 %d 基因；ALMP 覆盖：%s', nrow(expr),
  paste(sprintf('%s %d/%d', names(almp), sapply(almp, function(s) sum(s %in% rownames(expr))), sapply(almp, length)), collapse = ' | '))

## ================= 3. ★ 患者级折叠 =================
pat <- unique(clin$patient)
exprp <- sapply(pat, function(p) {
  cs <- clin$sample[clin$patient == p]; cs <- intersect(cs, colnames(expr))
  if (length(cs) == 1) expr[, cs] else rowMeans(expr[, cs, drop = FALSE], na.rm = TRUE)
})
rownames(exprp) <- rownames(expr)
w('[3] 患者级表达矩阵 %d 基因 x %d 患者', nrow(exprp), ncol(exprp))

clp <- clin[, .(n_cores = .N,
  os_w = suppressWarnings(as.numeric(os_w[!is.na(os_w)][1])),
  sex = sex[!is.na(sex)][1],
  age_dx = suppressWarnings(as.numeric(age_dx[!is.na(age_dx)][1])),
  age_death = suppressWarnings(as.numeric(age_death[!is.na(age_death)][1])),
  site_type = site_type[!is.na(site_type)][1], breslow = suppressWarnings(as.numeric(breslow[!is.na(breslow)][1])),
  elastosis = suppressWarnings(as.numeric(elastosis[!is.na(elastosis)][1])),
  mutation = mutation[!is.na(mutation)][1]), by = patient]

## ★★ 事件指示符推导（本队列无 vital status 字段）
##   规则：age at death 存在 → 认为已死亡（event=1）；不存在 → 删失（event=0）
clp[, event := fifelse(is.finite(age_death), 1L, 0L)]
w('[3] ★ 无 vital status 字段 → 以「age at death 是否存在」推导事件：事件 %d / 删失 %d',
  sum(clp$event == 1), sum(clp$event == 0))
## ★ 内部一致性核验：OS 周数 vs (age_death − age_dx) × 52
clp[, os_calc := (age_death - age_dx) * 52]
ok <- is.finite(clp$os_w) & is.finite(clp$os_calc)
if (sum(ok) >= 10) {
  ct <- suppressWarnings(cor.test(clp$os_calc[ok], clp$os_w[ok], method = 'spearman'))
  w('[3] ★ OS 一致性核验（n=%d）：Spearman rho=%.3f p=%.3g；中位差 %.0f 周',
    sum(ok), unname(ct$estimate), ct$p.value, median(clp$os_calc[ok] - clp$os_w[ok]))
}
setorder(clp, patient)
w('[3] 患者级临床 %d 行；OS 可用 %d；中位 OS %.0f 周', nrow(clp), sum(is.finite(clp$os_w)), median(clp$os_w, na.rm = TRUE))
print(clp[, .(n = .N, med_OS_wk = round(median(os_w, na.rm=TRUE)),
              n_primary = sum(site_type=='Primary', na.rm=TRUE),
              n_metastatic = sum(site_type=='Metastatic', na.rm=TRUE),
              n_event = sum(event==1))])
fwrite(clp, file.path(OUT, 'P2T11b_clinical_patient_level.tsv'), sep = '\t')

## ================= 4. 评分 + Cox =================
gs <- lapply(almp, function(s) intersect(s, rownames(exprp))); gs <- gs[sapply(gs, length) >= 3]
ss  <- tryCatch(gsva(ssgseaParam(exprp, gs, normalize = TRUE), verbose = FALSE), error=function(e){w('ssGSEA ERR %s', conditionMessage(e));NULL})
gv2 <- tryCatch(gsva(gsvaParam(exprp, gs, kcdf = 'Gaussian'), verbose = FALSE), error=function(e){w('GSVA ERR %s', conditionMessage(e));NULL})
w('[4] ssGSEA %s | GSVA %s', paste(dim(ss), collapse='x'), paste(dim(gv2), collapse='x'))

cox_one <- function(score, tt, ee, tag) {
  d <- data.frame(s = as.numeric(score), t = tt, e = ee)
  d <- d[is.finite(d$s) & is.finite(d$t) & d$t > 0, ]
  if (nrow(d) < 20 || isTRUE(sd(d$s) == 0)) return(NULL)   # ★ GSE59455 患者级仅 23 例，阈值放宽至 20（须在报告标注为功效不足）
  d$s <- as.numeric(scale(d$s))
  s1 <- summary(coxph(Surv(t, e) ~ s, data = d))
  a <- data.frame(endpoint = tag, model = 'continuous(per SD)', n = nrow(d),
                  HR = round(s1$conf.int[1,'exp(coef)'],4), ci_lo = round(s1$conf.int[1,'lower .95'],4),
                  ci_hi = round(s1$conf.int[1,'upper .95'],4), p = signif(s1$coefficients[1,'Pr(>|z|)'],4))
  dd <- d; dd$grp <- factor(ifelse(dd$s >= median(dd$s), 'high', 'low'), levels = c('low','high'))
  s2 <- summary(coxph(Surv(t, e) ~ grp, data = dd))
  b <- data.frame(endpoint = tag, model = 'median split(high vs low)', n = nrow(dd),
                  HR = round(s2$conf.int[1,'exp(coef)'],4), ci_lo = round(s2$conf.int[1,'lower .95'],4),
                  ci_hi = round(s2$conf.int[1,'upper .95'],4), p = signif(s2$coefficients[1,'Pr(>|z|)'],4))
  rbind(a, b)
}
## ★ 注意：GSE59455 只提供 OS（周），**无事件指示符** —— 该队列为「生存时间」型记录。
##   若全部视为事件（1），Cox 退化为加速失效/参数化假设，必须显式声明。
w('\n★★ 事件指示符缺失检查：字段列表 = %s', paste(names(clp), collapse=', '))
has_event <- any(grepl('status|event|dead|vital', names(clp), ignore.case = TRUE))
w('   event 字段存在 = %s → 若无，则以「死因为黑色素瘤且已死亡者」推断不可行，改用**时间型 Cox（全事件）**并在报告标注为局限', has_event)

res <- data.frame()
for (rd in c('ssGSEA','GSVA')) {
  sc <- if (rd=='ssGSEA') ss else gv2; if (is.null(sc)) next
  for (ft in rownames(sc)) {
    if (!ft %in% c(names(almp),'ALMP_ALL_308')) next
    r <- cox_one(sc[ft, ], clp$os_w, clp$event, 'OS')
    if (!is.null(r)) res <- rbind(res, cbind(readout = rd, feature = ft, r))
  }
}
if (nrow(res)) res$q <- signif(p.adjust(res$p,'BH'),4)
fwrite(res, file.path(OUT,'P2T11b_axis_OS.tsv'), sep='\t')
w('\n===== ALMP 轴 × OS（%d 条，患者级）=====', nrow(res))
if (nrow(res)) print(res[order(res$p), ], row.names = FALSE) else w('（无）')

w('\n===== ★ 主验证目标 ssGSEA SEG1_ADIPOKINE_CORE =====')
tv <- res[res$readout=='ssGSEA' & res$feature=='SEG1_ADIPOKINE_CORE', ]
if (nrow(tv)) print(tv, row.names = FALSE) else w('（未产出）')

if ('CDS2' %in% rownames(exprp)) {
  cds <- cox_one(exprp['CDS2', ], clp$os_w, clp$event, 'OS')
  if (!is.null(cds)) { fwrite(cds, file.path(OUT,'P2T11b_CDS2_OS.tsv'), sep='\t')
    w('\n===== CDS2 逐基因 OS ====='); print(cds, row.names = FALSE) } else w('\nCDS2 未产出（样本量不足）')
}

## ================= 5. 模型无关读出（OS 时间组间 Wilcoxon）=================
w('\n===== 模型无关读出：OS 周数 高 vs 低（中位切分，Wilcoxon）=====')
mw <- data.frame()
for (rd in c('ssGSEA','GSVA')) {
  sc <- if (rd=='ssGSEA') ss else gv2; if (is.null(sc)) next
  for (ft in rownames(sc)) {
    if (!ft %in% c(names(almp),'ALMP_ALL_308')) next
    x <- as.numeric(sc[ft, ]); ok <- is.finite(x) & is.finite(clp$os_w)
    if (sum(ok) < 20) next
    g <- x[ok] >= median(x[ok]); tt <- clp$os_w[ok]
    wt <- suppressWarnings(wilcox.test(tt[g], tt[!g]))
    mw <- rbind(mw, data.frame(readout=rd, feature=ft, n=sum(ok),
      med_high=round(median(tt[g]),1), med_low=round(median(tt[!g]),1),
      p=signif(wt$p.value,4)))
  }
}
if (nrow(mw)) { mw$q <- signif(p.adjust(mw$p,'BH'),4); fwrite(mw, file.path(OUT,'P2T11b_OS_wilcox.tsv'), sep='\t') }
print(mw[mw$feature %in% c('SEG1_ADIPOKINE_CORE','ALMP_ALL_308'), ], row.names = FALSE)

saveRDS(list(clin=clp, ss=ss, gv=gv2, exprp=exprp, res=res, mw=mw), file.path(OUT,'P2T11b_all.rds'))
writeLines(capture.output(sessionInfo()), file.path(OUT,'P2T11b_sessionInfo.txt'))
w('\nDONE')
