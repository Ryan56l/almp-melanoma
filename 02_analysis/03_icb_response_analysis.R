#!/usr/bin/env Rscript
# P2-T6: ICB 应答分层分析（GSE78220 / GSE91061）
# 四套免疫读出：ssGSEA / GSVA（ALMP 轴）+ MCP-counter + quanTIseq
# ★ 全部使用「全转录组」矩阵：quanTIseq 内部按列 TPM 归一化、MCP-counter 依赖全局归一化，
#   基因子集会引入样本特异偏倚
suppressMessages({library(GSVA); library(readxl); library(data.table); library(preprocessCore); library(quantiseqr)})
setwd('C:/Rtmp/deconv/icb')
w <- function(...) cat(sprintf(...), '\n', sep = '')
OUT <- 'C:/Rtmp/deconv/icb/out'; dir.create(OUT, showWarnings = FALSE)

## ================= 0. 基因集与标记 =================
almp <- list()
for (ln in readLines('C:/Rtmp/p2t4/ALMP_gene_set_v1.gmt')) {
  p <- strsplit(ln, '\t')[[1]]; almp[[p[1]]] <- unique(p[-c(1, 2)])
}
mcp <- fread('C:/Rtmp/deconv/MCPcounter_genes.txt', header = TRUE)
setnames(mcp, 1:3, c('symbol', 'population', 'entrez'))
mcp <- mcp[symbol != '' & population != '']
mcp_list <- split(mcp$symbol, mcp$population)
w('ALMP sets: %s', paste(sprintf('%s(%d)', names(almp), sapply(almp, length)), collapse = ' '))
w('MCP populations: %s', paste(sprintf('%s(%d)', names(mcp_list), sapply(mcp_list, length)), collapse = ' '))

# ★ 官方 MCP-counter 核心（ebecht/MCPcounter Source/R/MCPcounter.R）
appendSignatures <- function(xp, markers) as.data.frame(do.call(cbind, lapply(markers, function(x)
  apply(xp[intersect(rownames(xp), x), , drop = FALSE], 2, mean, na.rm = TRUE))))
mcp_score <- function(mat_linear) {
  lg <- log2(mat_linear + 1)
  lg <- preprocessCore::normalize.quantiles(lg); dimnames(lg) <- dimnames(mat_linear)
  as.data.frame(t(appendSignatures(lg, mcp_list)))   # populations x samples
}

score_cohort <- function(mat, tag) {
  gs <- lapply(almp, function(s) intersect(s, rownames(mat)))
  w('[%s] %d genes x %d samples | set coverage: %s', tag, nrow(mat), ncol(mat),
    paste(sprintf('%s %d/%d', names(almp), sapply(gs, length), sapply(almp, length)), collapse = ' | '))
  gs <- gs[sapply(gs, length) >= 3]
  ss <- tryCatch(gsva(ssgseaParam(mat, gs, normalize = TRUE), verbose = FALSE),
                 error = function(e) { message('ssGSEA ERR: ', conditionMessage(e)); NULL })
  gv <- tryCatch(gsva(gsvaParam(mat, gs, kcdf = 'Gaussian'), verbose = FALSE),
                 error = function(e) { message('GSVA ERR: ', conditionMessage(e)); NULL })
  mc <- mcp_score(mat)
  qt <- tryCatch({
    q <- run_quantiseq(mat, 'TIL10', FALSE, TRUE, TRUE, 'lsei', 'default')
    q <- as.data.frame(q)
    if (is.character(q[[1]])) { rownames(q) <- q[[1]]; q <- q[, -1, drop = FALSE] }
    ct <- c('B.cells','Macrophages.M1','Macrophages.M2','Monocytes','Neutrophils',
            'NK.cells','T.cells.CD4','T.cells.CD8','Tregs','Dendritic.cells')
    if (!all(ct %in% rownames(q)) && all(ct %in% colnames(q))) q <- as.data.frame(t(as.matrix(q)))
    q
  }, error = function(e) { w('quanTIseq ERR [%s]: %s', tag, conditionMessage(e)); NULL })
  if (!is.null(qt)) w('[%s] quanTIseq: %d cell types', tag, nrow(qt))
  # ★ quanTIseq 输出为成分数据（每样本和为 1），直接算相关会被"和恒定"约束衰减；
  #   补 CLR（中心化对数比）变换作为比例无关读出
  qtc <- NULL
  if (!is.null(qt)) {
    qm <- as.matrix(qt); qm[qm <= 0] <- min(qm[qm > 0]) / 2
    qtc <- as.data.frame(apply(qm, 2, function(x) log(x) - mean(log(x))))
  }
  list(ss = ss, gv = gv, mc = mc, qt = qt, qtclr = qtc)
}

## ================= 1-2. 评分（带缓存）=================
if (file.exists(file.path(OUT, 'icb_scores.rds'))) {
  L <- readRDS(file.path(OUT, 'icb_scores.rds'))
  clin78220 <- L$clin78220; clin91061 <- L$clin91061; s78220 <- L$s78220; s91061 <- L$s91061
  w('[cache] 复用 icb_scores.rds')
} else {
## ================= 1. GSE78220（symbol 全矩阵，FPKM）=================
d <- as.data.frame(read_excel('GSE78220_PatientFPKM.xlsx', sheet = 1), check.names = FALSE)
rownames(d) <- d[[1]]; d <- d[, -1, drop = FALSE]
mat78220 <- as.matrix(d); mat78220 <- mat78220[!duplicated(rownames(mat78220)), , drop = FALSE]

sm <- readLines(gzfile('GSE78220_series_matrix.txt.gz'))
grab <- function(key, pat = NULL) {
  l <- if (is.null(pat)) sm[grepl(paste0('^!Sample_', key), sm)][1] else sm[grepl(pat, sm)][1]
  gsub('"', '', strsplit(gsub('^[^\t]*\t', '', l), '\t')[[1]])
}
ttl <- grab('title'); rline <- grab(NULL, '^!Sample_characteristics_ch1.*anti-pd-1 response')
resp <- sub('.*anti-pd-1 response: *', '', rline)
clin78220 <- data.frame(col = colnames(mat78220), ttl = ttl, response = resp, stringsAsFactors = FALSE)
clin78220$patient <- sub('\\..*$', '', clin78220$col)
clin78220$visit <- sub('^[^.]*\\.', '', clin78220$col)
clin78220$bin <- ifelse(clin78220$response %in% c('Complete Response', 'Partial Response'), 'R',
                 ifelse(clin78220$response == 'Progressive Disease', 'NR', NA))
write.table(clin78220, file.path(OUT, 'icb_clin_GSE78220.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
w('\n[GSE78220] response:'); print(table(clin78220$response))
w('[GSE78220] pre-treatment analytic:'); print(table(clin78220$bin[clin78220$visit == 'baseline']))
s78220 <- score_cohort(mat78220, 'GSE78220')

## ================= 2. GSE91061（entrez→symbol 全矩阵，FPKM）=================
fp <- fread('GSE91061_BMS038109Sample.hg19KnownGene.fpkm.csv.gz', header = TRUE)
setnames(fp, 1, 'entrez')
e2f <- fread('C:/Rtmp/deconv/entrez2symbol_full.tsv')
w('\n[GSE91061] raw: %d entrez rows x %d samples; entrez2symbol: %d', nrow(fp), ncol(fp) - 1, nrow(e2f))
fp[, entrez := as.character(entrez)]
rn <- e2f$symbol[match(fp$entrez, as.character(e2f$entrez))]
mm <- as.matrix(fp[, -1, drop = FALSE])
ok <- !is.na(rn)
mm <- mm[ok, , drop = FALSE]; rr <- rn[ok]
mu <- rowMeans(mm); ord <- order(rr, -mu)
mm <- mm[ord, , drop = FALSE]; rr <- rr[ord]
keep <- !duplicated(rr)
mat91061 <- mm[keep, , drop = FALSE]; rownames(mat91061) <- rr[keep]
w('[GSE91061] symbol matrix: %d genes x %d samples', nrow(mat91061), ncol(mat91061))

sm9 <- readLines(gzfile('GSE91061_series_matrix.txt.gz'))
grab9 <- function(pat) gsub('"', '', strsplit(gsub('^[^\t]*\t', '', sm9[grepl(pat, sm9)][1]), '\t')[[1]])
t9 <- grab9('^!Sample_title'); v9 <- sub('.*: *', '', grab9('^!Sample_characteristics_ch1.*visit'))
r9 <- sub('.*: *', '', grab9('^!Sample_characteristics_ch1.*response'))
clin91061 <- data.frame(col = colnames(mat91061), patient = sub('_.*$', '', t9), visit = v9, response = r9, stringsAsFactors = FALSE)
stopifnot(identical(clin91061$col, colnames(mat91061)))
write.table(clin91061, file.path(OUT, 'icb_clin_GSE91061.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
w('[GSE91061] visit x response:'); print(table(clin91061$visit, clin91061$response))
s91061 <- score_cohort(mat91061, 'GSE91061')

saveRDS(list(clin78220 = clin78220, clin91061 = clin91061, s78220 = s78220, s91061 = s91061),
        file.path(OUT, 'icb_scores.rds'))
}

## ================= 3. 通用检验函数 =================
wtest <- function(score_mat, group, tag = '') {
  out <- data.frame()
  for (s in rownames(score_mat)) {
    x <- score_mat[s, group == 'R']; y <- score_mat[s, group == 'NR']
    x <- x[is.finite(x)]; y <- y[is.finite(y)]
    if (length(x) < 3 || length(y) < 3) next
    p <- suppressWarnings(wilcox.test(x, y)$p.value)
    dd <- (median(x) - median(y)) / sd(c(x, y))
    out <- rbind(out, data.frame(feature = s, n_R = length(x), n_NR = length(y),
      med_R = round(median(x), 3), med_NR = round(median(y), 3), delta = round(median(x) - median(y), 3),
      cohen_d = round(dd, 3), p = signif(p, 4)))
  }
  if (nrow(out)) out$q <- signif(p.adjust(out$p, 'BH'), 4)
  out
}

## ================= 4. 关联：ALMP 轴 + MCP + quanTIseq × 应答 =================
res_list <- list()
add <- function(tag, cohort, block) {
  for (nm in names(block)) {
    m <- block[[nm]]; if (is.null(m) || !nrow(m)) next
    r <- wtest(as.matrix(m), tg); if (nrow(r)) res_list[[paste(tag, nm)]] <<- cbind(cohort = cohort, readout = nm, r)
  }
}
# (a) GSE78220 预处理前
sub <- clin78220[clin78220$visit == 'baseline' & !is.na(clin78220$bin), ]; tg <- sub$bin
add('GSE78220_sample', 'GSE78220(pre, sample)', list(ssGSEA = s78220$ss, GSVA = s78220$gv, MCPcounter = s78220$mc, quanTIseq = s78220$qt, quanTIseq_CLR = s78220$qtclr))
# (a2) GSE78220 患者级（Pt27A/B 折叠）
subp <- sub[sub$patient != 'Pt27B', ]
two <- sub$col[sub$patient %in% c('Pt27A', 'Pt27B')]
fold <- function(m) { m2 <- m[, subp$col, drop = FALSE]; if (all(two %in% sub$col)) m2[, 'Pt27A.baseline'] <- rowMeans(m[, two, drop = FALSE]); m2 }
tg <- subp$bin
add('GSE78220_patientlevel', 'GSE78220(pre, patient-level)', list(ssGSEA = fold(s78220$ss), GSVA = fold(s78220$gv), MCPcounter = fold(as.matrix(s78220$mc)), quanTIseq = fold(as.matrix(s78220$qt)), quanTIseq_CLR = fold(as.matrix(s78220$qtclr))))
# (b) GSE91061 预处理前 PRCR vs PD
sub9 <- clin91061[clin91061$visit == 'Pre' & clin91061$response %in% c('PRCR', 'PD'), ]; tg <- ifelse(sub9$response == 'PRCR', 'R', 'NR')
add('GSE91061_pre_PRCRvsPD', 'GSE91061(pre, PRCR vs PD)', list(ssGSEA = s91061$ss, GSVA = s91061$gv, MCPcounter = s91061$mc, quanTIseq = s91061$qt, quanTIseq_CLR = s91061$qtclr))
# (c) GSE91061 预处理前 DCR（PRCR+SD vs PD）
sub9b <- clin91061[clin91061$visit == 'Pre' & clin91061$response %in% c('PRCR', 'SD', 'PD'), ]; tg <- ifelse(sub9b$response == 'PD', 'NR', 'R')
add('GSE91061_pre_DCR', 'GSE91061(pre, PRCR+SD vs PD)', list(ssGSEA = s91061$ss, GSVA = s91061$gv, MCPcounter = s91061$mc, quanTIseq = s91061$qt, quanTIseq_CLR = s91061$qtclr))
allres <- do.call(rbind, res_list)
write.table(allres, file.path(OUT, 'icb_response_assoc.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
w('\n===== 应答关联（全部 %d 条；q<0.1）=====', nrow(allres))
print(allres[allres$q < 0.1, ], row.names = FALSE)

## ================= 5. SEG1_ADIPOKINE_CORE × 免疫读出（跨算法）=================
coup <- data.frame()
for (nm in c('GSE78220', 'GSE91061')) {
  s <- if (nm == 'GSE78220') s78220 else s91061
  blocks <- list(ssGSEA = s$ss, MCPcounter = s$mc, quanTIseq = s$qt, quanTIseq_CLR = s$qtclr)
  for (ax in c('SEG1_ADIPOKINE_CORE', 'SEG2_LIPID_PATHWAY', 'SEG3_EXECUTOR')) {
    if (is.null(s$ss) || !(ax %in% rownames(s$ss))) next
    a <- s$ss[ax, ]
    for (rd in names(blocks)) {
      m <- blocks[[rd]]; if (is.null(m)) next
      for (ft in rownames(m)) {
        v <- as.numeric(m[ft, ])
        if (sum(is.finite(v)) < 5 || isTRUE(sd(v, na.rm = TRUE) == 0)) next
        ct <- suppressWarnings(cor.test(a, v, method = 'spearman'))
        coup <- rbind(coup, data.frame(cohort = nm, axis = ax, readout = rd, feature = ft,
                                       rho = round(unname(ct$estimate), 3), p = signif(ct$p.value, 4)))
      }
    }
  }
}
coup$q <- signif(p.adjust(coup$p, 'BH'), 4)
write.table(coup, file.path(OUT, 'icb_readout_axis_corr.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
w('\n===== SEG1_ADIPOKINE_CORE × 免疫读出（|rho|>=0.4, q<0.05）=====')
print(coup[coup$axis == 'SEG1_ADIPOKINE_CORE' & abs(coup$rho) >= 0.4 & coup$q < 0.05, ], row.names = FALSE)

## ================= 6. GSE91061 配对 On-Pre delta =================
paired <- clin91061[clin91061$response %in% c('PRCR', 'PD', 'SD'), ]
pat <- unique(paired$patient)
pat <- pat[sapply(pat, function(p)
  sum(paired$patient == p & paired$visit == 'Pre') == 1 && sum(paired$patient == p & paired$visit == 'On') == 1)]
delta_res <- data.frame()
for (rd in c('ssGSEA', 'GSVA')) {
  sc <- if (rd == 'ssGSEA') s91061$ss else s91061$gv
  if (is.null(sc)) next
  for (ax in rownames(sc)) {
    dv <- sapply(pat, function(p) sc[ax, paired$col[paired$patient == p & paired$visit == 'On']][1] -
                                     sc[ax, paired$col[paired$patient == p & paired$visit == 'Pre']][1])
    grp <- paired$response[match(pat, paired$patient)]
    x <- dv[grp == 'PRCR']; y <- dv[grp == 'PD']
    if (sum(is.finite(x)) >= 3 && sum(is.finite(y)) >= 3)
      delta_res <- rbind(delta_res, data.frame(readout = rd, set = ax, n_R = sum(is.finite(x)), n_NR = sum(is.finite(y)),
        med_delta_R = round(median(x, na.rm = TRUE), 3), med_delta_NR = round(median(y, na.rm = TRUE), 3),
        p = signif(suppressWarnings(wilcox.test(x, y)$p.value), 4)))
  }
}
if (nrow(delta_res)) delta_res$q <- signif(p.adjust(delta_res$p, 'BH'), 4)
write.table(delta_res, file.path(OUT, 'icb_gse91061_delta_assoc.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
w('\n===== GSE91061 On-Pre delta vs 应答（n_paired=%d）=====', length(pat))
print(delta_res, row.names = FALSE)

saveRDS(list(assoc = allres, coup = coup, delta = delta_res), file.path(OUT, 'icb_results.rds'))
writeLines(capture.output(sessionInfo()), file.path(OUT, 'icb_sessionInfo.txt'))
w('\nDONE')
