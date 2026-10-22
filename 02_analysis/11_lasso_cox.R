setwd("C:/Rtmp/p2t3")
suppressMessages({library(data.table);library(survival);library(glmnet);library(jsonlite)})
con<-file("p2t3_final.log","w"); w<-function(...){s<-paste0(...);cat(s,"\n");cat(s,"\n",file=con);flush(con)}

d<-fread("mrna_data.tsv"); samp<-fread("clinical_SAMPLE.tsv"); pat<-fread("clinical_PATIENT.tsv")
map<-fread("ALMP_ID_mapping_v1.csv"); sym<-setNames(map$gene_symbol, as.character(map$entrez))
M<-dcast(d,entrezGeneId~sampleId,value.var="value"); g<-as.character(M$entrezGeneId)
M<-as.matrix(M[,-1]); rownames(M)<-g
st<-samp[attr=="SAMPLE_TYPE",.(sampleId,sample_type=value)]
os<-dcast(pat[attr %in% c("OS_MONTHS","OS_STATUS")],patientId~attr,value.var="value")
os[,OS_MONTHS:=as.numeric(OS_MONTHS)][,event:=as.integer(sub(":.*","",OS_STATUS))][,OS_STATUS:=NULL]
s2p<-unique(d[,.(sampleId,patientId)]); ann<-merge(data.table(sampleId=colnames(M)),s2p,by="sampleId")
ann<-merge(ann,st,by="sampleId"); ann<-merge(ann,os,by="patientId")
ann2<-ann[!is.na(sample_type)&!is.na(OS_MONTHS)&OS_MONTHS>0]

SETS<-list(cohort="all", metastasis="Metastasis", primary="Primary")
THR<-0.05; NREP<-20

for(tag in names(SETS)){
  a<-if(tag=="cohort") ann2 else ann2[sample_type==SETS[[tag]]]
  X<-M[,a$sampleId]; lg<-log2(X+1); keep<-rowSums(X>=10)/ncol(X)>=0.5
  xall<-lg[keep,]; y<-Surv(a$OS_MONTHS,a$event)
  pv<-apply(xall,1,function(e){if(sd(e)==0) return(NA); summary(coxph(y~e))$coefficients[1,"Pr(>|z|)"]})
  fdr<-p.adjust(pv,"BH")
  uni<-data.table(entrez=rownames(xall),symbol=sym[rownames(xall)],
                  HR=apply(xall,1,function(e){if(sd(e)==0) return(NA); exp(summary(coxph(y~e))$coefficients[1,"coef"])}),
                  p=pv,FDR=fdr)
  uni<-uni[order(p)]
  fwrite(uni, sprintf("uni_%s.tsv",tag), sep="\t")

  sel<-names(pv)[!is.na(pv)&pv<THR]
  x<-scale(t(xall[sel,,drop=FALSE])); x[is.na(x)]<-0
  set.seed(20261012)
  fold<-sample(rep(1:10,length.out=nrow(x)))
  cv<-cv.glmnet(x,y,family="cox",alpha=1,foldid=fold,type.measure="deviance")
  imin<-cv$index[1]
  w(sprintf("===== %s n=%d kept=%d | uni p<0.05=%d FDR<0.05=%d | LASSO input=%d =====",
            tag,ncol(xall),sum(keep),sum(pv<0.05,na.rm=T),sum(fdr<0.05,na.rm=T),length(sel)))
  w(sprintf("  单次CV诊断: nlambda=%d  lambda.min=%.5f (idx %d, nzero=%d)  lambda.1se=%.5f (idx %d, nzero=%d)",
            length(cv$lambda), cv$lambda.min, imin, sum(as.numeric(coef(cv,s="lambda.min"))!=0),
            cv$lambda.1se, cv$index[2], sum(as.numeric(coef(cv,s="lambda.1se"))!=0)))
  w(sprintf("  cvm[1]=%.4f (null)  cvm[imin]=%.4f  cvsd[imin]=%.4f  →  1se 门槛=%.4f  %s",
            cv$cvm[1], cv$cvm[imin], cv$cvsd[imin], cv$cvm[imin]+cv$cvsd[imin],
            ifelse(cv$cvm[1] <= cv$cvm[imin]+cv$cvsd[imin], "⇒ 1SE 落入空模型【已确认】","⇒ 正常")))
  cm<-coef(cv,s="lambda.min"); v<-as.numeric(cm); nm<-rownames(cm)
  tab<-data.table(entrez=nm,symbol=sym[nm],coef=v)[v!=0][order(-abs(coef))]
  fwrite(tab, sprintf("lasso_min_%s.tsv",tag), sep="\t")
  w(sprintf("  → lambda.min 选中 %d 个基因:", nrow(tab)))
  w(sprintf("     %s", paste(tab$symbol, collapse=", ")))

  set.seed(20261012)
  cnt1<-list(); cntm<-list(); n1<-c(); nmv<-c(); l1<-c(); lm<-c()
  for(r in 1:NREP){
    foldid<-sample(rep(1:10,length.out=nrow(x)))
    c2<-cv.glmnet(x,y,family="cox",alpha=1,foldid=foldid,type.measure="deviance")
    c1<-coef(c2,s="lambda.1se"); c0<-coef(c2,s="lambda.min")
    a1<-rownames(c1)[as.numeric(c1)!=0]; a0<-rownames(c0)[as.numeric(c0)!=0]
    n1<-c(n1,length(a1)); nmv<-c(nmv,length(a0)); l1<-c(l1,c2$lambda.1se); lm<-c(lm,c2$lambda.min)
    for(z in a1) cnt1[[z]]<-c(cnt1[[z]],r)
    for(z in a0) cntm[[z]]<-c(cntm[[z]],r)
  }
  st1<-names(cnt1)[sapply(cnt1,length)>=10]; stm<-names(cntm)[sapply(cntm,length)>=10]
  w(sprintf("  重复CV(20x10, pre=%.2f): lambda.1se nzero 中位=%d (范围 %d-%d, 空集 %d/20)",
            THR,median(n1),min(n1),max(n1),sum(n1==0)))
  w(sprintf("  重复CV(20x10, pre=%.2f): lambda.min nzero 中位=%d (范围 %d-%d, 空集 %d/20)",
            THR,median(nmv),min(nmv),max(nmv),sum(nmv==0)))
  w(sprintf("  ★lambda.min 下 在>=10/20 次被选中 (%d): %s", length(stm), paste(sym[stm],collapse=", ")))
  w(sprintf("  ★lambda.1se 下 在>=10/20 次被选中 (%d): %s", length(st1), ifelse(length(st1)==0,"无",paste(sym[st1],collapse=", "))))
  json<-list(tag=tag,n=ncol(xall),kept=sum(keep),uni_p05=sum(pv<0.05,na.rm=T),
    uni_fdr05=sum(fdr<0.05,na.rm=T),lasso_input=length(sel),
    single_lambda_min=cv$lambda.min, single_lambda_1se=cv$lambda.1se,
    single_nzero_min=sum(as.numeric(coef(cv,s="lambda.min"))!=0),
    single_nzero_1se=sum(as.numeric(coef(cv,s="lambda.1se"))!=0),
    single_genes=tab$symbol, single_coef=tab$coef,
    rep_lambda_min_nzero_median=median(nmv), rep_lambda_min_nzero_range=c(min(nmv),max(nmv)),
    rep_lambda_1se_nzero_median=median(n1), rep_lambda_1se_nzero_range=c(min(n1),max(n1)),
    rep_lambda_1se_empty=sum(n1==0),
    stable_lambda_min=sym[stm], stable_lambda_min_n=sapply(cntm[stm],length),
    stable_lambda_1se=sym[st1])
  write(toJSON(json,auto_unbox=TRUE,pretty=TRUE), sprintf("p2t3_%s.json",tag))
}
close(con)
