

library(edgeR)
library(vsn)
library(ggplot2)
library(DESeq2)
library(Matrix)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyr)
library(tidyverse)
library(DEGreport)
library(reshape2)
library(RColorBrewer)
library(enrichplot)
library(msigdbr)
library(fgsea)
library(ggsci)


mcount<-read.csv("merge_count.csv")

colnames(mcount)<-c("Gene","0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")


######################
#####################

mcount_deg<-mcount
rownames(mcount_deg)<-mcount$Gene
mcount_deg<-mcount_deg[,-1]

group <- factor(c("0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")) 
d <- DGEList(counts = mcount_deg, group = group)
d <- calcNormFactors(d,method="TMM")
d$samples$lib.size
d$samples$norm.factors

for(i in 1:13){
  mcount_deg[,i]<-mcount_deg[,i]/d$samples$norm.factors[i]
}

genelis<-rownames(mcount_deg)

CTL_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H","2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")]
colnames(CTL_RNA_all)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")
CTL_RNA_all<-as.matrix(CTL_RNA_all)

INH_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H","2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")]
colnames(INH_RNA_all)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")
INH_RNA_all<-as.matrix(INH_RNA_all)

#apply(CTL_RNA_all,2,mean)
#apply(INH_RNA_all,2,mean)

CTL_RNA_all_sub<-CTL_RNA_all[apply(CTL_RNA_all,1,min)>1,]
INH_RNA_all_sub<-INH_RNA_all[apply(INH_RNA_all,1,min)>1,]

CTL_RNA_all_FC<-log2(CTL_RNA_all_sub/CTL_RNA_all_sub[,1])
INH_RNA_all_FC<-log2(INH_RNA_all_sub/INH_RNA_all_sub[,1])

allgene_FC_ctl<-apply(CTL_RNA_all_FC,2,mean)
allgene_FC_inh<-apply(INH_RNA_all_FC,2,mean)
apply(CTL_RNA_all_FC,2,mean)
apply(INH_RNA_all_FC,2,mean)
apply(INH_RNA_all_FC,2,mean)-apply(CTL_RNA_all_FC,2,mean)

CTLINH_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")]
colnames(CTLINH_RNA_all)<-c("0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")
CTLINH_RNA_all<-as.matrix(CTLINH_RNA_all)
CTLINH_RNA_all_sub<-CTLINH_RNA_all[apply(CTLINH_RNA_all,1,min)>1,]


write.csv(CTLINH_RNA_all_sub,"mergecount_TPMedgeR_min1.csv",quote=F)
######################
