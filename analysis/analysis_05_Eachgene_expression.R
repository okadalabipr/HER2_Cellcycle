setwd("D:/ichikawa-sensei")

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
library(dorothea)
library(rcompanion)


mcount<-read.csv("eachcount_edgeR.csv")

colnames(mcount)<-c("gene","0H_1","0H_2","0H_3","12H_1","12H_2","12H_3","16H_CDK4_1",
                    "16H_DMSO_1","16H_CDK4_2","16H_DMSO_2","16H_CDK4_3","16H_DMSO_3",
                    "20H_CDK4_1","20H_DMSO_1","20H_CDK4_2","20H_DMSO_2","20H_CDK4_3",
                    "20H_DMSO_3","24H_CDK4_1","24H_DMSO_1","24H_CDK4_2","24H_DMSO_2",
                    "24H_CDK4_3","24H_DMSO_3","28H_CDK4_1","28H_DMSO_1","28H_CDK4_2",
                    "28H_DMSO_2","28H_CDK4_3","28H_DMSO_3","2H_1","2H_2","2H_3",
                    "4H_1","4H_2","4H_3","8H_1","8H_2","8H_3")


agmx<-function(x){
  xdf<-as.data.frame(x)
  xdf<-t(xdf)
  return(argmax(xdf))
}  


##########

######################
#####################
#CTLINH#################
##########

mcount_deg<-mcount
rownames(mcount_deg)<-mcount$gene
mcount_deg<-mcount_deg[,-1]

genelis<-rownames(mcount_deg)

CTL_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3",
                                                          "4H_1","4H_2","4H_3","8H_1","8H_2","8H_3",
                                                          "12H_1","12H_2","12H_3",
                                                          "16H_DMSO_1","16H_DMSO_2","16H_DMSO_3",
                                                          "20H_DMSO_1","20H_DMSO_2",
                                                          "20H_DMSO_3","24H_DMSO_1","24H_DMSO_2",
                                                          "24H_DMSO_3","28H_DMSO_1",
                                                          "28H_DMSO_2","28H_DMSO_3")]
#colnames(CTL_RNA_all)<-c("0H","0H","0H","2H","2H","2H","4H","4H","4H",
#                         "8H","8H","8H","12H","12H","12H","16H","16H","16H",
#                         "20H","20H","20H","24H","24H","24H","28H","28H","28H")
colnames(CTL_RNA_all)<-c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3","4H_1","4H_2","4H_3",
                         "8H_1","8H_2","8H_3","12H_1","12H_2","12H_3","16H_1","16H_2","16H_3",
                         "20H_1","20H_2","20H_3","24H_1","24H_2","24H_3","28H_1","28H_2","28H_3")
CTL_RNA_all<-as.matrix(CTL_RNA_all)

INH_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3",
                                                          "4H_1","4H_2","4H_3","8H_1","8H_2","8H_3",
                                                          "12H_1","12H_2","12H_3",
                                                          "16H_CDK4_1","16H_CDK4_2","16H_CDK4_3",
                                                          "20H_CDK4_1","20H_CDK4_2",
                                                          "20H_CDK4_3","24H_CDK4_1","24H_CDK4_2",
                                                          "24H_CDK4_3","28H_CDK4_1",
                                                          "28H_CDK4_2","28H_CDK4_3")]
colnames(INH_RNA_all)<-c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3","4H_1","4H_2","4H_3",
                         "8H_1","8H_2","8H_3","12H_1","12H_2","12H_3","16H_1","16H_2","16H_3",
                         "20H_1","20H_2","20H_3","24H_1","24H_2","24H_3","28H_1","28H_2","28H_3")
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

CTLINH_RNA_all<-mcount_deg[rownames(mcount_deg)%in%genelis,c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3",
                                                             "4H_1","4H_2","4H_3","8H_1","8H_2","8H_3",
                                                             "12H_1","12H_2","12H_3",
                                                             "16H_DMSO_1","16H_DMSO_2","16H_DMSO_3","16H_CDK4_1","16H_CDK4_2","16H_CDK4_3",
                                                             "20H_DMSO_1","20H_DMSO_2","20H_DMSO_3",
                                                             "20H_CDK4_1","20H_CDK4_2","20H_CDK4_3",
                                                             "24H_DMSO_1","24H_DMSO_2","24H_DMSO_3",
                                                             "24H_CDK4_1","24H_CDK4_2","24H_CDK4_3",
                                                             "28H_DMSO_1","28H_DMSO_2","28H_DMSO_3","28H_CDK4_1",
                                                             "28H_CDK4_2","28H_CDK4_3")]
colnames(CTLINH_RNA_all)<-c("0H_1","0H_2","0H_3","2H_1","2H_2","2H_3",
                            "4H_1","4H_2","4H_3","8H_1","8H_2","8H_3",
                            "12H_1","12H_2","12H_3",
                            "16H_DMSO_1","16H_DMSO_2","16H_DMSO_3","16H_CDK4_1","16H_CDK4_2","16H_CDK4_3",
                            "20H_DMSO_1","20H_DMSO_2","20H_DMSO_3",
                            "20H_CDK4_1","20H_CDK4_2","20H_CDK4_3",
                            "24H_DMSO_1","24H_DMSO_2","24H_DMSO_3",
                            "24H_CDK4_1","24H_CDK4_2","24H_CDK4_3",
                            "28H_DMSO_1","28H_DMSO_2","28H_DMSO_3","28H_CDK4_1",
                            "28H_CDK4_2","28H_CDK4_3")
CTLINH_RNA_all<-as.matrix(CTLINH_RNA_all)


CTLINH_RNA_all_matome<-apply(CTLINH_RNA_all[,1:3],1,mean)
for(i in 1:12){
  stn<-i*3+1
  enn<-(i+1)*3
  tmpa<-apply(CTLINH_RNA_all[,stn:enn],1,mean)
  CTLINH_RNA_all_matome<-cbind(CTLINH_RNA_all_matome,tmpa)
}
colnames(CTLINH_RNA_all_matome)<-c("0H","2H","4H","8H","12H","16H_DMSO","16H_CDK4",
                                    "20H_DMSO","20H_CDK4","24H_DMSO","24H_CDK4","28H_DMSO","28H_CDK4")

CTLINH_RNA_all_matome_sub<-CTLINH_RNA_all_matome[apply(CTLINH_RNA_all_matome,1,min)>1,]
CTLINH_RNA_all_sub<-CTLINH_RNA_all[rownames(CTLINH_RNA_all) %in% rownames(CTLINH_RNA_all_matome_sub),]
CTLINH_RNA_all_sub<-CTLINH_RNA_all_sub[apply(CTLINH_RNA_all_sub,1,min)>0,]

#write.csv(CTLINH_RNA_all_sub,"mergecount_TPMedgeR_min1_each.csv",quote=F)

#############論文用　単品遺伝子


setgenelis="PCNA"
CTL_RNA_all_sub_set<-CTLINH_RNA_all_sub[rownames(CTLINH_RNA_all_sub) %in% setgenelis,c(1:15,16:18,22:24,28:30,34:36)]
INH_RNA_all_sub_set<-CTLINH_RNA_all_sub[rownames(CTLINH_RNA_all_sub) %in% setgenelis,c(1:15,19:21,25:27,31:33,37:39)]

CTL_RNA_all_sub_set<-log2(CTL_RNA_all_sub_set)
INH_RNA_all_sub_set<-log2(INH_RNA_all_sub_set)

CTL_RNA_all_FC_set_each<-CTL_RNA_all_sub_set
CTL_RNA_all_FC_set_each[c(1,4,7,10,13,16,19,22,25)]<-CTL_RNA_all_FC_set_each[c(1,4,7,10,13,16,19,22,25)]-CTL_RNA_all_FC_set_each[1]
CTL_RNA_all_FC_set_each[c(2,5,8,11,14,17,20,23,26)]<-CTL_RNA_all_FC_set_each[c(2,5,8,11,14,17,20,23,26)]-CTL_RNA_all_FC_set_each[2]
CTL_RNA_all_FC_set_each[c(3,6,9,12,15,18,21,24,27)]<-CTL_RNA_all_FC_set_each[c(3,6,9,12,15,18,21,24,27)]-CTL_RNA_all_FC_set_each[3]

INH_RNA_all_FC_set_each<-INH_RNA_all_sub_set
INH_RNA_all_FC_set_each[c(1,4,7,10,13,16,19,22,25)]<-INH_RNA_all_FC_set_each[c(1,4,7,10,13,16,19,22,25)]-INH_RNA_all_FC_set_each[1]
INH_RNA_all_FC_set_each[c(2,5,8,11,14,17,20,23,26)]<-INH_RNA_all_FC_set_each[c(2,5,8,11,14,17,20,23,26)]-INH_RNA_all_FC_set_each[2]
INH_RNA_all_FC_set_each[c(3,6,9,12,15,18,21,24,27)]<-INH_RNA_all_FC_set_each[c(3,6,9,12,15,18,21,24,27)]-INH_RNA_all_FC_set_each[3]

allgene_FC_set_ctl_each<-CTL_RNA_all_FC_set_each
allgene_FC_set_inh_each<-INH_RNA_all_FC_set_each

allgene_FC_set_ctl_each
allgene_FC_set_inh_each

CTL_RNA_all_FC_each_mean<-mean(allgene_FC_set_ctl_each[1:3])
for(i in 1:8){
  stn<-i*3+1
  enn<-(i+1)*3
  tmpa<-mean(allgene_FC_set_ctl_each[stn:enn])
  CTL_RNA_all_FC_each_mean<-cbind(CTL_RNA_all_FC_each_mean,tmpa)
}
colnames(CTL_RNA_all_FC_each_mean)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")


INH_RNA_all_FC_each_mean<-mean(allgene_FC_set_inh_each[1:3])
for(i in 1:8){
  stn<-i*3+1
  enn<-(i+1)*3
  tmpa<-mean(allgene_FC_set_inh_each[stn:enn])
  INH_RNA_all_FC_each_mean<-cbind(INH_RNA_all_FC_each_mean,tmpa)
}
colnames(INH_RNA_all_FC_each_mean)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")


CTL_RNA_all_FC_each_sd<-sd(allgene_FC_set_ctl_each[1:3])
for(i in 1:8){
  stn<-i*3+1
  enn<-(i+1)*3
  tmpa<-sd(allgene_FC_set_ctl_each[stn:enn])
  CTL_RNA_all_FC_each_sd<-cbind(CTL_RNA_all_FC_each_sd,tmpa)
}
colnames(CTL_RNA_all_FC_each_sd)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")


INH_RNA_all_FC_each_sd<-sd(allgene_FC_set_inh_each[1:3])
for(i in 1:8){
  stn<-i*3+1
  enn<-(i+1)*3
  tmpa<-sd(allgene_FC_set_inh_each[stn:enn])
  INH_RNA_all_FC_each_sd<-cbind(INH_RNA_all_FC_each_sd,tmpa)
}
colnames(INH_RNA_all_FC_each_sd)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")

CTL_RNA_all_FC_each_mean<-as.data.frame(t(CTL_RNA_all_FC_each_mean))
INH_RNA_all_FC_each_mean<-as.data.frame(t(INH_RNA_all_FC_each_mean))

CTL_RNA_all_FC_each_sd<-as.data.frame(t(CTL_RNA_all_FC_each_sd))
INH_RNA_all_FC_each_sd<-as.data.frame(t(INH_RNA_all_FC_each_sd))

df<-cbind(INH_RNA_all_FC_each_mean,CTL_RNA_all_FC_each_mean)
colnames(df)<-c("INH","CTL")
rownames(df)<-c(0,2,4,8,12,16,20,24,28)
df<-as.data.frame(cbind(rownames(df),df))
y <- melt(df)
colnames(y) <- c("hr", "treat", "log2FC")


dfsd<-cbind(INH_RNA_all_FC_each_sd,CTL_RNA_all_FC_each_sd)
colnames(dfsd)<-c("INH","CTL")
rownames(dfsd)<-c("0H","2H","4H","8H","12H","16H","20H","24H","28H")
dfsd<-as.data.frame(cbind(rownames(dfsd),dfsd))
ysd <- melt(dfsd)
colnames(ysd) <- c("hr", "treat", "sd")

y<-as.data.frame(cbind(y,ysd$sd))
colnames(y)<-c("hr", "treat", "log2FC","sd")
y$hr<-as.numeric(y$hr)

g <- ggplot(y, aes(x = hr, y = log2FC, color = treat))
g <- g + geom_line(size = 2)+geom_point(size = 5)
#g <- g + geom_errorbar(aes(ymin = log2FC - sd, ymax = log2FC + sd),width =.3)
g <- g + scale_color_nejm()
g <- g + ggtitle("PCNA")
g <- g + xlab("Time(hr)")
g <- g +theme_bw()
g <- g + ylab(expression(paste({Log[2]},"Fold Change")))
g <- g + scale_x_continuous(breaks=c(0,2,4,8,12,16,20,24,28))
g <- g +theme(text = element_text("Arial"),
              legend.position = "none",
              legend.title=element_text(size = 24),
              legend.text=element_text(size = 24),
              axis.title.x=element_text(size = 24),
              axis.title.y=element_text(size = 24),
              axis.text.x=element_text(size = 24),
              axis.text.y=element_text(size = 24),
              plot.title = element_text(size = 40,hjust = 0.5),
              panel.grid = element_blank())
g <- g + annotate("text", x=16,   y=0, label="Control", color="blue", size=10)
g <- g + annotate("text", x=16,   y=1.3, label="CDK4i", color="red", size=10)

plot(g)


png(paste0("",setgenelis,"_FC_eb.png"), width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off() 

