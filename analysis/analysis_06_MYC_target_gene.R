

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
library(gridExtra)
library(stats)

library(extrafont)
font_import(pattern="arial")
windowsFonts(sans="arial")
loadfonts(device="win")
names(windowsFonts())


mcount<-read.csv("merge_count.csv")

colnames(mcount)<-c("Gene","0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")


##########

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))
view(regulons)



######################
#####################
#CTLINH#################
##########

mcount_deg<-mcount
rownames(mcount_deg)<-mcount$Gene
mcount_deg<-mcount_deg[,-1]

group <- factor(c("0H","2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")) 
d <- DGEList(counts = mcount_deg, group = group)
d <- calcNormFactors(d,method="TMM")
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



#############ここを論文で使用fig3bcd
CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"16H_CDK4"]/CTLINH_RNA_all_sub[,"16H_DMSO"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$size[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]

pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)


sum((pd$stats$stat/pd$maxAbsStat)*(pd$spreadES/4)<0)

g <- with(pd,
     ggplot(data=curve) +
       geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
       geom_ribbon(data=stats,
                   mapping=aes(x=rank, ymin=0,
                               ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
       scale_fill_manual(values=c("blue", "red"), name="fill") +
       geom_segment(data=ticks,
                    mapping=aes(x=rank, y=-spreadES/16,
                                xend=rank, yend=spreadES/16),
                    size=0.2) +
       geom_hline(yintercept=posES, colour="red", linetype="dashed") +
       geom_hline(yintercept=negES, colour="red", linetype="dashed") +
       geom_hline(yintercept=0, colour="black") +
       ggtitle("HALLMARK_MYC_TARGETS_V2") +
       theme_bw() +
       labs(x="rank", y="enrichment score")+
       theme(text = element_text("Arial"),
             legend.position = "none",
             legend.title=element_text(size = 24),
             legend.text=element_text(size = 24),
             axis.title.x=element_text(size = 24),
             axis.title.y=element_text(size = 24),
             axis.text.x=element_text(size = 24),
             axis.text.y=element_text(size = 24),
             plot.title = element_text(size = 30,hjust = 0.5),
             panel.grid = element_blank()
       ))
plot(g)

png("MYCtarget_gsea_INH_CTL_16H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()


###

CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"20H_CDK4"]/CTLINH_RNA_all_sub[,"20H_DMSO"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]

pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)

g <- with(pd,
          ggplot(data=curve) +
            geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
            geom_ribbon(data=stats,
                        mapping=aes(x=rank, ymin=0,
                                    ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
            scale_fill_manual(values=c("blue", "red"), name="fill") +
            geom_segment(data=ticks,
                         mapping=aes(x=rank, y=-spreadES/16,
                                     xend=rank, yend=spreadES/16),
                         size=0.2) +
            geom_hline(yintercept=posES, colour="red", linetype="dashed") +
            geom_hline(yintercept=negES, colour="red", linetype="dashed") +
            geom_hline(yintercept=0, colour="black") +
            ggtitle("HALLMARK_MYC_TARGETS_V2") +
            theme_bw() +
            labs(x="rank", y="enrichment score")+
            theme(text = element_text("Arial"),
                  legend.position = "none",
                  legend.title=element_text(size = 24),
                  legend.text=element_text(size = 24),
                  axis.title.x=element_text(size = 24),
                  axis.title.y=element_text(size = 24),
                  axis.text.x=element_text(size = 24),
                  axis.text.y=element_text(size = 24),
                  plot.title = element_text(size = 30,hjust = 0.5),
                  panel.grid = element_blank()
            ))


png("MYCtarget_gsea_INH_CTL_20H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()


###

CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"24H_CDK4"]/CTLINH_RNA_all_sub[,"24H_DMSO"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]


pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)


g <- with(pd,
          ggplot(data=curve) +
            geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
            geom_ribbon(data=stats,
                        mapping=aes(x=rank, ymin=0,
                                    ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
            scale_fill_manual(values=c("blue", "red"), name="fill") +
            geom_segment(data=ticks,
                         mapping=aes(x=rank, y=-spreadES/16,
                                     xend=rank, yend=spreadES/16),
                         size=0.2) +
            geom_hline(yintercept=posES, colour="red", linetype="dashed") +
            geom_hline(yintercept=negES, colour="red", linetype="dashed") +
            geom_hline(yintercept=0, colour="black") +
            ggtitle("HALLMARK_MYC_TARGETS_V2") +
            theme_bw() +
            labs(x="rank", y="enrichment score")+
            theme(text = element_text("Arial"),
                  legend.position = "none",
                  legend.title=element_text(size = 24),
                  legend.text=element_text(size = 24),
                  axis.title.x=element_text(size = 24),
                  axis.title.y=element_text(size = 24),
                  axis.text.x=element_text(size = 24),
                  axis.text.y=element_text(size = 24),
                  plot.title = element_text(size = 30,hjust = 0.5),
                  panel.grid = element_blank()
            ))


png("MYCtarget_gsea_INH_CTL_24H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()

###

CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"28H_CDK4"]/CTLINH_RNA_all_sub[,"28H_DMSO"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$size[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]

pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)


g <- with(pd,
          ggplot(data=curve) +
            geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
            geom_ribbon(data=stats,
                        mapping=aes(x=rank, ymin=0,
                                    ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
            scale_fill_manual(values=c("blue", "red"), name="fill") +
            geom_segment(data=ticks,
                         mapping=aes(x=rank, y=-spreadES/16,
                                     xend=rank, yend=spreadES/16),
                         size=0.2) +
            geom_hline(yintercept=posES, colour="red", linetype="dashed") +
            geom_hline(yintercept=negES, colour="red", linetype="dashed") +
            geom_hline(yintercept=0, colour="black") +
            ggtitle("HALLMARK_MYC_TARGETS_V2") +
            theme_bw() +
            labs(x="rank", y="enrichment score")+
            theme(text = element_text("Arial"),
                  legend.position = "none",
                  legend.title=element_text(size = 24),
                  legend.text=element_text(size = 24),
                  axis.title.x=element_text(size = 24),
                  axis.title.y=element_text(size = 24),
                  axis.text.x=element_text(size = 24),
                  axis.text.y=element_text(size = 24),
                  plot.title = element_text(size = 30,hjust = 0.5),
                  panel.grid = element_blank()
            ))

png("MYCtarget_gsea_INH_CTL_28H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()



##########一部タイムコース

CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"16H_DMSO"]/CTLINH_RNA_all_sub[,"12H"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$size[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]

pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)


g <- with(pd,
          ggplot(data=curve) +
            geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
            geom_ribbon(data=stats,
                        mapping=aes(x=rank, ymin=0,
                                    ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
            scale_fill_manual(values=c("blue", "red"), name="fill") +
            geom_segment(data=ticks,
                         mapping=aes(x=rank, y=-spreadES/16,
                                     xend=rank, yend=spreadES/16),
                         size=0.2) +
            geom_hline(yintercept=posES, colour="red", linetype="dashed") +
            geom_hline(yintercept=negES, colour="red", linetype="dashed") +
            geom_hline(yintercept=0, colour="black") +
            ggtitle("HALLMARK_MYC_TARGETS_V2") +
            theme_bw() +
            labs(x="rank", y="enrichment score")+
            theme(text = element_text("Arial"),
                  legend.position = "none",
                  legend.title=element_text(size = 24),
                  legend.text=element_text(size = 24),
                  axis.title.x=element_text(size = 24),
                  axis.title.y=element_text(size = 24),
                  axis.text.x=element_text(size = 24),
                  axis.text.y=element_text(size = 24),
                  plot.title = element_text(size = 30,hjust = 0.5),
                  panel.grid = element_blank()
            ))
plot(g)
png("MYCtarget_gsea_12H16H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()


##############


CTL_RNA_all_FC<-log2(CTLINH_RNA_all_sub[,"2H"]/CTLINH_RNA_all_sub[,"0H"])
original_gene_list <- CTL_RNA_all_FC
names(original_gene_list) <- rownames(CTLINH_RNA_all_sub)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)

msigdbr_df = msigdbr(species = "human", category = "H")
msigdbr_list = split(x = msigdbr_df$gene_symbol, f = msigdbr_df$gs_name)
fgseaRes <- fgsea(msigdbr_list, gene_list, minSize=15, maxSize = 500, nperm=1000)
head(fgseaRes[order(padj, -abs(NES)), ], n=20)
fgseaRes[ES > 0][head(order(padj, -abs(NES)), n=30), pathway]
fgseaRes[ES < 0][head(order(padj, -abs(NES)), n=30), pathway]

fgseaRes$NES[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$padj[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$pval[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]
fgseaRes$size[fgseaRes$pathway=="HALLMARK_MYC_TARGETS_V2"]

pd <- plotEnrichmentData(
  pathway = msigdbr_list[["HALLMARK_MYC_TARGETS_V2"]],
  stats = gene_list
)


g <- with(pd,
          ggplot(data=curve) +
            geom_line(aes(x=rank, y=ES), color="green", linewidth=2) +
            geom_ribbon(data=stats,
                        mapping=aes(x=rank, ymin=0,
                                    ymax=stat/maxAbsStat*(spreadES/4), fill = stat/maxAbsStat*(spreadES/4) > 0)) +
            scale_fill_manual(values=c("blue", "red"), name="fill") +
            geom_segment(data=ticks,
                         mapping=aes(x=rank, y=-spreadES/16,
                                     xend=rank, yend=spreadES/16),
                         size=0.2) +
            geom_hline(yintercept=posES, colour="red", linetype="dashed") +
            geom_hline(yintercept=negES, colour="red", linetype="dashed") +
            geom_hline(yintercept=0, colour="black") +
            ggtitle("HALLMARK_MYC_TARGETS_V2") +
            theme_bw() +
            labs(x="rank", y="enrichment score")+
            theme(text = element_text("Arial"),
                  legend.position = "none",
                  legend.title=element_text(size = 24),
                  legend.text=element_text(size = 24),
                  axis.title.x=element_text(size = 24),
                  axis.title.y=element_text(size = 24),
                  axis.text.x=element_text(size = 24),
                  axis.text.y=element_text(size = 24),
                  plot.title = element_text(size = 30,hjust = 0.5),
                  panel.grid = element_blank()
            ))
plot(g)
png("MYCtarget_gsea_0H2H.png", width = 900, height = 700)  # 描画デバイスを開く
plot(g)
dev.off()



