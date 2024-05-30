library(edgeR)
library(vsn)
library(ggplot2)
library(Matrix)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)
library(tidyr)
library(tidyverse)
library(DEGreport)
library(ggrepel)


mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host="https://jan2020.archive.ensembl.org")
gene.annotations <- biomaRt::getBM(mart = mart, attributes=c("hgnc_symbol", "entrezgene_id"))

inputdf <- read.table("output.txt", header=T, stringsAsFactors=F,sep="\t")
rownames(inputdf) <- inputdf$Geneid
colnames(inputdf)
calcdf<-inputdf
genelength<-calcdf$Length
genecount<-calcdf[-1:-6]

TPMcount<-1000*genecount/genelength
readnumber<-apply(TPMcount, 2, sum)
TPMcount<-t(TPMcount)
TPMcount<-1000000*TPMcount/readnumber
TPMcount<-t(TPMcount)
TPMcount<-as.data.frame(TPMcount)

TPMsearch<-cbind(rownames(TPMcount),TPMcount)

#7,9 8 11 13 15 17,19
ched<-c(1:6,10,12,14,16,18,20:39)
degCheckFactors(genecount[,ched])
degCheckFactors(genecount)

colnames(TPMcount)<-c("0H_1","0H_2","0H_3","12H_1","12H_2","12H_3","16H_CDK4_1","16H_DMSO_1","16H_CDK4_2","16H_DMSO_2","16H_CDK4_3","16H_DMSO_3","20H_CDK4_1","20H_DMSO_1","20H_CDK4_2","20H_DMSO_2","20H_CDK4_3","20H_DMSO_3","24H_CDK4_1","24H_DMSO_1","24H_CDK4_2","24H_DMSO_2","24H_CDK4_3","24H_DMSO_3","28H_CDK4_1","28H_DMSO_1","28H_CDK4_2","28H_DMSO_2","28H_CDK4_3","28H_DMSO_3","2H_1","2H_2","2H_3","4H_1","4H_2","4H_3","8H_1","8H_2","8H_3")
write.csv(TPMcount,"TPMcount.csv")

#########################################each normalize

mcount_deg<-TPMcount

group <- factor(c("0H","0H","0H","12H","12H","12H","16H_CDK4","16H_DMSO","16H_CDK4","16H_DMSO","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","20H_CDK4","20H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","24H_CDK4","24H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO","28H_CDK4","28H_DMSO","28H_CDK4","28H_DMSO","2H","2H","2H","4H","4H","4H","8H","8H","8H")) 
d <- DGEList(counts = mcount_deg, group = group)
d <- calcNormFactors(d,method="TMM")
d$samples$lib.size
d$samples$norm.factors
normfc<-d$samples$lib.size*d$samples$norm.factors
normfc
mcount_deg_norm<-mcount_deg
for(i in 1:ncol(TPMcount)){
  mcount_deg_norm[,i]<-(mcount_deg[,i]*1000000)/normfc[i]
}

write.csv(mcount_deg_norm,"eachcount_edgeR.csv")


#########################################
#DEG,GO,KEGG compare to 0H

#cond<-c("12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO","2H","4H","8H")
cond<-c("12H")

TPMcount<-as.data.frame(TPMcount)
alldegtable<-c()
for (c in cond) {
  conn=c
  countp<-dplyr::select(.data=TPMcount,starts_with(conn),starts_with("0H"))
  groupp <- factor(c(conn, conn, conn,"0H", "0H", "0H"))
  designp <- model.matrix(~ groupp)

  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))

  degtable<-Table.DGE[Table.DGE$FDR < 0.01,]
  write.csv(degtable,paste0("DEG/",conn,"_0H_DEG.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01]
  
  png(paste0("DEG/",conn,"_0H_MA.png"), width = 700, height = 500)
  plotSmear(dfresult.DGE, de.tags = is.degs, cex=0.3)
  dev.off() 

  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"

  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")

  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)

  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_0H_GO_MF.csv"))

  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_0H_GO_CC.csv"))

  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_0H_GO_BP.csv"))

  kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  reskk <- summary(kk2)
  write.csv(reskk,paste0("GO_KEGG/",conn,"_0H_KEGG.csv"))
}

#########################################
#DEG,GO,KEGG compare to 12H

cond<-c("16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")

TPMcount<-as.data.frame(TPMcount)
alldegtable<-c()
for (c in cond) {
  conn=c
  countp<-dplyr::select(.data=TPMcount,starts_with(conn),starts_with("12H"))
  groupp <- factor(c(conn, conn, conn,"12H", "12H", "12H"))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01,]
  write.csv(degtable,paste0("DEG/",conn,"_12H_DEG.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01]
  
  png(paste0("DEG/",conn,"_12H_MA.png"), width = 700, height = 500)
  plotSmear(dfresult.DGE, de.tags = is.degs, cex=0.3)
  dev.off() 
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_12H_GO_MF.csv"))
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_12H_GO_CC.csv"))
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG/",conn,"_12H_GO_BP.csv"))
  
  kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  reskk <- summary(kk2)
  write.csv(reskk,paste0("GO_KEGG/",conn,"_12H_KEGG.csv"))
}

#######################
#DEG,GO,KEGG compare CDK-DMSO

cona<-c("16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
conb<-c("16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
alldegtable<-c()

for (c in 1:4) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01,]
  write.csv(degtable,paste0("DEG/",cona[c],"_",conb[c],"_DEG.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  #ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  #res <- summary(ego_result)
  #write.csv(res,paste0("GO_KEGG/",cona[c],"_",conb[c],"_GO_MF.csv"))
  
  #ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  #res <- summary(ego_result)
  #write.csv(res,paste0("GO_KEGG/",cona[c],"_",conb[c],"_GO_CC.csv"))
  
  #ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  #res <- summary(ego_result)
  #write.csv(res,paste0("GO_KEGG/",cona[c],"_",conb[c],"_GO_BP.csv"))
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG/",cona[c],"_",conb[c],"_KEGG.csv"))
}

##############

#DEG,GO,KEGG compare CDK-DMSO
cona<-c("16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
conb<-c("16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC>0,]
  write.csv(degtable,paste0("DEG_INH/",cona[c],"_",conb[c],"_DEG_2ndup.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC>0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_MF_2ndup.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_MF_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_CC_2ndup.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_CC_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_BP_2ndup.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_BP_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2ndup.csv"))
}

cona<-c("16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
conb<-c("16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")

alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC<0,]
  write.csv(degtable,paste0("DEG_INH/",cona[c],"_",conb[c],"_DEG_2nddown.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC<0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_MF_2nddown.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_MF_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_CC_2nddown.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_CC_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_BP_2nddown.csv"))
  
  png(paste0("GO_KEGG_INH/",cona[c],"_",conb[c],"_GO_BP_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2nddown.csv"))
}

#######################
#DEG,GO,KEGG compare next control

cona<-c("0H","2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO")
conb<-c("2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC>0,]
  write.csv(degtable,paste0("DEG_next/",cona[c],"_",conb[c],"_DEG_2ndup.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC>0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_MF_2ndup.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_MF_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_CC_2ndup.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_CC_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_BP_2ndup.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_BP_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2ndup.csv"))
}



#DEG,GO,KEGG compare next control

cona<-c("0H","2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO")
conb<-c("2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC<0,]
  write.csv(degtable,paste0("DEG_next/",cona[c],"_",conb[c],"_DEG_2nddown.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC<0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_MF_2nddown.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_MF_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_CC_2nddown.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_CC_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_BP_2nddown.csv"))
  
  png(paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_GO_BP_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2nddown.csv"))
}




cona<-c("0H","2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO")
conb<-c("2H","4H","8H","12H","16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  Table.DGE <- mutate(Table.DGE,mlogP=-1*log10(PValue))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01,]
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  png(paste0("DEG_next/",cona[c],"_",conb[c],"_MA.png"), width = 700, height = 500)
  plotSmear(dfresult.DGE, de.tags = is.degs, cex=0.3)
  dev.off() 
  
  Table.DGE$gene_name <- rep("",nrow(Table.DGE))
  Table.DGE$ideg <- rep("non_deg",nrow(Table.DGE))
  
  for(j in 1:nrow(Table.DGE)){
    if(rownames(Table.DGE)[j] %in% is.degs[1:50,1]){
      Table.DGE$gene_name[j]<-rownames(Table.DGE)[j]
    }
  }
  for(j in 1:nrow(Table.DGE)){
    if(rownames(Table.DGE)[j] %in% is.degs[,1]){
      Table.DGE$ideg[j]<-"deg"
    }
  }
  g <- ggplot(Table.DGE, aes(x = logFC, y = mlogP,label=gene_name,colour=ideg))
  g <- g + scale_color_manual(values = c("red", "black"))
  g <- g + geom_point()
  g <- g + geom_text_repel(na.rm = TRUE, size = 3.0, nudge_x = 1, nudge_y = -1, max.iter = 200,segment.alpha = 0.5)
  plot(g)
  ggsave(file = paste0("DEG_next/",cona[c],"_",conb[c],"_volcano.png"), plot = g)
  
}

#######################
#DEG,GO,KEGG compare next inhibitor

cona<-c("0H","2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4")
conb<-c("2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC>0,]
  #write.csv(degtable,paste0("DEG_next_inh/",cona[c],"_",conb[c],"_DEG_2ndup.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC>0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_MF_2ndup.csv"))

  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_MF_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()

  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_CC_2ndup.csv"))

  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_CC_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()

  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_BP_2ndup.csv"))

  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_BP_2ndup.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  reskk <- summary(kk2)
  write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2ndup.csv"))
}


#DEG,GO,KEGG compare next control

cona<-c("0H","2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4")
conb<-c("2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01&Table.DGE$logFC<0,]
  write.csv(degtable,paste0("DEG_next_inh/",cona[c],"_",conb[c],"_DEG_2nddown.csv"))
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01&Table.DGE$logFC<0]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_MF_2nddown.csv"))
  
  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_MF_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_CC_2nddown.csv"))
  
  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_CC_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  ego_result<-enrichGO(gene = is.degs$entrezgene_id, OrgDb=org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = FALSE)
  res <- summary(ego_result)
  write.csv(res,paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_BP_2nddown.csv"))
  
  png(paste0("GO_KEGG_next_inh/",cona[c],"_",conb[c],"_GO_BP_2nddown.png"), width = 900, height = 700)  # 描画デバイスを開く
  print(clusterProfiler::dotplot(ego_result))
  dev.off()
  
  #kk2<-enrichKEGG(is.degs$entrezgene_id, organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05, pAdjustMethod = "BH", all.genes$entrezgene_id, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.2, use_internal_data = FALSE)
  #reskk <- summary(kk2)
  #write.csv(reskk,paste0("GO_KEGG_next/",cona[c],"_",conb[c],"_KEGG_2nddown.csv"))
}




cona<-c("0H","2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4")
conb<-c("2H","4H","8H","12H","16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
alldegtable<-c()


for (c in 1:length(cona)) {
  countp<-dplyr::select(.data=TPMcount,starts_with(cona[c]),starts_with(conb[c]))
  groupp <- factor(c(cona[c], cona[c], cona[c],conb[c], conb[c], conb[c]))
  designp <- model.matrix(~ groupp)
  
  df.DGE <- DGEList(counts=countp,group=groupp)
  df.DGEo <- calcNormFactors(df.DGE)
  df.DGEt <- estimateCommonDisp(df.DGEo)
  df.DGEs <- estimateTagwiseDisp(df.DGEt)
  dfresult.DGE <- exactTest(df.DGEs)
  Table.DGE <- as.data.frame(topTags(dfresult.DGE, n = nrow(countp)))
  Table.DGE <- mutate(Table.DGE,mlogP=-1*log10(PValue))
  
  degtable<-Table.DGE[Table.DGE$FDR < 0.01,]
  alldegtable<-append(alldegtable,rownames(degtable))
  
  all.genes <- rownames(Table.DGE)
  is.degs <- all.genes[Table.DGE$FDR < 0.01]
  
  all.genes<-as.data.frame(all.genes)
  is.degs<-as.data.frame(is.degs)
  colnames(all.genes)<-"hgnc_symbol"
  colnames(is.degs)<-"hgnc_symbol"
  
  all.genes<-dplyr::inner_join(all.genes , gene.annotations ,by="hgnc_symbol")
  is.degs<-dplyr::inner_join(is.degs , gene.annotations ,by="hgnc_symbol")
  
  all.genes<-na.omit(all.genes)
  is.degs<-na.omit(is.degs)
  
  png(paste0("DEG_next_inh/",cona[c],"_",conb[c],"_MA.png"), width = 700, height = 500)
  plotSmear(dfresult.DGE, de.tags = is.degs, cex=0.3)
  dev.off() 
  
  Table.DGE$gene_name <- rep("",nrow(Table.DGE))
  Table.DGE$ideg <- rep("non_deg",nrow(Table.DGE))
  
  for(j in 1:nrow(Table.DGE)){
    if(rownames(Table.DGE)[j] %in% is.degs[1:50,1]){
      Table.DGE$gene_name[j]<-rownames(Table.DGE)[j]
    }
  }
  for(j in 1:nrow(Table.DGE)){
    if(rownames(Table.DGE)[j] %in% is.degs[,1]){
      Table.DGE$ideg[j]<-"deg"
    }
  }
  g <- ggplot(Table.DGE, aes(x = logFC, y = mlogP,label=gene_name,colour=ideg))
  g <- g + scale_color_manual(values = c("red", "black"))
  g <- g + geom_point()
  g <- g + geom_text_repel(na.rm = TRUE, size = 3.0, nudge_x = 1, nudge_y = -1, max.iter = 200,segment.alpha = 0.5)
  plot(g)
  ggsave(file = paste0("DEG_next_inh/",cona[c],"_",conb[c],"_volcano.png"), plot = g)
  
}



#######################
#FC calculation

tmp<-dplyr::select(.data=TPMcount,starts_with("0H"))
tmpm<-apply(tmp,1,mean)
tmpm<-as.data.frame(tmpm)
colnames(tmpm)<-"0H"
merge<-tmpm

cone<-c("2H","4H","8H","12H","16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")
for (c in cone) {
  tmp<-dplyr::select(.data=TPMcount,starts_with(c))
  tmpm<-apply(tmp,1,mean)
  tmpm<-as.data.frame(tmpm)
  colnames(tmpm)<-c
  merge<-cbind(merge,tmpm)
}
FC<-merge/merge[,1]

write.csv(merge,"merge_count.csv")
write.csv(FC,"fold_change.csv")

#######################
#FC calculation 12h

tmp<-dplyr::select(.data=TPMcount,starts_with("12H"))
tmpm<-apply(tmp,1,mean)
tmpm<-as.data.frame(tmpm)
colnames(tmpm)<-"12H"
merge<-tmpm

cone<-c("16H_CDK4","16H_DMSO","20H_CDK4","20H_DMSO","24H_CDK4","24H_DMSO","28H_CDK4","28H_DMSO")
for (c in cone) {
  tmp<-dplyr::select(.data=TPMcount,starts_with(c))
  tmpm<-apply(tmp,1,mean)
  tmpm<-as.data.frame(tmpm)
  colnames(tmpm)<-c
  merge<-cbind(merge,tmpm)
}
FC<-merge/merge[,1]

write.csv(merge,"merge_count_12.csv")
write.csv(FC,"fold_change_12.csv")

#######################
#FC calculation CDK/DMSO

cona<-c("16H_CDK4","20H_CDK4","24H_CDK4","28H_CDK4")
conb<-c("16H_DMSO","20H_DMSO","24H_DMSO","28H_DMSO")

tmp<-dplyr::select(.data=TPMcount,starts_with("12H"))
tmpm<-apply(tmp,1,mean)
tmpm<-as.data.frame(tmpm)
colnames(tmpm)<-"12H"
merge<-tmpm

for (i in 1:4) {
  tmp<-dplyr::select(.data=TPMcount,starts_with(conb[i]))
  tmpm<-apply(tmp,1,mean)
  tmpm<-as.data.frame(tmpm)
  colnames(tmpm)<-conb[i]
  
  tmp2<-dplyr::select(.data=TPMcount,starts_with(cona[i]))
  tmpm2<-apply(tmp2,1,mean)
  tmpm2<-as.data.frame(tmpm2)
  colnames(tmpm2)<-cona[i]
  
  tmp3<-tmpm2/tmpm
  
  if(i==1){
    FC<-tmp3
  }else{
    FC<-cbind(FC,tmp3)    
  }
}

colnames(FC)<-c("16H_CDK_INH","20H_CDK_INH","24H_CDK_INH","28H_CDK_INH")

#write.csv(merge,"merge_count_12.csv")
write.csv(FC,"fold_change_CDK_INH.csv")


