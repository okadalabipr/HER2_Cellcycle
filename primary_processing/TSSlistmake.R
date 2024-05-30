library(dplyr)
library(Seurat)
library(Matrix)
library(stringr)
library(biomaRt)

listEnsemblArchives(https = TRUE)
listMarts(host="https://jan2020.archive.ensembl.org")
db <- useMart("ensembl",host="https://jan2020.archive.ensembl.org")
listDatasets(db)

hg <- useDataset("hsapiens_gene_ensembl", mart = db)
res <- getBM(attributes = c("ensembl_gene_id", 'chromosome_name',"start_position", "end_position","external_gene_name","strand"),mart = hg)

chr<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y")
res_s<-res[res[,2] %in% chr,]

res_s_f<-res_s[res_s[,6]==1,]
res_s_r<-res_s[res_s[,6]==(-1),]

#######TSS2500

tss_f<-res_s_f[2:5]
tss_f[,2]<-res_s_f[,3]-250
tss_f[,3]<-res_s_f[,3]+250

tss_r<-res_s_r[2:5]
tss_r[,2]<-res_s_r[,4]-250
tss_r[,3]<-res_s_r[,4]+250

tss<-rbind(tss_f,tss_r)
tss[,1]<-paste0("chr",tss[,1])

#tss=apply(tss,c(1,2),as.character)

options(scipen = 999)
write.table(tss,"tss_250.bed",sep="\t",quote=F,col.names=F,row.names=F)
