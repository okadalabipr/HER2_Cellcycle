
library(dplyr)
library(Matrix)
library(edgeR)
library(tidyverse)
library(here)
library(OmnipathR)
library(dorothea)
library(decoupleR)
library(workflowr)
library(rmarkdown)
library(org.Hs.eg.db)
library(progeny)
library(dorothea)
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(readr)
library(gt)



data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>%
  dplyr::filter(confidence %in% c("A", "B"))

RNAcount<-read.table("../mergecount_TPMedgeR_min1.csv",sep=",")
colnames(RNAcount)<-RNAcount[1,]
rownames(RNAcount)<-RNAcount[,1]
RNAcount<-RNAcount[-1,]
RNAcount<-RNAcount[,-1]
RNAcount<-apply(RNAcount,c(1,2),as.numeric)

##########


RNAcount_b<-as.data.frame(cbind(rownames(RNAcount),RNAcount))
RNAcount_b[,2:14]<-apply(RNAcount_b[,2:14],c(1,2),as.numeric)
colnames(RNAcount_b)<-c("Gene.Name","RNA_0H","RNA_2H","RNA_4H","RNA_8H","RNA_12H","RNA_16H_INH","RNA_16H_CTL","RNA_20H_INH","RNA_20H_CTL","RNA_24H_INH","RNA_24H_CTL","RNA_28H_INH","RNA_28H_CTL")

            
######################

deglist_down<-read.csv("../DEG_INH/16H_DMSO_16H_CDK4_DEG_2nddown.csv")
deglist_up<-read.csv("../DEG_INH/16H_DMSO_16H_CDK4_DEG_2ndup.csv")
deglist<-rbind(deglist_down,deglist_up)

RNAcount_b_sub<-RNAcount_b[RNAcount_b$Gene.Name %in% deglist$X,]

lg2fc<-as.matrix(log2(RNAcount_b_sub$RNA_16H_INH/RNAcount_b_sub$RNA_16H_CTL))
colnames(lg2fc)<-"t"
rownames(lg2fc)<-RNAcount_b_sub$Gene.Name

lg2fc<-apply(lg2fc,c(1,2),as.numeric)
tf_activities_stat <- dorothea::run_viper(lg2fc, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

g<-ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")
g

png("16H_ctlinh_drothea_AB.png", width = 600, height = 400) 
print(g)
dev.off()  



######################


deglist_down<-read.csv("../DEG_INH/20H_DMSO_20H_CDK4_DEG_2nddown.csv")
deglist_up<-read.csv("../DEG_INH/20H_DMSO_20H_CDK4_DEG_2ndup.csv")
deglist<-rbind(deglist_down,deglist_up)



RNAcount_b_sub<-RNAcount_b[RNAcount_b$Gene.Name %in% deglist$X,]

lg2fc<-as.matrix(log2(RNAcount_b_sub$RNA_20H_INH/RNAcount_b_sub$RNA_20H_CTL))
colnames(lg2fc)<-"t"
rownames(lg2fc)<-RNAcount_b_sub$Gene.Name

lg2fc<-apply(lg2fc,c(1,2),as.numeric)
tf_activities_stat <- dorothea::run_viper(lg2fc, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE, p_value=TRUE))

as.list(args(dorothea::run_viper))

tidy = FALSE

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

g<-ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")
g

png("20H_ctlinh_drothea_AB.png", width = 600, height = 400) 
print(g)
dev.off()  

##############

deglist_down<-read.csv("../DEG_INH/24H_DMSO_24H_CDK4_DEG_2nddown.csv")
deglist_up<-read.csv("../DEG_INH/24H_DMSO_24H_CDK4_DEG_2ndup.csv")
deglist<-rbind(deglist_down,deglist_up)



RNAcount_b_sub<-RNAcount_b[RNAcount_b$Gene.Name %in% deglist$X,]

lg2fc<-as.matrix(log2(RNAcount_b_sub$RNA_16H_INH/RNAcount_b_sub$RNA_16H_CTL))
colnames(lg2fc)<-"t"
rownames(lg2fc)<-RNAcount_b_sub$Gene.Name

lg2fc<-apply(lg2fc,c(1,2),as.numeric)
tf_activities_stat <- dorothea::run_viper(lg2fc, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE))
tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

g<-ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")
g

png("24H_ctlinh_drothea_AB.png", width = 600, height = 400) 
print(g)
dev.off()  



######################


deglist_down<-read.csv("../DEG_INH/28H_DMSO_28H_CDK4_DEG_2nddown.csv")
deglist_up<-read.csv("../DEG_INH/28H_DMSO_28H_CDK4_DEG_2ndup.csv")
deglist<-rbind(deglist_down,deglist_up)

dim(deglist_down)
dim(deglist_up)

RNAcount_b_sub<-RNAcount_b[RNAcount_b$Gene.Name %in% deglist$X,]

lg2fc<-as.matrix(log2(RNAcount_b_sub$RNA_20H_INH/RNAcount_b_sub$RNA_20H_CTL))
colnames(lg2fc)<-"t"
rownames(lg2fc)<-RNAcount_b_sub$Gene.Name

lg2fc<-apply(lg2fc,c(1,2),as.numeric)
tf_activities_stat <- dorothea::run_viper(lg2fc, regulons,
                                          options =  list(minsize = 5, eset.filter = FALSE, 
                                                          cores = 1, verbose = FALSE, nes = TRUE, p_value=TRUE))

as.list(args(dorothea::run_viper))

tidy = FALSE

tf_activities_stat_top25 <- tf_activities_stat %>%
  as.data.frame() %>% 
  rownames_to_column(var = "GeneID") %>%
  dplyr::rename(NES = "t") %>%
  dplyr::top_n(25, wt = abs(NES)) %>%
  dplyr::arrange(NES) %>% 
  dplyr::mutate(GeneID = factor(GeneID))

g<-ggplot(tf_activities_stat_top25,aes(x = reorder(GeneID, NES), y = NES)) + 
  geom_bar(aes(fill = NES), stat = "identity") +
  scale_fill_gradient2(low = "darkblue", high = "indianred", 
                       mid = "whitesmoke", midpoint = 0) + 
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = 
          element_text(angle = 45, hjust = 1, size =10, face= "bold"),
        axis.text.y = element_text(size =10, face= "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Transcription Factors")
g

png("28H_ctlinh_drothea_AB.png", width = 600, height = 400) 
print(g)
dev.off()  

