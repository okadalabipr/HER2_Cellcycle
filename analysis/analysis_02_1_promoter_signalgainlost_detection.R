library(dplyr)
library(Matrix)
library(ggplot2)
library(pals)

pos<-c(2,3,4,16,20:21)
cnam<-c("Chr","Start","End","Gene","postsignal","targetsignal")

TE16HCTLINH_signal<-read.csv("H3K27Ac_homerhistone_16HCTLINH_promoter_nme1ts250_signal.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
TE16HCTLINH_signal<-TE16HCTLINH_signal[,pos]
colnames(TE16HCTLINH_signal)<-cnam

TE20HCTLINH_signal<-read.csv("H3K27Ac_homerhistone_20HCTLINH_promoter_nme1ts250_signal.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
TE20HCTLINH_signal<-TE20HCTLINH_signal[,pos]
colnames(TE20HCTLINH_signal)<-cnam

TE24HCTLINH_signal<-read.csv("H3K27Ac_homerhistone_24HCTLINH_promoter_nme1ts250_signal.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
TE24HCTLINH_signal<-TE24HCTLINH_signal[,pos]
colnames(TE24HCTLINH_signal)<-cnam

TE28HCTLINH_signal<-read.csv("H3K27Ac_homerhistone_28HCTLINH_promoter_nme1ts250_signal.txt", header = TRUE, sep = "\t", stringsAsFactors = F)
TE28HCTLINH_signal<-TE28HCTLINH_signal[,pos]
colnames(TE28HCTLINH_signal)<-cnam

#############################################################

TE16HCTLINH_signal<-dplyr::mutate(TE16HCTLINH_signal,FC=log2(targetsignal/postsignal))
TE_sg = dplyr::filter(TE16HCTLINH_signal, FC > 0.3) 
TE_su = dplyr::filter(TE16HCTLINH_signal, abs(FC) <= 0.3)
TE_sl = dplyr::filter(TE16HCTLINH_signal, FC < -0.3)

TE_sg$group = rep("H3K27hh_promoter_nme1ts250_sg", nrow(TE_sg))
TE_su$group = rep("H3K27hh_promoter_nme1ts250_su", nrow(TE_su))
TE_sl$group = rep("H3K27hh_promoter_nme1ts250_sl", nrow(TE_sl))

TE16HCTLINH_signal = dplyr::bind_rows(TE_sg,TE_su)
TE16HCTLINH_signal = dplyr::bind_rows(TE16HCTLINH_signal,TE_sl)

TE16HCTLINH_signal = dplyr::arrange(TE16HCTLINH_signal, desc(FC))

TE16HCTLINH_signal$N = seq(1,nrow(TE16HCTLINH_signal),1)
write.csv(TE16HCTLINH_signal,"H3K27hh16HCTLINH_promoter_nme1ts250_signal_3classified.csv", row.names = F)
write.table(TE_sg[,1:4],"H3K27hh16HCTLINH_promoter_nme1ts250_sg.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_su[,1:4],"H3K27hh16HCTLINH_promoter_nme1ts250_su.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_sl[,1:4],"H3K27hh16HCTLINH_promoter_nme1ts250_sl.bed",quote=F,sep="\t",row.names=F,col.names=F)


TE_g <- ggplot(TE16HCTLINH_signal,aes(x = N , y = FC, fill = -N))+
  geom_bar(stat = "identity")+labs(title="H3K27hh16H CTL to INH")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,2000,4000,6000,8000))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))+
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

png("H3K27hh16HCTLINH_promoter_nme1ts250_signal.png", width = 700, height = 700)
plot(TE_g)
dev.off()

######################################

TE20HCTLINH_signal<-dplyr::mutate(TE20HCTLINH_signal,FC=log2(targetsignal/postsignal))
TE_sg = dplyr::filter(TE20HCTLINH_signal, FC > 0.3) 
TE_su = dplyr::filter(TE20HCTLINH_signal, abs(FC) <= 0.3)
TE_sl = dplyr::filter(TE20HCTLINH_signal, FC < -0.3)

TE_sg$group = rep("H3K27hh_promoter_nme1ts250_sg", nrow(TE_sg))
TE_su$group = rep("H3K27hh_promoter_nme1ts250_su", nrow(TE_su))
TE_sl$group = rep("H3K27hh_promoter_nme1ts250_sl", nrow(TE_sl))

TE20HCTLINH_signal = dplyr::bind_rows(TE_sg,TE_su)
TE20HCTLINH_signal = dplyr::bind_rows(TE20HCTLINH_signal,TE_sl)

TE20HCTLINH_signal = dplyr::arrange(TE20HCTLINH_signal, desc(FC))

TE20HCTLINH_signal$N = seq(1,nrow(TE20HCTLINH_signal),1)
write.csv(TE20HCTLINH_signal,"H3K27hh20HCTLINH_promoter_nme1ts250_signal_3classified.csv", row.names = F)
write.table(TE_sg[,1:4],"H3K27hh20HCTLINH_promoter_nme1ts250_sg.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_su[,1:4],"H3K27hh20HCTLINH_promoter_nme1ts250_su.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_sl[,1:4],"H3K27hh20HCTLINH_promoter_nme1ts250_sl.bed",quote=F,sep="\t",row.names=F,col.names=F)

TE_g <- ggplot(TE20HCTLINH_signal,aes(x = N , y = FC, fill = -N))+
  geom_bar(stat = "identity")+labs(title="H3K27hh20H CTL to INH")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,2000,4000,6000,8000))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))+
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

png("H3K27hh20HCTLINH_promoter_nme1ts250_signal.png", width = 700, height = 700)
plot(TE_g)
dev.off()

######################################

TE24HCTLINH_signal<-dplyr::mutate(TE24HCTLINH_signal,FC=log2(targetsignal/postsignal))
TE_sg = dplyr::filter(TE24HCTLINH_signal, FC > 0.3) 
TE_su = dplyr::filter(TE24HCTLINH_signal, abs(FC) <= 0.3)
TE_sl = dplyr::filter(TE24HCTLINH_signal, FC < -0.3)

TE_sg$group = rep("H3K27hh_promoter_nme1ts250_sg", nrow(TE_sg))
TE_su$group = rep("H3K27hh_promoter_nme1ts250_su", nrow(TE_su))
TE_sl$group = rep("H3K27hh_promoter_nme1ts250_sl", nrow(TE_sl))

TE24HCTLINH_signal = dplyr::bind_rows(TE_sg,TE_su)
TE24HCTLINH_signal = dplyr::bind_rows(TE24HCTLINH_signal,TE_sl)

TE24HCTLINH_signal = dplyr::arrange(TE24HCTLINH_signal, desc(FC))

TE24HCTLINH_signal$N = seq(1,nrow(TE24HCTLINH_signal),1)
write.csv(TE24HCTLINH_signal,"H3K27hh24HCTLINH_promoter_nme1ts250_signal_3classified.csv", row.names = F)
write.table(TE_sg[,1:4],"H3K27hh24HCTLINH_promoter_nme1ts250_sg.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_su[,1:4],"H3K27hh24HCTLINH_promoter_nme1ts250_su.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_sl[,1:4],"H3K27hh24HCTLINH_promoter_nme1ts250_sl.bed",quote=F,sep="\t",row.names=F,col.names=F)

TE_g <- ggplot(TE24HCTLINH_signal,aes(x = N , y = FC, fill = -N))+
  geom_bar(stat = "identity")+labs(title="H3K27hh24H CTL to INH")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,2000,4000,6000,8000))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))+
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

png("H3K27hh24HCTLINH_promoter_nme1ts250_signal.png", width = 700, height = 700)
plot(TE_g)
dev.off()

######################################

TE28HCTLINH_signal<-dplyr::mutate(TE28HCTLINH_signal,FC=log2(targetsignal/postsignal))
TE_sg = dplyr::filter(TE28HCTLINH_signal, FC > 0.3) 
TE_su = dplyr::filter(TE28HCTLINH_signal, abs(FC) <= 0.3)
TE_sl = dplyr::filter(TE28HCTLINH_signal, FC < -0.3)

TE_sg$group = rep("H3K27hh_promoter_nme1ts250_sg", nrow(TE_sg))
TE_su$group = rep("H3K27hh_promoter_nme1ts250_su", nrow(TE_su))
TE_sl$group = rep("H3K27hh_promoter_nme1ts250_sl", nrow(TE_sl))

TE28HCTLINH_signal = dplyr::bind_rows(TE_sg,TE_su)
TE28HCTLINH_signal = dplyr::bind_rows(TE28HCTLINH_signal,TE_sl)

TE28HCTLINH_signal = dplyr::arrange(TE28HCTLINH_signal, desc(FC))

TE28HCTLINH_signal$N = seq(1,nrow(TE28HCTLINH_signal),1)
write.csv(TE28HCTLINH_signal,"H3K27hh28HCTLINH_promoter_nme1ts250_signal_3classified.csv", row.names = F)
write.table(TE_sg[,1:4],"H3K27hh28HCTLINH_promoter_nme1ts250_sg.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_su[,1:4],"H3K27hh28HCTLINH_promoter_nme1ts250_su.bed",quote=F,sep="\t",row.names=F,col.names=F)
write.table(TE_sl[,1:4],"H3K27hh28HCTLINH_promoter_nme1ts250_sl.bed",quote=F,sep="\t",row.names=F,col.names=F)

TE_g <- ggplot(TE28HCTLINH_signal,aes(x = N , y = FC, fill = -N))+
  geom_bar(stat = "identity")+ labs(title="H3K27hh28H CTL to INH")+
  theme_bw()+
  guides(fill = F)+
  coord_flip()+
  scale_fill_gradientn(colours=coolwarm(10))+
  scale_x_reverse(breaks=c(1,2000,4000,6000,8000))+
  xlab("Enhancers ranked by H3K27Ac Signal \n Fold Change (log2)")+
  ylab("H3K27Ac Signal Fold Change (log2)")+
  theme(axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))+
  theme(axis.title.x = element_text(size=15),axis.title.y = element_text(size=15))+
  theme(panel.grid = element_blank(), strip.background = element_blank())+
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size = 20))+
  scale_y_continuous(limits=c(-6,6), breaks=seq(-6,6,2))+
  geom_hline(yintercept = c(-0.2,0.2), linetype="dashed", size = 0.8)+
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 2.0))

png("H3K27hh28HCTLINH_promoter_nme1ts250_signal.png", width = 700, height = 700)
plot(TE_g)
dev.off()

######################################