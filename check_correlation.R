# Script to look at correlation between samples, and correlation between compounds
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(dplyr)
library(tidyr)
pal_cor<-colorRampPalette(brewer.pal(11, "RdYlGn")) #have one homogeneous color
nice_theme<-theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
###########################################################################
## Load data
All_media_pa<-read.csv("Clean_data/All_media_pa.csv")
All_media_ph<-read.csv("Clean_data/All_media_ph.csv")
All_pellet_pa<-read.csv("Clean_data/All_pellet_pa.csv")
All_pellet_ph<-read.csv("Clean_data/All_pellet_ph.csv")
All_pellet_pa_cnorm<-read.csv("Clean_data/All_pellet_pa_cnorm.csv")
All_pellet_ph_cnorm<-read.csv("Clean_data/All_pellet_ph_cnorm.csv")

###########################################################################
# Select data
#
# Data<-All_media_pa
# Name<-"Media_pa"
#
# Data<-All_media_ph
# Name<-"Media_ph"
#
# Data<-All_pellet_pa
# Name<-"All_pellet_pa"
#
Data<-All_pellet_ph
Name<-"All_pellet_ph"
#
# Data<-All_pellet_cnorm_pa
# Name<-"All_pellet_cnorm_pa"
#
# Data<-All_pellet_cnorm_ph
# Name<-"All_pellet_cnorm_ph"


###########################################################################
## Calculate correlation between samples
Cor_Data<-cor(t(log(Data[,-1],base=10)),method="pearson",use="pairwise.complete.obs")
## Plot these correlations as a heatmap
pheatmap(Cor_Data,cluster_cols=F,cluster_rows=F,breaks=c(seq(min(Cor_Data,na.rm=T),1,length.out=99)),color=pal_cor(99))
heatmap_name<-paste("Correlations/Full_heatmap_",Name,".pdf",sep="")
pheatmap(Cor_Data,cluster_cols=F,cluster_rows=F,breaks=c(seq(min(Cor_Data,na.rm=T),1,length.out=99)),color=pal_cor(99),width=12,height=9,filename=heatmap_name)
dev.off()
## Generate a boxplot of the correlations per sample to spot outliers
Df_cor<-data.frame(Sample_in=NA, Sample_out=NA, Corr=NA)[numeric(0), ]
n_line<-0
for (i in 1:nrow(Cor_Data)){
  for (j in 1:ncol(Cor_Data)){
    n_line<-n_line+1
    Df_cor[n_line,1]<-row.names(Cor_Data)[i]
    Df_cor[n_line,2]<-colnames(Cor_Data)[j]
    Df_cor[n_line,3]<-Cor_Data[i,j]
  }
}
ggplot(data=Df_cor) + geom_boxplot(aes(x=Sample_in,y=Corr)) + nice_theme + theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))
boxplot_name<-paste("Correlations/Full_boxplot_",Name,".pdf",sep="")
ggsave(file=boxplot_name,width=12,height=9)
