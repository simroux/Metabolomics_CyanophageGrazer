# Script to look at correlation between samples, and correlation between compounds
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(dplyr)
library(tidyr)
pal_cor<-colorRampPalette(brewer.pal(11, "RdYlGn")) #have one homogeneous color
nice_theme<-theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
heatmap_col_neg<-colorRampPalette(rev(brewer.pal(9,"YlGnBu")))
heatmap_col_pos<-colorRampPalette(brewer.pal(9,"YlOrRd"))
###########################################################################
#
# Name<-"Media_ph_CP"
# Code<-"Media_ph_s_CP_log10"
#
Name<-"Media_ph_CP_standards"
Code<-"Media_ph_s_CP_log10_std"
#
# Name<-"Media_ph_CPG"
# Code<-"Media_ph_s_CPG_log10"
#
# Name<-"Media_ph_CPG_standards"
# Code<-"Media_ph_s_CPG_log10_std"
#
# Name<-"Pellet_cnorm_ph_CP"
# Code<-"Pellet_cnorm_ph_s_CP_log10"
#
# Name<-"Pellet_cnorm_ph_CP_standards"
# Code<-"Pellet_cnorm_ph_s_CP_log10_std"
#
# Name<-"Pellet_cnorm_ph_CPG"
# Code<-"Pellet_cnorm_ph_s_CPG_log10"
#
# Name<-"Pellet_cnorm_ph_CPG_standards"
# Code<-"Pellet_cnorm_ph_s_CPG_log10_std"

###########################################################################
Data_pheatmap<-read.csv(paste("Heatmap_input/",Code,"_for_pheatmap.csv",sep=""))
Data_pheatmap_cluster<-read.csv(paste("Heatmap_input/",Code,"_for_clustering.csv",sep=""))
## Calculate clustering on all values (otherwise too many NAs)
tree_cluster<-hclust(dist(Data_pheatmap_cluster))
## Check the min and max (for adjusting color scale)
min_mat <- -1
if (min(Data_pheatmap,na.rm=T) < -1){min_mat<-min(Data_pheatmap,na.rm=T)}
max_mat <- 1
if (max(Data_pheatmap,na.rm=T) > 1){max_mat<-max(Data_pheatmap,na.rm=T)}
## And print the heatmap
pheatmap(Data_pheatmap,breaks=c(seq(min_mat,0,length.out=99),seq(0.000001,max_mat,length.out=100)),color=c(heatmap_col_neg(100),heatmap_col_pos(100)), cluster_cols = F, cluster_rows = tree_cluster)
## Now printing this heatmap in pdf
out_file<-paste("Heatmaps/",Name,".pdf",sep="")
pheatmap(Data_pheatmap,breaks=c(seq(min_mat,0,length.out=99),seq(0.000001,max_mat,length.out=100)),color=c(heatmap_col_neg(100),heatmap_col_pos(100)), cluster_cols = F, cluster_rows = tree_cluster, file=out_file)
dev.off()
