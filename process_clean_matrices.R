# Script to look at correlation between samples, and correlation between compounds
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(grDevices)
library(dplyr)
library(tidyr)
pal_cor<-colorRampPalette(brewer.pal(11, "RdYlGn")) #have one homogeneous color
nice_theme<-theme(axis.text=element_text(color="black",size=10),axis.ticks=element_line(color="black"))
log_base<-10 ## Pick one of the two
# log_base<-2
###########################################################################
## Load data
# All_media_pa_s<-read.csv("Selected_data/All_media_pa_selected.csv")
All_media_ph_s<-read.csv("Selected_data/All_media_ph_selected.csv")
# All_pellet_pa_cnorm_s<-read.csv("Selected_data/All_pellet_cnorm_pa_selected.csv")
All_pellet_ph_cnorm_s<-read.csv("Selected_data/All_pellet_cnorm_ph_selected.csv")

###########################################################################
# Select data
#
# Data<-All_media_pa_s
# Name<-"Media_pa_s"
#
Data<-All_media_ph_s
Name<-"Media_ph_s"
#
# Data<-All_pellet_pa_cnorm_s
# Name<-"Pellet_cnorm_pa_s"
#
Data<-All_pellet_ph_cnorm_s
Name<-"Pellet_cnorm_ph_s"

###########################################################################
## First, do a quick correlation between samples and a heatmap to verify that the data are clean
Cor_Data_clean<-cor(t(log(Data[,-1],base=log_base)),method="pearson",use="pairwise.complete.obs")
pheatmap(Cor_Data_clean,cluster_cols=F,cluster_rows=F,breaks=c(seq(min(Cor_Data_clean,na.rm=T),1,length.out=99)),color=pal_cor(99))
## Now, ready to calculate differential abundance between treatments (CG, CP, CPG) and control (C)
Out_dir<-"Diff_detections_2/"
# Chose which condition will be compared against which other
comp_tab<-c("CP","CG","CPG")
comp_1<-"C"
## This loop will calculate Fold Change (FC) and False-Discovery Rate (FDR) and print these in a csv file
for (k in 1:length(comp_tab)){
  comp_2<-comp_tab[k]
  # For each compound
  Diff_compound<-data.frame(Data=NA, Compound=NA, Timepoint=NA, Comparison=NA, logFC=NA, PValue=NA, FDR=NA)[numeric(0), ]
  n_comp<-0
  for (i in 2:ncol(Data)){
    for (t in 1:6){
      log_mat_1<-log(as.vector(na.omit(Data[which(Data[,1]==paste(comp_1,t,sep="")),i])),base=log_base)
      log_mat_2<-log(as.vector(na.omit(Data[which(Data[,1]==paste(comp_2,t,sep="")),i])),base=log_base)
      if (length(log_mat_1)>1 & length(log_mat_2)>1){ # We ask for at least two non-NA values in each vector
        fc<-mean(log_mat_2)-mean(log_mat_1)
        if (fc!=0){
          t_test<-t.test(log_mat_1,log_mat_2)
          n_comp<-n_comp+1
          Diff_compound[n_comp,1]<-Name
          Diff_compound[n_comp,2]<-colnames(Data)[i]
          Diff_compound[n_comp,3]<-t
          Diff_compound[n_comp,4]<-paste(comp_1,"_vs_",comp_2,sep="")
          Diff_compound[n_comp,5]<-fc
          Diff_compound[n_comp,6]<-t_test$p.value
          Diff_compound[n_comp,7]<-NA
        }
      }
    }
    # Now add the comparison without the time points
    log_mat_1<-log(as.vector(na.omit(Data[grepl(paste(comp_1,"\\d",sep=""),Data[,1],perl=T),i])),base=log_base)
    # if (Name!="Pellet_20161010_norm"){log_mat_1<-log_mat_1[log_mat_1>0]} # Remove all the zeroes
    # print(log_mat_1)
    log_mat_2<-log(as.vector(na.omit(Data[grepl(paste(comp_2,"\\d",sep=""),Data[,1],perl=T),i])),base=log_base)
    # if (Name!="Pellet_20161010_norm"){log_mat_2<-log_mat_2[log_mat_2>0]} # Remove all the zeroes
    # print(log_mat_2)
    if (length(log_mat_1)>1 & length(log_mat_2)>1){
      fc<-mean(log_mat_2)-mean(log_mat_1)
      if (fc!=0){
        t_test<-t.test(log_mat_1,log_mat_2)
        n_comp<-n_comp+1
        Diff_compound[n_comp,1]<-Name
        Diff_compound[n_comp,2]<-colnames(Data)[i]
        Diff_compound[n_comp,3]<-"All"
        Diff_compound[n_comp,4]<-paste(comp_1,"_vs_",comp_2,sep="")
        Diff_compound[n_comp,5]<-fc
        Diff_compound[n_comp,6]<-t_test$p.value
        Diff_compound[n_comp,7]<-NA
      }
    }
  }
  Diff_compound$FDR<-p.adjust(Diff_compound$PValue,method="BH",n=n_comp)
  # Look for Compounds with an abs(logFC)>=1 and FRD<=0.05
  if (log_base==10){df_sig<-Diff_compound[abs(Diff_compound$logFC)>=0.1 & Diff_compound$FDR<=0.05,]}
  else{df_sig<-Diff_compound[abs(Diff_compound$logFC)>=1 & Diff_compound$FDR<=0.05,]}
  write.csv(Diff_compound,file=paste(Out_dir,Name,"_",comp_2,"-vs-",comp_1,"_All_log",log_base,".csv",sep=""),row.names=FALSE)
  write.csv(df_sig,file=paste(Out_dir,Name,"_",comp_2,"-vs-",comp_1,"_Sig-only_log",log_base,".csv",sep=""),row.names=FALSE)
}
