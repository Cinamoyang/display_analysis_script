library(Rcpp)
library(readxl)
library(edgeR)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(RColorBrewer)
library(statmod)
library(ggsci)
library(dplyr)

#import the dataset pattern"MTTC.....B.....CSWD"
#the path should be changed to where the "allpeptides_only_MTTC_B_CSWD.txt" is in your PC
allpeptides <- read.table("C:/Users/user/OneDrive/Desktop/final_NGS_analysis/allpeptides_only_MTTC_B_CSWD.txt",check.names=FALSE, header = TRUE,sep = "\t")
allpeptides[is.na(allpeptides)] <- 0
allpeptides<-allpeptides[,c(2:31)]

#Convert to DGE list & filtering out the peptide sequence which show up less than 14 among 29 library
allpeptides_dge <-DGEList(counts = allpeptides[,2:30], genes = allpeptides[,1])
keep_allpeptides_dge <- rowSums((allpeptides_dge$counts)>=1)>= 14
allpeptides_dge_filtered<-allpeptides_dge[keep_allpeptides_dge,,keep.lib.sizes = FALSE]

#draw MDS(PCA),extract distance matrix&calculate eig &contribution
PCA_all_data<-plotMDS(allpeptides_dge_filtered,top = 200,gene.selection = "common",cex=0.5,plot=FALSE)
PCA_all_row_data<-cmdscale(as.dist(PCA_all_data$distance.matrix),eig = T)
all_eig_val<-PCA_all_row_data$eig
first_contrib<-all_eig_val[1]/sum(abs(all_eig_val))
second_contrib<-all_eig_val[2]/sum(abs(all_eig_val))

#extract all 7th round libraries & normrization of counts
all_7_dge<-allpeptides_dge_filtered[,c(10:25)]
all_7_dge<-calcNormFactors(all_7_dge,method="TMMwsp")

#save normrization of counts as count per million
a<-cpm(all_7_dge)
a<-as.data.frame(a)
a$name<-rownames(a)
write.table(a,"C:/Users/user/OneDrive/Desktop/final_NGS_analysis/cpm.txt",row.names=FALSE,col.names=TRUE,sep = "\t",fileEncoding="UTF-8")

#constract design matrix
Elution_method<-factor(design_group$Elution_method, levels=c("WBT","2PD1","HAC","Urea"))
immobilize_target<-factor(design_group$immobilize_target,levels=c("HAC","PDL1"))
design_matrix<-model.matrix(~Elution_method+immobilize_target+Elution_method:immobilize_target)

#calculate dispersion based on design_matrix then fit
y<-estimateGLMRobustDisp(all_7_dge,design_matrix) #common & trend dispersion
fit<-glmQLFit(y,design_matrix,robust=TRUE) #shrinked dispersion & fit for coeffiences

#(a) peptides have more elution against target of PD-L1 than HAC during 2-PD-1 elution
fit_a_HAC_vs_PDL1<- glmQLFTest(fit, contrast=c(0,0,0,0,1,1,0,0))
topTags_fit_a_HAC_vs_PDL1<-topTags(fit_a_HAC_vs_PDL1,n=2000,sort.by = "PValue")$table
keep_topTags_fit_a_HAC_vs_PDL1<-(topTags_fit_a_HAC_vs_PDL1[,6] < (0.1) & (topTags_fit_a_HAC_vs_PDL1[,2] > (1)))
topTags_fit_a_HAC_vs_PDL1_filter<-topTags_fit_a_HAC_vs_PDL1[keep_topTags_fit_a_HAC_vs_PDL1,]


#(b) peptides have more elution against target of PD-L1 than HAC during HAC elution
fit_b_HAC_vs_PDL1<-glmQLFTest(fit, contrast=c(0,0,0,0,1,0,1,0))
topTags_fit_b_HAC_vs_PDL1<-topTags(fit_b_HAC_vs_PDL1,n=2000,sort.by = "PValue")$table
keep_topTags_fit_b_HAC_vs_PDL1<-(topTags_fit_b_HAC_vs_PDL1[,6] < (0.1) & (topTags_fit_b_HAC_vs_PDL1[,2] > (1)))
topTags_fit_b_HAC_vs_PDL1_filter<-topTags_fit_b_HAC_vs_PDL1[keep_topTags_fit_b_HAC_vs_PDL1,]

#(c) peptides have averagely more elution with 2PD-1 and HAC than washing buffer when against target of PD-L1
fit_c_2PD1HAC_vs_WBT<-glmQLFTest(fit, contrast=c(0,0.5,0.5,0,0,0.5,0.5,0))
topTags_fit_c_2PD1HAC_vs_WBT<-topTags(fit_c_2PD1HAC_vs_WBT,n=2000,sort.by = "PValue")$table
keep_topTags_fit_c_2PD1HAC_vs_WBT<-(topTags_fit_c_2PD1HAC_vs_WBT[,6] < (0.1) & (topTags_fit_c_2PD1HAC_vs_WBT[,2] > (1)))
topTags_fit_c_2PD1HAC_vs_WBT_filter<-topTags_fit_c_2PD1HAC_vs_WBT[keep_topTags_fit_c_2PD1HAC_vs_WBT,]

#intersection of (a), (b) and (c) to give final candidate peptide
merge1<-inner_join(topTags_fit_a_HAC_vs_PDL1_filter,topTags_fit_b_HAC_vs_PDL1_filter, by = c("genes"="genes"),suffix=c(".a",".b"))
merge2<-inner_join(merge1,topTags_fit_c_2PD1HAC_vs_WBT_filter, by = c("genes"="genes"))
merge2[,"Final_FDR"]<-1-(1-merge2[,6])*(1-merge2[,11])*(1-merge2[,16])
colnames(merge2)[c(12:16)]<-c("logFC.c","logCPM.c","F.c","PValue.c","FDR.c")
colnames(merge2)[1]<-c("Peptides")
view(merge2)

