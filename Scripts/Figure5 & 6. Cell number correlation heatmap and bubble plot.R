library(reshape2)
library(Hmisc)
library(corrplot)
library(dplyr)
library(tibble)
library(Seurat)
library(ggplot2)
library(forcats)
library(ggcorrplot)
library(ggpubr)
library(RColorBrewer)
load("mydirectory/Allcell.RData")

#####1. Cell abundance correlation heatmap#####
All_meta=Allcell@meta.data
temp_all<-data.frame(table(All_meta$Sample_name,All_meta$Sub_cluster))
colnames(temp_all) <- c("Samples",'Cell_type',"Percentage of Cells")
temp2_all<-data.frame(table(All_meta$Sample_name))
colnames(temp2_all) <- c("Samples",'total')
temp_all=merge.data.frame(temp_all,temp2_all,all.x = T,by = 'Samples')
temp_all$`Percentage of Cells`=temp_all$`Percentage of Cells`/temp_all$total
temp_all_mat=dcast(data=temp_all,Cell_type ~Samples)
row.names(temp_all_mat)=temp_all_mat[,1]
temp_all_mat=temp_all_mat[,-1]
temp_all_mat=as.data.frame(t(temp_all_mat))
M=cor(temp_all_mat,method = "pearson")
p.mat<- cor_pmat(as.matrix(M[,1:8]),vars = NULL,method = "spearman",
                 alternative = "two.sided",conf.level = 0.95)
ggcorrplot(M, p.mat = p.mat,type = 'full',sig.level = 0.05,insig= "blank",show.diag = T,hc.order = T)+
  theme(axis.text.x = element_text(size = 8,face = "bold"),
        axis.text.y  = element_text(size = 8,face = "bold"))
ggsave(filename = 'mydirectory/Allcell_cor.pdf',width = 12,height = 12)

#####2. Cell abundance correlation with specific cell groups, bubble plot#####
load("mydirectory/Basalcell.RData")
load("mydirectory/immune.combined.RData")
## 1). Calculate the number of each cell group including the specific ones in each sample
Allcell$Sub_cluster[Allcell$Cellbarcode %in% Basal$Cellbarcode[Basal$branch=="Cell-Fate1"]]="Cell-fate1 Basal cells"
All_meta=Allcell@meta.data
temp_all<-data.frame(table(All_meta$Sample_name,All_meta$Sub_cluster))
colnames(temp_all) <- c("Samples",'Cell_type',"Percentage of Cells")
temp2_all<-data.frame(table(All_meta$Sample_name))
colnames(temp2_all) <- c("Samples",'total')
temp_all=merge.data.frame(temp_all,temp2_all,all.x = T,by = 'Samples')
temp_all$`Percentage of Cells`=temp_all$`Percentage of Cells`/temp_all$total
levels(temp_all$Cell_type)
needed=unique(temp_all$Cell_type)
needed=needed[!needed %in% c("Cell-fate1 Basal cells")]
temp_all_out=temp_all[temp_all$Cell_type %in% needed,]
temp_all_subcluster=temp_all[temp_all$Cell_type=="Cell-fate1 Basal cells",]

## 1). Calculate the correlation of specific cell groups with all other cell groups
temp_all_reorder=rbind(temp_all_subcluster,temp_all_out)
my_cell=as.vector(unique(temp_all_reorder$Cell_type))
my_cell
j=1
my_cell[j]
my_cell_tmp=my_cell[-j]
my_rev=c('PCC',"p.value",'cell')
for (i in 1:length(my_cell_tmp)) {
  my_target1=my_cell[j]
  my_target2=my_cell_tmp[i]
  temp_all1=temp_all_reorder[temp_all_reorder$Cell_type==my_target1,]
  temp_all2=temp_all_reorder[temp_all_reorder$Cell_type==my_target2,]
  temp_all1=temp_all1$`Percentage of Cells`
  temp_all2=temp_all2$`Percentage of Cells`
  my_cor=cor.test(temp_all1,temp_all2,method = "spearman")
  my_p=my_cor$p.value
  my_cor=my_cor$estimate
  my_cor=c(my_cor,my_p,my_target2)
  my_rev=rbind(my_rev,my_cor)
}
colnames(my_rev)=my_rev[1,]
my_rev=my_rev[-1,]
my_rev=as.data.frame(my_rev)
my_rev$PCC=as.numeric(my_rev$PCC)
my_rev$p.value=as.numeric(my_rev$p.value)
my_rev_cellfate1=my_rev
my_rev_cellfate1$group="Cell-fate1"
#Do the same for Cell-fate2

# Plotting
my_rev1=rbind(my_rev_cellfate2,my_rev_cellfate1)
my_rev1=subset(my_rev1,cell %in% unique(immune.combined$Sub_cluster))
my_rev1=my_rev1[order(my_rev1[,1],decreasing = T),]
order_name=my_rev1$cell
order_name=unique(order_name)
my_rev1$logP=-log(my_rev1$p.value)
ggplot(my_rev1, aes(x = group, y = cell, fill = PCC, size = logP)) +
  geom_point(shape = 21) +
  scale_size_continuous(range = c(1, 6)) +
  theme_classic()+
  scale_fill_gradient2(low = "blue",mid = "white",high = "red")+
  labs(x = "", y = "", fill = "Correlation", size = "-Log(p-value)",
       title = "Bubble Heatmap Plot")+ 
  scale_y_discrete(limits=rev(order_name))+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),
        axis.text.x = element_text(angle = 0, face="bold",size = 10),
        axis.text.y = element_text(size = 10, face="bold"))
ggsave(filename = paste0('mydirectory/Cell-fate correlation.pdf'),width = 6,height = 6)








