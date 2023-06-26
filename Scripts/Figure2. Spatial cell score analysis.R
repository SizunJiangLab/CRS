#####Spatial cell score assessment#####
library(GSVA)
library(dplyr)
library(patchwork)
library(tidyr)
library(ggpubr)
#Read in normalized GeoMX digital spatial profiling data
q3norm = read.csv("mydirectory/Q3 Normalization.csv")
rownames(q3norm) = q3norm$TargetName
q3norm = subset(q3norm, select = -c(TargetName))
GSVA.df.all=as.matrix(q3norm)
q3norm_cd45 = q3norm %>% select(matches("CD45"))
q3norm_PanCK = q3norm %>% select(matches("PanCK"))
#Read in signature genesets
genelist = read.csv("mydirectory/Supple.Table2 Genesets used in cell scoring analysis.csv")
EosinophilGene=genelist$Eosinophil
M2_PolarizationGene=genelist$M2_Polarization
Scoring_list=data.frame(Eosinophil=EosinophilGene,
                        M2_Polarization=M2_PolarizationGene)
#Take eosinophil signature score as an example for scoring
EosinophilScore=gsva(GSVA.df.all,Scoring_list,method="ssgsea",ssgsea.norm=TRUE,verbose=T)
EosinophilScore=t(EosinophilScore)
EosinophilScore=as.data.frame(EosinophilScore)
Eosinophil_score=data.frame(sample=rownames(EosinophilScore),value=EosinophilScore$Eosinophil)

#Match with grouping and region information
tmp_score=EosinophilScore
tmp_score$group=''
tmp_score$group[tmp_score$sample %in% colnames(q3norm)[1:12]|tmp_score$sample %in% colnames(q3norm)[21:32]]='CRSwNP'
tmp_score$group[tmp_score$group=='']='Ctrl'
tmp_score$region=''
tmp_score[grep("CD45",tmp_score$sample,fixed = T),'region']='CD45'
tmp_score[grep("PanCK",tmp_score$sample,fixed = T),'region']='PanCK'
tmp_score_region=subset(tmp_score,region=="CD45")

#Plot1: ggplot2
ggplot(tmp_score_region,aes(x=group,y=value,color=group))+
  scale_color_manual(values = subtype.colors)+
  stat_boxplot(geom = "errorbar",aes(color=group),lwd=0.8)+
  stat_compare_means(method = "wilcox.test",label = "..p.format..",label.x.npc = "center",hjust = 0.5,vjust=0.3)+ 
  geom_boxplot(outlier.fill = "white",lwd=0.8,outlier.color = "white",aes(color=group))+
  geom_jitter(aes(x=group,y=value,color=group),position=position_jitter(width = 0.1, height=0),alpha=0.5,size=1)+
  theme_classic()+
  theme(axis.text=element_text(face = "bold"),
        axis.text.x = element_text(vjust=0.9,hjust = 0.9,angle=45,size = 12,colour = "black"),
        axis.text.y = element_text(size = 12,colour = "black"),
        axis.title = element_text(size = 12,face = "bold",colour = "black")) +
  labs(x='',y="Eosinophil Spatial Score")

#Plot2: vilin plot
library("vioplot")
S1=tmp_score_region$value[tmp_score_region$group=='CRSwNP']
S2=tmp_score_region$value[tmp_score_region$group=='Ctrl']
vioplot(S1,S2,
        names=c('CRSwNP','Ctrl'),
        col=rev(subtype.colors),
        ylab = "Eosinophil Spatial Score",main = "CD45 region")+
  mtext(paste0("p=",wilcox.result$p.value,sep=''),side = 3)



