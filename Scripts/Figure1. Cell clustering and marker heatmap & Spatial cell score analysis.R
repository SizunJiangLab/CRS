#####1. Cell clustering#####
library(Seurat)
library(ggplot2)
library(SeuratData)
library(patchwork)
load("mydirectory/CRS.Immune.RData")

## 1). Integrate immune cells from each sample (listed in CRS.Immune.list) first
CRS.Immune.list <- lapply(X = CRS.Immune.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

Immune.features <- SelectIntegrationFeatures(object.list = CRS.Immune.list)
immune.anchors <- FindIntegrationAnchors(object.list = CRS.Immune.list, anchor.features = Immune.features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

## 2). Reduction and clustering on integrated dataset
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 50, verbose = FALSE)
ElbowPlot(immune.combined,ndims = 50)
#Select appropriate number of pca for next steps
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:40)
immune.combined <- FindClusters(immune.combined, resolution = 0.9)

## 3). Visualization
DimPlot(immune.combined, reduction = "umap",label=T)+scale_color_manual(values = my_color)


#####2. Marker heatmap#####
library(Seurat)
library(magrittr)
library(dplyr)
library(ggplot2)
library(pheatmap)
load("mydirectory/Allcell.RData")

## 1). Identify markers for each cell cluster

Allcell@active.ident=as.factor(Allcell$Sub_cluster)
Allcell.markers <- FindAllMarkers(Allcell,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10_Allcell <- Allcell.markers %>% group_by(cluster) %>% slice_max(n = 100, order_by = avg_log2FC)

## 2). Extract the normalized gene count matrix and calculate mean expression level for each cluster

Allcell_counts=as.matrix(Allcell@assays$RNA@data)
my_all_mat=t(Allcell_counts)
my_all_mat=as.data.frame(my_all_mat)
my_all_mat$Sub_cluster=Allcell$Sub_cluster
tmp=aggregate(my_all_mat[,-ncol(my_all_mat)],by=list(Sub_cluster=my_all_mat$Sub_cluster),FUN=mean)
row.names(tmp)=tmp$Sub_cluster
tmp=tmp[,-1]
my_order=unique(Allcell$Sub_cluster)
tmp=tmp[my_order,]
top10_Allcell$cluster=factor(top10_Allcell$cluster,levels = my_order)
top10_Allcell=top10_Allcell[order(top10_Allcell$cluster),]
top10_Allcell_gene=top10_Allcell$gene
top10_Allcell_gene=unique(top10_Allcell_gene)
tmp=tmp[,top10_Allcell_gene]
tmp=as.data.frame(t(tmp))

## 3). plotting
ht <- Heatmap(tmp, 
              col = colorRamp2(c(3, 0, -3),c("#D93223", "white", "#049DBF")),
              heatmap_legend_param = list(title= "Normalized expression",
                                          legend_height=unit(8,"cm"),
                                          legend_direction="vertical"),
              width = unit(13, "cm"), height = unit(13, "cm"),show_column_dend = F,
              cluster_rows = FALSE,column_dend_height = unit(1, "cm"),
              show_row_names = F, show_column_names = T,
              column_names_gp = gpar(fontsize = 9),
              cluster_columns = F,column_names_side = 'bottom',
              column_names_rot = 90
              )
pdf(file = paste0("scRNA_heatmap.pdf"),
    width = 16,height =18)
dev.off()


#####3.Spatial cell score assessment#####
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

#Plot2: violin plot
library("vioplot")
S1=tmp_score_region$value[tmp_score_region$group=='CRSwNP']
S2=tmp_score_region$value[tmp_score_region$group=='Ctrl']
vioplot(S1,S2,
        names=c('CRSwNP','Ctrl'),
        col=rev(subtype.colors),
        ylab = "Eosinophil Spatial Score",main = "CD45 region")+
  mtext(paste0("p=",wilcox.result$p.value,sep=''),side = 3)



