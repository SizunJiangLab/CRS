#####Differential gene expression and Pathway enrichment analysis in scRNA-seq data#####
library(Seurat)
library(ggplot2)
library(SeuratData)
library(patchwork)
library(clusterProfiler)
library(tidyverse)
library(ReactomePA)
library(Hmisc)

load("mydirectory/immune.combined.RData")
## 1). Differential gene expression  analysis, taking CD8+ T as an example
CD8=subset(immune.combined,cell_cluster=="CD8+T")
CD8@active.ident=factor(CD8$subtype)
CD8.subtype.DEG <- FindMarkers(CD8, ident.1 = "CRSwNP_NP",ident.2 = "CRSsNP_Eth",only.pos = F, min.pct = 0.1, logfc.threshold = 0)

## 2). GSEA
Enrichmentdata=CD8.subtype.DEG[,c(2,5)]
Enrichmentdata$gene_name=rownames(Enrichmentdata)
GSEAgene.df <- data.frame(logFC=Enrichmentdata$avg_log2FC,GENE=Enrichmentdata$gene_name)
GSEAgenelist <- GSEAgene.df$logFC
names(GSEAgenelist) <- GSEAgene.df$GENE
GSEAgenelist <- sort(GSEAgenelist,decreasing = T)
reactome_gmt <- read.gmt("mydirectory/c2.cp.reactome.v7.5.1.symbols.gmt")
CD8_GSEA_REACTOME<-GSEA(GSEAgenelist,TERM2GENE = reactome_gmt,pvalueCutoff = 1,eps = 0)

## 3). Visualization
up_draw=CD8_GSEA_REACTOME@result
up_draw=subset(up_draw,enrichmentScore>0)
up_draw=subset(up_draw,p.adjust<0.05)
up_draw$LogP=-log10(up_draw$p.adjust)
up_draw=up_draw[,c('Description','LogP')]
up_draw$type='CRSwNP_NP'

down_draw=CD8_GSEA_REACTOME@result
down_draw=subset(down_draw,enrichmentScore<0)
down_draw=subset(down_draw,p.adjust<0.05)
down_draw$LogP=-log10(down_draw$p.adjust)
down_draw=down_draw[,c('Description','LogP')]
down_draw$type='CRSsNP_Eth'
down_draw$LogP=-down_draw$LogP

CD8_pathway_draw=rbind(up_draw,down_draw)
pathway_draw=CD8_pathway_draw
pathway_draw=pathway_draw %>% arrange(desc(LogP))
pathway_draw$Description=gsub("REACTOME_","",pathway_draw$Description)
pathway_draw$Description=tolower(pathway_draw$Description)
pathway_draw$Description=capitalize(pathway_draw$Description)
pathway_draw$Description=as.factor(pathway_draw$Description)
ggplot(pathway_draw,aes(x=LogP,y=reorder(Description,LogP),fill=type))+geom_bar(stat = 'identity',position = position_dodge(1.5))+
  theme_classic()+
  scale_fill_manual(values = subtype.colors)+
  theme(axis.text.y = element_text(colour = "black"),axis.text.x = element_text(colour = "black",size = 10))+
  xlab('-log10(P-value)')+ylab('')
