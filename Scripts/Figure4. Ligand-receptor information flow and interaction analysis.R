library(CellChat)
library(ggplot2)
library(forcats)
library(RColorBrewer)

load("mydirectory/Allcell.RData")
Allcell_CRSwNP_NP=subset(Allcell,subtype=="CRSwNP_NP")
Allcell_CRSsNP_Eth=subset(Allcell,subtype=="CRSsNP_Eth")
## 1). Calculate Ligand-receptor signal separately with cellchat
cellchat_wNP_NP <-  createCellChat(object = Allcell_CRSwNP_NP, meta = Allcell_CRSwNP_NP@meta.data, group.by = "cell_cluster") 
cellchat_wNP_NP@DB <- CellChatDB.human
cellchat_wNP_NP <- subsetData(cellchat_wNP_NP)
cellchat_wNP_NP <- identifyOverExpressedGenes(cellchat_wNP_NP)
cellchat_wNP_NP <- identifyOverExpressedInteractions(cellchat_wNP_NP) 
cellchat_wNP_NP <- projectData(cellchat_wNP_NP, PPI.human)
cellchat_wNP_NP <- computeCommunProb(cellchat_wNP_NP,raw.use = FALSE,type="truncatedMean",trim=0.1)
cellchat_wNP_NP <- filterCommunication(cellchat_wNP_NP, min.cells = 10)
cellchat_wNP_NP <- computeCommunProbPathway(cellchat_wNP_NP)
cellchat_wNP_NP <- aggregateNet(cellchat_wNP_NP)
#Do the same for CRSsNP_Eth

## 2). Combine the results of two subtypes and compare information flow
object.list <- list(CRSsNP_Eth = cellchat_sNP_Eth, CRSwNP_NP = cellchat_wNP_NP)
cellchat_merge <- mergeCellChat(object.list, add.names = names(object.list))
rankNet(cellchat_merge,
        sources.use = 13,
        targets.use = 2,
        color.use = c("#CD96CD","#7CCD7C"),mode = 'comparison',comparison = c(2,1),
        cutoff.pvalue = 0.05,thresh = 0.05,stacked = T,do.stat = T)+
        scale_fill_manual(values = subtype.colors)

## 3). Ligand-receptor communication bubble plot
cellchat_merge@idents$joint=fct_relevel(cellchat_merge@idents$joint,unique(Allcell$Sub_cluster))
netVisual_bubble(cellchat_merge,
  sources.use  = 13,
  targets.use = 2,
  signaling = "IL4",
  thresh = 10,
  comparison = c(1,2),
  angle.x = 45,font.size = 10,
  remove.isolate = F,
  color.text = subtype.colors)+
  scale_color_distiller(palette = "Reds",direction = 1)+
  geom_point(shape = 21, color = "black") + 
  theme(axis.text.x = element_text(angle = 45,size = 8,vjust = 1),
        axis.text.y = element_text(size = 8))





