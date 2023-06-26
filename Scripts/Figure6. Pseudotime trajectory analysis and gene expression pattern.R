library(monocle) 
library(Seurat)
load("mydirectory/Basalcell.RData")

## 1). Constructing cell developmental trajectory
data <- as(as.matrix(Basal@assays$RNA@counts), 'sparseMatrix')
pDATA <-  Basal@meta.data
pDATA$cell.type=Basal$Sub_cluster
pd <- new('AnnotatedDataFrame', data =pDATA) 
fDATA <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
fd <- new('AnnotatedDataFrame', data = fDATA)
cds <- newCellDataSet(data,phenoData = pd,featureData = fd,expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds_ordering_genes <- Basal.markers$gene[Basal.markers$p_val_adj<=1e-30]
cds <- setOrderingFilter(cds, ordering_genes = cds_ordering_genes)
cds <- reduceDimension(cds,method = 'DDRTree',max_components = 2)
cds <- orderCells(cds)
plot_cell_trajectory(cds, color_by = "Pseudotime",cell_size = .1,show_branch_points = F)+viridis::scale_color_viridis(option="C")
plot_cell_trajectory(cds, color_by = "cell.type",cell_size = .6,show_branch_points = F)+scale_color_manual(values = c("#ffaaaa","#B8860B","#ff56ff"))
plot_cell_trajectory(cds, color_by = "subtype",cell_size = .6,show_branch_points = F)+facet_wrap(~subtype,nrow = 2)

## 2). Analyzing the gene expression pattern along trajectory
BEAM_res <- BEAM(cds, branch_point = 2, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
pdf("mydirectory/GenePattern.pdf",width = 6,height = 10)
plot_genes_branched_heatmap(cds[row.names(BEAM_res)[1:1000]], branch_point = 2, num_clusters = 6, 
                            cores=4, use_gene_short_name=TRUE, show_rownames=T)
dev.off()

## 3). Comparing the gene expression level/Cell score along trajectory between groups
Basal$pseudotime=cds$Pseudotime
Basal_sNP=subset(Basal,subtype=="CRSsNP_Eth")
Basal_sNP_data=as.data.frame(as.matrix(Basal_sNP@assays$RNA@counts))
dd1=data.frame(Pseudotime=Basal_sNP$pseudotime,exp=t(Basal_sNP_data['STAT3',]))
#dd1=data.frame(Pseudotime=Basal_sNP$pseudotime,exp=Basal_sNP$`IL4 and IL13 signaling`)
colnames(dd1)[2]="exp"

Basal_wNP=subset(Basal,subtype=="CRSwNP_NP")
Basal_wNP_data=as.data.frame(as.matrix(Basal_wNP@assays$RNA@counts))
dd2=data.frame(Pseudotime=Basal_wNP$pseudotime,exp=t(Basal_wNP_data['STAT3',]))
#dd2=data.frame(Pseudotime=Basal_wNP$pseudotime,exp=Basal_wNP$`IL4 and IL13 signaling`)
colnames(dd2)[2]="exp"
subtype.colors=c(CRSsNP_Eth="#7CCD7C",CRSwNP_NP="#CD96CD")  

ggplot() +
  geom_smooth(data = dd1,aes(x = Pseudotime,y = exp,colour = "CRSsNP_Eth",fill="CRSsNP_Eth"),size=1.5,se = T,alpha=0.2)+
  geom_smooth(data = dd2,aes(x = Pseudotime,y = exp,colour ="CRSwNP_NP",fill="CRSwNP_NP"),size=1.5,se = T,alpha=0.2)+
  scale_colour_manual("subtype",values = subtype.colors)+
  scale_fill_manual("subtype",values = subtype.colors)+
  theme_classic()+
  theme(
    axis.text.x=element_text(size = 10,vjust = 0.5, hjust = 0.5,colour = "black"),
    axis.text.y=element_text(size = 10,vjust = 0.5, hjust = 0.5,colour = "black"))+
  xlab("Pseudotime")+ylab("Relative expression")
ggsave(filename = paste0('mydirectory/pseudotime_STAT3.pdf'),width = 4,height = 5)





