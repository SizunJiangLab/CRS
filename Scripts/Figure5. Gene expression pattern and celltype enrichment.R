pacman::p_load(tidyverse, standR, SingleCellExperiment, Seurat, GSVA, msigdbr, ComplexHeatmap, ggpubr, limma, edgeR, fgsea, ggrepel, circlize)
load("mydirectory/GeoMX.RData")

## 1).  Gene expression pattern of each segment
# limma to separate the segment types
log2_counts <- log2_counts[, sampleAnnoFile$SegmentDisplayName] + 1 # reorder matrix columns, add pseudocount of 1 for negative values

groups <- factor(sampleAnnoFile$SegmentLabel)
design <- model.matrix(~0 + groups)
colnames(design) <- str_remove_all(colnames(design), 'groups')

dge <- DGEList(counts = log2_counts, genes = row.names(log2_counts))

v <- voom(dge, design, plot = F)

fit <- lmFit(v)

contr.matrix <- makeContrasts(
  EPI_v_OTH = EPI - (MAC + IMM)/2,
  IMM_v_OTH = IMM - (EPI + MAC)/2,
  MAC_v_OTH = MAC - (EPI + IMM)/2,
  levels = colnames(design)
)

fit_contrast <- contrasts.fit(fit, contrasts = contr.matrix)

efit <- eBayes(fit_contrast, robust = TRUE)

results_efit <- decideTests(efit, p.value = 0.05)
summary_efit <- summary(results_efit)

write.fit(efit, results = NULL, 'mydirectory/celltype_fit.csv', digits = NULL, adjust = 'BH', method = 'separate', sep = ',')

results <- read.csv('mydirectory/celltype_fit.csv', row.names = 1)

# HEATMAP of DEGs

# subset from deg analysis
epi_imm_genes_to_plot <- results %>%
  dplyr::select(Coef.EPI_v_OTH, Coef.IMM_v_OTH) %>%
  rownames_to_column('Target') %>%
  pivot_longer(starts_with('Coef.'), names_to = 'contrast', values_to = 'log2fc') %>%
  group_by(contrast) %>%
  slice_max(order_by = log2fc, n = 20) %>%
  dplyr::pull(Target)

# subset from sc macrophage clusters and intersect with highest DEGs in MAC_v_OTH
mac_genes <- readxl::read_xlsx('mydirectory/Macrophage markers.xlsx') %>%
  slice_max(order_by = avg_log2FC, n = 50) %>%
  dplyr::pull(gene)

mac_genes_to_plot <- results %>%
  dplyr::select(Coef.MAC_v_OTH) %>%
  rownames_to_column('Target') %>%
  dplyr::filter(Target %in% mac_genes) %>%
  slice_max(order_by = Coef.MAC_v_OTH, n = 20) %>%
  dplyr::pull(Target)

genes_to_plot <- c(epi_imm_genes_to_plot, mac_genes_to_plot)

filtered_counts <- log2_counts %>% t() %>% data.frame() %>%
  rownames_to_column(var = 'ROI') %>%
  separate(ROI, into = c('TMA', 'ROI', 'SegmentLabel'), sep = ' \\| ') %>%
  pivot_longer(cols = c(-TMA, -ROI, -SegmentLabel), names_to = 'Target', values_to = 'log2counts') %>%
  group_by(SegmentLabel, Target) %>%
  summarise(mean_log2counts = mean(log2counts), .groups = 'drop') %>%
  dplyr::filter(Target %in% genes_to_plot) %>%
  group_by(Target) %>%
  dplyr::mutate(z_mean_counts = scale(mean_log2counts)[,1]) %>%
  dplyr::select(-mean_log2counts) %>%
  pivot_wider(names_from = SegmentLabel, values_from = z_mean_counts) %>%
  column_to_rownames(var = 'Target') %>%
  as.matrix()

col_fun = colorRamp2(c(-2, 0, 2), c("#2188BB", "white", "#C00000"))

set.seed(2024)
pdf('mydirectory/GenePattern.pdf', width = 10, height = 8)
draw(Heatmap(filtered_counts, name = 'z-score',
             cluster_columns = F,
             cluster_rows = T,
             clustering_method_rows = 'centroid',
             row_names_gp = gpar(fontsize = 10),
             col = col_fun,
             heatmap_legend_param = list(legend_direction = "horizontal",
                                         legend_width = unit(4, "cm"))),
     padding = unit(c(10, 10, 10, 120), 'mm'),
     heatmap_legend_side = 'top')
dev.off()

## 2). Annotated cell subtype enrichment
# log2 fold enrichment
annots %>%
  dplyr::filter(!annotation == 'OTHER') %>%
  dplyr::filter(TissueType %in% c('CRSwNP', 'CRSsNP')) %>%
  group_by(annotation, TissueType) %>%
  dplyr::summarise(n = n(), .groups = 'drop') %>%
  group_by(TissueType) %>%
  dplyr::mutate(total_n = sum(n)) %>%
  ungroup() %>%
  group_by(annotation) %>%
  dplyr::mutate(diff = log2((lag(n)/lag(total_n)) / (n/total_n))) %>%
  drop_na() %>%
  ggplot(aes(x = reorder(annotation, diff), y = diff, fill = annotation)) +
  geom_col() +
  scale_fill_manual(values = c(EPI = '#CC6C25', IMM = '#008E47', MAC = '#036EB8')) +
  geom_rect(aes(xmin = 0, ymin = 0, xmax = 4, ymax = Inf), fill = '#7CCD7C', alpha = 0.05) +
  geom_rect(aes(xmin = 0, ymin = -Inf, xmax = 4, ymax = 0), fill = '#CD96CD', alpha = 0.05) +
  ylim(-0.35, 0.35) +
  annotate('text', x = 0.5, y = 0.30, label = 'CRSsNP', size = 8) +
  annotate('text', x = 3.5, y = -0.30, label = 'CRSwNP', size = 8) +
  labs(x = 'Cell Type', y = 'Log2 Fold Enrichment', title = 'Cell Type Enrichment (CRSsNP / CRSwNP)') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = 'none',
        text = element_text(size = 12))

ggsave('mydirectory/CellEnrichment.pdf', device = 'pdf', width = 6.8, height = 4)

# comparison of macrophage segments
annots %>%
   dplyr::filter(TissueType %in% c('CRSwNP', 'CRSsNP', 'Control')) %>%
   group_by(ROI_num, annotation, TissueType) %>%
   dplyr::summarise(n = n(), .groups = 'drop') %>%
   group_by(ROI_num) %>%
   dplyr::mutate(prop = n / sum(n)) %>%
   ungroup() %>%
   dplyr::filter(annotation == 'MAC') %>%
   dplyr::group_by(TissueType) %>%
   ggplot(aes(x = TissueType, y = prop)) +
     geom_boxplot(aes(fill = TissueType), outlier.alpha = 0) +
     scale_fill_manual(values = c(Control = 'FF7256', CRSsNP = '7CCD7C', CRSwNP = 'CD96CD')) +
     geom_jitter(alpha = 0.5, width = 0.1) +
     stat_compare_means(method = 'wilcox.test',
                        comparisons = list(c('Control', 'CRSsNP'), c('Control', 'CRSwNP'), c('CRSsNP', 'CRSwNP'))) +
     labs(title = 'ROI Proportion of Myeloid cells') +
     theme_bw()
 
ggsave('mydirectory/ROIProportion.pdf', device = 'pdf', width = 6, height = 6)




