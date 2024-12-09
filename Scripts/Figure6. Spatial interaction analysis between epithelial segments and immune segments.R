pacman::p_load(tidyverse, standR, SingleCellExperiment, Seurat, GSVA, msigdbr, ComplexHeatmap, ggpubr, limma, edgeR, fgsea, ggrepel, circlize)
load("mydirectory/GeoMX.RData")

## 1). Dumbbell plots for cell-cell interactions
source('mydirectory/dumbbell_scripts/triangulation_distances.R')

# round centroids to prevent overflow of XYcellID string causing string matching problems
rounded_annots <- annots %>%
  dplyr::mutate(Y_cent = round(Y_cent, 6), 
                X_cent = round(X_cent, 6))

cmap_celltypes <- c('IMM' = '#E41A1C',
                    'MAC' = '#984EA3',
                    'EPI' = '#A65628',
                    'OTHER' = '#999999')

ggplot(data = rounded_annots %>% dplyr::filter(ROI_num == 'TMA016001'),
       mapping = aes(x = X_cent, y = -Y_cent,
                     fill = annotation, color = annotation)) +
  geom_point(shape = 21, size = 1.5, alpha = 0.5) +
  theme_classic() +
  coord_fixed() +
  scale_fill_manual(values = cmap_celltypes) +
  scale_color_manual(values = cmap_celltypes) +
  guides(fill = guide_legend(override.aes = list(size = 3)),
         color = "none") +
  labs(fill = "Cell Annotation") +
  ggtitle("Cell Type Annotations")

triangulation_distances <- get_triangulation_distances(df = rounded_annots,
                                                       id = "cellLabel",
                                                       x_pos = "X_cent",
                                                       y_pos = "Y_cent",
                                                       cell_type = "annotation",
                                                       region = "ROI_num")

# head(triangulation_distances)

per_cell_summary <- triangulation_distances %>%
  group_by(celltype1_index, celltype1, celltype2, ROI_num) %>%
  summarize(per_cell_mean_dist = mean(distance)) %>%
  ungroup()

# head(per_cell_summary)

per_celltype_summary <- per_cell_summary %>%
  group_by(celltype1, celltype2, ROI_num) %>%
  summarize(mean_dist = mean(per_cell_mean_dist)) %>%
  ungroup()

# RUN PERMUTATIONS
iterated_distances_path <- 'mydirectory/iterated_triangulation_distances.rds'

if (file.exists(iterated_distances_path)){
  iterated_triangulation_distances <- readRDS(iterated_distances_path)
} else {
  iterated_triangulation_distances <- iterate_triangulation_distances(df = rounded_annots,
                                                                      id = "cellLabel",
                                                                      x_pos = "X_cent",
                                                                      y_pos = "Y_cent",
                                                                      cell_type = "annotation",
                                                                      region = "ROI_num",
                                                                      num_iterations = 1000,
                                                                      num_cores = 10)
  write_rds(iterated_triangulation_distances, iterated_distances_path)
}

# head(iterated_triangulation_distances)

# distance threshold for observed cell-cell interactions
distance_threshold = 250 # geomx --> 0.4 microns/pixel

observed_distances <- triangulation_distances %>%
  # Append metadata
  left_join(sampleAnnoFile[, c('ROILabel', 'group')] %>% distinct(),
            by = c('ROI_num' = 'ROILabel')) %>%
  dplyr::filter(!group %in% c('ImmuneTissue')) %>%
  dplyr::filter(distance <= distance_threshold) %>%
  # Calculate the average distance to every cell type for each cell
  group_by(celltype1_index, celltype1, celltype2, group, ROI_num) %>%
  summarize(mean_per_cell = mean(distance)) %>%
  ungroup() %>%
  # Calculate the average distance between cell type to cell type on a per group basis
  group_by(celltype1, celltype2, group) %>%
  summarize(observed = list(mean_per_cell),
            observed_mean = mean(unlist(observed), na.rm = TRUE)) %>%
  ungroup()

expected_distances <- iterated_triangulation_distances %>%
  left_join(sampleAnnoFile[, c('ROILabel', 'group')] %>% distinct(),
            by = c('ROI_num' = 'ROILabel')) %>%
  dplyr::filter(!group %in% c('ImmuneTissue')) %>%
  dplyr::filter(mean_dist <= distance_threshold) %>%
  group_by(celltype1, celltype2, group) %>%
  summarize(expected = list(mean_dist),
            expected_mean = mean(mean_dist, na.rm = TRUE)) %>%
  ungroup()

# Calculate pvalues and log fold differences
distance_pvals <- expected_distances %>%
  left_join(observed_distances,
            by = c('celltype1', 'celltype2', 'group')) %>%
  # Calculate wilcoxon test between observed and expected distances
  group_by(celltype1, celltype2, group) %>%
  mutate(pvalue = wilcox.test(unlist(expected), 
                              unlist(observed), 
                              exact = FALSE)$p.value) %>%
  ungroup() %>%
  select(-observed, -expected) %>%
  # Calculate log fold enrichment
  mutate(logfold_group = log2(observed_mean/expected_mean),
         interaction = paste0(celltype1, " --> ", celltype2))

# Get order of plot by magnitude of logfold differences between groups
ord <- (distance_pvals %>%
          select(interaction, group, logfold_group) %>%
          group_by(interaction) %>%
          mutate(largest_diff = max(logfold_group) - min(logfold_group)) %>% # calculate the largest diff between any two groups
          spread(key = group, value = logfold_group) %>%
          dplyr::filter(!is.na(largest_diff)) %>%
          arrange(`Control`))$interaction

distance_pvals$interaction <- factor(distance_pvals$interaction, levels = ord)

distance_pvals <- distance_pvals %>%
  dplyr::filter(celltype1 != 'OTHER' & celltype2 != 'OTHER') %>%
  mutate(sig = factor(ifelse(pvalue < 0.05, T, F), levels = c(T, F)))

ggplot(distance_pvals %>% dplyr::filter(!is.na(interaction)) %>% dplyr::mutate(group = ordered(group, levels = c('Control', 'CRSsNP', 'CRSwNP_NP', 'CRSwNP_UNC')))) +
  geom_vline(aes(xintercept = 0), linetype = 'dashed') +
  geom_line(aes(x = logfold_group, y = interaction), na.rm = T) +
  geom_point(aes(x = logfold_group, y = interaction, fill = group, shape = group, alpha = sig),
             size = 4, stroke = 0.5, na.rm = T) +
  scale_shape_manual(values = c('Control' = 22, 'CRSsNP' = 21, 'CRSwNP_NP' = 24, 'CRSwNP_UNC' = 25)) +
  scale_fill_manual(values = c(Control = '#FF7256', CRSsNP = '#7CCD7C', CRSwNP_NP = '#CD96CD', CRSwNP_UNC = '#AAD4FF')) +
  scale_alpha_manual(values = c('TRUE' = 1.0, 'FALSE' = 0.2), drop = F, name = 'p < 0.05') +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 14, angle = 0),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

ggsave('mydirectory/Dumbbell.pdf', device = 'pdf', width = 7, height = 5, bg = 'white')

## 2). Co-expressing analysis based on cNMF
# Prep data for cNMF
tt = as.data.frame(t(DSP_spe@assays@data$logcounts))
write.csv(tt, '../data/crs_logcount.csv', row.names = FALSE)

tt2 = as.data.frame(t(DSP_spe@assays@data$counts))
write.csv(tt2, '../data/crs_count.csv', row.names = FALSE)

df = as.data.frame(DSP_spe@colData)
write.csv(df, '../data/crs_meta.csv', row.names = FALSE)

# Run cNMF with cNMF_using_DSP_data.ipynb

# Load cNMF result files
cnmf_un = read.csv('../data/crs_cnmf35_usage.csv', row.names = 1)
meta = read.csv('../data/crs_meta.csv')
colnames(cnmf_un) = gsub('Usage', 'Program', colnames(cnmf_un))

cnmf_un$tissueid = meta$ROILabel
cnmf_imm_snp = subset(cnmf_un, meta$SegmentLabel == 'IMM' & meta$group == 'CRSsNP')
cnmf_imm_wnp = subset(cnmf_un, meta$SegmentLabel == 'IMM' & meta$group == 'CRSwNP')

cnmf_epi_snp = subset(cnmf_un, meta$SegmentLabel == 'EPI' & meta$group == 'CRSsNP')
cnmf_epi_wnp = subset(cnmf_un, meta$SegmentLabel == 'EPI' & meta$group == 'CRSwNP')

cnmf_epi_wnp = cnmf_epi_wnp[match(cnmf_imm_wnp$tissueid, cnmf_epi_wnp$tissueid),]

cnmf_imm_snp$tissueid = NULL
cnmf_imm_wnp$tissueid = NULL
cnmf_mac_snp$tissueid = NULL
cnmf_mac_wnp$tissueid = NULL
cnmf_epi_snp$tissueid = NULL
cnmf_epi_wnp$tissueid = NULL

library(RColorBrewer)
## make color break list
l = 0.8
r = -0.8
breaksList = seq(r, l, by = 0.1)
k = 35
## correlation
library(CCA)

# Identify co-expressing modules in CRSwNP
test = matcor(cnmf_imm_wnp, cnmf_epi_wnp)
plot_test = test$XYcor[c(1:k),c((k+1):(2*k))]
# remove rows or columns with all NAs
plot_test = plot_test[rowSums(is.na(plot_test)) != ncol(plot_test), ]
plot_test = plot_test[,colSums(is.na(plot_test)) != nrow(plot_test) ]
plot_test[plot_test> l] = l
plot_test[plot_test< r] = r

library(pheatmap)

pdf(file="../plots/imm-epi_wnp_35.pdf",height = 20, width = 20)
pheatmap(plot_test,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), 
         breaks = breaksList)
dev.off()

# Interpretate the process (gene pathways) involved in the interaction hotspots we identified

topg = read.csv('../data/top100_genes_cnmf35.csv')

library(enrichR)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")

# 2, 5, 9
genes = topg[c(1:20),c(2,5,9)]
enriched <- enrichr(genes = topg[c(1:20),c(2,5,9)], dbs)
a = enriched$GO_Biological_Process_2015
print(enriched$GO_Biological_Process_2015$Term[c(1:10)])
print(enriched$GO_Biological_Process_2015$Adjusted.P.value[c(1:10)])

## 3). Ligand-receptor communication analysis
# Run cellphoneDB with LR_analysis.ipynb
# Load the result files
base_dir <- 'mydirectory/run_cellphonedb'
pvals <- read.delim(file.path(base_dir, '/outputs/statistical_analysis_pvalues_03_28_2024_114344.txt'), check.names = F)
means <- read.delim(file.path(base_dir, 'outputs/statistical_analysis_means_03_28_2024_114344.txt'), check.names = F)
decon <- read.delim(file.path(base_dir, '/outputs/statistical_analysis_deconvoluted_03_28_2024_114344.txt'), check.names = F)
signif_means <- read.delim(file.path(base_dir, '/outputs/statistical_analysis_significant_means_03_28_2024_114344.txt'), check.names = F)

cell_interactions <- colnames(signif_means %>% dplyr::select(matches('\\|')))

mdata <- read_csv(file.path(base_dir, 'mdata.csv'))

counts <- read.csv(file.path(base_dir, 'log2_counts.csv')) %>%
  column_to_rownames(var = 'Target') %>%
  as.data.frame()


custom_genesets <- list(chemokines = grep('^CXC|CCL|CCR|CX3|XCL|XCR|^IL', signif_means$interacting_pair, value = T))

set <- 'ALL'
celltype1 <- 'ALL'
celltype2 <- 'ALL'

# for scaling mean expression
scale_means <- function(x) {
  (x - min(x))/(max(x) - min(x))
}

# subset and scale mean paired expression values
all_values <- means %>%
  {if(set != 'ALL') dplyr::filter(interacting_pair %in% custom_genesets[[set]]) else .} %>%
  {if(celltype1 != 'ALL') dplyr::select(., interacting_pair, any_of(grep(celltype1, cell_interactions, value = T, ignore.case = T))) else .} %>%
  {if(celltype2 != 'ALL') dplyr::select(., interacting_pair, any_of(grep(celltype2, cell_interactions, value = T, ignore.case = T))) else .} %>%
  group_by(interacting_pair) %>%
  dplyr::summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>% # take mean of duplicate pairs
  column_to_rownames(var = 'interacting_pair') %>%
  as.matrix()

all_values <- apply(all_values, 1, scale_means) %>% t()

# get significance from signif_means
signif <- signif_means %>%
  {if(set != 'ALL') dplyr::filter(interacting_pair %in% custom_genesets[[set]]) else .} %>%
  {if(celltype1 != 'ALL') dplyr::select(., interacting_pair, any_of(grep(celltype1, cell_interactions, value = T, ignore.case = T))) else .} %>%
  {if(celltype2 != 'ALL') dplyr::select(., interacting_pair, any_of(grep(celltype2, cell_interactions, value = T, ignore.case = T))) else .} %>%
  group_by(interacting_pair) %>%
  dplyr::summarise(across(where(is.numeric), ~mean(., na.rm = TRUE))) %>% # take mean of duplicate pairs
  column_to_rownames(var = 'interacting_pair') %>%
  as.matrix()

signif <- signif[rownames(all_values), colnames(all_values), drop = F] # match row/col order
signif <- ifelse(!is.na(signif), 1, NA) # set values to 1

only_signif <- all_values * signif

# write to csv
only_signif %>% as.data.frame() %>%
  dplyr::filter(!is.na(`IMM_CRSwNP_NP|EPI_CRSwNP_NP`) | !is.na(`EPI_CRSwNP_NP|IMM_CRSwNP_NP`) | !is.na(`IMM_CRSwNP_UNC|EPI_CRSwNP_UNC`) | !is.na(`IMM_CRSsNP|EPI_CRSsNP`) | !is.na(`EPI_CRSsNP|IMM_CRSsNP`)) %>%
  rownames_to_column('pair') %>%
  write_csv('mydirectory/run_cellphonedb/outputs/only_signif.csv')

# for scaling mean expression
scale_means <- function(x) {
  (x - min(x))/(max(x) - min(x))
}

# LR pairs for plotting
interactions <- c('IMM_CRSwNP_NP|EPI_CRSwNP_NP', 'EPI_CRSwNP_NP|IMM_CRSwNP_NP', 'IMM_CRSwNP_UNC|EPI_CRSwNP_UNC', 'EPI_CRSwNP_UNC|IMM_CRSwNP_UNC', 'IMM_CRSsNP|EPI_CRSsNP', 'EPI_CRSsNP|IMM_CRSsNP')
curated_pairs <- read_csv(file.path(base, 'run_cellphonedb/outputs/lr_curated.csv')) %>% dplyr::pull(pair)
curated_means <- means %>% dplyr::filter(interacting_pair %in% curated_pairs) %>% dplyr::select(interacting_pair, all_of(interactions)) %>%
  pivot_longer(-interacting_pair, names_to = 'interaction', values_to = 'mean') %>%
  mutate(scaled_mean = scale_means(mean)) %>%
  left_join(pvals %>% dplyr::filter(interacting_pair %in% curated_pairs) %>% 
              dplyr::select(interacting_pair, all_of(interactions)) %>% 
              pivot_longer(-interacting_pair, names_to = 'interaction', values_to = 'pval')) %>%
  mutate(signif = pval < 0.05) %>%
  mutate(interaction = ordered(str_replace(interaction, '\\|', ' -> '), levels = gsub('\\|', ' -> ', interactions)))

p <- curated_means %>%
  ggplot(aes(x = interaction, y = interacting_pair)) +
  geom_point(aes(fill = scaled_mean, size = scaled_mean, color = signif), shape = 21, size = 2.0) +
  scale_color_manual(values = c('FALSE' = 'white', 'TRUE' = 'black'), guide = 'none') +
  scale_fill_gradient2(low = 'white', high = '#bf0000', midpoint = 0.55) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  theme_bw()

ggsave('~/LR_analysis.pdf', p, device = 'pdf', width = 6, height = 6)


## 4). Geometric cell-cell interaction analysis
library(deldir)
library(igraph)
library(progress)
library(qs)
library(ggtern)
library(epitools)
Celltypepairs <- c("EPI","IMM")
source("mydirectory/WithinCrossDomainTableAndVisualization_preload_function.R")
figure_save_path <- file.path("mydirectory/GeometricInteraction/")
# setup your cell type color
cell_type_colors <- c("EPI" = "#f300ff",   # Green
                      "IMM" = "#59ff00", # Purple
                      "MAC" = "#66cbdb",      # Blue-Violet
                      "OTHER" = "#a87a66") # Brown

# Load the ROI annotation data
annotation <- read_csv("final_annotations.csv")
names(annotation)[names(annotation) == "annotation"] <- "Annotation"

drop.sample  <- annotation %>%
  group_by(ROI_num) %>%
  summarize(
    has_IMM = any(Annotation == "IMM"),
    has_EPI = any(Annotation == "EPI")
  ) %>%
  # Filter out ROI_num that do not have both annotations
  filter(!has_IMM | !has_EPI) %>%
  select(ROI_num)

filtered_annotation <- annotation %>%
  filter(!ROI_num %in% unlist(drop.sample))

stat.result <- list()
for (core in unique(filtered_annotation$ROI_num)) {
  print(paste("Working on", core, "..."))
  
# Subset data for the current core
  sub_df <- filtered_annotation %>%
    filter(ROI_num == core) 
  result <- generate_edge_data(sub_df,cell_types = c(Celltypepairs[1], Celltypepairs[2]),pixel_to_um = 0.5)
  p.print <- plot_cell_data_with_colors(edge_data = result$edge_data,
                                        cell_data = result$cell_data,
                                        cell_type_colors = cell_type_colors,
                                        cell_types = c(Celltypepairs[1], Celltypepairs[2]),size_domain = 2,size_other = 1
  )#+theme_void()
  save.name.fig <- paste0(core,"_",Celltypepairs[1],"_", Celltypepairs[2],"_Phenotype.tiff")
  ggsave(plot = p.print,filename = file.path(figure_save_path,save.name.fig),
         device = "tiff",dpi = 200,width = 8,height = 8)
  stat.result.tmp <- (calculate_odds(cell_data = result$cell_data,
                                     edge_data =result$edge_data ,
                                     ct1_name = Celltypepairs[1],
                                     ct2_name = Celltypepairs[2]))
  stat.result[[core]] <- stat.result.tmp
}

qsave(stat.result, file = file.path(figure_save_path,"stat.resultv2.qs"))

stat.result <- qread("stat.resultv2.qs")

# check gene expression
lr_curated <- read.csv("lr_curated.csv")
lr_curated <- lr_curated[,c(1,grep("CRSwNP",colnames(lr_curated)))]
sgtatsiticalSignaficance <- read.delim("statistical_analysis_significant_means_03_28_2024_114344.txt")

exp_count <- DSP_spe@assays@data$logcounts
meta_information <- DSP_spe@colData %>%
  as.data.frame() %>%
  select(Working_group,ROILabel, SlideName, SegmentLabel,)
identical(rownames(meta_information),colnames(exp_count))

# Define ligand-receptor pairs
LRP_pair <- list(
  "CCL26" = "CCR3",
  "CCL20" = "CCR6",
  "TSLP" = "IL7R",
  "IL33" = "IL1RL1"
)

# Extract expression values for ligands and receptors
ligands <- names(LRP_pair)
receptors <- unlist(LRP_pair)

# Create new columns for ligand and receptor expression
for (lig in ligands) {
  meta_information[, paste0(lig, "_expression")] <- as.numeric(exp_count[lig, rownames(meta_information)])
}

for (rec in receptors) {
  meta_information[, paste0(rec, "_expression")] <- as.numeric(exp_count[rec, rownames(meta_information)])
}

# Calculate ligand and receptor interaction score for EPI and IMM within each ROI

# Function to calculate the nth root of the product of all ligand-receptor interaction scores
calc_n_root_product <- function(values) {
  product_value <- prod(values, na.rm = TRUE)
  n <- length(values)
  return(product_value^(1/n))
}

# Loop through the ROIs and calculate interaction scores

meta_information$Overall_Geometric_Mean <- NA
meta_information$CCL26_CCR3_coexp<- NA
meta_information$CCL20_CCR6_coexp <- NA
meta_information$TSLP_IL7R_coexp <- NA
meta_information$IL33_IL1RL1_coexp <- NA

for (roi in unique(meta_information$ROILabel)) {
  # Subset for the current ROI
  epi_ligand <- meta_information$SegmentLabel == "EPI" & meta_information$ROILabel == roi
  imm_receptor <- meta_information$SegmentLabel == "IMM" & meta_information$ROILabel == roi
  
  # Initialize vectors to store interaction scores for each ligand-receptor pair
  interaction_scores <- c()
  
  # Iterate over each ligand-receptor pair
  for (lig in ligands) {
    rec <- LRP_pair[[lig]]
    
    # Extract ligand and receptor expression for EPI and IMM
    ligand_expression <- meta_information[epi_ligand, paste0(lig, "_expression")]
    receptor_expression <- meta_information[imm_receptor, paste0(rec, "_expression")]
    
    if (length(ligand_expression) > 0 && length(receptor_expression) > 0) {
      # Calculate the interaction score (product of ligand and receptor)
      interaction_score <- ligand_expression * receptor_expression/
        interaction_scores <- c(interaction_scores, interaction_score)
    }
  }
  if (length(interaction_scores) > 0) {
    meta_information[epi_ligand, c("CCL26_CCR3_coexp","CCL20_CCR6_coexp","TSLP_IL7R_coexp","IL33_IL1RL1_coexp")] <- interaction_scores
  }
  # Calculate the nth root of the product of interaction scores (arithmetic mean)
  if (length(interaction_scores) > 0) {
    meta_information[epi_ligand, "Overall_Geometric_Mean"] <- calc_n_root_product(interaction_scores)
    meta_information[,]
  }
  
  # Calculate the geometric mean as a product of pairwise scores
  
}

# View the updated meta_information with overall arithmetic and geometric means
head(meta_information)



# View the updated meta_information with interaction scores
head(meta_information)
epi_CCC <- meta_information %>%
  filter(SegmentLabel == "EPI")
head(epi_CCC)
# job ID 3915575
###
# Initialize an empty list to store the results for each TMA ID
data_list <- list()

# Loop through each TMA ID in stat.result
for (tma_id in names(stat.result)) {
  # Extract the values for each TMA ID and convert it to a data frame
  tma_data <- as.data.frame(t(unlist(stat.result[[tma_id]])))
  
  # Add the TMA ID as a column
  tma_data$TMA_ID <- tma_id
  
  # Append to the list
  data_list[[tma_id]] <- tma_data
}

# Combine the list into a single data frame
stat.result.df <- do.call(rbind, data_list)

# Reorder the columns to have TMA_ID as the first column
stat.result.df <- stat.result.df[, c("TMA_ID", setdiff(names(stat.result.df), "TMA_ID"))]

# Display the first few rows of the resulting data frame
head(stat.result.df)

stat.result.df$ROILabel <- rownames(stat.result.df)
##
stat.result.df <- left_join(stat.result.df, epi_CCC, by = "ROILabel")
table(DSP_spe$Working_group)
# Control       CRSsNP       CRSwNP

# extract condition

# Extract the columns for ternary plotting
ternary_data <- stat.result.df
ternary_data$GeometricOddsratio <- sqrt((ternary_data$odds_across/ternary_data$odds_within_ct1) * (ternary_data$odds_across/ternary_data$odds_within_ct2))
ternary_data$weighted_CCL26_CCR3_coexp <- ternary_data$CCL26_CCR3_coexp/ternary_data$GeometricOddsratio
ternary_data$weighted_CCL20_CCR6_coexp <- ternary_data$CCL20_CCR6_coexp/ternary_data$GeometricOddsratio
ternary_data$weighted_TSLP_IL7R_coexp <- ternary_data$TSLP_IL7R_coexp/ternary_data$GeometricOddsratio
ternary_data$weighted_IL33_IL1RL1_coexp <- ternary_data$IL33_IL1RL1_coexp/ternary_data$GeometricOddsratio
ternary_data$weighted_attraction_LRP2 <- log(apply(ternary_data[,c("weighted_CCL26_CCR3_coexp",
                                                                   "weighted_CCL20_CCR6_coexp",
                                                                   "weighted_TSLP_IL7R_coexp",
                                                                   "weighted_IL33_IL1RL1_coexp")],
                                                   1,calc_n_root_product))

write.csv(ternary_data,file = file.path(figure_save_path,paste0(Celltypepairs[1],"_", Celltypepairs[2],"ternary_data.csv")))
# Create the ternary plot
ternary_data <- ternary_data %>% filter(!is.na(Working_group))
colnames(ternary_data)
ggtern(data = ternary_data, aes(x = pp_within_ct1, y = pp_across, z = pp_within_ct2, color = weighted_attraction_LRP2)) +
  labs(x = "EPI-EPI",y="EPI-IMM",z="IMM-IMM")+
  geom_point(size =2) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)+
  ggtitle("Ternary Plot of pp_across, pp_within_ct1, pp_within_ct2") +
  theme_bw()+
  facet_wrap(~Working_group,ncol = 2)
# boxplot
ternary_data_sub <- ternary_data %>% filter(!is.na(Working_group)) %>%
  filter(!Working_group%in%c("CRSwNP_UNC","ImmuneTissue","Exclude") )

p <- ggboxplot(ternary_data_sub,
               x = "Working_group",
               y = "weighted_attraction_LRP2",
               color = "Working_group",
               palette = "jco",
               add = "jitter") + geom_pwc(method = "wilcoxon")
p.adj.barplot <-
  ggadjust_pvalue(p,
                  p.adjust.method = "bonferroni"
  )
ggsave("Interaction_EPI_IMM_plot.svg",device = "svg", plot = p.adj.barplot, width = 8, height = 6)



