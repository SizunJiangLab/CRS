# CRS Project Data Description

## Overview

This document describes the data available for the Chronic rhinosinusitis (CRS) project, including single-cell RNA sequencing (scRNA-seq) and Geographical Molecular Profiler (GeoMX) datasets. The purpose of this document is to provide a comprehensive understanding of the data formats, structures, and content to facilitate further analysis.

## Data Types

The CRS project contains three main types of data:

1. **Public scRNA-seq Data**: A publicly available dataset from a Nature Immunity 2022 publication.
2. **In-house scRNA-seq Data**: Two datasets representing epithelial cells and immune cells.
3. **In-house GeoMX Data**: Consisting of both sequencing data and imaging data.

## 1. Public scRNA-seq Dataset

### Data Location

- File: `CRS_Data/NatureImmunitydataPublic2022/crs.all_seurat_jan24_24.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

### Description

This is a publicly available scRNA-seq dataset from a Nature Immunity 2022 publication related to Chronic rhinosinusitis. The dataset has been processed and stored as a Seurat object.

### Detailed Analysis

- **Object class**: Seurat
- **Number of cells**: 91,212
- **Number of features (genes)**: 32,643
- **Available assays**: RNA, SCT
- **Dimensional reductions**: PCA, Harmony, UMAP (all available)
- **Disease groups**:
  - CRSsNP: 21,605 cells
  - Control: 27,550 cells
  - eCRSwNP: 14,240 cells
  - neCRSwNP: 27,817 cells
- **Cell types**: 24 distinct cell types identified
- **Batch information**: 12 different batches (n00, n05, n06, n07, n10, n11, n12, n94, n95, n96, n97, n99)

The public dataset is well-processed, with normalized data, batch correction using Harmony, and dimensionality reduction already performed. Cell clustering and annotation have also been completed, making this dataset ready for downstream analysis.

## 2. In-house scRNA-seq Datasets

### Dataset 1: Epithelial Cells

- File: `CRS_Data/InHouseScRNA/Epi_reannotated.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

#### Analysis

- **Object class**: Seurat
- **Number of cells**: 21,833
- **Number of features (genes)**: 24,407
- **Available assays**: RNA, ADT, integrated
- **Dimensional reductions**: PCA, UMAP (both available)
- **Sample types**:
  - Control_Eth: 7,034 cells
  - CRSsNP_Eth: 5,841 cells
  - CRSwNP_Eth: 1,457 cells
  - CRSwNP_NP: 7,501 cells
- **Cell subtypes**: 13 epithelial cell subtypes, including:
  - Basal cells: 2,242 cells
  - Ciliated cells: 4,709 cells
  - Cycling basal cells: 342 cells
  - Deuterosomal cells: 94 cells
  - FOXJ1low Ciliated cells: 812 cells
  - Goblet cells: 1,285 cells
  - Intermediate ciliated cells: 552 cells
  - Ionocytes: 383 cells
  - Mucous cell: 156 cells
  - Secretory cells: 6,544 cells
  - Serous cells: 2,187 cells
  - Suprabasal cells: 2,412 cells
  - Tuft cells: 115 cells

### Dataset 2: Immune Cells

- File: `CRS_Data/InHouseScRNA/ImmuneCell_processed.rds`
- Format: R Seurat object (.rds)
- Access Method: `readRDS()` function in R

#### Analysis

- **Object class**: Seurat
- **Number of cells**: 32,775
- **Number of features (genes)**: 24,407
- **Available assays**: RNA, ADT, integrated
- **Dimensional reductions**: PCA, UMAP (both available)
- **Sample types**:
  - Control_Eth: 6,633 cells
  - CRSsNP_Eth: 10,193 cells
  - CRSwNP_Eth: 4,510 cells
  - CRSwNP_NP: 11,439 cells
- **Cell types**: 13 immune cell types, including:
  - B cells: 2,545 cells
  - CD8 T cells: 6,017 cells
  - Macrophages: 1,502 cells
  - Mast cells: 985 cells
  - Mono/Macrophages: 2,890 cells
  - Naive T cells: 3,002 cells
  - Neutrophils: 2,154 cells
  - NK cells: 1,351 cells
  - NKT cells: 1,049 cells
  - pDCs: 832 cells
  - Plasma cells: 1,439 cells
  - Regulatory T cells: 1,467 cells
  - T cells: 5,540 cells

Both in-house scRNA-seq datasets are well-processed with integrated analysis, dimensionality reduction, and cell type annotation. The datasets include both RNA and Antibody-Derived Tag (ADT) assays, suggesting multi-modal data acquisition that combines transcriptomic and protein-level information.

## 3. GeoMX Datasets

### 3.1 GeoMX Sequencing Data

#### Data Files

- Raw Counts: `CRS_Data/GeoMX_sequencing_data/countFile_raw.csv`
- Processed Counts: `CRS_Data/GeoMX_sequencing_data/countFile_normalized_and_batch_effect_corrected.csv`
- Feature Annotations: `CRS_Data/GeoMX_sequencing_data/featureAnnoFile.csv`
- Sample Annotations: `CRS_Data/GeoMX_sequencing_data/sampleAnnoFile.csv`

#### Analysis

- **Raw counts matrix**: 18,677 genes × 545 samples
- **Processed counts matrix**: 18,676 genes × 537 samples
- **Feature annotation**: 18,815 genes × 10 attributes (including ProbeName, ProbeDisplayName, TargetName, HUGOSymbol, Accessions, GenomeBuild, GenomicPosition, etc.)
- **Sample annotation**: 545 samples × 37 attributes

#### Sample Metadata Structure

The sample metadata includes detailed information about each Area of Interest (AOI):

- **Slide and ROI information**: SlideName, ScanLabel, ROILabel, SegmentLabel
- **Tissue classification**: TissueType (CRSwNP, Spleen, LymphNode, Tonsil, CRSsNP), NewType (neCRSwNP, eCRSwNP, etc.)
- **Anatomical origin**: PathoType (NP, Spleen, LymphNode, Tonsil, UNC)
- **AOI metrics**: AOISurfaceArea, AOINucleiCount
- **ROI coordinates**: ROICoordinateX, ROICoordinateY
- **Sequencing metrics**: RawReads, AlignedReads, DeduplicatedReads, TrimmedReads, StitchedReads, SequencingSaturation

Each sample appears to represent a specific segment (MAC, IMM, EPI) within a tissue microarray (TMA) spot, suggesting that the data captures spatial information about different cellular compartments within the tissue.

### 3.2 GeoMX Imaging Data

#### Data Location

- Directory: `CRS_Data/GeoMX_image/`
- Number of TMA image directories: 183
- Sample TMA directories: TMA016001, TMA016002, TMA016003, etc.

#### Structure of Image Data

Each TMA directory contains a set of TIFF images representing different staining channels and segmentation information:

**Staining Channels**:

- CD45.tiff: Staining for CD45, a common leukocyte marker
- CD68.tiff: Staining for CD68, a macrophage marker
- PanCK.tiff: Pan-cytokeratin staining for epithelial cells
- SYTO13.tiff: Nuclear staining
- membrane.tiff: Membrane staining
- nuclear.tiff: Alternative nuclear staining

**Segmentation Files**:

- MESMER_mask.tiff: Cell segmentation mask generated by MESMER algorithm
- seg_outline.tiff: Segmentation outlines
- seg_overlay.tiff: Segmentation overlay on the original image

This imaging data provides the spatial context for the sequencing data, allowing for correlation between gene expression patterns and tissue morphology. The segmentation files indicate that cell-level resolution is available for spatial analysis.

## Conclusion

The CRS project data represents a comprehensive resource for studying Chronic rhinosinusitis through multiple modalities (scRNA-seq and spatial transcriptomics) and sources (public and in-house). The datasets are well-processed with detailed annotations, making them ready for advanced integrative analyses to understand the cellular and molecular basis of CRS pathology.
