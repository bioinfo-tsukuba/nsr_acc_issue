# 03_data_preprocessing_and_integration.R

# Seurat Basic Analysis Workflow

# Objective
# Preprocess quality-controlled scRNA-seq data to enable downstream analyses such as comparative analysis.

# Preprocessing Steps
#   - Normalize sequencing depth differences among cells
#   - Reduce dimensionality of high-dimensional data
#   - Integrate datasets measured under different conditions
#   - Classify cells into similar groups

# Seurat Basic Analysis Workflow
#   1. Data normalization
#   2. Identification of highly variable genes (HVG)
#   3. Data standardization
#   4. Principal Component Analysis (PCA)
#   5. Data integration
#   6. k-Nearest Neighbors (kNN) estimation
#   7. Clustering
#   8. UMAP visualization

# Load Required Libraries
# code:3-1
library(Seurat)
library(ggplot2)

# Load the Seurat object generated in the previous step.
# code:3-2
merged4norm_seurat_obj <- readRDS("./output/02_qc/obj/merged4norm_seurat_obj.RDS")

# Seurat Basic Analysis Workflow
# Seurat V5 enables simultaneous processing of multiple data layers stored in a single Seurat object.

# 1. Data Normalization
# The total number of reads varies between cells, requiring transformation into comparable values.
# By default, Seurat normalizes data by dividing by the total reads (nCounts) per cell,
# scales all detected reads across cells to a constant value (scale.factor: 10,000),
# and applies log transformation.
# code:3-3
norm_seurat_obj <- NormalizeData(merged4norm_seurat_obj)

# 2. Identification of Highly Variable Genes (HVG)
# Highly variable genes characterize cell properties, whereas low-variance genes
# (e.g., housekeeping genes) are consistently expressed across all cells.
# scRNA-seq analysis focuses on genes that define cellular differences.
# By default, Seurat selects 2,000 highly variable genes (common across layers).
# code:3-4
norm_seurat_obj <- FindVariableFeatures(norm_seurat_obj)

# Visualize the distribution of highly variable genes.
# code:3-5
VariableFeaturePlot(norm_seurat_obj)

# 3. Data Standardization
# Preparation for PCA
# Standardization enables gene-wise comparisons.
# By default, only HVGs are standardized to reduce computational cost.
# Data transformation to zero mean and unit variance per gene.
# code:3-6
norm_seurat_obj <- ScaleData(norm_seurat_obj)

# 4. Principal Component Analysis (PCA)
# Dimension reduction of high-dimensional scRNA-seq data (genes × cells).
# PCA is commonly used to extract key features while reducing dimensionality.
# Seurat adopts PCA as a standard approach.
# code:3-7
norm_seurat_obj <- RunPCA(norm_seurat_obj)

# 5. Data Integration
# Integrating datasets across experimental batches, donors, or conditions is a crucial step in scRNA-seq workflows.
# Integration enables matching shared cell types/states across datasets,
# facilitating cross-condition comparisons.
#
# Seurat V5 introduces `IntegrateLayers` for data integration.
# This function integrates multiple layers stored in a single Seurat object.
# Various integration methods, including CCA, are supported.

## Data Integration via Canonical Correlation Analysis (CCA)
# code:3-8
integ_norm_seurat_obj <- IntegrateLayers(
  object = norm_seurat_obj,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)

# 6. k-Nearest Neighbors (kNN) Estimation
# Preparation for the next step: clustering.
# `reduction`: Specifies the dimension reduction method used during integration.
# code:3-9
integ_norm_seurat_obj <- FindNeighbors(integ_norm_seurat_obj, reduction = "integrated.cca", dims = 1:30)

# 7. Clustering
# Groups cells with similar gene expression profiles.
# `resolution`: Parameter controlling cluster granularity (default: 0.8).
# Higher `resolution` → More clusters
# Lower `resolution` → Fewer clusters
# code:3-10
integ_norm_seurat_obj <- FindClusters(integ_norm_seurat_obj,
                                      cluster.name = "cca_clusters" # Metadata field for storing cluster labels
)

# 8. UMAP Visualization
# Reduces high-dimensional gene expression data to lower dimensions.
# Unlike PCA, UMAP is a nonlinear dimensionality reduction method used to visualize data structure and patterns.
# code:3-11
integ_norm_seurat_obj <- RunUMAP(integ_norm_seurat_obj, 
                                 reduction = "integrated.cca", # Specifies the reduction method used during integration
                                 dims = 1:30,                  # Number of dimensions to use from integration
                                 reduction.name = "umap.cca"   # Assigns "umap.cca" as the reduction name
)

# Visualize UMAP results.
# code:3-12
DimPlot(
  integ_norm_seurat_obj,
  reduction = "umap.cca",  # Specifies the dimensionality reduction result from RunUMAP
  group.by = c("condition", "cca_clusters") # Specifies metadata fields for color-coding
)

# Merge Layers into a Single Layer
# code:3-13
integ_norm_seurat_obj <- JoinLayers(integ_norm_seurat_obj)

# Save the integrated and clustered Seurat object.
# code:3-14
saveRDS(integ_norm_seurat_obj, "./output/03_norm2integ/obj/integ_norm_seurat_obj.RDS")
