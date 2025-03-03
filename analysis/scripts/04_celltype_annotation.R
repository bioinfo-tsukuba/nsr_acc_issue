# 04_celltype_annotation.R

# Cell Type Annotation

# Objective
# Assign cell type information to each cell using known marker genes and pre-annotated datasets.

# Load required libraries
# code:4-1
library(Seurat)
library(tidyverse)
source("./tools/viz_tools.R")

# Load data

# Load the Seurat object processed in 03_norm2integ
# code:4-2
integ_norm_seurat_obj <- readRDS("./output/03_norm2integ/obj/integ_norm_seurat_obj.RDS")

# code:4-3
integ_norm_seurat_obj$condition <- factor(integ_norm_seurat_obj$condition, levels = c("WT", "Fezf2KO"))

################################################################################

# 1. Assign cell types to clusters obtained through clustering

# 1.1 Cell type annotation using ScType
# ScType is used to automatically determine cell types.
# 
# ScType: A tool that performs automated cell type annotation using a list of marker genes.
#  Paper: <https://www.nature.com/articles/s41467-022-28803-w>
#  GitHub: <https://github.com/IanevskiAleksandr/sc-type>

# Load ScType R scripts
# ScType is not provided as a package but as a set of R scripts containing functions for cell type annotation.
# code:4-4
source("/home/app/sc-type/R/gene_sets_prepare.R")
source("/home/app/sc-type/R/sctype_score_.R")
library(HGNChelper)

# Specify the path to the ScType database file (list of marker genes)
# Modify `ScTypeDB_full.xlsx` if customization is needed.
# code:4-5
db_ = "/home/app/sc-type/ScTypeDB_full.xlsx"

# Specify the tissue type for the ScType database
# code:4-6
tissue = "Brain" # Other options: Immune system, Pancreas, Liver, Eye, etc.

# Prepare marker gene sets
# code:4-7
gs_list = gene_sets_prepare(db_, tissue)

# Extract data stored in the `scale.data` layer
# code:4-8
scaled_data <- GetAssayData(object = integ_norm_seurat_obj, assay = "RNA", layer = "scale.data")

# Compute cell-type enrichment scores
# code:4-9
es.max = sctype_score(scRNAseqData = scaled_data,
                      scaled = TRUE,
                      gs = gs_list$gs_positive,
                      gs2 = gs_list$gs_negative)


# Calculate ScType scores
# code:4-10
cL_results <- do.call("rbind", lapply(unique(integ_norm_seurat_obj[[]]$cca_clusters), function(cl){
  es.max.cl <- sort(rowSums(es.max[ ,rownames(integ_norm_seurat_obj[[]][integ_norm_seurat_obj[[]]$cca_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(integ_norm_seurat_obj[[]]$cca_clusters==cl)), 10)
}))

# Assign cell types to clusters based on ScType scores
# code:4-11
sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# Assign "Unknown" to clusters with low confidence (low ScType scores)
# code:4-12
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"

# Write the assigned cell types to the Seurat object metadata
# code:4-13
integ_norm_seurat_obj[[]]$celltype_sctype <- ""
for(j in unique(sctype_scores$cluster)){
  cl_type <- sctype_scores[sctype_scores$cluster==j,]
  integ_norm_seurat_obj[[]]$celltype_sctype[integ_norm_seurat_obj[[]]$cca_clusters == j] <- as.character(cl_type$type[1])
}

# Visualize the cell type annotation results using DimPlot
# code:4-14
DimPlot(
  integ_norm_seurat_obj,
  reduction = "umap.cca",
  group.by = c("cca_clusters", "celltype_sctype"),
  label = TRUE, # Display labels on the plot
  repel = TRUE  # Avoid label overlap
) + NoLegend() # Hide legend

# Store the annotated Seurat object separately
# code:4-15
annotated_integ_norm_seurat_obj <- integ_norm_seurat_obj

# Check the number of cells per cell type
# code:4-16
cell_number_pop_viz(annotated_integ_norm_seurat_obj, "condition", "celltype_sctype")


################################################################################

# 2. Assign cell types without using clustering results

# 2.1 Cell type annotation using Seurat FindTransferAnchors/DataTransfer
# Cell types are assigned by mapping the dataset to a pre-annotated scRNA-seq reference dataset.
# Seurat provides functions (FindTransferAnchors / DataTransfer) for this annotation approach.
# This method assumes that cells in the reference dataset are sufficiently similar to those in the target dataset.
# The annotation results depend on the reference dataset's cell types.
# If gene names in the reference dataset and target dataset are completely different, this method cannot be used.

# For demonstration, the reference dataset is the annotated wild-type P1 stage dataset from the original study.

# code:4-17
# Load the pre-annotated wild-type P1 dataset.
ref_seurat_obj <- readRDS("../data/reference/P1_ref_seurat_obj.RDS")

# Perform NormalizeData and FindVariableFeatures
# code:4-18
ref_seurat_obj <- NormalizeData(ref_seurat_obj)
ref_seurat_obj <- FindVariableFeatures(ref_seurat_obj)

# Find anchors (identify similar cells between query and reference datasets)
# query: The Seurat object to be annotated
# reference: The annotated wild-type P1 dataset
# code:4-19
anchors <- FindTransferAnchors(
  reference = ref_seurat_obj,               # Reference dataset
  query = annotated_integ_norm_seurat_obj,  # Target dataset
  normalization.method = "LogNormalize"     # Specify normalization method
)

# Transfer labels
# Predict cell types in the target dataset based on reference dataset annotations
# code:4-20
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_seurat_obj$New_cellType # Column name containing cell type annotations in reference data
)

# Add predicted cell types to the target dataset metadata
# code:4-21
annotated_integ_norm_seurat_obj <-
  AddMetaData(object = annotated_integ_norm_seurat_obj,
              metadata = predictions["predicted.id"])

# Visualize cell type annotation results using DimPlot
# code:4-22
DimPlot(
  annotated_integ_norm_seurat_obj,
  reduction = "umap.cca",
  group.by = c("cca_clusters", "predicted.id"),
)

# Check the number of cells per predicted cell type
# code:4-23
cell_number_pop_viz(annotated_integ_norm_seurat_obj, "condition", "predicted.id")


################################################################################

# #### 2.2 Cell Type Annotation Using Azimuth
# When a suitable reference dataset is available within Azimuth, annotation can be performed easily.
# The internal algorithm is the same as in Section 2.1.

# Load required libraries
# code:4-24
library(Azimuth)
library(SeuratData)

# Map the query dataset to Azimuth's reference data using the `RunAzimuth` function
# Azimuth mapping
# Mapping to Mouse Motor Cortex reference dataset
# code:4-25
tmp_seurat_obj <- RunAzimuth(query = annotated_integ_norm_seurat_obj, reference = "mousecortexref")

# Add the predicted cell type annotations to the metadata
# code:4-26
annotated_integ_norm_seurat_obj <- 
  AddMetaData(object = annotated_integ_norm_seurat_obj,
              metadata = tmp_seurat_obj[[]][c("predicted.class", "predicted.subclass")]
  )

# Visualize the cell type annotation results using DimPlot
# code:4-27
DimPlot(
  annotated_integ_norm_seurat_obj,
  reduction = "umap.cca",
  group.by = c("cca_clusters", "predicted.class", "predicted.subclass"),
)

# Check the number of cells in each predicted class
# code:4-28
cell_number_pop_viz(annotated_integ_norm_seurat_obj, "condition", "predicted.class")

# Check the number of cells in each predicted subclass
# code:4-29
cell_number_pop_viz(annotated_integ_norm_seurat_obj, "condition", "predicted.subclass")

################################################################################

# Comparison of ScType (Marker Gene-Based Annotation), Azimuth (Mouse Data Reference), and DataTransfer (Pre-prepared Reference Data)
# Visualize the correspondence between cell type annotation methods using a River Plot
# code:4-30

# Load required library
library(ggalluvial)

# Generate the River Plot
# code:4-31
ggplot(data = annotated_integ_norm_seurat_obj[[]],
       aes(axis1 = predicted.id,
           axis2 = celltype_sctype,
           axis3 = predicted.class,
           axis4 = predicted.subclass,
       )) +
  scale_x_discrete(limits = c("celltype_transfer\nfrom_reference",
                              "celltype_sctype",
                              "Azimuth_class",
                              "Azimuth_subclass"
  ),
  expand = c(.2, .05)) +
  xlab("Cell Type Annotation Methods") +
  geom_alluvium() +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Correspondence of Cell Types Across Annotation Methods")

# Save the annotated Seurat object
# code:4-32
saveRDS(annotated_integ_norm_seurat_obj, "./output/04_celltype_annotation/obj/annotated_integ_norm_seurat_obj.RDS")
