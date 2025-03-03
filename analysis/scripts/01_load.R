# 01_load.R

# Loading scRNA-seq data

# Objective
# - Load expression matrix data and convert it into a format compatible with Seurat
# - Add cell metadata to the Seurat object

# Load required libraries
# code:1-1
library(Seurat)
library(tidyverse)

# Load data

# Read metadata file containing paths and associated information for each dataset
# code:1-2
# Path to the text file containing metadata and expression matrix paths
data_sheet_path <- "./data_info/datasheet.tsv"

# code:1-3
# Read the metadata file into a dataframe
data_sheet_df <- read_tsv(data_sheet_path)

# Inspect the contents
# The dataframe contains metadata and file paths
# code:1-4
data_sheet_df

# Extract metadata associated with expression data as vectors
# code:1-5
age_vec <- data_sheet_df %>% pull(age)
area_vec <- data_sheet_df %>% pull(area)
condition_vec <- data_sheet_df %>% pull(condition)
sample_vec <- paste0(age_vec, "_", area_vec, "_", condition_vec)

names(age_vec) <- sample_vec
names(area_vec) <- sample_vec
names(condition_vec) <- sample_vec
names(sample_vec) <- sample_vec

# Extract expression matrix file paths as a vector
# code:1-6
path_vec <- data_sheet_df %>% pull(path)
names(path_vec) <- sample_vec

# Load expression matrix data
# code:1-7
orig_mtx_list <- list()
for(smn in names(path_vec)){
  print(smn)
  orig_mtx_list[[smn]] <- Read10X_h5(path_vec[smn])
}

# Save the loaded expression matrices
# code:1-8
saveRDS(orig_mtx_list, "./output/01_load/obj/orig_mtx_list.RDS")

# Create Seurat objects
# Convert the loaded expression matrices into Seurat objects
# `min.cells = 3`: Remove genes expressed in 3 or fewer cells
# code:1-9
# Initialize an empty list to store Seurat objects
created_seurat_obj_list <- list()
for(snm in names(orig_mtx_list)){ # Iterate over multiple files
  print(snm)
  created_seurat_obj_list[[snm]] <- CreateSeuratObject(orig_mtx_list[[snm]], min.cells = 3, project = snm)
}

# Inspect the created Seurat objects
# code:1-10
created_seurat_obj_list

# Check metadata in the Seurat objects
# code:1-11
head(created_seurat_obj_list[[1]][[]])
head(created_seurat_obj_list[[2]][[]])

# Add metadata to the Seurat objects
# code:1-12
for(smn in names(created_seurat_obj_list)){
  print(smn)
  created_seurat_obj_list[[smn]][[]]["age"]         <- age_vec[smn]
  created_seurat_obj_list[[smn]][[]]["area"]        <- area_vec[smn]
  created_seurat_obj_list[[smn]][[]]["condition"]   <- condition_vec[smn]
  created_seurat_obj_list[[smn]][[]]["sample_name"] <- sample_vec[smn]
}

# Verify that metadata has been successfully added
# code:1-13
head(created_seurat_obj_list[[1]][[]])
head(created_seurat_obj_list[[2]][[]])

# Save the Seurat objects before quality control
# code:1-14
saveRDS(created_seurat_obj_list, "./output/01_load/obj/created_seurat_obj_list.RDS")
