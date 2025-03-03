library(Seurat)
library(tidyverse)

# metadata (WT only)
metadata_df <- read_tsv("../data/SCP1290/metadata/metaData_scDevSC.txt")

NAMEs_vec <- metadata_df %>% pull(NAME)
metadata_df <- as.data.frame(metadata_df %>% dplyr::select(-NAME))
rownames(metadata_df) <- NAMEs_vec

# scRNA-seq data
path_vec <- c("../data/all/GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5",
              "../data/all/GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5")
names(path_vec) <- c("P1_S1", "P1_S2")

orig_mtx_list <- list()
for(smn in names(path_vec)){
  print(smn)
  orig_mtx_list[[smn]] <- Read10X_h5(path_vec[smn])
  colnames(orig_mtx_list[[smn]]) <- paste(smn, str_replace_all(colnames(orig_mtx_list[[smn]]), "-1", ""), sep = "_")
}

# Create Seurat objects
created_seurat_obj_list <- list()
for(snm in names(orig_mtx_list)){
  print(snm)
  created_seurat_obj_list[[snm]] <- CreateSeuratObject(orig_mtx_list[[snm]], min.cells = 3, project = snm)
  created_seurat_obj_list[[snm]] <- AddMetaData(created_seurat_obj_list[[snm]], metadata_df)
  }

merged4norm_seurat_obj <- created_seurat_obj_list[[1]]
for (snm in names(created_seurat_obj_list)[2:length(created_seurat_obj_list)]) {
  merged4norm_seurat_obj <- merge(merged4norm_seurat_obj, created_seurat_obj_list[[snm]])
}

ref_seurat_obj <- JoinLayers(merged4norm_seurat_obj)

dir.create("../data/reference/")
saveRDS(ref_seurat_obj, "../data/reference/P1_ref_seurat_obj.RDS")
