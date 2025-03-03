# 02_quality_control.R

# Quality Control

# Objective
# Perform Quality Control before executing the basic Seurat analysis workflow to remove low-quality cells.

# Quality Control Steps
#   1. Automated Quality Control using ddqc
#   2. Removal of doublet cells using DoubletFinder

# Load required libraries
# code:2-1
library(Seurat)
library(tidyverse)       # Includes dplyr, ggplot2, etc.
library(ddqcR)           # Automated Quality Control tool
library(DoubletFinder)   # Doublet removal tool

# Load objects required for Quality Control
# code:2-2
# Seurat object with added metadata from the previous step
created_seurat_obj_list <- readRDS("./output/01_load/obj/created_seurat_obj_list.RDS")

# Original expression matrix data without any processing
orig_mtx_list <- readRDS("./output/01_load/obj/orig_mtx_list.RDS")

# Calculate the proportion of mitochondrial and ribosomal genes
# Perform this calculation before executing Quality Control
# code:2-3
for(smn in names(created_seurat_obj_list)) {
  # Mitochondrial gene percentage
  created_seurat_obj_list[[smn]][["bfqc_percent.mt"]] <- 
    PercentageFeatureSet(created_seurat_obj_list[[smn]], 
                         pattern = "^mt-"  # Pattern differs for human ("^MT-") and mouse ("^mt-")
    )
  # Ribosomal gene percentage
  created_seurat_obj_list[[smn]][["bfqc_percent.rp"]] <- 
    PercentageFeatureSet(created_seurat_obj_list[[smn]],
                         pattern = "(?i)^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA"
    )
}

# 1. Automated Quality Control using ddqc
# code:2-4
ddqc_seurat_obj_list <- list()
ddqc_df_list <- list()

# Execute ddqc
for(smn in names(created_seurat_obj_list)){
  print(paste0("ddqc: ", smn))
  tmpdata <- initialQC(created_seurat_obj_list[[smn]])
  ddqc_df_list[[smn]] <- ddqc.metrics(tmpdata)
  ddqc_seurat_obj_list[[smn]] <- filterData(tmpdata, ddqc_df_list[[smn]])
}

# Check metadata: mitochondrial and ribosomal gene percentages should be added
# code:2-5
head(ddqc_seurat_obj_list[[1]][[]])

# Save Seurat object after ddqc execution
# code:2-6
saveRDS(ddqc_seurat_obj_list, "./output/02_qc/obj/ddqc_seurat_obj_list.RDS")
saveRDS(ddqc_df_list, "./output/02_qc/obj/ddqc_df_list.RDS")

# 2. Removal of doublet cells using DoubletFinder
# Preprocessing required for DoubletFinder
# DoubletFinder requires the basic Seurat analysis workflow beforehand
# code:2-7
for_doubletfinder_seurat_obj_list <- lapply(ddqc_seurat_obj_list, NormalizeData)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, FindVariableFeatures)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, ScaleData)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, RunPCA)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, RunUMAP, dims = 1:10)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, FindNeighbors)
for_doubletfinder_seurat_obj_list <- lapply(for_doubletfinder_seurat_obj_list, FindClusters)

# Detect doublets using DoubletFinder
# code:2-8
set.seed(1)
doublet_finder_res_df_list <- list()

for(nm in names(for_doubletfinder_seurat_obj_list)){
  print(paste0("DoubletFinder: ", nm))
  tmp_seurat_obj <- for_doubletfinder_seurat_obj_list[[nm]]
  
  ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
  sweep.res.list_pmbc <- paramSweep(tmp_seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats_pmbc <- summarizeSweep(sweep.res.list_pmbc, GT = FALSE)
  bcmvn_pmbc <- find.pK(sweep.stats_pmbc)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  homotypic.prop <- modelHomotypic(tmp_seurat_obj[[]]$seurat_clusters)           
  nExp_poi <- round(0.075*nrow(tmp_seurat_obj[[]]))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  pN <- 0.25
  pK <- bcmvn_pmbc %>% arrange(desc(BCmetric)) %>% pull(pK) %>% as.character %>% as.numeric %>% .[1] 
  
  base_char <- paste("DF.classifications", pN, pK, nExp_poi, sep = "_")
  adj_char  <- paste("DF.classifications", pN, pK, nExp_poi.adj, sep = "_")
  char_reuse_pANN <- paste("pANN", pN, pK, nExp_poi, sep = "_")
  
  tmp_seurat_obj <- doubletFinder(tmp_seurat_obj, PCs = 1:10, pN = pN, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  tmp_seurat_obj <- doubletFinder(tmp_seurat_obj, PCs = 1:10, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = char_reuse_pANN, sct = FALSE)
  tmp_doublet_finder_res_df <- tmp_seurat_obj[[]][c(base_char, adj_char)]
  colnames(tmp_doublet_finder_res_df) <- c("double_flg", "doublet_adj_flg")
  doublet_finder_res_df_list[[nm]] <- tmp_doublet_finder_res_df
}

# Remove doublet cells identified by DoubletFinder
# code:2-9
doubletfinder_ddqc_seurat_obj_list <- list()
for (nm in names(ddqc_seurat_obj_list)) {
  doubletfinder_ddqc_seurat_obj_list[[nm]] <- 
    AddMetaData(
      ddqc_seurat_obj_list[[nm]],
      doublet_finder_res_df_list[[nm]]
    )
  doubletfinder_ddqc_seurat_obj_list[[nm]] <- subset(doubletfinder_ddqc_seurat_obj_list[[nm]], subset = doublet_adj_flg %in% c("Singlet"))
}

# Save Seurat objects after doublet removal
# code:2-10
saveRDS(doubletfinder_ddqc_seurat_obj_list, "./output/02_qc/obj/doubletfinder_ddqc_seurat_obj_list.RDS")

################################################################################

# Aggregation Before and After Quality Control
# Collect metadata from each Seurat object for summary.
# code:2-11
get_metadata <- function(sample_name, obj) {
  return(obj[[sample_name]][[]])
}
before_qc_metadata_df <- do.call(bind_rows, lapply(names(created_seurat_obj_list), get_metadata, obj = created_seurat_obj_list))
colnames(before_qc_metadata_df)[which(colnames(before_qc_metadata_df) %in% c("bfqc_percent.mt", "bfqc_percent.rp"))] <- c("percent.mt", "percent.rb")
after_qc_metadata_df  <- do.call(bind_rows, lapply(names(ddqc_seurat_obj_list), get_metadata, obj = ddqc_seurat_obj_list))
after_qc_doubletf_metadata_df <- do.call(bind_rows, lapply(names(doubletfinder_ddqc_seurat_obj_list), get_metadata, obj = doubletfinder_ddqc_seurat_obj_list))

# Visualizing distribution changes of "nCount_RNA", "nFeature_RNA", "percent.mt", and "percent.rb" before and after Quality Control using ViolinPlot.

# Gather required information for ViolinPlot.
# code:2-12
plot_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.rb")
# Combine metadata information.
qc_result_metadata_df <- bind_rows(
  before_qc_metadata_df %>%
    dplyr::select(condition, {{ plot_vars }}) %>%
    mutate(qc="before_qc"),
  after_qc_metadata_df %>%
    dplyr::select(condition, {{ plot_vars }}) %>%
    mutate(qc="after_ddqc"),
  after_qc_doubletf_metadata_df %>%
    dplyr::select(condition, {{ plot_vars }}) %>%
    mutate(qc="after_ddqc_doubletfinder")
)
# Convert the 'qc' column into a factor and set category order: "before_qc" → "after_ddqc" → "after_ddqc_doubletfinder"
qc_result_metadata_df$qc <- factor(qc_result_metadata_df$qc, levels = c("before_qc", "after_ddqc", "after_ddqc_doubletfinder"))
qc_result_metadata_df$condition <- factor(qc_result_metadata_df$condition, levels = c("WT", "Fezf2KO"))

# Generate ViolinPlot
# code:2-13
violinplot_list <- list()
for(var in plot_vars) {
  print(var)
  target_var <- sym(var)  # Convert string to symbol
  violinplot_list[[var]] <- 
    ggplot(data = qc_result_metadata_df, mapping = aes(x = condition, y = {{ target_var }}, fill = qc)) +
    geom_violin(position = position_dodge(width = 0.9), width = 1.0, draw_quantiles = c(0.25, 0.5, 0.75)) + # Violin plot
    ggtitle(paste("QC result: ", var)) + # Add title
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Display ViolinPlots
# code:2-14
violinplot_list

# Visualizing changes in the number of genes and cells before and after Quality Control using BarPlot.
# Extract gene/cell count information from Seurat objects.
# code:2-15
data_vec <- c("rawdata", "before_qc", "after_ddqc", "after_ddqc_doubletfinder")
# Collect gene/cell count information from Seurat objects.
gene_cell_pop_list <- list(
  as.data.frame(t(as.data.frame(lapply(orig_mtx_list, dim)))),
  as.data.frame(t(as.data.frame(lapply(created_seurat_obj_list, dim)))),
  as.data.frame(t(as.data.frame(lapply(ddqc_seurat_obj_list, dim)))),
  as.data.frame(t(as.data.frame(lapply(doubletfinder_ddqc_seurat_obj_list, dim))))
)

# Prepare data for BarPlot visualization.
# code:2-16
names(gene_cell_pop_list) <- data_vec
for(dtn in data_vec) {
  colnames(gene_cell_pop_list[[dtn]]) <- c("gene_number", "cell_number")
  gene_cell_pop_list[[dtn]] <-
    gene_cell_pop_list[[dtn]] %>%
    mutate(data_name = dtn)
  gene_cell_pop_list[[dtn]]["sample"] <- rownames(gene_cell_pop_list[[dtn]])
  row.names(gene_cell_pop_list[[dtn]]) <- NULL
}
gene_cell_pop_df <- 
  do.call(bind_rows, gene_cell_pop_list) %>%
  mutate(data_name = factor(data_name, levels = c("rawdata", "before_qc", "after_ddqc", "after_ddqc_doubletfinder"))) %>%
  separate(col = sample, into = c("age", "area", "condition"), sep = "_") %>%
  mutate(condition = factor(condition, levels = c("WT", "Fezf2KO")))

# Generate BarPlots for gene/cell counts.
# code:2-17
plot_vars <- c("gene_number", "cell_number")
barplot_list <- list()
for(var in plot_vars) {
  target_var <- sym(var)  # Convert string to symbol
  barplot_list[[var]] <-
    ggplot(data = gene_cell_pop_df, mapping = aes(x = condition, y = {{ target_var }}, fill = data_name)) +
    geom_col(position = position_dodge()) +
    ggtitle(paste("QC result: ", var)) + # Add title
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Display BarPlots
# code:2-18
barplot_list

# Merge Seurat objects
# Previously, each scRNA-seq dataset was handled separately, but now they are merged into a single Seurat object.
# Each scRNA-seq dataset is stored in Layers.
# code:2-19
merged4norm_seurat_obj <- doubletfinder_ddqc_seurat_obj_list[[1]]
for (sample_name in names(doubletfinder_ddqc_seurat_obj_list)[2:length(doubletfinder_ddqc_seurat_obj_list)]) {
  merged4norm_seurat_obj <- merge(merged4norm_seurat_obj, doubletfinder_ddqc_seurat_obj_list[[sample_name]])
}

# Save the merged Seurat object
# code:2-20
saveRDS(merged4norm_seurat_obj, "./output/02_qc/obj/merged4norm_seurat_obj.RDS")
