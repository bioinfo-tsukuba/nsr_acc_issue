# 05_comparative_analysis.R

# Cell-level DEG Analysis

# Purpose:
# Using preprocessed data, identify differentially expressed genes (DEGs) between conditions in each cell type.
# This script explores a quantitative approach utilizing Seurat's differential gene expression analysis functions.
# The goal is to detect genes that exhibit differential expression between conditions within each cell type population.

# *Differential Gene Expression (DEG) Analysis*
# A method used to identify differences in gene expression between different experimental conditions or biological states 
# (e.g., disease vs. healthy, pre- and post-drug treatment).

# Load required libraries
# code:5-1
library(Seurat)
library(tidyverse)
library(patchwork)
source("./tools/viz_tools.R") 

# Load the annotated Seurat object with cell type annotations
# code:5-2
annotated_integ_norm_seurat_obj <- readRDS("./output/04_celltype_annotation/obj/annotated_integ_norm_seurat_obj.RDS")

# Add a new metadata column "population" to represent distinct cell populations
# code:5-3
annotated_integ_norm_seurat_obj$population <- 
  paste(annotated_integ_norm_seurat_obj$condition,
        annotated_integ_norm_seurat_obj$celltype_sctype,
        sep = "_")

# Identify differentially expressed genes between wild-type and Fezf2 mutant cells within a specific cell type population.

# As an example, compare "Mature neurons" between wild-type (WT) and Fezf2 knockout (Fezf2KO).
# "WT_Mature neurons" vs. "Fezf2KO_Mature neurons"

# Set "population" as the active identity
# code:5-4
Idents(annotated_integ_norm_seurat_obj) <- "population"

# Detect differentially expressed genes (DEGs) in each population using `FindMarkers`.
# Specify two populations to compare within the metadata field set as the active identity.
# Example: If the "celltype" metadata contains {A, B, C}, we can compare population A with population B.
# Idents(seurat_obj) <- "celltype"
# FindMarkers(object = seurat_obj, ident.1 = "A", ident.2 = "B")

# code:5-5
de.markers <- FindMarkers(annotated_integ_norm_seurat_obj,
                          ident.1 = "WT_Mature neurons", # Group 1 for comparison
                          ident.2 = "Fezf2KO_Mature neurons" # Group 2 for comparison
)

# Examine the results of FindMarkers
# Each row corresponds to a gene with the following test results:
# p_val: p-value of the test
# avg_log2FC: Log2 fold-change in expression
#   avg_log2FC: log2(mean expression in ident.1 / mean expression in ident.2)
#   avg_log2FC > 0: Higher in WT, avg_log2FC < 0: Higher in Fezf2KO
# pct.1: Proportion of cells expressing the gene in Group 1
# pct.2: Proportion of cells expressing the gene in Group 2
# p_val_adj: Adjusted p-value
# code:5-6
head(de.markers)

# Visualization of results using a Volcano plot
# X-axis: avg_log2FC (Log2 fold-change in gene expression)
# Y-axis: -Log(adjusted p-value), where larger values indicate smaller p-values

# code:5-7
# Define threshold values for p_val_adj and avg_log2FC
P_THR <- 0.05 # p_val_adj threshold
FC_THR <- 2.0 # avg_log2FC threshold

# Generate the Volcano plot
# code:5-8
volcano_plot(de.markers, P_THR, FC_THR)

# Extract genes that exceed the defined thresholds
# code:5-9
de.markers %>%
  filter(p_val_adj <= P_THR) %>%
  filter(abs(avg_log2FC) > FC_THR)

# Identify differentially expressed genes between wild-type and Fezf2 mutant cells across all cell types.

# Instead of analyzing a single cell type, perform comparisons across all cell types.
# The same process is repeated for all cell types.

# code:5-10
condition_vec <- unique(annotated_integ_norm_seurat_obj$condition)
celltype_vec <- unique(annotated_integ_norm_seurat_obj$celltype_sctype)

# Generate a list of comparisons
compare_list <- lapply(as.data.frame(t(matrix(c(paste(condition_vec[1], celltype_vec, sep="_"), 
                                                paste(condition_vec[2], celltype_vec, sep="_")), ncol=2))), function(x){x})
names(compare_list) <- celltype_vec
compare_list

# Set "population" as the active identity
# code:5-11
Idents(annotated_integ_norm_seurat_obj) <- "population"

# code:5-12
# Initialize a list to store results
sc_de_res_list <- list()

# Perform DEG analysis for each cell type population by running `FindMarkers`
# The process is repeated for each cell type.
for(comp in names(compare_list)){
  print(comp)
  print(paste0(compare_list[[comp]][1], " VS ", compare_list[[comp]][2]))
  
  min_cell_num <-
    annotated_integ_norm_seurat_obj[[]] %>%
    filter(population %in% compare_list[[comp]]) %>%
    group_by(condition) %>%
    summarise(cell_number = n(), .groups = "drop") %>%
    pull(cell_number) %>%
    min()
  if(min_cell_num < 3){next} 
  
  sc_de_res_list[[comp]] <- 
    FindMarkers(annotated_integ_norm_seurat_obj,
                ident.1 = compare_list[[comp]][1],
                ident.2 = compare_list[[comp]][2]
    )
  sc_de_res_list[[comp]]$celltype <- comp
  sc_de_res_list[[comp]]$gene <- row.names(sc_de_res_list[[comp]])
}

# Compile results from all cell types into a single data frame
# code:5-13
sc_de_all_res_df <- do.call(bind_rows, sc_de_res_list)

# Examine the results
# code:5-14
head(sc_de_all_res_df)


# Summary and Visualization of Differential Expression Test Results
# Adding additional information for easier aggregation.
# code:5-15
sc_de_all_res_df_with_flg <-
  sc_de_all_res_df %>%
  mutate(test_sig_flg = ifelse(p_val_adj < P_THR, 1, 0)) %>%
  mutate(test_res_significant = case_when(
    test_sig_flg == 1 ~ "Significant",
    test_sig_flg == 0 ~ "NOT_significant",
    TRUE ~ "test_NA"
  )) %>%
  mutate(wt_high_flg = ifelse(avg_log2FC > FC_THR, 1, 0)) %>%    # positive (WT high)
  mutate(ko_high_flg = ifelse(avg_log2FC < -FC_THR, 1, 0)) %>%   # negative (Fezf2KO high)
  mutate(condition_diff = case_when(
    wt_high_flg == 1 & ko_high_flg == 0 ~ "WT_high",
    wt_high_flg == 0 & ko_high_flg == 1 ~ "Fezf2KO_high",
    TRUE ~ "No_difference"
  ))
sc_de_all_res_df_with_flg$condition_diff <- factor(sc_de_all_res_df_with_flg$condition_diff, levels = c("WT_high", "Fezf2KO_high", "No_difference"))

# Aggregation of gene counts with/without significant differential expression, considering cell types
# code:5-16
sc_de_all_res_df_with_flg %>%
  group_by(celltype, test_res_significant, condition_diff) %>%
  summarise(gene_number = n(), .groups = "drop")

# Bar plot visualization of gene counts with/without significant differential expression, considering cell types
# code:5-17
sc_de_all_res_df_with_flg %>%
  group_by(celltype, test_res_significant, condition_diff) %>%
  summarise(gene_number = n(), .groups = "drop") %>%
  filter(test_res_significant != "NOT_significant") %>% # Exclude "NOT_significant" to improve plot clarity
  filter(condition_diff != "No_difference") %>% # 
  ggplot(., aes(x = test_res_significant, y = gene_number, fill = condition_diff)) +
  geom_col(position = position_dodge2(preserve = "single")) +
  facet_wrap(. ~ celltype, scales = "free_y") +
  ggtitle("Comparison of Gene Counts Based on Significance of WT/Fezf2KO Expression Differences (Single-cell level)") +
  labs(x = "Significance Test Result", y = "Number of Genes", fill = "WT/Fezf2KO") 

# Accessing an example result
# Checking the difference between wild-type and Fezf2KO in "Microglial cells"
# code:5-18
example_df <- 
  sc_de_all_res_df %>%
  filter(celltype == "Microglial cells")

# Volcano plot visualization
# code:5-19
volcano_plot(example_df, P_THR, FC_THR)

# Identifying genes exceeding threshold
# code:5-20
example_df %>%
  filter(p_val_adj <= P_THR) %>%
  filter(abs(avg_log2FC) > FC_THR)


# Visualization of DEG expression differences using ViolinPlot

# Creating lists of DEGs for each cell type
# code:5-21
sc_deg_positive_fc_list <- list() # positive (WT high)
sc_deg_negative_fc_list <- list() # negative (Fezf2KO high)
for(ct in celltype_vec){
  sc_deg_positive_fc_list[[ct]] <-
    sc_de_all_res_df %>%
    filter(p_val_adj < P_THR) %>%
    filter(avg_log2FC > FC_THR) %>%      # positive (WT high)
    filter(celltype == ct) %>%
    pull(gene)
  
  sc_deg_negative_fc_list[[ct]] <-
    sc_de_all_res_df %>%
    filter(p_val_adj < P_THR) %>%
    filter(avg_log2FC < -FC_THR) %>%      # negative (Fezf2KO high)
    filter(celltype == ct) %>%
    pull(gene)
}

# code:5-22
# Preparing lists to store results
sc_vln_plot_positive_list <- list()
sc_vln_plot_negative_list <- list()

for(tmp_celltype in celltype_vec){ # Iterate over cell types
  # Extract cell type-specific data from the single-cell Seurat object
  tmp_obj_sc <- subset(annotated_integ_norm_seurat_obj, subset = celltype_sctype %in% c(tmp_celltype))
  
  # avgLog2FC > 0 (positive) DEG
  if(length(sc_deg_positive_fc_list[[tmp_celltype]]) !=0){
    # Violin plot visualization of single-cell level DEG expression
    sc_vln_plot_positive_list[[tmp_celltype]] <- 
      VlnPlot(object = tmp_obj_sc,
              features = sc_deg_positive_fc_list[[tmp_celltype]],
              group.by = 'celltype_sctype',
              split.by = 'condition',
              combine = F)
  }else{sc_vln_plot_positive_list[[tmp_celltype]] <- NULL} # Assign NULL if no DEGs exist for the cell type
  
  # avgLog2FC < 0 (negative) DEG
  if(length(sc_deg_negative_fc_list[[tmp_celltype]]) != 0){
    # Violin plot visualization of single-cell level DEG expression
    sc_vln_plot_negative_list[[tmp_celltype]] <- 
      VlnPlot(object = tmp_obj_sc,
              features = sc_deg_negative_fc_list[[tmp_celltype]],
              group.by = 'celltype_sctype',
              split.by = 'condition',
              combine = F)
  }else{sc_vln_plot_negative_list[[tmp_celltype]] <- NULL} # Assign NULL if no DEGs exist for the cell type
}

# Checking an example Violin Plot output
# Limiting the number of displayed plots when too many are available.
# Viewing DEG expression differences in "Immature neurons" as an example.

# Violin Plot for genes highly expressed in wild-type
# code:5-23
sc_vln_plot_positive_list[["Immature neurons"]][1:3]

# Violin Plot for genes highly expressed in Fezf2KO
# code:5-24
sc_vln_plot_negative_list[["Immature neurons"]][1:3]


################################################################################
# GO Enrichment Analysis
# This analysis aims to identify characteristic patterns and trends among the differentially expressed genes (DEGs) and understand how they relate to specific biological processes (BP), molecular functions (MF), and cellular components (CC).

# Load required libraries
# code:5-25
library(clusterProfiler)
library(org.Mm.eg.db)

# Filter DEG analysis results based on adjusted p-value and avg_log2FC.
# DEGs upregulated in WT
# code:5-26
fc_positive_pfiltered_sc_de_all_res_df <- 
  sc_de_all_res_df %>%
  filter(p_val_adj < P_THR) %>%
  filter(avg_log2FC > FC_THR) # WT DEGs

# DEGs upregulated in Fezf2KO
# code:5-27
fc_negative_pfiltered_sc_de_all_res_df <- 
  sc_de_all_res_df %>%
  filter(p_val_adj < P_THR) %>%
  filter(avg_log2FC < -FC_THR) # Fezf2KO DEGs

# Perform GO enrichment analysis using clusterProfiler.

# Prepare a vector for ontology types: Biological Process (BP), Molecular Function (MF), Cellular Component (CC)
# code:5-28
GO_type_vec <- c("CC", "BP", "MF") # Vector for iterative processing

# Create lists to store GO enrichment results
positive_go_res_list <- list()
negative_go_res_list <- list()

# Perform GO enrichment analysis for each ontology type
for (ont_type in GO_type_vec) {
  print(ont_type)
  
  # GO analysis for DEGs upregulated in WT
  print("fc positive")
  positive_go_res_list[[ont_type]] <- 
    # The `compareCluster` function allows comparison of results across different cell types
    compareCluster(gene~celltype, # Specify gene name (SYMBOL) and comparison group (celltype) from DEG results
                   data = fc_positive_pfiltered_sc_de_all_res_df, # Specify DEG results dataframe
                   fun = 'enrichGO',
                   OrgDb = org.Mm.eg.db, # Specify OrgDB
                   keyType = 'SYMBOL',
                   ont = ont_type, # Select ontology: {BP, MF, CC}
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
  
  # GO analysis for DEGs upregulated in Fezf2KO
  print("fc negative")
  negative_go_res_list[[ont_type]] <- 
    compareCluster(gene~celltype,
                   data = fc_negative_pfiltered_sc_de_all_res_df,
                   fun = 'enrichGO',
                   OrgDb = org.Mm.eg.db, # Specify OrgDB
                   keyType = 'SYMBOL',
                   ont = ont_type, # Select ontology: {BP, MF, CC}
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
}

# Visualize GO enrichment analysis results.
# Select from {BP, MF, CC}.
## X-axis: Cell type (numeric) - Number of genes used in GO analysis
## Y-axis: Gene Ontology (GO) terms
### GeneRatio: Proportion of genes belonging to a GO term among the genes used in GO analysis
### p.adjust: Adjusted p-value

# DotPlot for GO enrichment analysis of DEGs in WT
# code:5-29
dotplot(positive_go_res_list[["MF"]])
# code:5-30
dotplot(positive_go_res_list[["CC"]])
# code:5-31
dotplot(positive_go_res_list[["BP"]])

# DotPlot for GO enrichment analysis of DEGs in Fezf2KO
# code:5-32
dotplot(negative_go_res_list[["MF"]])
# code:5-34
dotplot(negative_go_res_list[["CC"]])
# code:5-35
dotplot(negative_go_res_list[["BP"]])
