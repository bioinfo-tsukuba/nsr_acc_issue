library(ggplot2)
library(dplyr)
library(patchwork)

cell_number_pop_viz <- function(seurat_obj, comp, celltype){
  
  # Aggregate cell counts by cell type
  celltype_df <- 
    seurat_obj[[]] %>%
    group_by(.data[[comp]], .data[[celltype]]) %>%
    summarise(cell_number = n(), .groups = "drop")
  
  # Bar plot of cell counts for each cell type
  p1 <- ggplot(data = celltype_df, mapping = aes(x = .data[[comp]], y = cell_number)) +
    geom_col(position = position_dodge()) +
    facet_wrap(vars(.data[[celltype]]), scales = "free_y") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Stacked bar plot of cell counts for each cell type
  p2 <- ggplot(data = celltype_df, aes(x = .data[[comp]], y = cell_number, fill = .data[[celltype]])) +
    geom_col() + 
    coord_flip()
  
  # Proportional bar plot  
  # Create a horizontal stacked bar plot showing the proportion of each cell type
  p3 <- ggplot(data = celltype_df %>%
                 group_by(.data[[comp]]) %>%
                 mutate(
                   total_cell_num = sum(cell_number),
                   ratio = cell_number / total_cell_num,
                   label = ifelse(round(ratio, 2) == 0, "", sprintf("%.2f", ratio))
                 ),
               aes(x = .data[[comp]], y = ratio, fill = .data[[celltype]])) +
    geom_col() +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5)) +
    coord_flip()
  
  return(p1 + (p2 / p3))
}



volcano_plot <- function(fmark_res, p_thr, fc_thr){
  tmp_df <- 
    fmark_res %>%
    mutate(label = case_when(
      avg_log2FC > fc_thr & p_val_adj < p_thr ~ "ident1_high",
      avg_log2FC < -fc_thr & p_val_adj < p_thr ~ "ident2_high",    
      TRUE ~ NA
    )
    )
  
  volcano_p <- ggplot(tmp_df, aes(x = avg_log2FC, y = -log10(p_val_adj), color=label)) +
    geom_point() +
    geom_hline(yintercept = -log10(p_thr), linetype = "dashed", color = "blue") +
    geom_vline(xintercept = fc_thr, linetype = "dashed", color = "red") +
    geom_vline(xintercept = -fc_thr, linetype = "dashed", color = "red") +
    theme_classic() +
    labs(x = "Average log2 Fold Change",
         y = "-log(Adjusted p-value)",
         title = "Volcano Plot of Markers")
  
  return(volcano_p)
}



