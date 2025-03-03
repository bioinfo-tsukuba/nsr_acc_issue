# 必要なライブラリをロード
library(Seurat)
library(ggplot2)

custom_dotplot <- function(obj, features, group1, group2){
  exp_assay_data <- GetAssayData(object = obj, assay = "RNA", layer = "data")
  count_assay_data <- GetAssayData(object = obj, assay = "RNA", layer = "counts")
  exp_mat <- exp_assay_data[features, , drop = FALSE]
  count_mat <- count_assay_data[features, , drop = FALSE]
  meta <- obj@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  
  # get the average expression       
  exp_df <- as.matrix(exp_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "expression") %>%
    left_join(meta) %>%
    group_by(gene, !!sym(group1), !!sym(group2)) %>%
    summarise(average_expression = mean(expression)) 
  
  # get percentage of positive cells
  percent_df <- as.matrix(count_mat) %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column(var="gene") %>%
    tidyr::pivot_longer(!gene, names_to = "cell", values_to = "count") %>%
    left_join(meta) %>%
    group_by(gene, !!sym(group1), !!sym(group2)) %>%
    summarise(percentage = mean(count > 0)) 
  
  summarised_df <-
    exp_df %>%
    full_join(percent_df, by=c("gene", group1, group2))
  
  # プロット作成
  ggplot(summarised_df, aes(x = !!sym(group2), y = !!sym(group1))) +
    geom_point(aes(size = percentage, color = average_expression), alpha = 0.7) +
    scale_size_continuous(name = "Size (percentage)") + # サイズの凡例
    scale_color_gradient(name = "Color (average_expression)", low = "gray", high = "red") + # 色のグラデーション
    facet_grid(.~gene, scales = "free", space = "free") +
    labs(
      title = "Average gene expression",
    ) +
    theme(
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1), # X軸ラベルを回転
      panel.grid.major = element_blank(), # 主な背景線を消す
      panel.grid.minor = element_blank(), # 補助的な背景線を消す
      axis.line = element_line(color = "black"), # 軸線を表示
      axis.ticks = element_line(color = "black") # 目盛りを表示
    )
}
