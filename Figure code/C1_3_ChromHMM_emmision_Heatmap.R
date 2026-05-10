library(TheBestColors)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
source("/home/HD1/EP300/0_code/C1_2_ChromHMM_anno.R")

output_dir <- "/home/HD1/EP300/C_ChromHMM/1_discovery_output/plot"
ChromHMM_15 <- read.table("/home/HD1/EP300/C_ChromHMM/1_discovery_output/15_model_without_BRD4_EP300/emissions_15.txt", header = T, sep = "\t", row.names = 1)

rownames(ChromHMM_15) <- state_reorder

ChromHMM_15 <- ChromHMM_15 %>%
  mutate(anno = state_anno) %>%
  arrange(rownames(.))

#开始绘图  
# 设置颜色
colors <- colorRampPalette(Best100(48) )(256)

# 创建热图
data_matrix <- as.matrix(ChromHMM_15[,c(6,3,4,5,2,7,1)])
 
plot_emission <- Heatmap(data_matrix, 
                         
        col = colors,
        row_title = "Chromatin States",
        row_title_side = "left",
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_labels = rownames(data_matrix),
        row_names_side = "left",
        column_labels = colnames(data_matrix),
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        show_heatmap_legend = F)


# 添加注释
row_annotation <- rowAnnotation(
  text = anno_text(ChromHMM_15$anno), 
                   gp = gpar(col = "black"))

draw(plot_emission+row_annotation) 

