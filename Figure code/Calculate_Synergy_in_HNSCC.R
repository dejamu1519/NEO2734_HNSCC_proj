#BiocManager::install("synergyfinder")
#ainstall.packages("rgdal")
#install.packages("sf")
#直接删除原有的sf包，然后conda安装r-sf包和r-stars包

library(synergyfinder)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(data.table)
library(patchwork)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

source("/home/HD1/EP300/0_code/adj_function/Synergyfinder_Plot_adj.R")
data("mathews_screening_data")
data("ONEIL_screening_data")


combo_plate_id <- fread("/home/HD1/EP300/Drug_synergy/combo_id.csv") %>% as.data.frame()
combo_plate_set <- fread("/home/HD1/EP300/Drug_synergy/combo_set.csv") %>% as.data.frame()

input_dir <- "/home/HD1/EP300/Drug_synergy/2_splitdata/"
reshap_dir <- "/home/HD1/EP300/Drug_synergy/3_reshapedata/"
output_dir <- "/home/HD1/EP300/Drug_synergy/4_output/"

# Load the data
#panel_cell_names <- c("HN4","HN30","Cal27","HN6","HSC3","HSC4","Cal33","HOK")
panel_cell_names <- c("HN4","HN30","Cal27","Cal33","HN6" ,"Fadu")
#panel_cell_names <- c("HN6","HSC4","HSC3","Fadu")
#panel_cell_names <- c("HN4","HN30","Cal27","Cal33")
#panel_cell_names <- c("HN6","HSC4")
drug_r <- "SGC-CBP30"
drug_c <- "Dasatinib"
rep <- "re1"

block_id <- combo_plate_id$block_id[which(combo_plate_id$combo == paste0(drug_r,"_plus_", drug_c  ) )] 
print(paste0("Block ID: Plate", block_id," for ", paste0(drug_r,"_plus_", drug_c  ) ))

no_df <- data.frame() 

for (cn in panel_cell_names) {
  cell_name <- cn
  write.csv(no_df, file = paste0(input_dir,"/tmp/",block_id ,"_",drug_r,"_plus_", drug_c, "_", cell_name, "_", rep, ".csv"), row.names = FALSE)
}
################### 贴入数据准备代码，手动拉到外面 ###################

#数据准备

#single_combo_df <- read.csv(paste0(input_dir,block_id ,"_", drug_r,"_plus_", drug_c, "_", cell_name, "_", rep, ".csv") , header = TRUE,check.names = F, row.names = 1, sep = ",") %>%
#  as.matrix() %>%
#  as.data.frame() %>%
#  rownames_to_column("drug_r") %>%
# gather(key = "drug_c", value = "response", -drug_r) %>%
#  mutate(drug_r = as.numeric(drug_r),
#         drug_c = as.numeric(drug_c),
#         response = as.numeric(response))

#write.csv(single_combo_df, file = paste0(reshap_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", cell_name, "_", rep, "_Reshaped.csv"), row.names = FALSE)
combo_data_list <- list()
single_combo_list <- list()
res <- list()
res_syn <- list()
pp1 <- list()
pp1_3D <- list()
pp2 <- list()
pp2_3D <- list()
pp3 <- list()
ht <- list()
ht_list <- list()
combin_list_2D <- list()
combin_list_3D <- list()

for (i in 1:length(panel_cell_names)) {
 cell_name <-panel_cell_names[i]
 
 single_combo_df <- read.csv(paste0(input_dir,block_id ,"_", drug_r,"_plus_", drug_c, "_", cell_name, "_", rep, ".csv") ,header = F)%>%
  reshape2::melt(value.name = "response") 
 
 title_for_plot <- paste0(cell_name, " with ", drug_r, " + ", drug_c)
 
#需要为长数据 
combo_data_df <- data.frame(
  block_id = block_id,
  drug_row = drug_r, 
  drug_col = drug_c,
  conc_r = combo_plate_set$conc_r[which(combo_plate_set$block_id == block_id)],
  conc_c = combo_plate_set$conc_c[which(combo_plate_set$block_id == block_id)],
  response = single_combo_df$response,
  conc_r_unit = combo_plate_set$conc_r_unit[which(combo_plate_set$block_id == block_id)],
  conc_c_unit = combo_plate_set$conc_c_unit[which(combo_plate_set$block_id == block_id)]
)

combo_data_list[[i]] <- as.data.frame(combo_data_df)
##开始画抑制热图

single_combo_plot_df <- combo_data_df %>% reshape2::dcast(conc_r ~ conc_c, value.var = "response") %>%
  column_to_rownames("conc_r") %>%
  mutate(across(everything(), ~ if_else(. < 0, 0, .))) %>%
  as.matrix() %>% t()

draw_cell_fun <- function(j, i, x, y, width, height, fill) {   #自定义函数添加边框
  grid.rect(x = x, y = y, 
            width = unit(1.5, "cm"), 
            height = unit(1.5, "cm"),
            gp = gpar(col = "white", fill = fill, lwd = 2.5))
}

cor_fun = colorRamp2(c(0,10,20,30,40,50,60,70,80,90,100), 
                     c("navy", "#1C6FB3","#72ABEA","#BBD7FD", "#D8FEC1",
                       "#FEDAA0",
                       "#F3B79B", "#EF857D", "#ff7855", "#E92334", "darkred"))

#开始绘图
ht[[i]] <-    Heatmap(single_combo_plot_df, name = "Response", 
                      col = cor_fun, 
                      rect_gp = gpar(type = "none"),
                      cell_fun = draw_cell_fun,
                      border = "black",
                      border_gp = gpar(lwd = 3),
                      show_row_names = T, 
                      show_column_names = T,
                      cluster_rows = F,
                      cluster_columns = F,
                      column_title = paste0(drug_r, " (uM)"),
                      row_title = paste0(drug_c, " (uM)"),
                      column_title_gp = gpar(fontsize = 14,
                                             fontface= 1),
                      row_title_gp = gpar(fontsize = 14,
                                          fontface= 1),
                      column_title_side = "bottom",
                      column_names_side = "bottom",
                      column_names_centered = T,
                      column_names_rot = 0,
                      row_title_side = "left",
                      row_names_side = "left",
                      row_names_centered = T,
                      
                      row_order = nrow(single_combo_plot_df):1,
                      heatmap_legend_param = list(title = "Response", 
                                                  at = seq(0, 100, by = 20), 
                                                  labels = c("0", "20",  "40",  "60", "80",  "100"),
                                                  legend_gp = gpar(col = cor_fun(c(0,10,20,30,40,50,60,70,80,90,100))  ),
                                                  #direction = "horizontal",
                                                  #ncol = 11,
                                                  #title_position = "topcenter",
                                                  title_gp = gpar(fontsize = 12),
                                                  grid_height = unit(0.5, "cm"),
                                                  grid_width = unit(0.5, "cm"),
                                                  #tick_length = unit(0.2, "cm"),
                                                  legend_height = unit(7, "cm")
                      ),
                      
                      width = unit(ncol(single_combo_plot_df)*1.5, "cm"), 
                      height = unit(nrow(single_combo_plot_df)*1.5, "cm")) 



##########################################
#### 2. Reshaping and pre-processing   ###
##########################################
res <- ReshapeData(
  data = combo_data_list[[i]],
  data_type = "inhibition",     #或者为 "inhibition" or "viability",
  impute = F,
  impute_method = NULL,
  noise = F,
  iteration = 100,
  seed = 1)

#str(res)


#############################################
#### 3. Synergy and sensitivity analysis  ###
#############################################
res_syn <- CalculateSynergy(
  data = res,
  method = c("ZIP", "HSA", "Bliss", "Loewe"),
  Emin = NA,
  Emax = NA,
  correct_baseline = "non")

#str(res_syn[i]$synergy_scores)

#############################################
#### 3. Synergy and sensitivity analysis  ###
#############################################
pp1[[i]] <- Plot2DrugContour(
  data = res_syn,
  plot_block = block_id,
  drugs = c(1, 2),
  plot_value = "Bliss_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean"),
  text_size_scale = 1.1
)

pp1_3D[[i]] <- Plot2DrugSurface(
  data = res_syn,
  plot_block = block_id,
  drugs = c(1, 2),
  plot_value = "Bliss_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean"),
  text_size_scale = 1.1
)

#ZIP
pp2[[i]] <- Plot2DrugContour(
  data = res_syn,
  plot_block = block_id,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean"),
  text_size_scale = 1.1
)

pp2_3D[[i]] <- Plot2DrugSurface(
  data = res_syn,
  plot_block = block_id,
  drugs = c(1, 2),
  plot_value = "ZIP_synergy",
  dynamic = FALSE,
  summary_statistic = c("mean"),
  text_size_scale = 1.1
)


pp3[[i]] <- PlotMultiDrugBar_adj(
  data = res_syn,
  plot_block = block_id,
  plot_value = c("response", "ZIP_synergy", "Loewe_synergy", "HSA_synergy", "Bliss_synergy"),
  sort_by = "Bliss_synergy",
  #highlight_row = c(9.7656, 50),
  highlight_label_size = 8,
  data_table = F,
  cell_name = cell_name,
  drug_r = drug_r,
  drug_c = drug_c)

#汇总图
ht_list[[i]] <-  grid.grabExpr(draw(ht[[i]]))

combined_3D <- plot_grid(ht_list[[i]], pp1_3D[[i]], pp2_3D[[i]], nrow = 1, align = "hv", rel_widths = c(0.9, 1, 1))
combined_2D <- plot_grid(ht_list[[i]], pp1[[i]], pp2[[i]], nrow = 1, align = "hv", rel_widths = c(0.9, 1, 1))

title <- ggdraw() + draw_label(paste(title_for_plot), fontface = 'bold', x = 0.5, hjust = 0.5,size = 22)

combin_list_3D[[i]] <- plot_grid(title, combined_3D, ncol = 1, rel_heights = c(0.1, 1))
combin_list_2D[[i]] <- plot_grid(title, combined_2D, ncol = 1, rel_heights = c(0.1, 1))
}


#ht_list <- lapply(ht, function(x) grid.grabExpr(draw(x)))
#htunion <- wrap_plots(ht_list, ncol = 2)
#ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_Response_matrix.pdf"), plot = htunion, width = 14, height = 24)
#pp1union <- plot_grid(plotlist = pp1, ncol = 2)
#ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_Bliss_2D.pdf"), plot = pp1union, width = 14, height = 24)
#pp1_3Dunion <- plot_grid(plotlist = pp1_3D, ncol = 2)
#ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_Bliss_3D.pdf"), plot = pp1_3Dunion, width = 14, height = 24)

pp3union <- plot_grid(plotlist = pp3, ncol = 2)
ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_All_methods_Bar.pdf"), plot = pp3union, width = 20, height = 8*3)

combin_list_union_3D <- plot_grid(plotlist = combin_list_3D, ncol = 1, align = "hv")
ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_Two_methods_3D.pdf"), plot = combin_list_union_3D, width = 18.5, height = 5.5*length(panel_cell_names) ,limitsize = FALSE)

combin_list_union_2D <- plot_grid(plotlist = combin_list_2D, ncol = 1, align = "hv")
ggsave( filename = paste0(output_dir, block_id ,"_", drug_r,"_plus_", drug_c, "_", rep, "_Two_methods_2D.pdf"), plot = combin_list_union_2D, width = 18.5, height = 5.5*length(panel_cell_names) ,limitsize = FALSE)

