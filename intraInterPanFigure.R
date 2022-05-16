## #############################################################
## This file integrates the intra_fraction values at different 
## top MI thresholds for all cancer and normal tissues into a 
## single heatmap
## Figure 1 in article. 
################################################################

library(readr)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(stringr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

mi_data <- lapply(tissues, function(tissue) {
  bind_rows(
    read_tsv(paste0(tissue, "/distance_analysis/normal-intra-inter-count-log-bins.tsv")),
    read_tsv(paste0(tissue, "/distance_analysis/cancer-intra-inter-count-log-bins.tsv"))
  ) %>% mutate(tissue = tissue)
})
mi_data <- bind_rows(mi_data) 
mi_data$cond <- factor(mi_data$cond, levels = c("normal", "cancer"), labels = c("Normal", "Cancer"))
mi_data$tissue <- factor(mi_data$tissue, levels = tissues, labels = str_to_title(tissues))
mi_data <- mi_data %>% filter(bin <= 5e7)

mi_data_matrix <-
  mi_data %>% mutate(id = paste(tissue, cond, sep = "-")) %>% 
    pivot_wider(id_cols = bin, names_from = id, values_from = intra_fraction)

rnames <- mi_data_matrix %>% pull(bin)
mi_data_matrix <- mi_data_matrix %>% select(-bin) %>% as.matrix()
cnames <- colnames(mi_data_matrix)

col_fun_fraction = colorRamp2(c(0, 1), c("snow2","slateblue4"))
annotations <- unlist(lapply(strsplit(cnames, "-"), "[[", 2))
color_pal <-  list(condition = c("Normal" = "#e3a098", "Cancer" = "#a32e27"))
cnames <-  unlist(lapply(strsplit(cnames, "-"), "[[", 1))

mi_data_matrix <- t(mi_data_matrix)

colnames(mi_data_matrix) <- rnames
rownames(mi_data_matrix) <- cnames
cnames[seq(from = 2, to = length(cnames), by = 2)] <- ""
rnames[seq(from = 2, to = length(rnames), by = 2)] <- ""
row_ha = rowAnnotation(condition = annotations, 
                              col = color_pal,
                              annotation_legend_param =  list(title = "Condition", 
                                                              title_gp = gpar(fontsize = 32),
                                                              labels_gp = gpar(fontsize = 28)),
                              show_annotation_name = F)

hm <- Heatmap(mi_data_matrix,  cluster_rows = F,  cluster_columns = F, 
        row_labels = cnames, column_labels = rnames, right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 35),
        column_names_gp = gpar(fontsize = 25),
        col = col_fun_fraction, na_col = "white", 
        heatmap_legend_param = list(title = "Intra- fraction", 
                                    legend_height = unit(8, "cm"), 
                                    title_gp = gpar(fontsize = 32), 
                                    labels_gp = gpar(fontsize = 28), 
                                    direction = "vertical",
                                    title_position = "leftcenter-rot"))

png(paste0("pan-loss/heatmaps/intra_inter.png"), width = 1200, height = 1000)
draw(hm, padding = unit(10, "mm"), annotation_legend_side = "left")
dev.off()
