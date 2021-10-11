library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)

saveHeatmap <- function(type, n) {
  barcodes <- read_tsv(paste0("pan-loss/", type, "-interaction-barcodes-", n, ".tsv"))
  tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
               "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
               "testis", "thyroid","skin", "uterus")
  
  inter_fracc <- lapply(tissues, function(tissue) {
    stats <- read_tsv(paste0(tissue, "/network_aracne/cancer-network-stats-", n,".tsv")) %>%
      filter(statistic == "inter_fraction") %>% mutate(tissue = tissue)
  })
  inter_fracc <- bind_rows(inter_fracc)
  inter_fracc_v <- inter_fracc %>% pull(value)
  names(inter_fracc_v) <- inter_fracc %>% pull(tissue)
  
  bcs <- barcodes %>% group_by(barcode) %>% tally() %>%
    filter(str_count(barcode, "1") > 1)%>%
    separate(barcode, into = tissues, sep = 1:length(tissues), convert = T) %>%
    mutate_at(all_of(tissues), ~.*n) %>% arrange(desc(n)) 
  
  bcs_75 <- bcs %>% filter(n > 75) %>%
    select(all_of(tissues)) %>% as.matrix()
  
  bcs <- bcs %>%
    select(all_of(tissues)) %>% as.matrix()
  
  substring(colnames(bcs_75), 1, 1) <- toupper(substring(colnames(bcs_75), 1, 1))
  substring(colnames(bcs), 1, 1) <- toupper(substring(colnames(bcs), 1, 1))
  substring(names(inter_fracc_v), 1, 1) <- toupper(substring(names(inter_fracc_v), 1, 1))
  
  col_fun = colorRamp2(c(0, 10,  100, 200, 400, max(bcs)), 
                       c("white", "#fff1bd", "#ffe587",
                         "#c199de", "#5f84cf", "#a32e27"))
  col_fun_fraction = colorRamp2(c(0, 1), c("white","#a32e27"))
  
  column_ha = HeatmapAnnotation(inter_fraction = inter_fracc_v, 
                                col = list(inter_fraction = col_fun_fraction), 
                                annotation_legend_param =  list(title = "Inter- fraction", 
                                                                legend_height = unit(6, "cm"), 
                                                                title_gp = gpar(fontsize = 13),
                                                                labels_gp = gpar(fontsize = 14),
                                                                direction = "horizontal",
                                                                title_position = "topleft"),
                                show_annotation_name = F)
  ## Get Heatmap for entire set
  ht <- Heatmap(bcs, cluster_rows = F, 
                col = col_fun, show_row_dend = F, use_raster = F, 
                column_names_gp =  gpar(fontsize = 16),
                heatmap_legend_param = list(title = "Shared Interactions", 
                                            legend_height = unit(6, "cm"), 
                                            title_gp = gpar(fontsize = 13),
                                            title_position = "topleft",
                                            labels_gp = gpar(fontsize = 14),
                                            direction = "horizontal"),
                bottom_annotation = column_ha)
  ht_opt("legend_gap" =unit(1, "cm"))
  hd <- draw(ht,  merge_legend = TRUE, heatmap_legend_side = "bottom", 
             annotation_legend_side = "bottom" )
  cd <- column_dend(hd)
  co <- column_order(hd)
  
  png(paste0("pan-loss/heatmaps/", type, "_heatmap_complete_", n, ".png"), width = 300, height = 1000)
  print(hd)
  dev.off()
  
  column_ha = HeatmapAnnotation(inter_fraction = inter_fracc_v, 
                                col = list(inter_fraction = col_fun_fraction), 
                                annotation_legend_param =  list(title = "Inter- fraction", 
                                                                legend_height = unit(6, "cm"), 
                                                                title_gp = gpar(fontsize = 23),
                                                                labels_gp = gpar(fontsize = 20),
                                                                direction = "vertical",
                                                                title_position = "leftcenter-rot"),
                                show_annotation_name = F)
  
  ht <- Heatmap(bcs_75, column_order =  co, use_raster = F, 
                col = col_fun, cluster_rows = F, show_column_dend = T, 
                cluster_columns  = cd, column_dend_height = unit(2.5, "cm"), 
                column_names_gp =  gpar(fontsize = 24),
                heatmap_legend_param = list(title = "Shared Interactions", 
                                            legend_height = unit(6, "cm"), 
                                            title_gp = gpar(fontsize = 23), 
                                            labels_gp = gpar(fontsize = 20), 
                                            direction = "vertical",
                                            title_position = "leftcenter-rot"),
                bottom_annotation = column_ha)
  ht_opt("legend_gap" =unit(2, "cm"))
  
  png(paste0("pan-loss/heatmaps/", type, "_heatmap_up75_", n, ".png"), width = 800, height = 1000)
  draw(ht,  merge_legend = TRUE, padding = unit(1, "cm"))
  dev.off()
}

saveHeatmap("cancer", "10000")
saveHeatmap("cancer", "100000")

saveHeatmap("cancer_only", "10000")
saveHeatmap("cancer_only", "100000")

saveHeatmap("normal", "10000")
saveHeatmap("normal", "100000")
