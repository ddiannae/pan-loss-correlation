## #############################################################
## This file integrates the p values from the Kolmogorov-Smirnov
## tests to evaluate the differences between intra fraction counts
## from cancer and normal tissues at different top MI thresholds.
## Data companion of Figure 1 in article. 
################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

mi_data <- lapply(tissues, function(tissue) {
    read_tsv(paste0(tissue, "/distance_analysis/intra-inter-ks-log-bins.tsv")) %>% 
      mutate(tissue = tissue)
})
mi_data <- bind_rows(mi_data) 

mi_data_matrix <- mi_data %>%
  pivot_wider(id_cols = from_bin, names_from = tissue, values_from = p.value)

bins <- mi_data_matrix %>% pull(from_bin)

t_mi_data_matrix <- t(mi_data_matrix %>% select(-from_bin)) %>%
  as_tibble(colnames = bins, .name_repair = "minimal")

colnames(t_mi_data_matrix) <- bins
pval_stats <- tibble(mean = colMeans(t_mi_data_matrix), sd = sapply(t_mi_data_matrix, sd), 
       from_bin = bins)
mi_data_matrix %>% 
  inner_join(pval_stats, by = "from_bin") %>%
  write_tsv("pan-loss/distance_data/ks_pvalues.tsv")

mi_data_matrix %>% 
  inner_join(pval_stats, by = "from_bin") %>%
  mutate(across(all_of(tissues), ~format(., digits=3))) %>%
  write_tsv("pan-loss/distance_data/ks_pvalues_rounded.tsv")

htcolors <- circlize::colorRamp2(breaks = c(0, 0.05, 1),
                                 colors = c("#528B8B", "violetred4", "lightyellow1"))

mdm <- mi_data_matrix %>% 
  select(-from_bin) %>%
  as.matrix()

rnames <- mi_data_matrix %>% 
  pull(from_bin)
  
rnames[seq(from = 2, to = length(rnames), by = 2)] <- ""

Heatmap(mdm, cluster_rows = FALSE, cluster_columns =  FALSE, show_row_names = T, 
        show_column_names = TRUE, row_labels = rnames, col = htcolors, name = "p-val", 
        heatmap_legend_param = list(legend_height = unit(8, "cm"), 
                                    title_gp = gpar(fontsize = 30), 
                                    labels_gp = gpar(fontsize = 30), 
                                    direction = "vertical"),
        column_title_gp = grid::gpar(fontsize = 30), na_col = "lightyellow1")
