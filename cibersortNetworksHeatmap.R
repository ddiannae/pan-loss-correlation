library(readr)
library(tidyr)
library(dplyr)
library(stringr)
library(circlize)
library(ComplexHeatmap)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","uterus")

t <- "bladder"
cell <- "CD31"
ninter <- "100000"
cond <- "normal"
cell_types <- c("CD10", "CD31", "CD45", "EPCAM")
intras <- lapply(cell_types, function(cell) {
  intras <- lapply(c("normal", "cancer"), function(cond) {
    intras <- lapply(tissues, function(t) {
        interactions <- read_tsv(paste0(t, "/cibersort_networks_tr4/", cond, "_", cell, "_", ninter, "_interactions.tsv"))
        inter <- interactions %>% 
          filter(interaction_type == "Inter") %>%
          nrow()
        intra <-  interactions %>% 
          filter(interaction_type == "Intra") %>%
          nrow() 
        tibble(tissue = t, cond = cond, intra = (intra/(intra+inter)), cell = cell)
    })
    bind_rows(intras)
  })
  bind_rows(intras)
})
intras <- bind_rows(intras)
intras <- intras %>%
  mutate(name = paste0(toupper(substr(tissue,1, 1)), substr(tissue,2, nchar(tissue)), 
                       " ", toupper(substr(cond,1, 1)), substr(cond,2, nchar(cond))))

intras <- intras %>% 
  pivot_wider(id_cols = c(name, cond), names_from = cell, values_from = intra)
conds <- intras %>% pull(cond, name = name)
substr(conds,1, 1) <- substr(toupper(conds),1, 1)

intras %>% 
  write_tsv("pan-loss/immune/intra_fractions.tsv")

mintras <- intras %>%
  select(all_of(cell_types)) %>% 
  as.matrix()
rownames(mintras) <- intras$name

color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")
names(color_pal) <- labels
col_fun = colorRamp2(c(0, 1), c("white", "slateblue4"))

ht <- Heatmap(mintras, cluster_rows = TRUE, cluster_columns =  TRUE, show_row_names = TRUE, 
              show_column_names = TRUE, name = "Score", col = col_fun,
              column_labels = cell_types, 
              column_names_rot = 0,
              column_names_gp =  gpar(fontsize = 20),
              row_names_gp =  gpar(fontsize = 20),
              column_names_centered = T,
              left_annotation = rowAnnotation(Condition = conds, 
                                              col = list(Condition = color_pal), 
                                              show_annotation_name = F, 
                                              annotation_legend_param = list(title_gp = gpar(fontsize = 20), 
                                                                             labels_gp = gpar(fontsize = 18))), 
              heatmap_legend_param = list(legend_height = unit(3, "cm"), 
                                          title_gp = gpar(fontsize = 20), 
                                          labels_gp = gpar(fontsize = 18), 
                                          title_position = "leftcenter-rot")) 

png("pan-loss/immune/cibersort_networks.png", width = 700 , height = 1000)
draw(ht,   padding = unit(0.5, "cm"))
dev.off()
                         