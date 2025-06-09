library(vroom)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)

tissues <- c("brain","colon", "esophagus", 
             "liver","ovary", "pancreas", "prostate", 
             "testis", "skin")

subids <- lapply(tissues, function(tissue) {
  file <- dir(path = paste0(tissue, "/results"), pattern = "*_*_*_si-arsyn_normal.tsv", full.names = TRUE)
  samples <- vroom(paste0(tissue, "/results/", tissue, "-normal-samples.tsv"))
  in_matrix <- vroom(file = file, n_max = 1, col_names = T)
  in_matrix <- colnames(in_matrix)[-1]
  samples %>% 
    filter(id %in% in_matrix) %>%
    mutate(tissue = tissue,
           subid = map_chr(file, .f = ~ paste(unlist(strsplit(.x, "-"))[1:2], collapse = "-"))) %>%
    select(subid, tissue)
  
})

subids <- bind_rows(subids)
grouped_subids <- subids %>% 
  group_by(subid, tissue) %>%
  summarize(tissue = tissue,
            n = n()) %>% ungroup() %>%
  distinct() %>%
  pivot_wider(id_cols = subid, names_from = tissue, values_from = n, values_fill = 0) %>%
  arrange(esophagus, skin) %>%
  rowwise() %>%
  mutate(nts = sum(across(all_of(tissues))))

grouped_subids %>% 
  filter(nts == 1) %>%
  select(all_of(tissues)) %>%
  colSums()
 
mgs <- grouped_subids %>%
  select(all_of(tissues)) %>%
  as.matrix() 

Heatmap(mgs)
  
subids %>% 
  group_by(subid) %>%
  tally() %>%
  ggplot(aes(x = n)) +
  geom_histogram(bins = 8)

grouped_subids %>%
  filter(colon == 1 & liver == 1)

grouped_subids %>%
  filter(testis == 1 & prostate == 1)
  