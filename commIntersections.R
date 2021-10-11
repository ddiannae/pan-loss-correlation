library(readr)
library(dplyr)
library(tidyr)
library(igraph)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")

jacc_idxs <- read_tsv("pan-loss/network_intersections/all_jacc_idx.tsv", 
                      col_names = c("c1", "c2", "jacc_idx", "t1", "t2"))

all_comms <- lapply(tissues, function(tissue) {
  communities <- read_tsv(paste0(tissue, "/network_aracne/cancer-comm-info-100000.tsv"))
  communities$tissue <- tissue
  return(communities)
})
all_comms <- bind_rows(all_comms)

all_counts <- jacc_idxs %>% group_by(t1, t2) %>% tally()

jacc_idxs <- jacc_idxs %>% filter(!(t1 == "brain" & t2 =="bladder"))

jacc_idxs <- jacc_idxs %>% 
  unite(col = "source", t1, c1) %>% unite(col = "target", t2, c2) %>%
  mutate(jacc_idx = round(jacc_idx, 2))

all_comms <- all_comms %>% 
  unite(col = "id_comm", tissue, com_id, remove = F)

all_comms <- all_comms %>% 
  filter(id_comm %in% union(jacc_idxs$source, jacc_idxs$target))

jacc_idxs %>% write_tsv("pan-loss/network_intersections/jaccard_indexes_interactions.tsv")
all_comms %>% write_tsv("pan-loss/network_intersections/jaccard_indexes_vertices.tsv")
