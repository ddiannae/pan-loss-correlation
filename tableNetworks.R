library(readr)
library(dplyr)
library(stringr)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

getNetworkAttributes <- function(ts, cond) {
  networks_attr <- lapply(ts, function(tissue) {
    nattr <- read_tsv(paste0(tissue, "/network_aracne/", cond, "-network-stats-100000.tsv"))
    nattr$tissue <- tissue
    return(nattr)
  })
  networks_attr <- bind_rows(networks_attr)
  networks_attr$cond <- cond
  return(networks_attr)
}

normal_networks <- getNetworkAttributes(tissues, "normal")
cancer_networks <- getNetworkAttributes(tissues, "cancer")
attrs <- c("vertices", "no_components", "inter_fraction")

networks <- bind_rows(normal_networks, cancer_networks) %>%
  filter(statistic %in% attrs) %>% 
  mutate(tissue = str_to_title(tissue), cond = str_to_title(cond)) %>%
  arrange(tissue, desc(cond)) %>% 
  pivot_wider(id_cols = c(tissue, cond), 
              names_from = statistic, values_from = value) %>%
  select(tissue, cond, vertices, no_components, inter_fraction) %>%
  rename("Tissue" = "tissue", "Condition" = "cond", 
         "Genes" = "vertices", "Fraction of inter-chromosomal\ninteractions" = "inter_fraction",
         "Components" = "no_components")

write_tsv(networks, file = "pan-loss/network_article_table.tsv")
 