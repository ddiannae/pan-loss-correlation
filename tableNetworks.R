library(readr)
library(dplyr)
library(stringr)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
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

getNetworkCommunities <- function(ts, cond) {
  networks_attr <- lapply(tissues, function(tissue) {
    nattr <- read_tsv(paste0(tissue, "/network_aracne/communities/", cond, "-comm-info-all-louvain-100000.tsv"))
    return(tibble(statistic = "coms", value = max(nattr$com_id), tissue = tissue))
  })
  networks_attr <- bind_rows(networks_attr)
  networks_attr$cond <- cond
  return(networks_attr)
}

getDEG <- function(ts) {
  networks_attr <- lapply(tissues, function(tissue) {
    deg_file <- dir(path = paste0(tissue, "/deg"), pattern = "*_*_*_si-arsyn_deg_results.tsv", full.names = TRUE)[1]
    deg_genes <- read_tsv(paste0(tissue, "/network_aracne/cancer-vertices-100000.tsv")) %>%
      left_join(read_tsv(deg_file, col_types = cols_only("ensembl_id" = col_character(), "log_fc" = col_double())), 
                by = c("ensembl" = "ensembl_id"))
    return(tibble(statistic = c("up", "down"), value = c(sum(deg_genes$log_fc > 0), sum(deg_genes$log_fc < 0)),
                  tissue = c(tissue, tissue)))
  })
  networks_attr <- bind_rows(networks_attr)
  networks_attr$cond <- "cancer"
  return(networks_attr)
}

normal_networks <- getNetworkAttributes(tissues, "normal") %>%
  bind_rows(getNetworkCommunities(tissues, "normal"))

cancer_networks <- getNetworkAttributes(tissues, "cancer") %>%
  bind_rows(getNetworkCommunities(tissues, "cancer")) %>%
  bind_rows(getDEG(tissues))

attrs <- c("vertices", "inter_fraction",  "coms", "up", "down")

networks <- bind_rows(normal_networks, cancer_networks) %>%
  filter(statistic %in% attrs) %>% 
  mutate(tissue = str_to_title(tissue), cond = str_to_title(cond)) %>%
  arrange(tissue, desc(cond)) %>% 
  pivot_wider(id_cols = c(tissue, cond), 
              names_from = statistic, values_from = value) %>%
  select(tissue, cond, vertices, inter_fraction, coms, up, down) %>%
  mutate(inter_fraction = 1 - inter_fraction) %>%
  rename("Tissue" = "tissue", "Condition" = "cond", 
         "Genes" = "vertices", 
         "Intra- fraction" = "inter_fraction", 
         "Communities" = "coms")

write_tsv(networks, file = "pan-loss/network_article_table.tsv")
 