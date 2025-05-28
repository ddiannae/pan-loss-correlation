#######################################################################
## Script to get bipartite communities for enrichments network.
#######################################################################

library(readr)
library(dplyr)
library(igraph)
library(stringr)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")
ni <- "100000"

enrichments_normal <- lapply(tissues, function(tissue) {
  enrch <- read_tsv(paste0(tissue, "/network_aracne/enrichments/go-normal-comm-all-", ni, ".tsv"))
  enrch$tissue <- tissue
  return(enrch)
})

enrichments_normal <- bind_rows(enrichments_normal)

enrichments_normal <- enrichments_normal %>% filter(p_adjust < 1e-10)
enrichments_normal %>% mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
  select(id, comm_tissue) %>% write_tsv("pan-loss/enrichments/comm-enrich-normal-interactions_10.tsv")

enrichments_normal %>% mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
    write_tsv("pan-loss/enrichments/comm-enrich-normal_10.tsv")

enrichments_normal %>% select(id, description) %>%
  mutate(description = stringr::str_to_title(description), type = "GO") %>%
  rename("label" = "description") %>% bind_rows(
  enrichments_normal %>% select(commun, tissue) %>%
    mutate(label = paste( stringr::str_to_title(tissue), commun, sep = "_"), id = str_to_lower(label), 
           tissue =  stringr::str_to_title(tissue), type = "COM") %>%
    select(-commun)
) %>%  write_tsv("pan-loss/enrichments/comm-enrich-normal-vertices_10.tsv")

enrichments_normal %>% select(id, description, tissue) %>% distinct() %>%
  group_by(id, description) %>% tally() %>% filter(n == 1)  %>% 
  write_tsv("pan-loss/enrichments/comm-enrich-normal-unicos_10.tsv")

enrichments_cancer <- lapply(tissues, function(tissue) {
  enrch <- read_tsv(paste0(tissue, "/network_aracne/enrichments/go-cancer-comm-all-", ni, ".tsv"))
  enrch$tissue <- tissue
  return(enrch)
})

enrichments_cancer <- bind_rows(enrichments_cancer)
enrichments_cancer <- enrichments_cancer %>% filter(p_adjust < 1e-10)

enrichments_cancer %>% mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
  select(id, comm_tissue) %>% write_tsv("pan-loss/enrichments/comm-enrich-cancer-interactions_10.tsv")

enrichments_cancer %>% mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
   write_tsv("pan-loss/enrichments/comm-enrich-cancer_10.tsv")


enrichments_cancer %>% select(id, description) %>% 
  mutate(description = stringr::str_to_title(description), type = "GO") %>%
  rename("label" = "description") %>% bind_rows(
    enrichments_cancer %>% select(commun, tissue) %>%
      mutate(label = paste( stringr::str_to_title(tissue), commun, sep = "_"), id = str_to_lower(label), 
             tissue =  stringr::str_to_title(tissue), type = "COM") %>%
      select(-commun)
  ) %>%  write_tsv("pan-loss/enrichments/comm-enrich-cancer-vertices_10.tsv")

enrichments_cancer %>% select(id, description, tissue) %>% distinct() %>%
  group_by(id, description) %>% tally() %>% filter(n == 1)  %>% 
  write_tsv("pan-loss/enrichments/comm-enrich-cancer-unicos_10.tsv")

go_normal <- enrichments_normal %>%
  pull(id) %>% unique()

go_cancer <- enrichments_cancer %>%
  pull(id) %>% unique()

enrichments_normal %>% 
  filter(!id %in% go_cancer) %>%
  mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
  select(id, comm_tissue) %>% 
  write_tsv("pan-loss/enrichments/comm-enrich-normal-only-interactions_10.tsv")

enrichments_normal %>% 
  filter(!id %in% go_cancer) %>%
  select(id, description) %>%
  mutate(description = stringr::str_to_title(description), type = "GO") %>%
  rename("label" = "description") %>% bind_rows(
    enrichments_normal %>% select(commun, tissue) %>%
      mutate(label = paste( stringr::str_to_title(tissue), commun, sep = "_"), id = str_to_lower(label), 
             tissue =  stringr::str_to_lower(tissue), type = "COM") %>%
      select(-commun)
  ) %>%  write_tsv("pan-loss/enrichments/comm-enrich-normal-only-vertices_10.tsv")


enrichments_cancer %>% 
  filter(!id %in% go_normal) %>%
  mutate(comm_tissue = paste(tissue, commun, sep = "_")) %>%
  select(id, comm_tissue) %>% 
  write_tsv("pan-loss/enrichments/comm-enrich-cancer-only-interactions_10.tsv")


enrichments_cancer %>% 
  filter(!id %in% go_normal) %>%
  select(id, description) %>% 
  mutate(description = stringr::str_to_title(description), type = "GO") %>%
  rename("label" = "description") %>% bind_rows(
    enrichments_cancer %>% select(commun, tissue) %>%
      mutate(label = paste( stringr::str_to_title(tissue), commun, sep = "_"), id = str_to_lower(label), 
             tissue =  stringr::str_to_lower(tissue), type = "COM") %>%
      select(-commun)
  ) %>%  write_tsv("pan-loss/enrichments/comm-enrich-cancer-only-vertices_10.tsv")
