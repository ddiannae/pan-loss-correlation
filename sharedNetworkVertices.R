library(readr)
library(dplyr)

ninter <- 10000
tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

vertices <- lapply(tissues, function(tissue) {
  file_vert <- paste0(tissue, "/network_aracne/cancer-vertices-", ninter, ".tsv")
  tissue_ver <- read_tsv(file=file_vert) %>% select(ensembl, chr, symbol) 
  tissue_ver$tissue <- tissue
  return(tissue_ver)
})
vertices <- bind_rows(vertices)


normal_vertices <- lapply(tissues, function(tissue) {
  file_vert <- paste0(tissue, "/network_aracne/normal-vertices-", ninter, ".tsv")
  tissue_ver <- read_tsv(file=file_vert) %>% select(ensembl, chr, symbol) 
  tissue_ver$tissue <- tissue
  return(tissue_ver)
})
normal_vertices <- bind_rows(normal_vertices)

rep_vertices <- vertices %>% group_by(ensembl, chr, symbol) %>% 
  tally() %>%filter(n>1) %>% arrange(desc(n))

rep_normal_vertices <- normal_vertices %>% group_by(ensembl, chr, symbol) %>% 
  tally() %>%filter(n>1) %>% arrange(desc(n)) %>% rename("n_sanos" = "n")

cancer_normal <- rep_vertices %>% inner_join(rep_normal_vertices, by=c("ensembl", "chr", "symbol"))


normal_repetidos 


table(rep_normal_vertices$n)

table(rep_vertices$n)
