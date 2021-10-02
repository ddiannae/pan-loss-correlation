library(readr)
library(dplyr)
library(tidyr)

load("pan-loss/rdata/annot.RData")

getExclusiveNetworks <- function(cond, ts, ni) {
  
  barcodes <- read_tsv(paste0("pan-loss/", cond, "-interaction-barcodes-", ni, ".tsv"))
  barcodes <- barcodes %>% 
    separate(inter, into =c("source_ensembl", "target_ensembl")) %>%
    separate(barcode, into = tissues, sep = 1:length(tissues))
  
  l <- lapply(tissues, function(tissue) { 
    interactions <- barcodes %>% filter(nts == 1 & get(tissue) == "1") %>% 
      select(source_ensembl, target_ensembl, interaction_type)
    
    all_annot %>%
      filter(ensembl_id %in% union(interactions$source_ensembl, interactions$target_ensembl)) %>%
      write_tsv(paste0("pan-loss/exclusive_networks/", tissue, "-", cond, "-vertices-", ni, ".tsv"))
    
    interactions %>%
      write_tsv(paste0("pan-loss/exclusive_networks/", tissue, "-", cond, "-interactions-", ni, ".tsv"))

  })
}

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")
ninter <- "10000"

getExclusiveNetworks("cancer_only", tissues, ninter)
getExclusiveNetworks("normal_only", tissues, ninter)
getExclusiveNetworks("shared", tissues, ninter)
