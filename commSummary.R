########################################################################
## Script to integrate information from communities per condition. 
## The output file contains expr and chr assortativity, intra-fraction, 
## enriched terms, and tissue per community.
######################################################################

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")
inter <- "100000"

getSummaries <- function(cond, tss) {
  all_summaries <- lapply(tss, function(ts) {
    read_tsv(paste0(ts, "/network_aracne/assortativity/", cond, "-comm-summary-100000.tsv")) %>%
      mutate(tissue = ts, id = paste(tissue, community_id, sep="_")) %>%
      select(id, tissue, community_id, everything())
  })  
  return(bind_rows(all_summaries))
}

getSummaries("normal", tissues) %>% 
  write_tsv("pan-loss/network_aracne/normal-communities-summary-100000.tsv")

getSummaries("cancer", tissues) %>% 
  write_tsv("pan-loss/network_aracne/cancer-communities-summary-100000.tsv")
