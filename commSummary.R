library(readr)
library(dplyr)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
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

getSummaries("cancer", tissues) %>% 
  write_tsv("pan-loss/network_aracne/cancer-communities-summary-100000.tsv")
getSummaries("normal", tissues) %>% 
  write_tsv("pan-loss/network_aracne/normal-communities-summary-100000.tsv")

