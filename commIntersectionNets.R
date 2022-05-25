########################################################################
## Script to get a network of communities with overlap indexes 
## or jaccard indexes > 0.5
######################################################################

library(readr)
library(dplyr)
library(tidyr)
library(igraph)

getScoresNetwork <- function(cond, type) {
  idxs <- read_tsv(paste0("pan-loss/community_intersections/", cond, "_", type, "_scores.tsv"))
  idxs <- idxs %>%
    filter(score >= 0.5)
  
  all_comms <- read_tsv(paste0("pan-loss/network_aracne/", cond, "-communities-summary-100000.tsv"))
  
  idxs <- idxs %>% 
    unite(col = "source", tissue1, comm1) %>% unite(col = "target", tissue2, comm2) %>%
    mutate(score = round(score, 3))
  
  all_comms <- all_comms %>% 
    filter(id %in% union(idxs$source, idxs$target))
  
  idxs %>% write_tsv(paste0("pan-loss/community_intersections/", cond, "_", type, "_interactions.tsv"))
  all_comms %>% write_tsv(paste0("pan-loss/community_intersections/", cond, "_", type, "_vertices.tsv"))
}

getScoresNetwork("cancer", "overlap")
getScoresNetwork("normal", "overlap")
getScoresNetwork("cancer", "jaccard")
getScoresNetwork("normal", "jaccard")
