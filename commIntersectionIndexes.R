#######################################################################
## Script to get overlap indexes and jaccard indexes between 
## communities from all tissues and from cancer and normal 
######################################################################

library(readr)
library(tidyr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
            "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
            "testis", "thyroid","skin", "uterus")

setwd("/datos/ot/diana/regulacion-trans/")
tissue_combs <- expand_grid(t1 = tissues, t2 = tissues) %>% 
  filter(t1 != t2) %>%
  mutate(id = paste0(pmin(t1, t2), "-", pmax(t1, t2))) %>% 
  distinct(id) %>% separate(id, into = c("t1", "t2"))

jacc_index <- function(net1, net2, comm1, comm2) {
  links1 <- net1 %>% filter(source_ensembl %in% comm1 & 
                              target_ensembl %in% comm1) %>% pull(id)
  links2 <- net2 %>% filter(source_ensembl %in% comm2 & 
                              target_ensembl %in% comm2) %>% pull(id)
  return(length(intersect(links1, links2))/length(union(links1, links2)))
}

get_scores <- function(cond, type) {
  cat("Working with condition ", cond, "\n")
  all_interactions <- lapply(tissues, function(tissue) {
    interactions <- read_tsv(paste0(tissue, 
                                    "/network_aracne/", cond, "-interactions-100000.tsv"))
    interactions$tissue <- tissue
    return(interactions)
  })
  all_interactions <- bind_rows(all_interactions)
  all_interactions <- all_interactions %>% 
    mutate(id = paste0(pmin(source_ensembl, target_ensembl), "-", 
                       pmax(source_ensembl, target_ensembl)))
  
  all_comms <- lapply(tissues, function(tissue) {
    communities <- read_tsv(paste0(tissue, "/network_aracne/communities/", cond, "-comm-all-louvain-100000.tsv"))
    communities$tissue <- tissue
    return(communities)
  })
  all_comms <- bind_rows(all_comms)
  max_comms <- all_comms %>% group_by(tissue) %>% 
    summarise(max = max(community)) %>% pull(max, name = tissue)
  
  cat("Getting parallel scores")
  all_scores <- parallel::mclapply(X = 1:nrow(tissue_combs), mc.cores = 77, FUN = function(i) {
    tissue1 <-  unlist(tissue_combs[i, "t1"])
    tissue2 <-  unlist(tissue_combs[i, "t2"])
    cat("Getting ", tissue1, " ", tissue2, " intersections\n")
    j_scores <- lapply(1:max_comms[tissue1], function(j) {
      set1 <- all_comms %>% filter(tissue == tissue1 & community == j) %>% pull(ensembl)
      links1 <- all_interactions %>% 
        filter(tissue == tissue1) %>% 
        filter(source_ensembl %in% set1 & 
                 target_ensembl %in% set1) %>% pull(id)
      k_scores <- lapply(j:max_comms[tissue2], function(k) {
        set2 <- all_comms %>% filter(tissue == tissue2 & community == k) %>% pull(ensembl)
        links2 <- all_interactions %>% 
          filter(tissue == tissue2) %>% 
          filter(source_ensembl %in% set2 & 
                   target_ensembl %in% set2) %>% pull(id)
        
        intersection <- length(intersect(links1, links2))
        if(type == "overlap") {
          return(intersection/min(length(links1), length(links2)))
        } else {
          return(intersection/length(union(links1, links2)))
        }
      })
      tibble(comm1 = j, comm2 = j:max_comms[tissue2], score = unlist(k_scores))
    })
    
    bind_rows(j_scores)  %>%
      filter(score > 0) %>%
      mutate(tissue1 = tissue1, tissue2 = tissue2)
  })
  
  bind_rows(all_scores) %>%
    write_tsv(paste0("pan-loss/community_intersections/", cond, "_", type, "_scores.tsv"))
  
}

get_scores("normal", "overlap")
get_scores("cancer", "overlap")
get_scores("normal", "jaccard")
get_scores("cancer", "jaccard")
