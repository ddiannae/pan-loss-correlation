library(readr)
library(tidyr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")
ni <- "100000"

saveExprAssort <- function(ts, cond, ni) {
  expr_assort <- lapply(ts, function(tissue) {
    eat <- read_tsv(paste0(tissue, "/network_aracne/assortativity/", cond, "-comm-assort-enrich-",
                           ni, ".tsv"))
    eat$tissue <- tissue
    return(eat)
  })
  expr_assort <- bind_rows(expr_assort)
  expr_assort <- expr_assort %>% rename("symbol" = "name")  %>% 
    unite(col = "name", tissue, community_id)
  
  jacc_idxs <- read_tsv(paste0("pan-loss/community_intersections/", cond, "_all_jacc_idxs.tsv"),
                        col_names = c("c1", "c2", "jacc_idx", "t1", "t2"))
  
  jacc_idxs <- jacc_idxs %>% unite(col = "source", t1, c1) %>%
    unite(col = "target", t2, c2)
  
  expr_assort <- expr_assort %>% 
    filter(name %in% union(jacc_idxs$source, jacc_idxs$target)) %>%
    separate(col = "name", into = c("tissue", "comm"), remove = F)
  
  jacc_idxs <- jacc_idxs %>% filter(source %in% expr_assort$name &
                                      target %in%  expr_assort$name)
  
  jacc_idxs %>% write_tsv(file = paste0("pan-loss/community_intersections/", cond, "-comms-enrich-assort-", 
                                        ni, "-interactions.tsv"))
  expr_assort %>% write_tsv(file = paste0("pan-loss/community_intersections/", cond, "-comms-enrich-assort-", 
                                          ni, "-vertices.tsv"))
}

saveExprAssort(tissues, "cancer", ni )
saveExprAssort(tissues, "normal", ni )
