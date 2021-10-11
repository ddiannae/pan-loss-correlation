library(readr)
library(tidyr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
            "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
            "testis", "thyroid","skin", "uterus")

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

all_interactions <- lapply(tissues, function(tissue) {
  interactions <- read_tsv(paste0("/datos/ot/diana/regulacion-trans/", tissue, 
                                  "/network_aracne/cancer-interactions-100000.tsv"))
  interactions$tissue <- tissue
  return(interactions)
})
all_interactions <- bind_rows(all_interactions)
all_interactions <- all_interactions %>% 
  mutate(id = paste0(pmin(source_ensembl, target_ensembl), "-", 
                     pmax(source_ensembl, target_ensembl)))

all_comms <- lapply(tissues, function(tissue) {
  communities <- read_tsv(paste0(tissue, "/network_aracne/cancer-comm-100000.tsv"))
  communities$tissue <- tissue
  return(communities)
})
all_comms <- bind_rows(all_comms)
max_comms <- all_comms %>% group_by(tissue) %>% summarise(max = max(community))

all_jis <- apply(tissue_combs, 1, function(row) {
  tissue1 =  unlist(row["t1"])
  tissue2 = unlist(row["t2"])
  cat("Getting ", tissue1, " ", tissue2, " intersections\n")
  m1 <- max_comms %>% filter(tissue == tissue1) %>% pull(max)
  m2 <- max_comms %>% filter(tissue == tissue2) %>% pull(max)
  
  high_jis <- parallel::mclapply(X = 1:m1, mc.cores = 70, FUN = function(c1) {
    set1 <- all_comms %>% filter(tissue == tissue1 & community == c1) %>% pull(ensembl)
    jis <- lapply(1:m2, function(c2) {
      set2 <- all_comms %>% filter(tissue == tissue2 & community == c2) %>% pull(ensembl)
      ji <- jacc_index(all_interactions %>% filter(tissue == tissue1), 
                       all_interactions %>% filter(tissue == tissue2), 
                       set1, set2)
      return(list(c1 = c1, c2 = c2, ji = ji))
    })
    jis <- bind_rows(jis) %>% filter(ji >= 0.5)
  })
  high_jis <- bind_rows(high_jis)
  high_jis$tissue1 <- tissue1
  high_jis$tissue2 <- tissue2
  write_tsv(high_jis, file = paste0("/datos/ot/diana/regulacion-trans/pan-loss/network_intersections/", 
                                    tissue1, "-", tissue2, "_jacc_idx.tsv"))
  return(high_jis)
})
