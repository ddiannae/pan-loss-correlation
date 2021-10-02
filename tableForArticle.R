library(readr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

interactions <- lapply(tissues, function(tissue) {
  lines_matrix <- as.numeric(strsplit(system2(command = "wc", 
          args = c("-l", paste0(tissue, "/results/*_*_*_si-arsyn_cancer.tsv")),
          stdout = TRUE), split = " ")[[1]][1]) -1
  lines_intra <- as.numeric(strsplit(system2(command = "wc", 
          args = c("-l", paste0(tissue, "/distance_analysis/cancer-all-distance-mi.tsv")),
          stdout = TRUE), split = " ")[[1]][1]) -1
  
  return(list(tissue = tissue, total=(lines_matrix*(lines_matrix - 1))/2, intra = lines_intra))
})
interactions <- bind_rows(interactions)
interactions <- interactions %>% mutate(fraction_intra = intra/total)

interactions_normal <- lapply(tissues, function(tissue) {
  inter_count <- read_tsv(paste0(tissue, "/distance_analysis/normal-intra-inter-count-onek-bins.tsv"))
  return(list(tissue = tissue, inter = sum(inter_count$inter), intra = sum(inter_count$intra)))
})
interactions_normal <- bind_rows(interactions_normal)
interactions_normal <- interactions_normal %>% mutate(total = inter + intra, fraction = intra/total)


interactions2 <- bind_rows(interactions2)
interactions2 <- interactions2 %>% mutate(total = inter + intra, fraction = intra/total)

interactions$total == interactions2$total
interactions$intra == interactions2$intra
interactions$fraction_intra == interactions2$fraction

mean(interactions2$fraction)

interactions2 <- interactions2 %>% arrange(desc(total))
interactions_normal <- interactions_normal %>% arrange(desc(total))
