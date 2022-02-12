library(readr)
library(dplyr)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")

algorithms <- c("fgreedy", "infomap", "louvain", "leadeigen")

mvals_cancer <- lapply(tissues, function(tissue) {
  moduls <- sapply(algorithms, function(alg) {
    file <-  paste0(tissue, "/network_aracne/log/cancer_communities_all_", 
                    alg, "_100000.log")
    mod <- as.numeric(trimws(strsplit(grep("modularity",readLines(file), value = TRUE), "  ")[[1]][3]))
    return(mod)
  })
  moduls <- tibble(val = moduls, algorithm = names(moduls))
  moduls$tissue <- tissue
  return(moduls)
})
mvals_cancer <- bind_rows(mvals_cancer)
mvals_cancer <- mvals_cancer %>% pivot_wider(id_cols = tissue, 
                                                   names_from = algorithm, 
                                                   values_from = val)
mvals_cancer$condition <- "cancer"

mvals_normal <- lapply(tissues, function(tissue) {
  moduls <- sapply(algorithms, function(alg) {
    file <-  paste0(tissue, "/network_aracne/log/normal_communities_all_", 
                    alg, "_100000.log")
    mod <- as.numeric(trimws(strsplit(grep("modularity",readLines(file), value = TRUE), "  ")[[1]][3]))
    return(mod)
  })
  moduls <- tibble(val = moduls, algorithm = names(moduls))
  moduls$tissue <- tissue
  return(moduls)
})
mvals_normal <- bind_rows(mvals_normal)
mvals_normal <- mvals_normal %>% pivot_wider(id_cols = tissue, 
                                             names_from = algorithm, 
                                             values_from = val)
mvals_normal$condition <- "normal"

mvals <- bind_rows(mvals_normal, mvals_cancer) %>% 
  arrange(tissue) %>% select(tissue, condition, louvain, fgreedy, infomap, leadeigen)

mvals %>% write_tsv("pan-loss/network_aracne/modularity_vals.tsv")
