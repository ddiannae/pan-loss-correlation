library(readr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")

ni <- "100000"

comm_cytobands <- lapply(tissues, function(tissue) {
  membership <- read_tsv(paste0(tissue, "/network_aracne/communities/cancer-comm-all-",
                         ni, ".tsv"))
  vertices <-  read_tsv(paste0(tissue, "/network_aracne/cancer-vertices-",
                              ni, ".tsv"))
  all_cytobands <- lapply(unique(membership$community), function(comm) {
    cytobands_in_comm <- membership %>% filter(community == comm) %>%
      inner_join(vertices, by = "ensembl") %>% 
      group_by(chr, cytoband) %>% tally(name = "vertices")
    
    cytobands_in_comm$community <- comm
    return(cytobands_in_comm)
  })
  all_cytobands <- bind_rows(all_cytobands)
  all_cytobands$tissue <- tissue
  return(all_cytobands)
})
comm_cytobands <- bind_rows(comm_cytobands)
barcodes_cancer <- read_tsv("pan-loss/cytobands/barcodes_fraction_cancer.tsv")
barcodes_cancer <- barcodes_cancer %>% 
  pivot_longer(cols = all_of(tissues), names_to = "tissue", values_to = "fraction") %>%
  filter(fraction > 0)

barcodes_cancer <- barcodes_cancer %>% 
  inner_join(comm_cytobands, by = c("tissue", "chr", "cytoband"))

enrich_by_comm <- read_tsv("pan-loss/enrichments/go_by_communities.tsv")
enrich_by_comm <- enrich_by_comm %>% rename("community" = "commun")

barcodes_cancer <- barcodes_cancer %>% 
  left_join(enrich_by_comm, by = c("tissue", "community"))


barcodes_cancer %>% write_tsv("pan-loss/cytobands/enriched_communities_cytbands.tsv")

barcodes_cancer <- barcodes_cancer %>% arrange(tissue, community)
