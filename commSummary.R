library(readr)
library(dplyr)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")
inter <- "100000"
load("pan-loss/rdata/annot.RData")
saveCommSummary <- function(ts, cond, ni) {
  
  enrich_by_comm <- lapply(ts, function(tissue) {
    go_enrich <- read_tsv(paste0(tissue, "/network_aracne/enrichments/go-", cond, "-comm-all-",
                                 ni, ".tsv"))
    go_enrich <- go_enrich %>% filter(p_adjust < 1e-10) %>%
      group_by(commun) %>% tally(name = "terms")
    go_enrich$tissue <- tissue
    return(go_enrich)
  })
  enrich_by_comm <- bind_rows(enrich_by_comm)
  
  chr_assort <- lapply(ts, function(tissue) {
    chra <- read_tsv(paste0(tissue,
                            "/network_aracne/assortativity/", cond, "-chr-assortativity-",
                            ni, ".tsv"))
    chra <- chra %>% mutate(tissue = tissue) %>%
      select(community_id, diffraction, tissue)
    chra$type <- "chr_assort"
    return(chra)
  })
  chr_assort <- bind_rows(chr_assort)
  
  if(cond == "cancer") {
    expr_assort <- lapply(ts, function(tissue) {
      ea <- read_tsv(paste0(tissue,
                            "/network_aracne/assortativity/cancer-expr-assortativity-100000.tsv"))
      ea <- ea %>% mutate(tissue = tissue) %>%
        select(community_id, diffraction, tissue)
      ea$type <- "expression_assort"
      return(ea)
    })
    expr_assort <- bind_rows(expr_assort)
    
    assort <- bind_rows(chr_assort, expr_assort)
  } else {
    assort <- chr_assort
  }
  assort <- assort %>% pivot_wider(id_cols = c("tissue", "community_id"),
                                   names_from = "type", 
                                   values_from = "diffraction")
  
  all_infos <- lapply(ts, function(tissue) {
    cinfo <- read_tsv(paste0(tissue,
                          "/network_aracne/communities/", cond, "-comm-info-all-louvain-100000.tsv"),
                      col_types = cols(chr = col_character()))
    cinfo$tissue <- tissue
    return(cinfo)
  })
  all_infos <- bind_rows(all_infos)
  
all_infos <- all_infos %>%
    left_join(assort, by=c("com_id" = "community_id", "tissue" = "tissue")) %>%
    left_join(enrich_by_comm, by=c("com_id" = "commun", "tissue" = "tissue")) %>%
    unite(col = "id", tissue, com_id, remove = FALSE) %>%
    mutate(terms = ifelse(is.na(terms), 0, terms)) %>%
    select(id, tissue, everything()) %>%
    left_join(all_annot %>% select(ensembl_id, gene_name), 
              by = c("pg_gene" = "ensembl_id"))
  write_tsv(all_infos, paste0("pan-loss/network_aracne/", cond, "-communities-summary-", ni,".tsv"))
}
saveCommSummary(tissues, "cancer", inter)
saveCommSummary(tissues, "normal", inter)
