library(readr)
library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(ggthemes)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus",
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
             "testis", "thyroid","skin", "uterus")
ni <- "100000"

enrichments <- lapply(tissues, function(tissue) {
  expr_assortativity <- read_tsv(paste0(tissue, "/network_aracne/assortativity/cancer-expr-assortativity-100000.tsv")) %>%
    select(community_id, diffraction, mean_log_fc)
  expr_assortativity$type <- "expr_assort"
  
  chr_assortativity <- read_tsv(paste0(tissue, "/network_aracne/assortativity/cancer-chr-assortativity-100000.tsv")) %>%
    select(community_id, diffraction)
  chr_assortativity$type <- "chr_assort"
  assort <- bind_rows(chr_assortativity, expr_assortativity)
  
  assort <- assort %>% pivot_wider(id_cols ="community_id",
                                   names_from = "type", 
                                   values_from = "diffraction")
  
  
  comm_info <- read_tsv(paste0(tissue, "/network_aracne/communities/cancer-comm-info-all-louvain-100000.tsv"), 
                        col_types = cols_only(com_id = col_double(), order = col_double(), 
                                              pg_gene = col_character())) %>%
    left_join(read_tsv(paste0(tissue, "/network_aracne/cancer-vertices-100000.tsv"), 
                       col_types = cols_only(ensembl = col_character(), symbol = col_character())),
              by = c("pg_gene" = "ensembl"))
  
  comm_enrich <- read_tsv(paste0(tissue, "/network_aracne/enrichments/go-cancer-comm-all-100000.tsv"),
                          col_types = cols( p_adjust = col_double())) %>%
    filter(p_adjust < 1e-10) %>%
    select(id, commun) %>% group_by(commun) %>% 
    tally(name = "nterms") 
  
  plot <- assort %>% 
    inner_join(expr_assortativity %>% select(community_id, mean_log_fc),
               by = c("community_id" = "community_id")) %>%
    inner_join(comm_info, by = c("community_id" = "com_id")) %>%
    left_join(comm_enrich, by = c("community_id" = "commun")) %>%
    filter(nterms > 0)
  
  plot$tissue <- tissue
  return(plot)
})
enrichments <- bind_rows(enrichments)
enrichments <- enrichments %>% mutate(community_id = paste0(tissue, "_", community_id))
enrichments %>% write_tsv("pan-loss/enrichments/cancer_expr_assortativity_enrichments.tsv")

p <- ggplot(enrichments, aes(x = chr_assort , y = expr_assort, color = tissue,
                        size = nterms, label = symbol)) +
 geom_label(check_overlap = TRUE) +
  theme_base() +
  ylim(-1, 1) +
  xlim(-1, 1) + 
  facet_wrap(~tissue, nrow = 3)

png(paste0("pan-loss/enrichments/expr_assortativity.png"), width = 3000, height = 1000)
print(p)
dev.off()
