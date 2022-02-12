library(readr)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

cat("Reading files\n")

load("pan-loss/rdata/annot.RData") 

getEnrichments <- function(cond) {
  membership <- read_tsv(paste0("pan-loss/network_intersections/", cond, "-shared-communities-100000.tsv"))
  
  all_annot <- all_annot %>% pull(ensembl_id)
  
  all_enrichments <- lapply(X = unique(membership$community),
                            FUN = function(com){
                              
                              cat("Working with community: ", com, "\n")
                              
                              gene_list <- membership %>% filter(community == com) %>%
                                pull(ensembl)
                              
                              if(length(gene_list) >= 5) {
                                
                                terms <- enrichGO(gene          = gene_list,
                                                  universe      = all_annot,
                                                  OrgDb         = org.Hs.eg.db,
                                                  keyType       = "ENSEMBL",
                                                  ont           = "BP",
                                                  pAdjustMethod = "BH",
                                                  pvalueCutoff  = 0.005,
                                                  qvalueCutoff  = 0.01,
                                                  minGSSize     = 10,
                                                  readable      = FALSE)
                                
                                tryCatch({
                                  bp <- pairwise_termsim(terms)
                                  bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
                                  simple_results <- bp2@result
                                  if(nrow(simple_results) > 0) {
                                    simple_results$commun <- com 
                                  }
                                  return(simple_results)
                                }, error = function(cond) {
                                  return(NULL)
                                })
                                
                              }
                              return(NULL)
                            })
  
  bind_rows(all_enrichments) %>% 
    janitor::clean_names() 
  
  all_enrichments <- bind_rows(all_enrichments) %>% 
    janitor::clean_names()
  
  all_enrichments %>%
    write_tsv(paste0("pan-loss/network_intersections/", cond, "-shared-enrichments-100000.tsv"))
  
  all_enrichments %>% filter(p_adjust <= 1e-5)  %>%
    write_tsv(paste0("pan-loss/network_intersections/", cond, "-shared-enrichments-e_5-100000.tsv"))
  
}

getEnrichments("normal")
getEnrichments("cancer")