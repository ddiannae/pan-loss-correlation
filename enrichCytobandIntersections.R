library(readr)
library(tidyr)

load("pan-loss/rdata/annot.RData")
intersections <- read_tsv("pan-loss/cytobands/cancer_interactions_in_intersections.tsv") %>% 
  unite(col = "region", chr, cytoband, sep = "")
isumm <- read_tsv("pan-loss/cytobands/intersections_summary.tsv") %>% 
  unite(col = "region", chr, cytoband, sep = "")
isumm <- isumm %>% filter(n_intersections >= 3) 

intersections <- intersections %>% filter(region %in% isumm$region)

enrichments <- lapply(isumm %>% pull(region),  function(reg) {      
  genes_in_region <- intersections %>% filter(region == reg) %>%
    dplyr::select(source_ensembl, target_ensembl) %>% unique() 
  
  genes_in_region <- union(genes_in_region$source_ensembl, 
                           genes_in_region$target_ensembl)
  
  terms <- enrichGO(gene = genes_in_region, 
                    ont = "BP",
                    keyType = 'ENSEMBL',
                    universe = all_annot$ensembl_id,
                    OrgDb  = org.Hs.eg.db,
                    pAdjustMethod = "BH",
                    readable = TRUE,
                    minGSSize = 10,
                    maxGSSize = 500)
  if(!is.null(terms)) {
    results <- terms@result %>% filter(pvalue < 0.05 & qvalue < 0.1)
    if(nrow(results) > 0) {
      results$region <- reg
    }
  } else {
    results <- NULL
  }
  tryCatch({
    bp <- pairwise_termsim(terms)
    bp2 <- simplify(bp, cutoff=0.7, by="p.adjust", select_fun=min)
    simple_results <- bp2@result
    if(nrow(simple_results) > 0) {
      simple_results$region <- reg
    }
    return(list(all_results = results, simple_results = simple_results))
  }, error = function(cond) {
    return(list(all_results = results, simple_results = NULL))
  })
})
all_results <- bind_rows(lapply(enrichments, "[[", "all_results"))
simple_results <- bind_rows(lapply(enrichments, "[[", "simple_results"))


simple_results  %>% write_tsv("pan-loss/cytobands/enrichments_intersections_simple.tsv")
all_results  %>% write_tsv("pan-loss/cytobands/enrichments_intersections.tsv")
