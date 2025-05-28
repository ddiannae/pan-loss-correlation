library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(ggthemes)

ninter <- "100000"
tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

getRepeatedInteractions <- function(cond, ts, ni) {
  
  ### Get interactions and add symbols and chromosomes
  interactions <- lapply(ts, function(tissue) {
    file_inter <- paste0(tissue, "/network_aracne/", cond, "-interactions-", ni, ".tsv")
    file_vert <- paste0(tissue, "/network_aracne/", cond, "-vertices-", ni, ".tsv")
    tissue_inter <- read_tsv(file=file_inter)
    tissue_ver <- read_tsv(file=file_vert,  col_types = cols(chr = col_character())) %>% 
      select(ensembl, chr, symbol) 
    colnames(tissue_ver) <- c("source_ensembl", "source_chr", "source_symbol")
    tissue_inter <- tissue_inter %>% inner_join(tissue_ver)
    colnames(tissue_ver) <- c("target_ensembl", "target_chr", "target_symbol")
    tissue_inter <- tissue_inter %>% inner_join(tissue_ver)
    tissue_inter$tissue <- tissue
    return(tissue_inter)
  })
  interactions <- bind_rows(interactions) 
  
  ### Setup id interactions
  interactions <- interactions %>% 
    mutate(inter = paste0(pmin(source_ensembl, target_ensembl),"-",
                          pmax(source_ensembl, target_ensembl)))
  
  ### Are interactions repeated on different tissues?
  reps <- interactions %>% group_by(inter, interaction_type) %>% 
    summarize(nts = n_distinct(tissue)) %>%
    filter(nts > 1) %>% arrange(desc(nts)) %>% 
    separate(inter, into=c("source_ensembl", "target_ensembl"))
  
  ### Get the vertices for repeated interactions
  sources <- interactions %>% select(source_ensembl, source_chr, source_symbol)
  colnames(sources) <- c("ensembl", "chr", "symbol")
  targets <- interactions %>% select(target_ensembl, target_chr, target_symbol) 
  colnames(targets) <- c("ensembl", "chr", "symbol")
  
  repsv <- sources %>% bind_rows(targets) %>% distinct() %>% 
    filter(ensembl %in% union(reps$source_ensembl, reps$target_ensembl))
  
  return(list(all_interactions=interactions, rep_interactions = reps, rep_vertices = repsv))
}

cancer <- getRepeatedInteractions("cancer", tissues, ninter)
normal <- getRepeatedInteractions("normal", tissues, ninter)

cancer$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer-shared-interactions-", ninter, ".tsv"))
cancer$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer-shared-vertices-ininter-", ninter, ".tsv"))

normal$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/normal-shared-interactions-", ninter, ".tsv"))
normal$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/normal-shared-vertices-ininter-", ninter, ".tsv"))

getBarcodes <-  function(cond, interactions, ts, ninter) {
  
  rep_groups <- interactions %>% group_by(inter, interaction_type) %>% 
    summarize(nts = n_distinct(tissue), tissue=tissue, mi=mi) %>%
    arrange(desc(nts)) %>% ungroup() %>%
    pivot_wider(id_cols = c(inter, interaction_type, nts), names_from = tissue, values_from = mi) %>%
    mutate_at(all_of(ts), function(x) (ifelse(is.na(x), FALSE, TRUE))) 
  
  return(rep_groups)
}

barcodes_cancer <- getBarcodes("cancer", cancer$all_interactions, tissues, ninter)

barcodes_cancer %>% mutate_at(all_of(tissues), as.numeric) %>% 
  unite(col = "barcode", all_of(tissues), sep="") %>%
  write_tsv(file = paste0("pan-loss/network_intersections/cancer-interaction-barcodes-", ninter, ".tsv"))

barcodes_normal <- getBarcodes("normal", normal$all_interactions, tissues, ninter)

barcodes_normal %>% mutate_at(all_of(tissues), as.numeric) %>% 
  unite(col = "barcode", all_of(tissues), sep="") %>%
  write_tsv(file = paste0("pan-loss/network_intersections/normal-interaction-barcodes-", ninter, ".tsv"))

savePlot <- function(barcodes, ts, cond, ni, color) {
  csums <- colSums(barcodes %>% select(all_of(ts)))
  
  barcodes <- barcodes %>% filter(nts > 1) 
  
  substring(colnames(barcodes), 1, 1) <- toupper(substring(colnames(barcodes), 1, 1))
  substring(ts, 1, 1) <- toupper(substring(ts, 1, 1))
  substring(names(csums), 1, 1) <-  toupper(substring(names(csums), 1, 1))
  inter_colors <- c("#C85200", "#1696AC")
  
  p <- upset(barcodes, ts, name = "", min_size=100, width_ratio = 0.15, 
             height_ratio = 0.70, stripes='white', sort_sets = "descending", n_intersections = 20,
             base_annotations = list('size' = 
                                       intersection_size(counts=FALSE, mode='exclusive_intersection',
                                                         mapping=aes(fill=Interaction_type)) +
                                       scale_fill_manual(name="Interaction type", values = inter_colors,
                                                         labels=c("Inter-chromosomal", "Intra-chromosomal"))
             ),
             themes = upset_modify_themes(
               list(default = theme(axis.text.x=element_blank(),
                                    panel.grid.major=element_blank(), 
                                    panel.grid.minor=element_blank(),
                                    text=element_text(size=30)),
                    'intersections_matrix'=theme(text=element_text(size=30))
               )),  
             set_sizes = (
               upset_set_size(geom=geom_bar(stat='count', color=color, fill=color),
                              mapping=aes(y=csums[x]-..count..),  position='right') +
               theme(axis.text.x=element_text(angle=90), text=element_text(size=30)) + 
               ylab("Unique interactions") +
               scale_y_continuous(limits = c(0, as.numeric(ni)), 
                                  breaks = seq(0, as.numeric(ni), by = as.numeric(ni)/4)) 
             ),  
             matrix = intersection_matrix(geom = geom_point(size = 8)),
             guides='over', wrap=TRUE
  )
  png(paste0("pan-loss/network_intersections/", cond, "-shared-interactions-", ni, ".png"), width = 2000, height = 1000)
  print(p)
  dev.off() 
}
### Generate plots
savePlot(barcodes_cancer, tissues, "cancer", ninter, "#a32e27")
savePlot(barcodes_normal, tissues, "normal", ninter, "#e3a098")
