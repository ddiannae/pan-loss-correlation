library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexUpset)
library(ggthemes)

ninter <- "100000"
tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
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
cancer_only <- getRepeatedInteractions("cancer_only", tissues, ninter)
normal_only <- getRepeatedInteractions("normal_only", tissues, ninter)
shared <- getRepeatedInteractions("shared", tissues, ninter)

cancer$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer-shared-interactions-", ninter, ".tsv"))
cancer$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer-shared-vertices-ininter-", ninter, ".tsv"))

normal$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/normal-shared-interactions-", ninter, ".tsv"))
normal$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/normal-shared-vertices-ininter-", ninter, ".tsv"))

cancer_only$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer_only-shared-interactions-", ninter, ".tsv"))
cancer_only$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/cancer_only-shared-vertices-ininter-", ninter, ".tsv"))

normal_only$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/normal_only-shared-interactions-", ninter, ".tsv"))
normal_only$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/normal_only-shared-vertices-ininter-", ninter, ".tsv"))

shared$rep_interactions %>% write_tsv(file = paste0("pan-loss/network_intersections/shared-shared-interactions-", ninter, ".tsv"))
shared$rep_vertices %>% write_tsv(file = paste0("pan-loss/network_intersections/shared-shared-vertices-ininter-", ninter, ".tsv"))


getBarcodes <-  function(cond, interactions, ts, ninter) {
  
  rep_groups <- interactions %>% group_by(inter, interaction_type) %>% 
    summarize(nts = n_distinct(tissue), tissue=tissue, mi=mi) %>%
    arrange(desc(nts)) %>% ungroup() %>%
    pivot_wider(id_cols = c(inter, interaction_type, nts), names_from = tissue, values_from = mi) %>%
    mutate_at(all_of(ts), function(x) (ifelse(is.na(x), FALSE, TRUE))) 
  
  csums <- colSums(rep_groups %>% select(all_of(ts)))
  
  rep_groups <- rep_groups %>% mutate_at(all_of(ts), as.numeric) %>% 
    unite(col = "barcode", all_of(ts), sep="") 
  
  rep_groups %>% write_tsv(file = paste0("pan-loss/network_intersections/", cond, "-interaction-barcodes-", ninter, ".tsv"))
  
  return(rep_groups)
}

barcodes_cancer <- getBarcodes("cancer", cancer$all_interactions, tissues, ninter)
barcodes_normal <- getBarcodes("normal", normal$all_interactions, tissues, ninter)
barcodes_cancer_only <- getBarcodes("cancer_only", cancer_only$all_interactions, tissues, ninter)
barcodes_normal_only <- getBarcodes("normal_only", normal_only$all_interactions, tissues, ninter)

shared$all_interactions$mi <- 1
barcodes_shared <- getBarcodes("shared", shared$all_interactions, tissues, 0)



### Generate plots
rep_cancer_groups <- rep_cancer_groups %>% filter(nts > 1) 

substring(colnames(rep_cancer_groups), 1, 1) <- toupper(substring(colnames(rep_cancer_groups), 1, 1))
substring(ts, 1, 1) <- toupper(substring(ts, 1, 1))
substring(names(csums), 1, 1) <-  toupper(substring(names(csums), 1, 1))
inter_colors <- c("#C85200", "#1696AC")

p <- upset(rep_cancer_groups, ts, name = "", min_size=100, width_ratio = 0.15, 
           height_ratio = 0.75, stripes='white', sort_sets = F, n_intersections = 60,
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
                  'intersections_matrix'=theme(text=element_text(size=25))
             )),  
           set_sizes = (
             upset_set_size(geom=geom_bar(stat='count', color=color, fill=color),
                            mapping=aes(y=csums[x]-..count..),  position='right')
             + theme(axis.text.x=element_text(angle=90), text=element_text(size=30)) 
             + ylab("Unique interactions")
           ),  
           matrix = intersection_matrix(geom = geom_point(size = 5)),
           guides='over', wrap=TRUE
)
png(paste0("pan-loss/", cond, "-shared-interactions-", ninter, ".png"), width = 2000, height = 800)
print(p)
dev.off()  


