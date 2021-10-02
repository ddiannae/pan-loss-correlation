library(igraph)
library(readr)
library(dplyr)
library(tidyr)

load("pan-loss/rdata/annot.RData")
getInteractionsByChr <- function(cond, ts, ni) {

  nets <- lapply(ts, function(tissue) {
    
    interactions <- read_tsv(paste0(tissue, "/network_aracne/", cond,  "-interactions-", ni, ".tsv")) 
    
    vertices <- read_tsv(paste0(tissue, "/network_aracne/", cond,  "-vertices-", ni, ".tsv"), 
                         col_types = cols(chr=col_character())) %>% 
      select(ensembl, chr, symbol) 
    chr_annot <- all_annot %>% group_by(chr) %>% tally(name = "total")
    
    colnames(vertices) <- c("source_ensembl", "source_chr", "source_symbol")
    interactions <- interactions %>% inner_join(vertices, by="source_ensembl")
    colnames(vertices) <- c("target_ensembl", "target_chr", "target_symbol")
    interactions <- interactions %>% inner_join(vertices, by="target_ensembl")
    
    interactions <- interactions %>% filter(interaction_type == "Inter")
    
    colnames(vertices) <- c("ensembl", "chr", "symbol")
    chr_vertices <- vertices %>% 
      filter(ensembl %in% union(interactions$source_ensembl, 
                                       interactions$target_ensembl)) %>%
      group_by(chr) %>% tally(name = "in_inter") %>% inner_join(
        vertices %>% group_by(chr) %>% tally(name = "in_net"), 
        by = "chr") %>% 
      inner_join(chr_annot, by = "chr") 
    
    interactions <- interactions %>% mutate(chrs = paste0(pmin(source_chr, target_chr), "-", 
                                          pmax(source_chr, target_chr))) 
    
    chr_interactions <- interactions %>% group_by(chrs) %>% 
      summarize(n_interactions = n()) %>% 
      separate(chrs, into = c("source_chr", "target_chr"), sep = "-")
    
    colnames(chr_vertices) <- c("source_chr", "source_in_inter", "source_in_net", "source_total")
    chr_interactions <- chr_interactions %>% inner_join(chr_vertices, by="source_chr")
    colnames(chr_vertices) <- c("target_chr", "target_in_inter", "target_in_net", "target_total")
    chr_interactions <- chr_interactions %>% inner_join(chr_vertices, by="target_chr")
    
    chr_interactions <- chr_interactions %>% 
      mutate(total_interactions = if_else(source_chr == target_chr, 
                                          as.integer((source_in_net*(source_in_net-1))/2), 
                                          source_in_net*target_in_net),
             fraction_interactions = n_interactions/total_interactions) %>% 
      select(source_chr, target_chr, n_interactions, total_interactions, fraction_interactions)
    
    colnames(chr_vertices) <- c("chr", "in_inter",  "in_net","total")
    chr_vertices <- chr_vertices %>% mutate(fraction_vertices = in_inter/in_net)
    
    chr_interactions %>% 
      write_tsv(paste0(tissue, "/network_aracne/", cond,  "-chr-interactions-", ni, ".tsv"))
    
    chr_vertices %>% 
      write_tsv(paste0(tissue, "/network_aracne/", cond,  "-chr-vertices-", ni, ".tsv"))
    
    net <- igraph::graph_from_data_frame(chr_interactions, directed=FALSE, vertices = chr_vertices)
    return(net)
  })
}

readNets <- function(cond, ts, ni) {
  nets <- lapply(ts, function(tissue) {
    chr_interactions <- read_tsv(paste0(tissue, "/network_aracne/", cond,  "-chr-interactions-", ni, ".tsv"))
    chr_interactions <- chr_interactions %>% filter(source_chr != target_chr)
    net <- igraph::graph_from_data_frame(chr_interactions, directed=FALSE, vertices = NULL)
    return(net)
  })
}

saveNetworkPlots <- function(nets, cond, ts, ni) {
  
  mapply(function(net, tissue) {
    ttitle <- tissue
    substring(ttitle, 1, 1) <- toupper(substring(ttitle, 1, 1))
    
   
    E(net)$width <- E(net)$fraction_interactions*5000
    V(net)$size <- V(net)$fraction_vertices*100
    l <- layout_with_fr(net)
    
    png(filename = paste0(tissue, "/network_aracne_plots/", cond,
                          "-chr-inter-interactions-", ni, ".png"), 
        width=1000, height=1000)
    plot.igraph(net, layout=l, vertex.color = rgb(233/255, 178/255, 0/255, alpha = 0.5),
                vertex.label.family = "Helvetica", vertex.frame.color = NA,
                vertex.label.color = "grey0", edge.color = "grey80",
                vertex.label.cex = 2)
    title(ttitle, cex.main = 3)
    dev.off()
    
    cutoff <- mean(E(net)$fraction_interactions)
    net_sp <- delete_edges(net, E(net)[fraction_interactions < cutoff])
    net_sp <- delete_vertices(net_sp, degree(net_sp)==0)
    l <- layout_with_fr(net_sp)
    
    png(filename = paste0(tissue, "/network_aracne_plots/", cond,
                          "-chr-inter-interactions-sparse-", ni, ".png"),
        width=1000, height=1000)
    plot.igraph(net_sp, layout= l, vertex.color = rgb(233/255, 178/255, 0/255, alpha = 0.5),
                vertex.label.family = "Helvetica", vertex.frame.color = NA,
                vertex.label.color = "grey0", edge.color = "grey80", 
                vertex.label.cex=2)
    title(ttitle, cex.main = 3)
    dev.off()
  }, nets, ts)
}


tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")
ninteractions <- "10000"
cancer_nets <- getInteractionsByChr("cancer", tissues, ninteractions)
#cancer_nets <- readNets("cancer", tissues, ninteractions)
names(cancer_nets) <- tissues
saveNetworkPlots(cancer_nets, "cancer", tissues, ninteractions)

normal_nets <- getInteractionsByChr("normal", tissues, ninteractions)
names(normal_nets) <- tissues
saveNetworkPlots(normal_nets, "normal", tissues, ninteractions)

cancer_only_nets <- getInteractionsByChr("cancer_only", tissues, ninteractions)
names(cancer_only_nets) <- tissues
saveNetworkPlots(cancer_only_nets, "cancer_only", tissues, ninteractions)

normal_only_nets <- getInteractionsByChr("normal_only", tissues, ninteractions)
names(normal_only_nets) <- tissues
saveNetworkPlots(normal_only_nets, "normal_only", tissues, ninteractions)

shared_nets <- getInteractionsByChr("shared", tissues, ninteractions)
names(shared_nets) <- tissues
saveNetworkPlots(shared_nets, "shared", tissues, ninteractions)
