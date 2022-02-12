library(readr)
library(dplyr)
library(igraph)


getComInfo <- function(cmembership, network){
  comp_info <- lapply(unique(cmembership), function(idc){
    mem <- names(cmembership[cmembership == idc])
    com <- induced.subgraph(network, mem)
    intras <- igraph::as_data_frame(com, what = "edges") %>% 
      filter(interaction_type == "Intra") %>% nrow()
    prs <- page.rank(com)
    chrs <- table(V(com)$chr)
    edges <- length(E(com))
    return(data.frame(com_id = idc,
                      pg_gene = names(which.max(prs$vector))[1], 
                      chr = names(which.max(chrs))[1], order = length(V(com)), 
                      size = edges,
                      intra_fraction = intras/edges))
  })
  comp_info <- plyr::compact(comp_info)
  comp_info <- plyr::ldply(comp_info)
  comp_info <- comp_info %>% arrange(desc(size))
  return(comp_info)
}

getCommunities <- function(cond) {
  interactions <- read_tsv(paste0("pan-loss/network_intersections/", cond, "-shared-interactions-100000.tsv"))
  vertices <- read_tsv(paste0("pan-loss/network_intersections/", cond, "-shared-vertices-ininter-100000.tsv"))
  
  interactions <- interactions %>% filter(nts >= 10) %>%
    dplyr::rename("from" = "source_ensembl", "to" = "target_ensembl")
  
  vertices <- vertices %>% 
    filter(ensembl %in% union(interactions$from,
                              interactions$to))  %>% 
    dplyr::rename("name" = "ensembl")
  
  net <- graph_from_data_frame(interactions, 
                               directed=F, vertices = vertices)
  
  comm <- cluster_louvain(graph = net)
  names(comm$membership) <- comm$names
  
  df_comm <- data.frame(comm$names, comm$membership)
  colnames(df_comm) <- c("ensembl", "community")
  
  comm_info <- getComInfo(comm$membership, net)
  
  cat("Saving files\n")
  write_tsv(df_comm, paste0("pan-loss/network_intersections/", cond, "-shared-communities-100000.tsv"))
  write_tsv(comm_info, paste0("pan-loss/network_intersections/", cond, "-shared-communities-info-100000.tsv"))
  
  
}

getCommunities("normal")
getCommunities("cancer")