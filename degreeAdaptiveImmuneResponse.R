########################################################################
## Script to generate degree distribution of the  comms to comms 
## network with a highlight of the GO:0002250 communities. 
######################################################################

library(dplyr)
library(readr)
library(igraph)
library(ggplot2)
library(ggthemes)
library(scales)

getDegreesDistribution <-  function(cond) {
  tissues <- c("bladder", "brain", "breast", "colon", "esophagus",
               "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
               "testis", "thyroid","skin", "uterus")
  
  all_degrees <- lapply(tissues, function(tissue) {
    interactions <- read_tsv(paste0(tissue, "/network_aracne/communities/", cond, "-comm-network-all-louvain-100000.tsv")) %>%
      rename("weight" = "n")
    net <- igraph::graph_from_data_frame(interactions, directed = FALSE)
    net <- igraph::simplify(net)
    degree <- igraph::degree(net)
    tibble(community = names(degree), degree = unname(degree), tissue = tissue)
  })
  
  all_degrees <- bind_rows(all_degrees)
  all_degrees <- all_degrees %>% 
    filter(degree > 0)
  
  adaptive_enrich <- read_tsv(paste0("pan-loss/enrichments/comm-enrich-", cond, "_10.tsv")) %>%
    filter(id == "GO:0002250") %>% 
    select(id = comm_tissue)
  all_degrees %>%
    mutate(adaptive_degree = if_else(community %in% adaptive_enrich$id, degree, NA_real_))
}

degrees_cancer <- getDegreesDistribution("cancer")
degrees_normal <- getDegreesDistribution("normal")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

p <- ggplot(degrees_cancer) +
  geom_boxplot(aes(y = degree, x = tissue)) +
  geom_hline(aes(yintercept = adaptive_degree, linetype = "GO:0002250. Adaptive immune response"),
             color = "orange", size = 1) +
  facet_wrap(~tissue, scales="free", labeller = labeller(tissue = capitalize), nrow = 3) +
  theme_base(base_size = 25) +
  ylab("Community Degree") +
  xlab("") +
  theme(axis.text.x = element_blank(), legend.title=element_blank(), 
        axis.ticks.x=element_blank(),  legend.position="bottom",
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) 

png(paste0("pan-loss/network_aracne_plots/cancer_comms_comms_adaptive_immune.png"), 
    width = 850, height = 650)
print(p)
dev.off()

p <- ggplot(degrees_normal) +
  geom_boxplot(aes(y = degree, x = tissue)) +
  geom_hline(aes(yintercept = adaptive_degree, linetype = "GO:0002250. Adaptive immune response"),
             color = "orange", size = 1) +
  facet_wrap(~tissue, scales="free", labeller = labeller(tissue = capitalize), nrow = 3) +
  theme_base(base_size = 25) +
  ylab("Community Degree") +
  xlab("") +
  theme(axis.text.x = element_blank(), legend.title=element_blank(), 
        axis.ticks.x=element_blank(),  legend.position="bottom",
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) 

png(paste0("pan-loss/network_aracne_plots/normal_comms_comms_adaptive_immune.png"), 
    width = 850, height = 650)
print(p)
dev.off()
