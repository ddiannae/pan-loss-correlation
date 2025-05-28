library(igraph)
library(vroom)
library(ggplot2)
library(ggthemes)

getCommKnessDistribution <-  function(cond) {
  tissues <- c("bladder", "brain", "breast", "colon", "esophagus",
               "kidney", "liver", "lung", "ovary", "pancreas", "prostate",
               "testis", "thyroid","skin", "uterus")
  
  all_comm_kness <- lapply(tissues, function(tissue) {
    interactions <- vroom::vroom(paste0(tissue, "/network_aracne/", cond, "-interactions-100000.tsv"))
    vertices <- vroom::vroom(paste0(tissue, "/network_aracne/", cond, "-vertices-100000.tsv"), 
                             col_types = cols_only("ensembl" = col_character())) %>%
      dplyr::inner_join(vroom::vroom(paste0(tissue, 
                                            "/network_aracne/communities/", cond, "-comm-all-louvain-100000.tsv")),
                        by = "ensembl")
    comms_info <- vroom::vroom(paste0(tissue, "/network_aracne/communities/", cond, "-comm-info-all-louvain-100000.tsv"),
                               col_types = cols_only("com_id" = col_double(), "order" = col_double()))
    
    net <- igraph::graph_from_data_frame(interactions, directed = F, vertices = vertices)
    
    vertices <- igraph::coreness(net) %>%
      stack() %>%
      dplyr::rename(kness = values, ensembl = ind) %>%
      dplyr::inner_join(vertices, by = "ensembl")
    
    kness_comm <- vertices %>% 
      dplyr::group_by(community) %>%
      dplyr::summarise(mean_kness = mean(kness)) %>%
      dplyr::inner_join(comms_info, by = c("community" = "com_id")) %>%
      dplyr::filter(order > 5) %>%
      dplyr::mutate(community = paste(tissue, community, sep = "_" ), 
                    tissue = tissue)
  
  })
  
  all_comm_kness <- bind_rows(all_comm_kness)
  
  adaptive_enrich <- read_tsv(paste0("pan-loss/enrichments/comm-enrich-", cond, "_10.tsv")) %>%
    filter(id == "GO:0002250") %>% 
    select(id = comm_tissue)
  
  all_comm_kness %>%
    mutate(adaptive_kness = if_else(community %in% adaptive_enrich$id, mean_kness, NA_real_))
}


kness_cancer <- getCommKnessDistribution("cancer")
kness_normal <- getCommKnessDistribution("normal")

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

p <- ggplot(kness_cancer) +
  geom_boxplot(aes(y = mean_kness, x = tissue)) +
  geom_hline(aes(yintercept = adaptive_kness, linetype = "GO:0002250. Adaptive immune response"),
             color = "orange", size = 1) +
  facet_wrap(~tissue, scales="free", labeller = labeller(tissue = capitalize), nrow = 3) +
  theme_base(base_size = 25) +
  ylab("Mean Community Coreness") +
  xlab("") +
  theme(axis.text.x = element_blank(), legend.title=element_blank(), 
        axis.ticks.x=element_blank(),  legend.position="bottom",
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) 

png(paste0("pan-loss/network_aracne_plots/cancer_kness_adaptive_immune.png"), 
    width = 850, height = 650)
print(p)
dev.off()

p <- ggplot(kness_normal) +
  geom_boxplot(aes(y = mean_kness, x = tissue)) +
  geom_hline(aes(yintercept = adaptive_kness, linetype = "GO:0002250. Adaptive immune response"),
             color = "orange", size = 1) +
  facet_wrap(~tissue, scales="free", labeller = labeller(tissue = capitalize), nrow = 3) +
  theme_base(base_size = 25) +
  ylab("Mean Community Coreness") +
  xlab("") +
  theme(axis.text.x = element_blank(), legend.title=element_blank(), 
        axis.ticks.x=element_blank(),  legend.position="bottom",
        plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1)) 

png(paste0("pan-loss/network_aracne_plots/normal_kness_adaptive_immune.png"), 
    width = 850, height = 650)
print(p)
dev.off()
