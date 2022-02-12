library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(tidyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "skin", "testis", "thyroid","uterus")
names(tissues) <- tissues
substr(names(tissues),1, 1) <- toupper(substr(names(tissues),1, 1))
  
# colors <- c("yellow", "grey", "pink", "navyblue", "#ccccff",
#             "orange", "#50C878", "white", "darkcyan", "purple", "lightblue",
#             "black", "orchid", "blue",  "peachpuff")
# names(colors) <- tissues

color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")
names(color_pal) <- labels

chr_cancer_assort <- lapply(tissues, function(tissue) {
  chra <- read_tsv(paste0(tissue,
                        "/network_aracne/assortativity/cancer-chr-assortativity-100000.tsv"))
  chra <- chra %>% mutate(tissue = tissue) %>%
    select(community_id, diffraction, tissue)
  chra$cond <- "cancer"
  return(chra)
})

chr_cancer_assort <- bind_rows(chr_cancer_assort)

chr_normal_assort <- lapply(tissues, function(tissue) {
  chra <- read_tsv(paste0(tissue,
                          "/network_aracne/assortativity/normal-chr-assortativity-100000.tsv"))
  chra <- chra %>% mutate(tissue = tissue) %>%
    select(community_id, diffraction, tissue)
  chra$cond <- "normal"
  return(chra)
})
chr_normal_assort <- bind_rows(chr_normal_assort)
chr_assort <- bind_rows(chr_normal_assort, chr_cancer_assort)

chr_assort$cond <- factor(chr_assort$cond, levels = c("normal", "cancer"), labels = labels)

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

p <- ggplot(chr_assort , aes(x = tissue, y = diffraction, fill = cond,
                          color = cond)) +
  geom_violin() +
  scale_fill_manual(values = color_pal) +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  ylab("Chromosomal\nAssortativity") +
  xlab("") +
  theme(legend.position = "none", axis.text.x = element_blank(),
        plot.background = element_blank(), 
        strip.text.x = element_text(size = 22),
        strip.text.y = element_text(size = 20)) + 
  facet_grid(cond~tissue, scale= "free_x", labeller = labeller(tissue = capitalize))

png(paste0("pan-loss/network_aracne_plots/assortativity/chr_assortativity.png"), 
    width = 2000, height = 350)
print(p)
dev.off()


expr_assort <- lapply(tissues, function(tissue) {
  ea <- read_tsv(paste0(tissue,
                        "/network_aracne/assortativity/cancer-expr-assortativity-100000.tsv"))
  ea <- ea %>% mutate(tissue = tissue) %>%
    select(community_id, diffraction, tissue)
  ea$cond <- "cancer"
  return(ea)
})

expr_assort <- bind_rows(expr_assort)
expr_assort$cond <- factor(expr_assort$cond, levels = c("normal", "cancer"), labels = labels)

p <- ggplot(expr_assort , aes(x = tissue, y = diffraction, fill = cond,
                             color = cond)) +
  geom_violin() +
  scale_fill_manual(values = color_pal) +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  ylab("Expression\nAssortativity") +
  xlab("") +
  theme(legend.position = "none", axis.text.x = element_blank(),
        plot.background = element_blank(), 
        strip.text.x = element_text(size = 22),
        strip.text.y = element_text(size = 20)) + 
  facet_grid(cond~tissue, scale= "free_x", labeller = labeller(tissue = capitalize))

png(paste0("pan-loss/network_aracne_plots/assortativity/expr_assortativity.png"), 
    width = 2000, height = 200)
print(p)
dev.off()
