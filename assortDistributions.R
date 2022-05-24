########################################################################
## Script to generate assortativity distribution plots. 
## Chromosomal assortativity per tissue and condition, expression 
## assortativity in cancer and enriched vs not enriched communities 
## in cancer with chromosomal and expression assortativity.
######################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
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

getSummaries <- function(cond, tss) {
  all_summaries <- lapply(tss, function(ts) {
    read_tsv(paste0(ts, "/network_aracne/assortativity/", cond, "-comm-summary-100000.tsv")) %>%
      mutate(tissue = ts, id = paste(tissue, community_id, sep="_")) %>%
      select(id, tissue, community_id, everything())
  })  
  return(bind_rows(all_summaries))
}

normal <- getSummaries("normal", tissues)
normal$cond <- "normal"
cancer <- getSummaries("cancer", tissues)
cancer$cond <- "cancer"
all <- bind_rows(normal, cancer)
all$cond <- factor(all$cond, levels = c("normal", "cancer"), labels = labels)

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

p <- ggplot(all , aes(x = tissue, y = chr_assortativity, fill = cond,
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


### Only Cancer and compare enriched vs not enriched communities
cancer$cond <- "Cancer"

p <- ggplot(cancer , aes(x = tissue, y = expr_assortativity, fill = cond,
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


cancer_enriched <- cancer %>%
  filter(size >= 5) %>%
  mutate(isEnriched = if_else(enriched_terms > 0, "go", "no_go" ),
         isEnriched = factor(isEnriched, levels = c("no_go", "go"), labels = c("No GO", "GO"), 
                             ordered = TRUE))

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", ""))

p <- ggplot(cancer_enriched,  aes(x=isEnriched, y=expr_assortativity, color=cond)) +
  geom_boxplot(size=1) + 
  ylab("Expression\nAssortativity") +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  xlab("") +
  theme(legend.position = "none") +
  stat_compare_means(aes(label = ..p.signif..),  label.x = 1.5, label.y = 1.1, symnum.args = symnum.args)+
  ylim(c(-1, 1.15)) +
  facet_grid(~tissue, scale= "free_x", labeller = labeller(tissue = capitalize))

png(paste0("pan-loss/network_aracne_plots/assortativity/expr_assort_enrichement.png"), 
    width = 1800, height = 250)
print(p)
dev.off()

p <- ggplot(cancer_enriched,  aes(x=isEnriched, y=chr_assortativity, color=cond)) +
  geom_boxplot(size=1) + 
  ylab("Chromosomal\nAssortativity") +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  xlab("")  +
  theme(legend.position = "none") +
  stat_compare_means(aes(label = ..p.signif..),  label.x = 1.5, label.y = 1.1, symnum.args = symnum.args)+
  ylim(c(-1, 1.15)) +
  facet_grid(~tissue, scale= "free_x", labeller = labeller(tissue = capitalize))

png(paste0("pan-loss/network_aracne_plots/assortativity/chr_assort_enrichement.png"), 
    width = 1800, height = 250)
print(p)
dev.off()
