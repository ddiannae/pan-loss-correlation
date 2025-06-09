library(vroom)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggthemes)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")
names(color_pal) <- labels

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

ciber_cond <- lapply(tissues, function(t) {
  nciber <- vroom(paste0(t, "/immune/normal_tr4_cibersort.tsv"))
  cciber <- vroom(paste0(t, "/immune/cancer_tr4_cibersort.tsv"))
  
  nciber$type <- "Normal"
  cciber$type <- "Cancer"
  
  ciber <- bind_rows(nciber, cciber)
  ciber %>%
    select(-rmse, -correlation) %>%
    pivot_longer(cols = c("cd10", "cd31", "cd45", "epcam", "p_value"), 
                 names_to = "cell", values_to = "score") %>%
    mutate(tissue = t)
  
})

ciber_cond <- bind_rows(ciber_cond)
ciber_cond$type <- factor(ciber_cond$type, levels = c("Normal", "Cancer"))

p <- ggplot(ciber_cond) +
  geom_boxplot(aes(x = cell, y = score, color = type)) + 
  scale_color_manual(name = "Condition", values = color_pal) +
  theme_base(base_size = 20) +
  facet_wrap(~tissue, scales = "free_x", nrow = 3, labeller = labeller(tissue = capitalize)) +
  xlab("") +
  ylab("Score") + 
  ggtitle("Cibersort results") + 
  scale_x_discrete(breaks=c("cd10", "cd31", "cd45", "epcam", "p_value"),
                   labels=c("CD10", "CD31", "CD45", "EPCAM", "P-val")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.text.x = element_text(size = 16))
  


png(paste0("pan-loss/immune/cibersort_scores_by_condition.png"), 
    width = 1200, height = 1000)
print(p)
dev.off()

ciber <- lapply(tissues, function(t) {
  ciber_all <- vroom(paste0(t, "/immune/all_tr4_cibersort.tsv"))
  
  ciber_all <- ciber_all %>%
    mutate(type = if_else(grepl("cancer", sample), "Cancer", "Normal")) %>%
    select(-rmse, -correlation) %>%
    pivot_longer(cols = c("cd10", "cd31", "cd45", "epcam", "p_value"), 
                 names_to = "cell", values_to = "score") %>%
    mutate(tissue = t)
  
})

ciber <- bind_rows(ciber)
ciber$type <- factor(ciber$type, levels = c("Normal", "Cancer"))

p <- ggplot(ciber) +
  geom_boxplot(aes(x = cell, y = score, color = type)) + 
  scale_color_manual(name = "Condition", values = color_pal) +
  theme_base(base_size = 20) +
  facet_wrap(~tissue, scales = "free_x", nrow = 3, labeller = labeller(tissue = capitalize)) +
  xlab("") +
  ylab("Score") + 
  ggtitle("Cibersort results") + 
  scale_x_discrete(breaks=c("cd10", "cd31", "cd45", "epcam", "p_value"),
                   labels=c("CD10", "CD31", "CD45", "EPCAM", "P-val")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
        strip.text.x = element_text(size = 16))

png(paste0("pan-loss/immune/cibersort_scores.png"), 
    width = 1200, height = 1000)
print(p)
dev.off()


xcell <- lapply(tissues, function(t) {
  
  xcell_all <- vroom(paste0(t, "/immune/all_tr4_xcell.tsv"))
  
  xcell_all %>%
    select(sample, immune_score, stroma_score, microenvironment_score) %>%
    pivot_longer(cols = c("immune_score", "stroma_score", "microenvironment_score"), 
                 names_to = "cell", values_to = "score") %>%
    mutate(tissue = t) %>%
    mutate(type = if_else(grepl("cancer", sample), "Cancer", "Normal")) 

})

xcell <- bind_rows(xcell)
xcell$type <- factor(xcell$type, levels = c("Normal", "Cancer"))

p <- ggplot(xcell) +
  geom_boxplot(aes(x = cell, y = score, color = type)) + 
  scale_color_manual(name = "Condition", values = color_pal) +
  theme_base(base_size = 20) +
  facet_wrap(~tissue, nrow = 3, labeller = labeller(tissue = capitalize)) +
  xlab("") +
  ylab("Score") + 
  ggtitle("xCell results") + 
  scale_x_discrete(breaks=c("immune_score", "microenvironment_score",  "stroma_score"),
                   labels=c("Immune", "Microenvironment", "Stroma")) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1), 
        strip.text.x = element_text(size = 16))

png(paste0("pan-loss/immune/xcell_scores.png"), 
    width = 1200, height = 1000)
print(p)
dev.off()


