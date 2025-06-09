library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)

tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "skin", "testis", "thyroid","uterus")

getTotals <- function(data_frame, type) {
  data_frame %>%
    rowwise() %>%
    mutate(total = sum(c_across(where(is.double)))) %>%
    ungroup() %>%
    arrange(desc(total)) %>% 
    select(sample, total) %>%
    mutate(cond = type)
}

capitalize <- function(string) {
  substr(string, 1, 1) <- toupper(substr(string, 1, 1))
  string
}

all <- lapply(tissues, function(tissue) {
  cancer_cibersort_vals <- read_tsv(paste0(tissue, "/cancer_immune/cibersort.tsv"))
  normal_cibersort_vals <- read_tsv(paste0(tissue, "/normal_immune/cibersort.tsv"))
  
  cancer_totals <- getTotals(cancer_cibersort_vals, "cancer") %>%
    mutate(tissue = tissue)
  normal_totals <- getTotals(normal_cibersort_vals, "normal") %>%
    mutate(tissue = tissue)
  
  return(bind_rows(cancer_totals, normal_totals))
})
all <- bind_rows(all)

color_pal <- c("#e3a098", "#a32e27")
labels <- c( "Normal", "Cancer")

all$cond <- factor(all$cond, levels = c("normal", "cancer"), labels = labels)

p <- ggplot(all, aes(x = cond, y = total, color = cond, fill = cond)) +
  geom_violin() +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), label.x = 1.5, label.y = 0.2) +
  scale_fill_manual(name = "Condition", values = color_pal) +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  guides(color="none") +
  ylab("Total fraction") +
  xlab("") +
  facet_wrap(~ tissue, labeller = labeller(tissue = capitalize), ncol = 3)

png(paste0("pan-loss/immune/cibersorts_fraction_lm22.png"), 
    width = 800, height = 1200)
print(p)
dev.off()

all <- lapply(tissues, function(tissue) {
  cancer_cibersort_vals <- read_tsv(paste0(tissue, "/cibersort_tr4/cancer_fractions.tsv")) %>%
    select(sample, cd10, cd31, cd45) 
  normal_cibersort_vals <- read_tsv(paste0(tissue, "/cibersort_tr4/normal_fractions.tsv"))%>%
    select(sample, cd10, cd31, cd45) 
  
  cancer_totals <- getTotals(cancer_cibersort_vals, "cancer") %>%
    mutate(tissue = tissue)
  normal_totals <- getTotals(normal_cibersort_vals, "normal") %>%
    mutate(tissue = tissue)
  
  return(bind_rows(cancer_totals, normal_totals))
})
all <- bind_rows(all)
all$cond <- factor(all$cond, levels = c("normal", "cancer"), labels = labels)

ggplot(all, aes(x = cond, y = total, color = cond, fill = cond)) +
  geom_violin() +
 # stat_compare_means(aes(label = paste0("p = ", after_stat(p.format))), label.x = 1.5, label.y = 0.2) +
  scale_fill_manual(name = "Condition", values = color_pal) +
  scale_color_manual(values = color_pal) +
  theme_base(base_size = 20) +
  guides(color="none") +
  ylab("Total fraction") +
  xlab("") +
  facet_wrap(~ tissue, labeller = labeller(tissue = capitalize), ncol = 3)
