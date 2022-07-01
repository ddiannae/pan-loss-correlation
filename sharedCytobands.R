########################################################################
## Script to get hot-spot cytobands (cytobands with more than half of
## their total interactions present in the top 100 MI interactions),
## the p-values from the null model, and their associated figure. 
######################################################################

library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggthemes)

getSharedCytobands <- function(ts, cond) {
  cytobands_count <- lapply(ts, function(tissue) {
    file_intra_count <- paste0(tissue, "/distance_analysis/", cond, "-intra-interactions-by_cytoband-count.tsv")
    intra_count <- read_tsv(file=file_intra_count)
    intra_count <- intra_count %>% filter(top <= 1e5) %>% 
      group_by(chr, cytoband, total_inter) %>% 
      summarise(inter = sum(n_inter), fraction = sum(fraction)) %>%
      filter(fraction >= 0.5)
    intra_count$tissue <- tissue
    return(intra_count)
  })
  cytobands_count <- bind_rows(cytobands_count)
  return(cytobands_count)
}

getPValues <- function(ts, cond) {
  all_pvals <- lapply(ts, function(tissue) {
  
        t_pvals <- read_tsv(paste0(tissue, "/distance_analysis/", cond, "-intra-interactions-by_cytoband-null_model-100000.tsv"),
                          col_types = cols(chr = col_character()))
        if(nrow(t_pvals > 0)) {
          t_pvals %>%
            mutate(tissue = tissue)
        } else {
          NULL
        }
  })
  all_pvals <- bind_rows(all_pvals)
  return(all_pvals)
}

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

cyto_normal <- getSharedCytobands(tissues, "normal")
cyto_cancer <- getSharedCytobands(tissues, "cancer") 
pvals_normal <- getPValues(tissues, "normal")
pvals_cancer <- getPValues(tissues, "cancer")

counts_normal <- cyto_normal %>% group_by(chr, cytoband) %>% 
  tally(name = "ts", sort = T)
counts_cancer <- cyto_cancer %>% group_by(chr, cytoband) %>% 
  tally(name = "ts", sort = T)

counts_normal %>% filter(ts == 1) %>% 
  inner_join(cyto_normal, by = c("chr", "cytoband")) %>% select(-ts) %>%
  write_tsv("pan-loss/cytobands/unique_normal.tsv")

cyto_cancer_unique <- counts_cancer %>% filter(ts == 1) %>%
  inner_join(cyto_cancer, by = c("chr", "cytoband")) %>% select(-ts) %>%
  write_tsv("pan-loss/cytobands/unique_cancer.tsv")

cyto_normal_mas5 <- cyto_normal %>% 
  semi_join(counts_normal %>% filter(ts >= 5),by = c("chr", "cytoband")) %>%
  filter(total_inter > 1) %>% 
  inner_join(pvals_normal, by = c("chr", "cytoband", "tissue"))

# A tibble: 20 × 6
# Groups:   chr, cytoband [3]
# chr   cytoband total_inter inter fraction tissue    
# <chr> <chr>          <dbl> <dbl>    <dbl> <chr>     
# 1 Y     q11.221            3     3    1     brain     
# 2 2     p25.2              1     1    1     breast    
# 3 2     p25.2              1     1    1     colorectal
# 4 2     p25.2              1     1    1     esophagus 
# 5 Y     q11.221            3     2    0.667 esophagus 
# 6 Y     q11.223            1     1    1     esophagus 
# 7 2     p25.2              1     1    1     kidney    
# 8 Y     q11.221           10     5    0.5   kidney    
# 9 Y     q11.223            1     1    1     kidney    
# 10 Y     q11.221            6     6    1     liver     
# 11 Y     q11.223            1     1    1     liver     
# 12 Y     q11.221            3     2    0.667 lung      
# 13 Y     q11.221            3     2    0.667 pancreas  
# 14 Y     q11.223            1     1    1     pancreas  
# 15 2     p25.2              1     1    1     prostate  
# 16 2     p25.2              1     1    1     thyroid   
# 17 2     p25.2              1     1    1     skin      
# 18 Y     q11.221            3     2    0.667 skin      
# 19 Y     q11.223            1     1    1     skin      
# 20 2     p25.2              1     1    1     uterus

cyto_cancer_mas5 <- cyto_cancer %>%
  semi_join(counts_cancer %>% filter(ts >= 5), by = c("chr", "cytoband")) %>%
  filter(total_inter > 1) %>% 
  arrange(chr, cytoband) %>% 
  inner_join(pvals_cancer, by = c("chr", "cytoband", "tissue"))

chrs <- c(as.character(1:22), "X", "Y")
labels_chr <- paste0("Chr ", chrs)
cytobands <- unique(cyto_cancer_mas5$cytoband)
tissues_labels <- str_to_title(tissues)

cyto_cancer_mas5 <- cyto_cancer_mas5 %>% 
  mutate(chr = factor(chr, levels=chrs, labels = labels_chr), 
        cytoband = factor(cytoband, levels = c(sort(cytobands[grep("p", cytobands)], decreasing = TRUE),
                                                sort(cytobands[grep("q", cytobands)]))),
        tissue = factor(str_to_title(tissue), levels = rev(tissues_labels)))


g <- ggplot(cyto_cancer_mas5, aes(x=cytoband, y= tissue)) +
  geom_point(aes(size = total_inter, fill = fraction),stroke = 0.5, 
             shape = 21, color = "grey45") + 
  theme_base(base_size = 20) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5),
        strip.text.x = element_text(face = "bold"), plot.background=element_blank())  +
  facet_wrap(~chr, scales = "free_x", ncol = 5) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("") +
  xlab("") +
  scale_fill_distiller(palette = "Oranges", direction = 1, 
                        name = "Fraction of \ninteractions \nin top 100k") +
  scale_size(name = "Number of \ninteractions")

png("pan-loss/cytobands/shared_cytobands_cancer.png", width = 1200, height = 1100)
print(g)
dev.off()

cytobands <- unique(cyto_normal_mas5$cytoband)

cyto_normal_mas5 <- cyto_normal_mas5 %>% 
  mutate(chr = factor(chr, levels=chrs, labels=labels_chr), 
         cytoband = factor(cytoband, levels = c(sort(cytobands[grep("p", cytobands)], decreasing = TRUE),
                                                sort(cytobands[grep("q", cytobands)]))),
         tissue = factor(str_to_title(tissue), levels = rev(tissues_labels))) 

g <- ggplot(cyto_normal_mas5, aes(x=cytoband, y= tissue)) +
  geom_point(aes(size = total_inter, fill = fraction),stroke = 0.5, 
             shape = 21, color = "grey45") + 
  theme_base(base_size = 20) +
  theme(axis.text.x = element_text(vjust = 0.5, hjust=0.5),
        strip.text.x = element_text(face = "bold"), plot.background=element_blank())  +
  facet_wrap(~chr, scales = "free_x", ncol = 5) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("") +
  xlab("") +
  scale_fill_distiller(palette = "Oranges", direction = 1, 
                       name = "Fraction of \ninteractions \nin top 100k") +
  scale_size(name = "Number of \ninteractions")

png("pan-loss/cytobands/shared_cytobands_normal.png", width = 450, height = 250)
print(g)
dev.off()

