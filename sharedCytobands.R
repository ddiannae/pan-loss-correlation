library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

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

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

cyto_normal <- getSharedCytobands(tissues, "normal")
cyto_cancer <- getSharedCytobands(tissues, "cancer")  

counts_normal <- cyto_normal %>% group_by(chr, cytoband) %>% tally(name = "ts", sort = T)
counts_cancer <- cyto_cancer %>% group_by(chr, cytoband) %>% tally(name = "ts", sort = T)

cyto_normal_mas5 <- cyto_normal %>% semi_join(counts_normal %>% filter(ts >= 5), 
                                              by = c("chr", "cytoband"))

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

counts_cancer_mas5 <- cyto_cancer %>% semi_join(counts_cancer %>% filter(ts >= 5), 
                                                by = c("chr", "cytoband"))
chrs <- c(as.character(1:22), "X", "Y")
cytobands <- unique(counts_cancer_mas5$cytoband)
tissues_labels <- str_to_title(tissues)

counts_cancer_mas5 <- counts_cancer_mas5 %>% 
  mutate(chr = factor(chr, levels=chrs), 
        cytoband = factor(cytoband, levels = c(sort(cytobands[grep("p", cytobands)], decreasing = TRUE),
                                                sort(cytobands[grep("q", cytobands)]))),
        tissue = factor(str_to_title(tissue), levels = rev(tissues_labels)))


g <- ggplot(counts_cancer_mas5, aes(x=cytoband, y= tissue)) +
  geom_point(aes(size = total_inter, fill = fraction),stroke = 0.5, 
             shape = 21, color = "grey45") + 
  theme_base(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_wrap(~chr, scales = "free_x", ncol = 3) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("") +
  xlab("") +
  scale_fill_distiller(palette = "Oranges", direction = 1, 
                        name = "Fraction of \ninteractions in top 100k") +
  scale_size(name = "Number of interactions")

png("pan-loss/cytobands/shared_cytobands.png", width = 1000, height = 2000)
print(g)
dev.off()
