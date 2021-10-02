library(readr)
library(dplyr)
library(tidyr)
library(igraph)
library(stringr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

barcodes_cancer <- read_tsv("pan-loss/cytobands/barcodes_cancer.tsv") %>%
  filter(str_count(barcode, "1") > 1) %>%
  unite(col = "region", chr, cytoband, sep = "-")
  
cis <- read_tsv("pan-loss/cytobands/cancer_interactions_in_intersections.tsv")
cis <- cis %>% select(source_ensembl, chr, cytoband, tissue) %>%
  rename("ensembl_id" = "source_ensembl") %>% bind_rows(
    cis %>% select(target_ensembl, chr, cytoband, tissue) %>%
      rename("ensembl_id" = "target_ensembl")) %>% distinct() %>%
  unite(col = "region", chr, cytoband, sep = "-")

deg_cis <- lapply(tissues, function(tissue) {
  files <- dir(path =  paste0(tissue, "/deg"), pattern ="*.tsv", 
              full.names = T)
  deg <- read_tsv(files[1]) %>% select(-chr)
  deg$tissue <- tissue
  cis %>% inner_join(deg, by = c("ensembl_id", "tissue"))
})
deg_cis <- bind_rows(deg_cis)

deg_cis <- deg_cis %>% group_by(region, tissue) %>% 
  summarise(mean_lfc = mean(log_fc), mean_adj_p = mean(adj_p_val))

cytobands_count <- lapply(tissues, function(tissue) {
  file_intra_count <- paste0(tissue, "/distance_analysis/cancer-intra-interactions-by_cytoband-count.tsv")
  count <- read_tsv(file_intra_count)
  count <- count %>% filter(top <= 1e5)
  count$tissue <- tissue
  return(count)
})
cytobands_count <- bind_rows(cytobands_count) %>%
  unite(col = "region", chr, cytoband, sep = "-") %>%
  group_by(region, tissue) %>%
  summarise(n_inter = sum(n_inter), n_genes = first(n_genes), 
            total_inter = first(total_inter), top = max(top)) %>%
  mutate(fraction = n_inter/total_inter) 

deg_cis <- deg_cis %>% 
  inner_join(cytobands_count, by=c("region","tissue")) %>%
  inner_join(barcodes_cancer, by="region") 

deg_cis <- deg_cis %>% separate(col = region, into=c("chr", "cytoband"), sep = "-")
chrs <- c(as.character(1:22), "X", "Y")
cytobands <- unique(deg_cis$cytoband)
deg_cis <- deg_cis %>% mutate(chr = factor(chr, levels=chrs), 
                                            cytoband = factor(cytoband, levels = c(sort(cytobands[grep("p", cytobands)], decreasing = TRUE),

                                                                                                                                                                     sort(cytobands[grep("q", cytobands)]))))

deg_cis_min <- deg_cis %>% filter(chr %in% c("1", "2"))
library(ggplot2)
library(ggthemes)

g <- ggplot(deg_cis, aes(x=cytoband, y= tissue)) +
  geom_point(aes(size = fraction, fill = mean_lfc), stroke = 0.5, 
             shape = 21, color = "grey45") + 
  theme_base(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_wrap(~chr, scales = "free_x", ncol = 3) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("") +
  xlab("") +
  scale_fill_distiller(palette = "RdBu", direction = -1, name = "LFC") +
  scale_size(name = "Fraction of interactions", breaks = c(0.5, 0.75, 1))

png("pan-loss/cytobands/shared_cytobands_lfc.png", width = 1000, height = 3000)
print(g)
dev.off()
