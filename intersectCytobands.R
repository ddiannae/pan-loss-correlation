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
  
  barcodes <- cytobands_count %>%
    mutate(present = 1) %>% distinct() %>%
    pivot_wider(id_cols = c(cytoband, chr), names_from = tissue, 
                values_from = present, values_fill = 0)
    
    for(t in tissues[!tissues %in% names(barcodes)]){
      barcodes <- barcodes %>% mutate("{t}" := 0)
    }
  
  barcodes %>%
    select(cytoband, chr, all_of(tissues)) %>%
    unite(col = "barcode", all_of(tissues), sep="") %>% 
    mutate(chr = factor(chr, levels = c(as.character(1:22), "X", "Y"))) %>%
    arrange(desc(str_count(barcode, "1")))   
  
}

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

barcodes_normal <- getSharedCytobands(tissues, "normal")
barcodes_cancer <- getSharedCytobands(tissues, "cancer")  
barcodes_cancer %>% write_tsv("pan-loss/cytobands/barcodes_cancer.tsv")
barcodes_normal %>% write_tsv("pan-loss/cytobands/barcodes_normal")

cytoband_annot <- read_tsv("pan-loss/input/Biomart_Ensembl80_GRCh38_p2_regions.tsv")
colnames(cytoband_annot) <- c("ensembl", "cytoband")

barcodes_normal <- barcodes_normal %>% 
    separate(barcode, into = tissues, sep = 1:length(tissues), convert = T)

barcodes_cancer <- barcodes_cancer %>%
    separate(barcode, into = tissues, sep = 1:length(tissues), convert = T)
  
getInteractionsInRegion <- function(is, vs, rw) {
  vertices_inr <- vs %>% 
    filter(cytoband == rw %>% pull(cytoband) & 
             chr == rw %>% pull(chr))
  
  interactions_inr <- is %>%
    filter(source_ensembl %in% (vertices_inr %>% pull(ensembl)) &
             target_ensembl %in% (vertices_inr %>% pull(ensembl)))
  
  return(interactions_inr)
}


all_regions <- lapply(tissues, function(tissue) {
  cancer_inter <- read_tsv(paste0(tissue, "/network_aracne/cancer-interactions-100000.tsv"))
  cancer_vert <- read_tsv(paste0(tissue, "/network_aracne/cancer-vertices-100000.tsv"))
  cancer_vert <- cancer_vert %>% inner_join(cytoband_annot)
  
  normal_inter <- read_tsv(paste0(tissue, "/network_aracne/normal-interactions-100000.tsv"))
  normal_vert <- read_tsv(paste0(tissue, "/network_aracne/normal-vertices-100000.tsv"))
  normal_vert <- normal_vert %>% inner_join(cytoband_annot)
  
  tbc <- barcodes_cancer %>%  filter(get(tissue) == 1)
  
  in_region <- parallel::mclapply(X = 1:nrow(tbc), mc.cores = 20, FUN = function(i) {
    row <- tbc[i, ]
    c_in_region <- getInteractionsInRegion(cancer_inter, cancer_vert, row)
    n_in_region <- getInteractionsInRegion(normal_inter, normal_vert, row)
    
    region_in_normal <- barcodes_normal %>% filter(chr == row %>% pull(chr) &
                                                     cytoband == row %>% pull(cytoband) &
                                                     get(tissue) == 1)
    if(nrow(region_in_normal) > 0) {
      c_in_region$in_normal <- TRUE
    } else {
      c_in_region$in_normal <- FALSE
    }
    c_in_region$chr = row %>% pull(chr)
    c_in_region$cytoband = row %>% pull(cytoband)
    n_in_region$chr = row %>% pull(chr)
    n_in_region$cytoband = row %>% pull(cytoband)
    return(list(cancer = c_in_region, normal = n_in_region))
  })
  cancer_regions <- bind_rows(lapply(in_region, "[[", "cancer"))
  normal_regions <- bind_rows(lapply(in_region, "[[", "normal"))
  cancer_regions$tissue <- tissue
  normal_regions$tissue <- tissue
  return(list(cancer = cancer_regions, normal = normal_regions))
})

cancer_regions <- bind_rows(lapply(all_regions, "[[", "cancer"))
normal_regions <- bind_rows(lapply(all_regions, "[[", "normal"))

cancer_regions %>% write_tsv("pan-loss/cytobands/cancer_interactions_in_cytobands.tsv")
normal_regions %>% write_tsv("pan-loss/cytobands/normal_interactions_in_cytobands.tsv")

count_cancer <- cancer_regions %>% 
  group_by(tissue,chr, cytoband) %>% 
  tally() %>% arrange(chr, cytoband, n) %>% ungroup()

all_intersections <-  apply(count_cancer %>% select(chr, cytoband) %>% unique(), 1,
                            function(region) {
  interactions_in_region <- cancer_regions %>% filter(chr == region[1] & 
                                                        cytoband == region[2]) %>%
    mutate(id_interaction = paste0(pmin(source_ensembl, target_ensembl), "_",
                                   pmax(source_ensembl, target_ensembl)))
  
  n_tissues <- length(unique(interactions_in_region %>% pull(tissue)))
  inter_intersection <- interactions_in_region %>% group_by(id_interaction) %>% 
    tally() %>% filter(n == n_tissues)
  
  interactions_in_region %>% 
    filter(id_interaction %in% (inter_intersection %>% pull(id_interaction)))
  
})
all_intersections <- bind_rows(all_intersections)
all_intersections %>% write_tsv("pan-loss/cytobands/cancer_interactions_in_intersections.tsv")

cytobands_count <- lapply(tissues, function(tissue) {
  file_intra_count <- paste0(tissue, "/distance_analysis/cancer-intra-interactions-by_cytoband-count.tsv")
  count <- read_tsv(file_intra_count)
  ## Only 1e5 because I only need total_inter and they're the same on every top
  count <- count %>% filter(top == 1e5)
  count$tissue <- tissue
  return(count)
})
cytobands_count <- bind_rows(cytobands_count)
cytobands_count %>% inner_join(count_cancer, by= c("cytoband", "chr", "tissue"))

cytobands_means <- cytobands_count %>% group_by(chr, cytoband) %>% 
  summarise(mean_inter = mean(total_inter), max_inter = max(total_inter))

inter_for_plot <- all_intersections %>% group_by(chr, cytoband, tissue) %>% tally() %>% 
  group_by(chr, cytoband) %>% summarise(n_tissues = n(), n_intersections = first(n)) %>%
  arrange(desc(n_tissues), desc(n_intersections)) %>% 
  inner_join(cytobands_means, by = c("chr", "cytoband")) %>%
  mutate(fraction = n_intersections/mean_inter) %>% filter(n_tissues > 1)

chrs <- c(as.character(1:22), "X", "Y")
cytobands <- unique(inter_for_plot$cytoband)
inter_for_plot <- inter_for_plot %>% mutate(chr = factor(chr, levels=chrs), 
                                          cytoband = factor(cytoband, levels = c(sort(cytobands[grep("p", cytobands)], decreasing = TRUE),
                                                                                 sort(cytobands[grep("q", cytobands)]))))
inter_for_plot %>% write_tsv("pan-loss/cytobands/intersections_summary.tsv")
cat("Building plot\n")
g <- ggplot(inter_for_plot, aes(x=cytoband, y=fraction, fill = n_tissues/15)) +
  geom_bar(stat="identity", width = 0.5) + 
  theme_base(base_size = 20) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  +
  facet_wrap(~chr, scales = "free_x", ncol = 2) +
  theme(panel.spacing.y = unit(1, "lines")) +
  ylab("Fraction of total interacions in region") +
  xlab("") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, 
                       name = "Fraction of tissues") 

png("pan-loss/cytobands/shared_cytobands.png", width = 1000, height = 2000)
print(g)
dev.off()
