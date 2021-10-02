library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggthemes)

getSharedRegions <- function(ts, cond) {
  windows <- lapply(ts, function(tissue) {
    file_intra_windows <- paste0(tissue, "/distance_analysis/", cond, "-intra-interactions-count.tsv")
    intra_windows <- read_tsv(file=file_intra_windows)
    intra_windows <- intra_windows %>% filter(top <= 1e5) %>% 
      group_by(chr, start, end, inter_total) %>% 
      summarise(inter = sum(inter), fraction = sum(fraction)) %>%
      filter(fraction >= 0.5)
    intra_windows$tissue <- tissue
    return(intra_windows)
  })
  windows <- bind_rows(windows)
  barcodes <- windows %>% group_by(start, end, chr) %>% 
    summarize(nts = n_distinct(tissue), tissue=tissue) %>%
    arrange(desc(nts)) %>% ungroup() %>% mutate(present=1) %>% distinct() %>%
    pivot_wider(id_cols = c(start, end, chr), names_from = tissue, 
                values_from = present, values_fill = 0) %>%
    unite(col = "barcode", all_of(tissues), sep="") %>% 
    mutate(chr = factor(chr, levels = c(as.character(1:22), "X", "Y"))) %>%
    arrange(desc(str_count(barcode, "1"))) 
  
  barcodes %>% write_tsv(paste0("pan-loss/", cond, "-intra-regions-barcodes.tsv"))
  return(barcodes)
}
 
readSharedRegions <- function(tissues, cond) {
  barcodes <-  read_tsv(paste0("pan-loss/", cond, "-intra-regions-barcodes.tsv")) %>%
    separate(barcode, into = tissues, sep = 1:length(tissues), convert = T) %>%
    unite(col = "id",  chr, start, end, sep = "_", remove = F)
}

getInteractionsInRegion <- function(is, vs, rw) {
  vertices_inr <- vs %>% 
    filter(start >= rw %>% pull(start) & 
             start <= rw %>% pull(end) & 
             chr == rw %>% pull(chr))
  
  interactions_inr <- is %>%
    filter(source_ensembl %in% (vertices_inr %>% pull(ensembl)) &
             target_ensembl %in% (vertices_inr %>% pull(ensembl)))
  
  return(interactions_inr)
}

### First time
# barcodes_normal <- getSharedRegions(tissues, "normal")
# barcodes_cancer <-  getSharedRegions(tissues, "cancer")
# 
# barcodes_normal <- barcodes_normal %>% 
#   separate(barcode, into = tissues, sep = 1:length(tissues), convert = T) %>%
#   unite(col = "id",  chr, start, end, sep = "_", remove = F)
# 
# barcodes_cancer <- barcodes_cancer %>% 
#   separate(barcode, into = tissues, sep = 1:length(tissues), convert = T) %>%
#   unite(col = "id",  chr, start, end, sep = "_", remove = F)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")

barcodes_normal <- readSharedRegions(tissues, "normal")
barcodes_cancer <- readSharedRegions(tissues, "cancer")

all_regions <- lapply(tissues, function(tissue) {
  cancer_inter <- read_tsv(paste0(tissue, "/network_aracne/cancer-interactions-100000.tsv"))
  cancer_vert <- read_tsv(paste0(tissue, "/network_aracne/cancer-vertices-100000.tsv"))
  
  normal_inter <- read_tsv(paste0(tissue, "/network_aracne/normal-interactions-100000.tsv"))
  normal_vert <- read_tsv(paste0(tissue, "/network_aracne/normal-vertices-100000.tsv"))
  
  tbc <- barcodes_cancer %>%  filter(get(tissue) == 1)
  
  in_region <- parallel::mclapply(X = 1:nrow(tbc), mc.cores = 20, FUN = function(i) {
    row <- tbc[i, ]
    c_in_region <- getInteractionsInRegion(cancer_inter, cancer_vert, row)
    n_in_region <- getInteractionsInRegion(normal_inter, normal_vert, row)
    
    region_in_normal <- barcodes_normal %>% filter(id == row %>% pull(id) &
                                 get(tissue) == 1)
    if(nrow(region_in_normal) > 0) {
      c_in_region$in_normal <- TRUE
    } else {
      c_in_region$in_normal <- FALSE
    }
    c_in_region$id <- row %>% pull(id)
    n_in_region$id <- row %>% pull(id)
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

cancer_regions_in_normal <- cancer_regions %>% filter(in_normal == TRUE)

count_cancer <- cancer_regions %>% 
  group_by(tissue, id) %>% 
  tally() %>% arrange(id, n) %>% ungroup()

count_normal <- normal_regions %>% 
  group_by(tissue, id) %>% 
  tally() %>% arrange(id, n) %>% ungroup()

cancer_mas3 <- 
  count_cancer %>% group_by(id) %>% filter(n() > 1 & all(n > 3))

sum(cancer_regions %>% filter(id %in% (cancer_mas3 %>%pull(id))) %>% pull(in_normal))

length(unique(cancer_mas3$id))

normal_mas3 <- 
  count_normal %>% ungroup() %>% group_by(id) %>% filter(n() > 1 & all(n > 3))

length(unique(normal_mas3$id))

cancer_mas3 %>% group_by(id) %>% tally(sort = T)

all_intersections <-  lapply(cancer_mas3 %>% pull(id) %>% unique(), function(region) {
 
   interactions_in_region <- cancer_regions %>% filter(id == region) %>%
    mutate(id_interaction = paste0(pmin(source_ensembl, target_ensembl), "_",
                                    pmax(source_ensembl, target_ensembl)))
    
  n_tissues <- length(unique(interactions_in_region %>% pull(tissue)))
  inter_intersection <- interactions_in_region %>% group_by(id_interaction) %>% 
    tally() %>% filter(n == n_tissues)
  
  interactions_in_region %>% 
    filter(id_interaction %in% (inter_intersection %>% pull(id_interaction)))
  
  })
all_intersections <- bind_rows(all_intersections)

inter_for_plot <- all_intersections %>% group_by(id, tissue) %>% tally() %>% 
  group_by(id) %>% summarise(n_tissues = n(), n_intersections = first(n)) %>%
  arrange(desc(n_tissues), desc(n_intersections)) %>% 
  separate(col="id", into = c("chr", "start", "end"), sep = "_") %>%
  mutate(start = as.numeric(start), end = as.numeric(end), 
         chr = factor(chr, levels = rev(c(as.character(1:22), "X", "Y"))),
         frac_tissues = n_tissues/15)
inter_top <- inter_for_plot %>% top_n(16)

ggplot(inter_top) + 
  geom_segment(aes(x = chr, 
                   xend = chr, 
                   y = start, 
                   yend= end, 
                   size = n_intersections,
                   color = frac_tissues)) +
  coord_flip() + 
  theme_base()
  





reg18q112 <- all_intersections %>% filter(id == "18_20800000_21800000")
reg18q112 %>% group_by(tissue) %>% tally()
reg18q112 %>% filter(tissue == "bladder")
reg8q243 <- all_intersections %>% filter(id == "8_1.44e+08_1.45e+08")
reg8q243 %>% group_by(tissue) %>% tally()


uniq_cancermas3 <- unique(cancer_mas3$id)
uniq_intersections <- all_intersections %>% pull(id) %>% unique()

which(!uniq_cancermas3 %in% uniq_intersections )
uniq_cancermas3[201]
reg7 <- cancer_regions %>% filter(id == "7_26400000_27400000")
reg7 %>% group_by(tissue) %>% tally()

reg7 %>% filter(tissue == "ovary")
regions_by_tissue <- all_intersections %>% group_by(id) %>% select(tissue) %>% distinct()%>% 
  group_by(id) %>% tally(sort = T)
