library(readr)
library(dplyr)

tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
             "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
             "testis", "thyroid","skin", "uterus")


cancer_barcodes <- read_tsv("pan-loss/cancer-interaction-barcodes-10000.tsv")
normal_barcodes <- read_tsv("pan-loss/normal-interaction-barcodes-10000.tsv")

cancer_barcodes <- cancer_barcodes %>% filter(interaction_type == "Inter")
normal_barcodes <- normal_barcodes %>% filter(interaction_type == "Inter")

tbc <- cancer_barcodes %>% group_by(barcode) %>% tally() %>% arrange(desc(n))


### Prostate Pancreas
pros_pan_cancer <- cancer_barcodes %>% filter(barcode == "000000000110000")
pros_pan_normal <-  normal_barcodes %>% filter(barcode == "000000000110000")

pros_pan_shared <- pros_pan_cancer %>% filter(inter %in% pros_pan_normal$inter)
pros_pan_only_cancer <-  pros_pan_cancer %>% filter(!inter %in% pros_pan_normal$inter)

pros_pan_shared <- pros_pan_shared %>% separate(col=inter, into=c("source", "target"), sep="-")
pros_pan_shared_genes <- union(pros_pan_shared$source, pros_pan_shared$target)

pros_pan_only_cancer <- pros_pan_only_cancer %>% separate(col=inter, into=c("source", "target"), sep="-")
pros_pan_only_cancer_genes <- union(pros_pan_only_cancer$source, pros_pan_only_cancer$target)

pros_pan_only_cancer_genes %in% pros_pan_shared_genes


### Thiroid Kidney
### Pancreas Prostate
### Prostate Brain
### Pancreas Prostate Brain
### Pancreas Prostate Kidney
### Prostate Thyroid
### Pancreas Prostate Kidney Brain
### Prostate Kidney
### Thyroid Brain
## Thyroid Prostate

all_annot <- lapply(tissues, function(tissue){
  load(paste0(tissue, "/rdata/annot.RData"))
  return(annot)
})
all_annot <- bind_rows(all_annot)
all_annot <- all_annot %>% arrange(ensembl_id) %>%
  filter(!is.na(strand)) %>% distinct() %>%
  select(ensembl_id, chr, start, end, length, gene_name, gene_type)

save(all_annot, file="pan-loss/rdata/annot.RData", compress="xz")
