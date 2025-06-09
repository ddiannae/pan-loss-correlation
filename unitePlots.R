########################################################################
## Script to get one single plot from plots of all tissues 
#######################################################################

library(grid)
library(png)
library(gridExtra)
library(ggplot2)

unitePlots <- function(plotname) {
  tissues <- c("bladder", "brain", "colon", "esophagus", 
                "liver", "ovary", "pancreas", "prostate", 
               "testis", "thyroid")
  
  # tissues <- c("bladder", "brain", "breast", "colon", "esophagus", 
  #              "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
  #              "testis", "thyroid","skin", "uterus")
  # 
  plot_tissues <- paste0(tissues, plotname)
  
  thePlots <- lapply (plot_tissues, function(pngFile) {
    rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
  })
  uniSize <- dim(thePlots[[1]]$raster)
  
  ggsave(paste0("pan-loss/", plotname), 
         arrangeGrob(grobs=thePlots, nrow=8, ncol=2), 
         height=uniSize[1]*8, width=uniSize[2]*2, units = "px")
}

unitePlots("/distance_plots_tpm/intra-inter-count-onek-bins.png")
unitePlots("/distance_plots_tpm/intra-inter-count-log-bins.png")
unitePlots("/distance_plots_tpm/intra-inter-count-onek-chunks.png")
unitePlots("/distance_plots_tpm/intra-inter-count-log-chunks.png")

unitePlots("/distance_plots_tpm/bin-size-1000-mean.png")
unitePlots("/distance_plots_tpm/bin-size-1000-median.png")
unitePlots("/distance_plots_tpm/bin-distance-100000-mean.png")
unitePlots("/distance_plots_tpm/bin-size-1000-mean_fitted.png")
unitePlots("/distance_plots_tpm/bin-distance-100000-mean_fitted.png")
# unitePlots("/distance_plots_tpm/heatmap-bins-size-all-ttests-1000.png")

unitePlots("/networks_tpm_plots/degree-distribution-100000.png")
unitePlots("/networks_tpm_plots/degree-distribution-cumulative-100000.png")
unitePlots("/networks_tpm_plots/mi-density-network-100000.png")
unitePlots("/networks_tpm_plots/mi-boxplot-network-100000.png")

unitePlots("/networks_tpm_plots/assortativity/cancer-comm-assort-enrich-100000.png")
unitePlots("/networks_tpm_plots/assortativity/normal-comm-assort-enrich-100000.png")
unitePlots("/networks_deseq2_plots/assortativity/cancer-comm-assort-enrich-10000.png")
unitePlots("/networks_deseq2_plots/assortativity/normal-comm-assort-enrich-10000.png")

