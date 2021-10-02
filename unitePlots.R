library(grid)
library(png)
library(gridExtra)
library(ggplot2)

unitePlots <- function(plotname) {
  tissues <- c("bladder", "brain", "breast", "colorectal", "esophagus", 
               "kidney", "liver", "lung", "ovary", "pancreas", "prostate", 
               "testis", "thyroid","skin", "uterus")
  
  plot_tissues <- paste0(tissues, plotname)
  
  thePlots <- lapply (plot_tissues, function(pngFile) {
    rasterGrob(readPNG(pngFile, native = FALSE), interpolate = FALSE)
  })
  uniSize <- dim(thePlots[[1]]$raster)
  
  ggsave(paste0("pan-loss/", plotname), 
         arrangeGrob(grobs=thePlots, nrow=8, ncol=2), 
         height=uniSize[1]*8, width=uniSize[2]*2, units = "px")
}

unitePlots("/distance_plots/intra-inter-count-onek-bins.png")
unitePlots("/distance_plots/intra-inter-count-log-bins.png")
unitePlots("/distance_plots/bin-size-50.png")
unitePlots("/distance_plots/bin-distance-100000.png")

unitePlots("/network_aracne_plots/comm-diameter-boxplot-network-10000.png")
unitePlots("/network_aracne_plots/comm-diameter-histogram-network-10000.png")
unitePlots("/network_aracne_plots/comm-meandistance-boxplot-network-10000.png")
unitePlots("/network_aracne_plots/comm-meandistance-histogram-network-10000.png")
unitePlots("/network_aracne_plots/comm-density-boxplot-network-10000.png")
unitePlots("/network_aracne_plots/comm-density-histogram-network-10000.png")
unitePlots("/network_aracne_plots/comm-order-boxplot-network-10000.png")
unitePlots("/network_aracne_plots/comm-order-histogram-network-10000.png")
unitePlots("/network_aracne_plots/comm-size-boxplot-network-10000.png")
unitePlots("/network_aracne_plots/comm-size-histogram-network-10000.png")

unitePlots("/network_aracne_plots/degree-distribution-10000.png")
unitePlots("/network_aracne_plots/degree-distribution-cumulative-10000.png")
unitePlots("/network_aracne_plots/mi-density-network-100000.png")
unitePlots("/network_aracne_plots/mi-boxplot-network-100000.png")

unitePlots("/network_aracne_plots/normal_only-chr-inter-interactions-sparse-10000.png")
unitePlots("/network_aracne_plots/normal_only-chr-inter-interactions-10000.png")
unitePlots("/network_aracne_plots/cancer_only-chr-inter-interactions-sparse-10000.png")
unitePlots("/network_aracne_plots/cancer_only-chr-inter-interactions-10000.png")
unitePlots("/network_aracne_plots/shared-chr-inter-interactions-sparse-10000.png")
unitePlots("/network_aracne_plots/shared-chr-inter-interactions-10000.png")

unitePlots("/distance_plots_no_arsyn/intra-inter-count-onek-bins.png")
unitePlots("/distance_plots_no_arsyn/intra-inter-count-log-bins.png")
unitePlots("/distance_plots_no_arsyn/bin-size-50.png")
unitePlots("/distance_plots_no_arsyn/bin-distance-100000.png")



