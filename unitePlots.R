########################################################################
## Script to get one single plot from plots of all tissues 
#######################################################################

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
unitePlots("/distance_plots/intra-inter-count-onek-chunks.png")
unitePlots("/distance_plots/intra-inter-count-log-chunks.png")

unitePlots("/distance_plots/bin-size-1000-mean.png")
unitePlots("/distance_plots/bin-distance-100000-mean.png")
unitePlots("/distance_plots/bin-size-1000-mean_fitted.png")
unitePlots("/distance_plots/bin-distance-100000-mean_fitted.png")
unitePlots("/distance_plots/heatmap-bins-size-all-ttests-1000.png")

##### Changed to /network_aracne_plots/communities
unitePlots("/network_aracne_plots/communities/comm-diameter-boxplot-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-diameter-histogram-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-meandistance-boxplot-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-meandistance-histogram-network-intra-100000.png")

unitePlots("/network_aracne_plots/communities/comm-density-boxplot-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-density-histogram-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-order-boxplot-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-order-histogram-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-size-boxplot-network-intra-100000.png")
unitePlots("/network_aracne_plots/communities/comm-size-histogram-network-intra-100000.png")

unitePlots("/network_aracne_plots/degree-distribution-100000.png")
unitePlots("/network_aracne_plots/degree-distribution-cumulative-100000.png")
unitePlots("/network_aracne_plots/mi-density-network-100000.png")
unitePlots("/network_aracne_plots/mi-boxplot-network-100000.png")

# unitePlots("/network_aracne_plots/normal_only-chr-inter-interactions-sparse-10000.png")
# unitePlots("/network_aracne_plots/normal_only-chr-inter-interactions-10000.png")
# unitePlots("/network_aracne_plots/cancer_only-chr-inter-interactions-sparse-10000.png")
# unitePlots("/network_aracne_plots/cancer_only-chr-inter-interactions-10000.png")
# unitePlots("/network_aracne_plots/shared-chr-inter-interactions-sparse-10000.png")
# unitePlots("/network_aracne_plots/shared-chr-inter-interactions-10000.png")

unitePlots("/distance_plots_no_arsyn/intra-inter-count-onek-bins.png")
unitePlots("/distance_plots_no_arsyn/intra-inter-count-log-bins.png")
unitePlots("/distance_plots_no_arsyn/bin-size-50.png")
unitePlots("/distance_plots_no_arsyn/bin-distance-100000.png")

unitePlots("/network_aracne_plots/assortativity/cancer-comm-assort-enrich-100000.png")
unitePlots("/network_aracne_plots/assortativity/normal-comm-assort-enrich-100000.png")


