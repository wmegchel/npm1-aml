##################################################################################
# 
# Wout Megchelenbrink
# May 08, 2025
#
# UMAP of the scRNAseq for paient 01 colored by stage and cluster
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scRNA_uncorrected.rds")
sox <- subset(so, patient == "NPM1-01")

sox <- NormalizeData(sox) %>%
      FindVariableFeatures() %>%
      ScaleData() %>%
      RunPCA() %>%
      RunUMAP(dims=1:20) %>%
      FindNeighbors() %>%
      FindClusters(resolution = .3)

## Rasterized UMAP colored by stage (DX, RL)
DimPlot(sox, raster = T, raster.dpi = c(320,320), pt.size = 2, 
        group.by = "stage", cols = c("#4955F3", "#82B300")) + 
theme_void()

ggsave("scRNAseq/img/SupFig5a_P01_UMAP_by_stage.pdf", width = 5, height = 4)

## Rasterized UMAP colored by cluster
DimPlot(sox, raster = T, raster.dpi = c(320,320), pt.size = 2) +
theme_void()

ggsave("scRNAseq/img/SupFig5b_P01_UMAP_by_cluster.pdf", width = 5, height = 4)

