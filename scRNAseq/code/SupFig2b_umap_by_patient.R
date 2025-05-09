#############################################################################################
#
# Wout Megchelenbrink
# 
#############################################################################################

# Clear workspace
rm(list=ls())

# Load scRNA
so <- readRDS("processed_data/scRNA_uncorrected.rds")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 1, group.by = "patient") + 
theme_void()

# Save as PDF
ggsave("scRNAseq/img/SupFig4b_UMAP_by_patient.pdf", width = 4, height = 4)
