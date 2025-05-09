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
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 1, group.by = "stage") + 
scale_color_manual(values = c("#4955F3", '#82B300')) +
theme_void()

# Save as PDF
ggsave("scRNAseq/img/SupFig4c_UMAP_by_stage.pdf", width = 4, height = 4)
