##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# UMAP of the scTapestri ADT colored by stage (DX, RL)
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 0, group.by = "stage") + 
scale_color_manual(values = c("#4955F3", '#82B300')) +
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig3a_UMAP_by_stage.pdf", width = 4, height = 4)


