##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# UMAP of the scTapestri ADT colored by main clone (clone.label)
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 0, group.by = "clone.label") + 
scale_color_manual(values = c('#4daf4a',  '#935503', '#888888', '#8c21ea',  '#377eb8', '#f2c406', '#e41a1c')) +
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig3e_UMAP_by_clone_label.pdf", width = 4, height = 4)
