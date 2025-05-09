##################################################################################
# 
# Wout Megchelenbrink
# May 08, 2025
#
# UMAP of the scTapestri ADT colored by patient
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 0, group.by = "patient") + 
theme_void()

# Save as PDF
ggsave("scTapestri/img/SupFig4e_UMAP_by_patient.pdf", width = 5, height = 4)


