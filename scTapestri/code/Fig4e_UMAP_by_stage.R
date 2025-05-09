##################################################################################
# 
# Wout Megchelenbrink
# May 08, 2025
#
# UMAP of the scTapestri ADT for patient 01 colored by stage (DX, RL)
##################################################################################

# Clear workspace
rm(list=ls())

# Select patient 01. You can choose other patients from the Seurat list object
donor <- "P01"
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 2, group.by = "sample_name") + 
scale_color_manual(values = c("#4955F3", '#82B300')) +
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig4e_UMAP_by_stage.pdf", width = 4, height = 4)
