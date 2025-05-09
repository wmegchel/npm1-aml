##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# UMAP of the scTapestri ADT colored by cell type
##################################################################################

# Clear workspace
rm(list = ls())

# Select patient 01. You can choose other patients from the Seurat list object
donor <- "P01"
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Rasterized UMAP with colors for each clone
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 2, group.by = "seurat_clusters") + 
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig4g_UMAP_by_celltype.pdf", width = 4, height = 4)
