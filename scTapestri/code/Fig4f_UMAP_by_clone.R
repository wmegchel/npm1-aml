##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# UMAP of the scTapestri ADT colored by somatic (sub)clone
##################################################################################

# Clear workspace
rm(list = ls())

# Select patient 01. You can choose other patients from the Seurat list object
donor <- "P01"
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Rasterized UMAP with colors for each clone
clone.colors <- c("#666666", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#BCB815", "#A65628")
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 2, group.by = "clone", cols =  clone.colors) + 
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig4f_UMAP_by_clone.pdf", width = 4, height = 4)
