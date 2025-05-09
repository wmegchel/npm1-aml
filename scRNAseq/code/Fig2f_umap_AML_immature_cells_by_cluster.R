######################################################################################
#
# Wout Megchelenbrink
# March 22, 2025
#
# HSPC-like cells batch integrated using RPCA, reclustered
######################################################################################

# Clear workspace
rm(list=ls())

# read seurat object
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")

# Plot the UMAP
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 0) + 
theme_void()

# Save PDF
ggsave("scRNAseq/img/Fig2f_UMAP_HSPC_like_cells_RPCA_reclustered.pdf", width = 6, height = 4)
