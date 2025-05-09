#############################################################################################
#
# Wout Megchelenbrink
# 
# May 05, 2025 

#############################################################################################

# Clear workspace
rm(list=ls())

# Load scRNA
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 01, group.by = "Phase") + 
theme_void()

# Save as PDF
ggsave("scRNAseq/img/SupFig4e_UMAP_by_cell_cycle_phase.pdf", width = 4, height = 4)
