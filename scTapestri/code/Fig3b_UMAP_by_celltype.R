##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# UMAP of the scTapestri ADT colored by cell type
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

# Swap a few colors to better discriminate adjacent cell types
cols <- c("#F8766D", "#E68613", "#CD9600",
          "#ABA300", "#7CAE00","#0CB702",
          "#00BE67", "#00C19A", "#00BFC4",
          "#00B8E7", "#00A9FF", "#C77CFF",
          "#8494FF", "#ED68ED", "#FF61CC", "#FF68A1")

# Rasterized UMAP with colored by stage (DX, RL)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 0, group.by = "cell.type") +
scale_color_manual(values = cols) +
theme_void()

# Save as PDF
ggsave("scTapestri/img/Fig3b_UMAP_by_celltype.pdf", width = 4, height = 4)


