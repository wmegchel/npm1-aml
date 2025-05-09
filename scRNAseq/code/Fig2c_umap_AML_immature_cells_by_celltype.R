######################################################################################
#
# Wout Megchelenbrink
# May 08, 2025
#
# scRNA-seq UMAP of the stem- and progenitor population, reintegrated
# using RPCA integration, that better separates cell states compared to MNN
#
######################################################################################

# Clear worksapce
rm(list=ls())

# Read Seurat object
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")

# Get metadata, rename and factorize
meta <- as.data.table(so@meta.data, keep.rownames = "BC")
meta[, cell_type:=predicted.cell_type_2]
meta[cell_type == "early-Neu", cell_type:="GMP"]
meta[, cell_type:=factor(cell_type, levels=c("HSC", "MPP", "GMP", "MEP", "Myeloid intermediate"))]

so <- AddMetaData(so, meta$cell_type, 'cell_type')

# Set colors
cols <- c("#D5377E", "#008AD0", "#428900", "#B79F00", "#BD81FF")

# Plot the UMAP
DimPlot(so, group.by = "cell_type", cols = cols,
        raster = T, raster.dpi = c(320,320), pt.size = 0) + 
theme_void()

# Save PDF
ggsave("scRNAseq/img/Fig2c_HSPC_like_AML_cells.pdf", width = 6, height = 4)
