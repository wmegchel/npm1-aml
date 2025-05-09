#############################################################################################
#
# Wout Megchelenbrink
# May 08, 2025
#
# UMAP of the CITE-seq reference by Zhang et al. Nature Immunology (2024)
#############################################################################################

# Clear workspace
rm(list=ls())

# Load MNN integrated reference data
ref <- readRDS("processed_data/Zhang_et_al_hemato_progenitors_CITEseq_MNN_integrated_seurat_object.rds")

# Uncolor lymphopid cell types
meta <- as.data.table(ref@meta.data, keep.rownames = "BC")
meta[cell_type_2 %in% c("CLP", "pre-B", "Transitional-B", "B cell", "CD4 T-cell", "T/NK"), cell_type_2:=NA]

# Set colors such that adjacent clusters can be easily discriminated
cols <- c("HSC"="#D5377E",
          "MPP"="#008AD0",
          "MultiLin"= "#9456E0",
          "GMP"="#9F7100",
          "early-Neu"="#428900",
          "Ba/Ma/Eo"="#CB4D32",
          "MEP"="#B79F00",
          "MKP"="#00BDD2",
          "early-Erythroid"="#1DB700",
          "late-Erythroid"="#EC69EF",
          "Myeloid intermediate"= "#BD81FF",
          "early-Mono"="#86AC00",
          "Monocyte"="#00BE6D",
          "Mac"="#00B1F5",
          "pre-DC"="#00C1A5",
          "ASDC"="#619CFF",
          "cDC"="#F27C55",
          "MDP"="#FF61C8",
          "pDC"="#DA8F00",
          "Stroma"="#af0921",
          "Plasma Cell"="#a309af")

# Add cell type metadata
ref <- AddMetaData(ref, meta$cell_type_2, 'cell_type_2')
ref$cell_type_2 <- factor(ref$cell_type_2, levels = names(cols))

# Plot the UMAP
DimPlot(ref, group.by = "cell_type_2", cols = cols, raster = T, raster.dpi = c(320,320), pt.size = 0, na.value = "#DDDDDD") + 
theme_void()

# Save as PDF
ggsave("scRNAseq/img/Fig2a_CITEseq_reference_annotation_level_2_umap.pdf", width = 8, height = 6)

