############################################
#
# Wout Megchelenbrink 
# May 08, 2025
#
#
############################################

# Clear workspace
rm(list=ls())

# Load unfiltered scTapestri data
so <- readRDS("processed_data/scTapestri_90K_cells_unfiltered.rds")

### No need to rerun, but these settings were used for dim reduc.
# VariableFeatures(so) <- rownames(so)
# so <- ScaleData(so)
# so <- RunPCA(so, npcs = 15)
# so <- RunUMAP(so, dims=1:15)
# so <- FindNeighbors(so, dims = 1:15)
# so <- FindClusters(so, resolution = .4)

meta <- as.data.table(so@meta.data)
meta[, status:="Singlet"]
meta[sample.genotype == "mixed genotype", status:="Mixed genotype"]
meta[cell.adt.profile == "doublet", status:="Mixed protein"]

so <- AddMetaData(so, meta$status, col.name = 'doublet_status')
DimPlot(so, group.by = "doublet_status")

types <- c("Singlet", "Mixed genotype", "Mixed protein")

# Get percentages for barplot
DTX <- meta[, .(pct=.N/nrow(meta) * 100), by=status]
DTX[, type:=factor(status, levels=rev(types) )]

# Colors
cols <- c("#DDDDDD", "#D62528", "#06A1EF")

# Plot percentages
ggplot(DTX, aes(x=1, y=pct,  fill=status)) +
geom_bar(stat="identity", color="white") + 
scale_fill_manual(values = cols, breaks=types ) + 
coord_flip() +
theme_void()

# Save PDF
ggsave("scTapestri/img/SupFig3f_barplot_mixed_genotype_and_protein.pdf", width = 5, height = 1)

# Plot UMAP
so$doublet_status <- factor(so$doublet_status, levels = types)
DimPlot(so, raster = T, raster.dpi = c(320,320), pt.size = 1, group.by = "doublet_status") +
scale_color_manual(values = cols) +
theme_minimal()

# Save PDF
ggsave("scTapestri/img/SupFig3f_umap_mixed_genotype_and_protein.pdf", width = 5, height = 4)





