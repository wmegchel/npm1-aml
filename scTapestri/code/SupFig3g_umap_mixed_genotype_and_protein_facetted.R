#########################################################################
#
# Wout Megchelenbrink
# May 05, 2025
#
#
#########################################################################

# Clear workspace
rm(list=ls())

# Read Seurat object
so <- readRDS("processed_data/scTapestri_90K_cells_unfiltered.rds")

# Define doublet types
types <- c("Singlet", "Mixed genotype", "Mixed protein")

# Get and annotate metadata
meta <- as.data.table(so@meta.data, keep.rownames = "BC")
meta[, type:="Singlet"]
meta[sample.genotype == "mixed genotype", type:="Mixed genotype"]
meta[cell.adt.profile == "doublet", type:="Mixed protein"]
meta[ , type:=factor(type, levels = types)]

# Get UMAP embedddings
umap <- as.data.table(so@reductions$umap@cell.embeddings, keep.rownames = "BC")

# Merge UMAP and doublet types
DT <- merge(umap, meta[, .(BC, library, type)], by="BC")

# Plot rasterized UMAP, facetted by "doublet type"
ggplot(DT, aes(x=umap_1, y=umap_2, color=library)) + 
geom_scattermore(pointsize = 2, pixels = c(320,320)) +
facet_wrap(~type) +
theme_minimal()

# Save PDF
ggsave("scTapestri/img/SupFig3g_umap_mixed_genotype_and_protein_facetted.pdf", width = 10, height = 5)  
