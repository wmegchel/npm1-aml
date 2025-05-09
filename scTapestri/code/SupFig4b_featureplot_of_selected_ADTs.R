##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Customized featureplot of selected ADTs, facetted by ADT
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

# Selected antibody derived tags (ADTs)
adts <- c("CD14", "CD16",  "CD11b", "CD11c", "CD64") 

# Get the UMAP cell embeddings from the Seurat object
umap <- as.data.table(so@reductions$umap@cell.embeddings, keep.rownames = "BC")

# Get the scaled data from the Seurat object
scaled.data <- as.data.table(t(GetAssayData(so)[adts, ]), keep.rownames = "BC")
scaled.data <- melt.data.table(scaled.data, id.vars = "BC", variable.name = "ADT", value.name = "value")

# Merge UMAP coords and scaled values into one data.table
DT <- merge(umap, scaled.data, by="BC")

# Make the featureplot; facet_wrap by ADT
ggplot(DT, aes(x=umap_1, y=umap_2, color=value)) + 
geom_scattermore(  pixels = c(320,320),  pointsize = 0) +
scale_color_viridis_c(limits = c(-4,4), oob = squish) + # set all data to range [-4, 4] and set out of boundary (oob) values to these limits
facet_wrap(~ADT, nrow = 1) + 
theme_void()

# Save as PDF
ggsave("scTapestri/img/SupFig4b_featureplots_selected_ADTs.pdf", width = 10, height = 4)

