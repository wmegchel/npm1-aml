##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Barplots that show per clone the relative abundance of DX/RL cells and cell types
##################################################################################

library(Seurat)
library(scales)
library(forcats)

# Clear workspace
rm(list = ls())

# Load data
donor <- "P01"
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Get the UMAP cell embeddings from the Seurat object
meta <- as.data.table(so@meta.data, keep.rownames = "BC")

DT <- merge(meta[, .N, by=.(clone, cluster=seurat_clusters)],
            meta[, .N, by=.(clone)], 
            by="clone",
            suffixes = c("", ".tot"))
DT[, pct:=N/N.tot]

# Reverse order of clones and clusters, because we will flip the bars horizontally
clusters <- sort(unique(DT$cluster))

DT[, clone:=fct_rev(factor(clone))]
DT[, cluster:=fct_rev(cluster)]

# Plot
ggplot(DT, aes(x=clone, y=pct, fill = cluster)) + 
geom_bar(stat="identity", color="white") + 
coord_flip() + 
scale_fill_manual(values = hue_pal()(length(clusters)), breaks=clusters) +
scale_y_continuous(labels=percent) + 
theme_minimal()  

# Save as PDF
ggsave("scTapestri/img/Fig4i_Barplots_relative_nr_of_cells_per_celltype_by_clone.pdf", width = 6, height = 4)

