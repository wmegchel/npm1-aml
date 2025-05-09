##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Barplots that show per clone the relative abundance of DX/RL cells and cell types
##################################################################################

library(Seurat)
library(scales)

# Clear workspace
rm(list = ls())

# Load data
donor <- "P01"
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Get the UMAP cell embeddings from the Seurat object
meta <- as.data.table(so@meta.data, keep.rownames = "BC")

DT <- merge(meta[, .N, by=.(clone, stage=sample_name)],
            meta[, .N, by=.(clone)], 
            by="clone",
            suffixes = c("", ".tot"))

DT[, pct:=N/N.tot]


# Reverse order of clones, because we will flip the bars horizontally
DT[, clone:=factor(clone, levels = 8:1)]

# Create barplot
ggplot(DT, aes(x=clone, y=pct, fill=stage)) + 
geom_bar(stat = "identity") + 
scale_fill_manual(values = c("#4955F3", '#82B300')) +
scale_y_continuous(label=percent) +
coord_flip() + 
theme_minimal()

# Save as PDF
ggsave("scTapestri/img/Fig4i_Barplots_relative_nr_of_cells_per_stage_by_clone.pdf", width = 6, height = 4)

