##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Heatmap of the SNV, CNV and ADT for a selected patient (here patient P08)
##################################################################################

rm(list=ls())
source("scTapestri/code/lib/process_tapestri_data.R")

# Select patient 02. You can choose other patients from the Seurat list object
donor <- "P08"

selected.proteins <- c("CD71", "CD44", "CD123", "CD117", "CD69", "CD141", "CD38",
                       "CD33", "CD45", "CD4", "CD45RA", "CD62L", "CD13", "CD83", "CD10",
                       "CD138", "CD30", "CD25", "CD163", "CD22", "CD303", "CD34",
                       "CD11c", "CD2", "HLA-DR", "CD3", "CD90", "CD19", "CD8")

so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Make heatmaps of the SNVs, CNVs and ADTs
hm.snv <- get.snv.heatmap(so)
#hm.cnv <- get.cnv.heatmap(so)
hm.protein <- get.protein.heatmap(so, selected.proteins)

# Save PDF
pdf(file = sprintf("scTapestri/img/SupFig6de_heatmap_%s.pdf", donor), width = 12, height = 14)
draw(hm.snv %v% hm.protein, gap=unit(5, "mm"))
dev.off()
