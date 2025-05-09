##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Heatmap of the SNV, CNV and ADT for a selected patient (here patient P01)
##################################################################################

library(ggrastr)
library(ggpubr)
library(scales)

rm(list=ls())
source("scTapestri/code/lib/process_tapestri_data.R")

set.seed(42)

# Select patient 01. You can choose other patients from the Seurat list object
donor <- "P01"

selected.proteins <- c("CD3", "CD4", "CD8",  "CD19", "CD22",
                       "CD83", "CD138", "CD30", "CD10",
                       "CD11c", "CD44",
                       "CD163", "CD13", "HLA-DR", "CD141", "CD38", "CD71", "CD69", "CD123",
                       "CD33", "CD45", "CD25", "CD303", "CD117", "CD45RA", "CD34", "CD62L", "CD2", "CD90")

so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]

# Make heatmaps of the SNVs, CNVs and ADTs
hm.snv <- get.snv.heatmap(so)
hm.cnv <- get.cnv.heatmap(so)
hm.protein <- get.protein.heatmap(so, selected.proteins)

# Save PDF
pdf(file = sprintf("scTapestri/img/Fig4bcd_heatmap_%s.pdf", donor), width = 12, height = 14)
draw(hm.snv %v% hm.cnv %v% hm.protein, gap=unit(5, "mm"))
dev.off()
