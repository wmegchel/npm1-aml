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

selected.proteins <- c("CD69", "CD123", "CD117", "CD45RA", "CD38", "CD25", "CD34")
so <- readRDS("processed_data/scTapestri_filtered_and_annotated_Seurat_object_list.rds")[[donor]]
head(so@meta.data)

# Clone 4 in the paper is clone 5 here, because clone 1 = WT cells
sox <- subset(so, clone == 5)

# Make heatmaps of the SNVs, CNVs and ADTs
hm.protein <- get.protein.heatmap(sox, selected.proteins)

# Save PDF
pdf(file = sprintf("scTapestri/img/Fig4j_heatmap_%s_clone4.pdf", donor), width = 10, height = 5)
draw(hm.protein, gap=unit(5, "mm"))
dev.off()