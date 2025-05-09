##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Piechart of lineage commitment by at DX and RL
##################################################################################

library(Seurat)

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

meta <- as.data.table(so@meta.data, keep.rownames = "BC")
DT <- meta[, .N, by=.(commitment, stage)]
DT[, pct:=N/sum(N), by=stage]

DT[, commitment:=factor(commitment, levels = c("immature", "myeloid", "lymphoid"))]

DT[, sum(N), by=stage]
# Basic piechart
ggplot(DT, aes(x="", y=pct, fill=commitment)) +
geom_bar(stat="identity", width=1, color="white") +
coord_polar("y", start=0, direction = -1) +
facet_wrap(~stage) + 
scale_fill_manual(values = c("#CCCCCC", "#E0860E","#AE09B7" )) + 
theme_void()


# Save as PDF
ggsave("scTapestri/img/Fig3c_Piechart_commited_cells.pdf", width = 6, height = 4)

