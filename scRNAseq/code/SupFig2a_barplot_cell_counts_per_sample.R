#############################################################################################
#
# Wout Megchelenbrink
# 
# 03 May 2025
# Bar plots of cell counts per patient sample (DX, RL)
#
#############################################################################################

# Clear workspace
rm(list=ls())

# Load scRNA
so <- readRDS("processed_data/scRNA_uncorrected.rds")
meta <- as.data.table(so@meta.data, keep.rownames = "BC")

DT <- meta[, .N, by=.(patient, stage)]
DT[, patient:=factor(patient, levels=rev(unique(patient)))]

# Plot
ggplot(DT, aes(x=patient, y=N, fill=stage)) + 
geom_bar(stat = "identity", position = "dodge") + 
scale_fill_manual(values = c("#4955F3", "#82B300")) +
coord_flip() + 
theme_minimal()

# Save PDF
ggsave("scRNAseq/img/SupFig2a_barplot_cell_counts_per_sample.pdf", width = 5, height = 3)
