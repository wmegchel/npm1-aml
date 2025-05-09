######################################################################################
#
# Wout Megchelenbrink
# March 22, 2025
# Boxplot of the scRNA-seq inferred cycle scores for GMP- and MPP-like cells
#
######################################################################################

library(scales)
library(data.table)
library(ggplot2)

# clear workspace
rm(list=ls())

# read data for progenitors
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")
meta <- as.data.table(so@meta.data, keep.rownames = "BC")

# rename early neturophils and GMPs to GMP for simplicity
meta[predicted.cell_type_2 %in% c("early-Neu", "GMP"), predicted.cell_type_2 := "GMP"]

# define cell cycle phases as "G1" or "G2M/S" for simplicity
meta[Phase == "G1", phase:="G1"]
meta[Phase != "G1", phase:="G2M/S"]

# discard cell types that are not GMP or MPP
meta <- meta[predicted.cell_type_2 %in% c("GMP", "MPP")]

# count fractions
DT <- merge(meta[, .N, by=.(cell.type=predicted.cell_type_2, phase, stage, donor)], 
            meta[, .N, by=.(cell.type=predicted.cell_type_2, stage, donor)], 
            by=c("cell.type", "stage", "donor"), suffixes = c("",".tot"))
DT[, pct:=N/N.tot]

# create boxplot
ggplot(DT, aes(x=cell.type, y=pct, fill=phase)) + 
geom_boxplot() + 
scale_fill_manual(values = c("#c77cff", "#00be67")) +
scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, by=.2), labels=percent) + 
theme_minimal()  

# save as PDF
ggsave("scRNAseq/img/Fig2e_boxplot_inferred_cell_cycle_scores.pdf", width = 4, height = 3)
