##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Barplot of the cell composition per patient per clone at DX vs RL
##################################################################################

# Clear workspace
rm(list=ls())

# Read data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")
meta <- as.data.table(so@meta.data, keep.rownames = "BC")
meta <- meta[!patient %in% c("P06D", "P09D")]

# Get the immature cell types
immature.cts <- meta[commitment == "immature", sort(unique(cell.type))]
meta[!cell.type %in% immature.cts, cell.type:="committed"]

# Compute percentage per clone
DTX <- merge(meta[, .N, by=.(patient, stage, patient.clone, cell.type)], meta[, .N, by=.(patient, stage)], by=c("patient", "stage"), suffixes = c("", ".tot"))
DTX[, pct:=N/N.tot * 100]
DTX[, patient:=factor(patient, levels=c("P01", "P02", "P03", "P05", "P08","P10"))]
DTX[stage == "Rl", pct:=pct * -1]

# Since clone 3 was removed (unreliable), renumber clones for patient 10
DTX[patient == "P10" & patient.clone > 2, patient.clone:=patient.clone -1]

# Set colors
DTX[, cell.type:=factor(cell.type, levels=c("committed", "early GMP", "late GMP", "CD56+ pre-neu.", "CD10+ LMPP", "CD25+ LMPP"))]
cols <- c("#DDDDDD", hue_pal()(5))

# Make baplot
ggplot(DTX, aes(x=patient.clone, y=pct, fill=cell.type)) +
geom_bar(stat="identity") +
scale_fill_manual(values = cols) +
facet_wrap(~patient, nrow = 1) 

#+ theme_void() 

# Save PDF
ggsave("scTapestri/img/SupFig6a_Barplot_cell_compositon_per_clone_dx_rl.pdf", width = 12, height = 6)
