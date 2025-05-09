#############################################################################################
#
# Wout Megchelenbrink
# 
# May 05, 2025 
#
# Boxplot of average number of protein reads per cell per library
# The averages are the column means of the protein read matrix in the Tapestri H5 file
#############################################################################################

# Clear workspace
rm(list=ls())

# Load the proteins reads counts per ADT
DT <- fread("processed_data/scTapestri_mean_protein_reads_per_cell.tsv")

# Factorize
libs <- DT[, sort(unique(library))]
DT[, library:=factor(library, levels=rev(sort(unique(library))))]

# Make boxplot 
ggplot(DT, aes(x=library, y=reads, fill=library)) + 
geom_boxplot() + 
geom_hline(yintercept = DT[, median(reads)], linetype = "dashed", color = "#555555") + 
scale_fill_manual(values = hue_pal()(length(libs)), breaks = libs) + 
coord_flip(ylim = c(0,500)) + 
theme_minimal()

# Save PDF
ggsave("scTapestri/img/SupFig3b_boxplot_protein_read_depth_per_cell.pdf", width = 6, height = 2.5)
