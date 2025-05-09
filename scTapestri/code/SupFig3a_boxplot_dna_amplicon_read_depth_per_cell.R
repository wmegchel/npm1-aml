#############################################################################################
#
# Wout Megchelenbrink
# 
# May 05, 2025 
#
# Boxplot of depths per amplicon per scTapestri library
# The depth per amplicon is obtained from the scTapestri H5 files before any filtering
#############################################################################################

# Clear workspace
rm(list=ls())

# Load the depths per amplicon
DT <- fread("processed_data/scTapestri_depth_per_amplicon_per_library.tsv")

# Factor libraries in reverse order (flipped plot)
libs <- DT[, sort(unique(library))]
DT[, library:=factor(library, levels=rev(libs))]

# Make boxplots
ggplot(DT, aes(x=library, y=DP, fill=library)) + 
geom_boxplot() + 
geom_hline(yintercept = DT[, median(DP)], linetype = "dashed", color = "#555555") + 
scale_fill_manual(values = hue_pal()(length(libs)), breaks = libs) + 
coord_flip(ylim = c(0,500)) + 
theme_minimal()

# Save PDF
ggsave("sctapestri/img/boxplot_amplicon_read_depth_per_cell.pdf", width = 6, height = 2.5)
