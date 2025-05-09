#############################################################################################
#
# Wout Megchelenbrink
# 
# May 05, 2025 
#
# Barplots of cell counts for sample1, sample2 and doublets per library
#############################################################################################

# Clear workspace
rm(list=ls())

# Read data
demux <- fread("processed_data/scTapestri_90K_cell_barcodes_demultiplexed.tsv")

# Get percentages of sample1, sample2 and doublets per library
DT <- merge(demux[, .N, by=.(library, sample.1, sample.2, demultiplexed.sample)], 
            demux[, .N, by=.(library, sample.1, sample.2)],
            by=c("library", "sample.1", "sample.2"), 
            suffixes = c("", ".tot"))
DT[, pct:=N/N.tot * 100]

# Rename stuff
DT[demultiplexed.sample == sample.1, demultiplexed.sample := "sample 1"]
DT[demultiplexed.sample == sample.2, demultiplexed.sample := "sample 2"]
DT[demultiplexed.sample == "mixed", demultiplexed.sample := "doublet"]
DT[str_detect(demultiplexed.sample, "alloSCT"), demultiplexed.sample := "allo-SCT donor"]

# Factorize
DT[, demultiplexed.sample:=factor(demultiplexed.sample, levels=rev(c("sample 1", "sample 2", "allo-SCT donor", "doublet")))]
DT[, library:=factor(library, levels=rev(sort(unique(library))))]

# Make barplots
ggplot(DT, aes(x=library, y=N, fill=demultiplexed.sample)) + 
geom_bar(stat="identity", color = "white")  +
scale_fill_manual(values = c("#1B77B5", "#F07E1B", "#0A770A", "#D62528"), breaks = c("sample 1", "sample 2", "allo-SCT donor", "doublet")) +
scale_y_continuous(labels = comma) + 
coord_flip() + 
theme_minimal()

# Save PDF
ggsave("scTapestri/img/SupFig3d_barplot_cell_counts_and_doublets_per_library.pdf", width = 5, height = 4)
