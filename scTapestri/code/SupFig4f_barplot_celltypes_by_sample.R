##################################################################################
# 
# Wout Megchelenbrink
# May 08, 2025
#
# Barplot of 
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

DT <- as.data.table(so@meta.data)

DTX <- merge(DT[, .N, by=.(patient, stage, cell.type)],
             DT[, .N, by=.(patient, stage)], 
             by = c("patient", "stage"),
             suffixes = c("", ".tot"))

DTX[, pct:=N/N.tot]

# Plot
DTX[, patient:=factor(patient, levels=rev(sort(unique(patient))))]

ggplot(DTX, aes(x=patient, y=pct, fill=cell.type)) +
geom_bar(stat = "identity") +
facet_wrap(~stage) +
coord_flip() +
theme_minimal()

# Save PDF
ggsave("scTapestri/img/SupFig4f_barplot_cell_types_per_sample.pdf", width = 10, height = 5)






