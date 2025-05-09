#############################################################################################
#
# Wout Megchelenbrink
# May 05 2025 
#
#
#############################################################################################

# Clear workspace
rm(list=ls())

# Define cell types
cell_types <- c("HSC", "MPP", "GMP", "MEP", "Myeloid intermediate", "Lineage committed")

# Load projected scRNA
so <- readRDS("processed_data/scRNA_projected_on_Zhang_et_al_hemato_progenitors_seurat_object.rds")

# Get meta data
meta <- as.data.table(so@meta.data, keep.rownames = "BC")

# Simplify cell types a bit
meta[predicted.cell_type_2 == "early-Neu", predicted.cell_type_2:="GMP"]
meta[predicted.cell_type_2 == "MKP", predicted.cell_type_2:="MEP"]
meta[!predicted.cell_type_2 %in% cell_types, predicted.cell_type_2:= "Lineage committed"]

# Compute percentage per cell type
DT <- merge(meta[, .N, by=.(cell_type=predicted.cell_type_2, donor, stage)],
            meta[, .N, by=.(donor, stage)], by=c("donor", "stage"),
            suffixes = c("", ".tot"))
DT[, pct:=N/N.tot]

# Factorize and order
DT[, stage:=factor(stage, levels=c("DX", "RE"))]
DT[, cell_type:=factor(cell_type, levels=rev(cell_types))]

# Plot
ggplot(DT, aes(x=stage, y=pct, fill = cell_type)) + 
geom_bar(stat = "identity", color = "white") + 
scale_y_continuous(labels=percent) + 
facet_wrap(~donor, nrow = 1) + 
scale_fill_manual(values = c("#D5377E","#008AD0", "#428900", "#B79F00", "#BD81FF", "#DDDDDD"),
                  breaks = cell_types) + 
theme_minimal() 
  
# Save PDF
ggsave("scRNAseq/img/SupFig4d_Barplots_cell_types_per_sample.pdf", width = 8, height = 4)
