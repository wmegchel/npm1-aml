#############################################################################################
#
# Wout Megchelenbrink
# April 18, 2025
#
# Pie chart of the relative abundance of HSPC-like AML cells
# per cell type at DX and RL
#############################################################################################

library(data.table)
library(stringr)
library(ggplot2)
library(scales)
library(stringr)

rm(list=ls())

# read data
qry <- readRDS("processed_data/scRNA_projected_on_Zhang_et_al_hemato_progenitors_seurat_object.rds")

# Simplify cell types
aml.DT <- as.data.table(qry@meta.data, keep.rownames = "BC")[, .(BC, donor, stage, ct2=predicted.cell_type_2)] 
aml.DT[, cell_type:="other"]
aml.DT[ct2 %in% c("HSC", "MPP", "MEP"), cell_type:=ct2]
aml.DT[ct2 %in% c("early-Neu", "GMP"), cell_type:="GMP"]
aml.DT[ct2 %in% c("Myeloid intermediate"), cell_type:="Myeloid intermediate"]

# calculate relative abundance per cell type
DTX <- merge(aml.DT[, .N, by=.(cell_type, stage)], aml.DT[, .N, by=.(stage)], by=c("stage"), suffixes=c("", ".tot"))
DTX[, pct:=N/N.tot * 100]

# define cell type colors
cols <- c("HSC"="#D5377E", "MPP"="#008AD0", "GMP"="#428900", "MEP"="#B79F00",
          "Myeloid intermediate"= "#BD81FF", "other"="#DDDDDD")

# factorize and relevel factors
DT.col <- data.table(cell.type=names(cols), col=cols)
DT.col[, cell.type:=factor(names(cols), levels=names(cols))]
DT.col[, idx:=as.numeric(cell.type)]
DTX[, stage:=str_replace(stage, "RE", "RL")]
DTX[, stage:=factor(stage, levels=c("DX", "RL"))]
DTX[, cell_type:=factor(cell_type, levels=DT.col$cell.type)]

# plot pie charts and split into DX and RL
ggplot(DTX, aes(x="", y=pct, fill=cell_type)) +
geom_bar(stat="identity", width=1, color="white") +
coord_polar("y", start=0, direction = -1) +
scale_fill_manual(values = cols) +
facet_wrap(~stage) + 
theme_void() 

# save PDF
ggsave("scRNAseq/img/Fig2d_piecharts_HSPC_like_cells_DX_RL.pdf", width = 8, height = 4)
