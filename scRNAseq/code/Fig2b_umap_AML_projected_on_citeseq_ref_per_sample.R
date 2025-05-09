#############################################################################################
#
# Wout Megchelenbrink
# May 08, 2025
#
# NPM1 data projected on the CITE-seq hematopoietic progenitor dataset
#############################################################################################

# Clear workspace
rm(list=ls())

# Read reference and query into Seurat
ref <- readRDS("processed_data/Zhang_et_al_hemato_progenitors_CITEseq_MNN_integrated_seurat_object.rds")
qry <- readRDS("processed_data/scRNA_projected_on_Zhang_et_al_hemato_progenitors_seurat_object.rds")

# Reference (gray background)
ref.DT <- as.data.table(ref@reductions$umap@cell.embeddings, keep.rownames="BC")
setnames(ref.DT, c("BC", "UMAP_1", "UMAP_2"))

# Query (DX=blue, RL=green)
aml.DT <- as.data.table(qry@reductions$ref.umap@cell.embeddings, keep.rownames="BC")
setnames(aml.DT, c("BC", "UMAP_1", "UMAP_2"))
meta <- as.data.table(qry@meta.data, keep.rownames = "BC")[, .(BC, donor, stage, ct1=predicted.cell_type_1, ct2=predicted.cell_type_2, ct3=predicted.cell_type_3, score2=predicted.cell_type_2.score)]
aml.DT <- merge(aml.DT, meta, by="BC")

# UMAP of projected NPM1 scRNA data with facets for each of the 10 donors
ggplot(aml.DT, aes(x=UMAP_1, y=UMAP_2, color = stage, fill=stage)) +
geom_scattermore(data=ref.DT, mapping = aes(x=UMAP_1, y=UMAP_2), pixels = c(320, 320), pointsize = 0, colour="#DDDDDD", fill="#DDDDDD")  +
geom_scattermore(pixels = c(320, 320), pointsize = 1) +
scale_color_manual(values = c("#4955F3", "#82B300")) +
facet_wrap(~donor, nrow = 2) + 
theme_void() 

# Save PDF
ggsave("scRNAseq/img/Fig2b_umap_aml_projected_on_citeseq_all_patients.pdf", width = 10, height = 8)


