##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Full annotation of the single cell Tapestri data
# - cluster
# - cell type
# - clone (WT, pre-leukemic, NPM1, NRAS, FLT3-TKD, FLT3-ITD, FLT3-LOH)
# - stage (DX, RL)
#
##################################################################################

library(ggrastr)
library(ggpubr)
library(scales)

rm(list=ls())

source("scTapestri/code/lib/process_tapestri_data.R")
set.seed(42)

so <- readRDS("processed_data/scTapestri_filtered_and_annotated_seurat_object.rds")
so[['ADT_snn_res.0.4']] <- NULL
so[['seurat_clusters']] <- NULL
so <- cluster.protein(so, res=.4, n.dims = 12)

# Mono-Mac
so <- FindSubCluster(so, 4, resolution = .15, graph.name = "ADT_snn", subcluster.name = "Mono_Mac")

# Split Erythrocytes and erythroblasts
so <- FindSubCluster(so, 7, resolution = .1, graph.name = "ADT_snn", subcluster.name = "Erythrocytes_Erythroblasts")

# Split CD8-T and NK cells
so <- FindSubCluster(so, 10, resolution = .1, graph.name = "ADT_snn", subcluster.name = "CD8_NK")

meta <- as.data.table(so@meta.data, keep.rownames = "BC")
meta[seurat_clusters == 4, seurat_clusters:=Mono_Mac]
meta[seurat_clusters == 7, seurat_clusters:=Erythrocytes_Erythroblasts]
meta[seurat_clusters == 10, seurat_clusters:=CD8_NK]
meta[, seurat_clusters:=factor(seurat_clusters,
                   levels=c("0", "1", "2", "3", "4_0", "4_1", "4_2", "5", "6", "7_0", "7_1", "8", "9", "10_0", "10_1", "11", "12"))]

so <- AddMetaData(so, meta$seurat_clusters, "seurat_clusters")
so <- SetIdent(so,value = so$seurat_clusters)
DimPlot(so, group.by = 'seurat_clusters')




############ T ools to annotate cell types ############
### Marker proteins for all clusters 
#so <- SetIdent(so,value = so$sub.cluster)
# m <- FindAllMarkers(so)
# # head(so@meta.data)
# markers <- as.data.table(m)
# markers[cluster == "11"][order(-avg_log2FC)]
# # markers[cluster == "1_1"][order(-avg_log2FC)]
# # markers[cluster == "1_2"][order(-avg_log2FC)]
# # 
# # DimPlot(so)
# m <- FindMarkers(so, ident.1 = "11", ident.2 = "12")
# 
# FeaturePlot(so, features = c("CD117", "CD34", "CD25", "CD10", "CD22", "CD19", "CD1c"), min.cutoff = "q5", max.cutoff = "q95", ncol = 3)

## Annotated cell types
meta <- as.data.table(so@meta.data, keep.rownames = "barcode")
meta[seurat_clusters == "0", cell.type := "early GMP"] 
meta[seurat_clusters == "1", cell.type := "late GMP"] 
meta[seurat_clusters == "2", cell.type := "CD25+ LMPP"]
meta[seurat_clusters == "3", cell.type := "CD10+ LMPP"]
meta[seurat_clusters == "4_0", cell.type := "CD14+ monocyte"]
meta[seurat_clusters == "4_1", cell.type := "Early monocyte"] 
meta[seurat_clusters == "4_2", cell.type := "Macrophage"] 
meta[seurat_clusters == "5", cell.type := "CD56+ pre-neu."]
meta[seurat_clusters == "6", cell.type := "pre-DC"]
meta[seurat_clusters == "7_0", cell.type := "Erythroblasts"]
meta[seurat_clusters == "7_1", cell.type := "Erythrocytes"]
meta[seurat_clusters == "8", cell.type := "CD16+ (pre-)neu."]
meta[seurat_clusters == "9", cell.type := "CD4-T"]
meta[seurat_clusters == "10_0", cell.type := "CD8-T"]
meta[seurat_clusters == "10_1", cell.type := "NK"]
meta[seurat_clusters == "11", cell.type := "B-cell"]
meta[seurat_clusters == "12", cell.type := "late GMP"]


## Annotate cell commitment into 'myeloid', 'lymphoid' and 'immature'
meta[, commitment:="myeloid"]
meta[cell.type %in% c("CD4-T", "CD8-T", "NK", "B-cell"), commitment:="lymphoid"]
meta[cell.type %in% c("early GMP", "late GMP", "CD56+ pre-neu.", "MDP", "CD25+ LMPP", "CD10+ LMPP"), commitment:="immature"]

# Annotate the NRAS clone for patient P05
meta[patient == "P05" & clone == 3, label:="NRAS"]

## Remove obsolete metadata
so[["SNV"]] <- NULL
so[["CNV"]] <- NULL
so[['Mono_Mac']] <- NULL
so[['Erythrocytes_Erythroblasts']] <- NULL
so[['CD8_NK']] <- NULL
so[['orig.ident']] <- NULL
so[['nCount_SNV']] <- NULL 
so[['nFeature_SNV']] <- NULL
so[['nCount_CNV']] <- NULL
so[['nFeature_CNV']] <- NULL
so[['original_sample_name']] <- NULL
so[['sample_name']] <- NULL
so[['label']] <- NULL
so[['clone']] <- NULL
so[['seurat_clusters']] <- NULL

## Annotate required metadata
so <- AddMetaData(so, meta$original_sample_name, 'sample')
so <- AddMetaData(so, meta$sample_name, 'stage')
so <- AddMetaData(so, meta$label, 'clone.label')
so <- AddMetaData(so, meta$clone, 'patient.clone')
so <- AddMetaData(so, meta$seurat_clusters, 'cluster')
so <- AddMetaData(so, meta$commitment, 'commitment')
so <- AddMetaData(so, meta$cell.type, 'cell.type')

## UMAPs by feature
DimPlot(so, group.by = 'cluster')
DimPlot(so, group.by = 'commitment')
DimPlot(so, group.by = 'cell.type')
DimPlot(so, group.by = 'clone.label')
DimPlot(so, group.by = 'stage')

#meta <- as.data.table(so@meta.data, keep.rownames = "BC")

so$clone.label <- factor(so$clone.label, 
                         levels = c("WT", "pre-leukemic", "NPM1", "NRAS", "FLT3-TKD", "FLT3-ITD", "FLT3-LOH"))


## Check FLT3-LOH cells
# DT <- as.data.table(so@meta.data, keep.rownames = "BC")
# DTX <- merge(DT[clone.label == "FLT3-LOH", .N, by=cluster], 
#              DT[, .N, by=cluster], by="cluster",suffixes = c("", ".tot"))
# DTX[, pct:=N/N.tot * 100]
# DTX[order(-pct)]

saveRDS(so, file = "processed_data/sctapestri_integrated.rds")













################# ALL BELOW CAN GO INTO SEPARATE FILES #####################



## SUp

## by commitment
ggplot(DT, aes(x=umap_1, y=umap_2, color=commitment)) + 
geom_point_rast(scale = .1, raster.dpi = 50) + 
scale_color_manual(values = c('#CCCCCC', '#00b6eb', '#f8766d')) +
theme_void()
ggsave("sctapestri/img/UMAP_by_commitment.pdf", width = 4, height =4)




ggsave("sctapestri/img/barplot_mutation_counts.pdf", width = 8, height = 1)


# by patient
ggplot(DT, aes(x=umap_1, y=umap_2,  color=patient)) +
geom_point_rast(scale = .1, raster.dpi = 50) +
theme_void()
ggsave("sctapestri/img/UMAP_by_patient.pdf", width = 4, height = 4)



DTX <- merge(DT[, .N, by=.(sample_name, commitment)], DT[, .N, by=.(sample_name)], by="sample_name", suffixes = c("", ".tot"))
DTX[, pct:=N/N.tot * 100]





# barplots

DTX <- merge(DT[, .N, by=.(sample_name, patient, cell.type)], DT[, .N, by=.(sample_name, patient)], by=c("sample_name", "patient"), suffixes = c("", ".tot"))
DTX[, pct:=N/N.tot * 100]
DTX[, patient:=factor(patient, levels=rev(c("P01", "P02", "P03", "P05", "P06D", "P08", "P09D", "P10")))]

ggplot(DTX, aes(x=patient, y=pct, fill=cell.type)) +
geom_bar(stat="identity", color="white") +
coord_flip() + 
facet_wrap(~sample_name) +
theme_minimal()

ggsave("sctapestri/img/barplot_celltype_per_patient.pdf", width = 8, height = 5)


ggsave("sctapestri/img/piechart_commitment_per_stage.pdf", width = 6, height = 4)





DTX <- merge(DT[, .N, by=.(clone, sample_name, patient)], DT[, .N, by=.(sample_name, patient)], by=c("sample_name", "patient"), suffixes = c("", ".tot"))
DTX[, pct:=N/N.tot]
DTX[patient == "P10" & clone >= 4, clone:=clone -1]

DTX <- DTX[!patient %in% c("P06D", "P09D")]
DTX[patient == "P05"]

DTX[sample_name == "Rl", pct:=pct * -1]

ggplot(DTX, aes(x=clone, y=pct, fill=sample_name)) + 
geom_bar(stat = "identity") + 
facet_wrap(~patient, ncol = 6, scales = "free_x") +
scale_y_continuous(labels=percent) +
theme_minimal() + 
scale_fill_manual(values = c("#4955F3", '#82B300'))

ggsave("sctapestri/img/clone_percentage_dx_rl.pdf", width = 10, height = 5)




DT[clone >= 4, clone:=999]
DT[, clone:=factor(clone)]

ggplot(DT, aes(x=umap_1, y=umap_2, color=clone)) + 
geom_point_rast(scale = .1, raster.dpi = 50) + 
theme_void()


sox





markers <- FindAllMarkers(sox)
m <- as.data.table(markers)
m[cluster == 0]
m[cluster == 1]
m[cluster == 2]
m[cluster == 3]
m[cluster == 4]
m[cluster == 5][order(-avg_log2FC)]

m[gene == "CD11"]

# by cluster
DT[, cluster:=factor(seurat_clusters)]
ggplot(DT, aes(x=umap_1, y=umap_2,  color=cluster)) +
geom_point_rast(scale = .1, raster.dpi = 50) +
theme_void()

FeaturePlot(sox, features = c("CD56", "CD16", "CD19", "CD22", "CD14", "HLA-DR"), min.cutoff = "q5", max.cutoff = "q95")

ggsave("sctapestri/img/UMAP_by_cluster.pdf", width = 4, height = 4)




# sox <- SetIdent(sox, value=sox$cluster)
# markers <- FindAllMarkers(sox)
# m <- as.data.table(markers)
# 
# m[cluster == "7_0"][order(-avg_log2FC)]
# FindMarkers(sox, ident.1 = "7_0", ident.2 = "7_1")
# m[cluster == "7_0"][order(-avg_log2FC)]
# DimPlot(sox, group.by = 'cluster')
