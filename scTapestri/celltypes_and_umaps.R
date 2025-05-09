

library(ggrastr)
library(ggpubr)
library(scales)

rm(list=ls())

source("sctapestri/code/process_tapestri_data.R")
set.seed(42)

so <- readRDS("../d2/sctapestri/sctapestri.Rds")
doublets <- fread("sctapestri/doublets_annotated.tsv")

meta1 <- as.data.table(so@meta.data, keep.rownames = "barcode")
meta1 <- merge(meta1, doublets, by=c("barcode", "clone", "original_sample_name", "sample_name", "label"), sort=F)

so <- AddMetaData(so, meta1$is.doublet, "is_doublet")
sox <- subset(so, is_doublet==FALSE)

sox <- cluster.protein(sox, res=.3)
sox <- Seurat::FindSubCluster(sox, 7, resolution = .1, graph.name = "ADT_snn", subcluster.name = "c7")
sox <- Seurat::FindSubCluster(sox, 10, resolution = .1, graph.name = "ADT_snn", subcluster.name = "c10")
DimPlot(sox)

FeaturePlot(sox, features = c("CD303", "CD56"), min.cutoff = "q5", max.cutoff = "q95")


meta <- as.data.table(sox@meta.data, keep.rownames = "barcode")
meta[, cluster:=seurat_clusters]
meta[seurat_clusters == 7, cluster:=c7]
meta[seurat_clusters == 10, cluster:=c10]
sox <- AddMetaData(sox, meta$cluster, 'cluster')

meta <- as.data.table(sox@meta.data, keep.rownames = "BC")
meta[cluster == "0", cell.type := "GMP"]
meta[cluster == "1", cell.type := "CD25+ LMPP"]
meta[cluster == "2", cell.type := "CD10+ LMPP"]
meta[cluster == "3", cell.type := "CD56+ pre-neu."]
meta[cluster == "4", cell.type := "Mono"]
meta[cluster == "5", cell.type := "MDP"]
meta[cluster == "6", cell.type := "pre-DC"]
meta[cluster == "7_0", cell.type := "Erythroid"]
meta[cluster == "7_1", cell.type := "Erythroid"]
meta[cluster == "8", cell.type := "Mac"]
meta[cluster == "9", cell.type := "CD4-T"]
meta[cluster == "10_0", cell.type := "CD8-T"]
meta[cluster == "10_1", cell.type := "NK"]
meta[cluster == "11", cell.type := "B-cell"]

meta[, commitment:="myeloid"]
meta[cell.type %in% c("CD4-T", "CD8-T", "NK", "B-cell"), commitment:="lymphoid"]
meta[cell.type %in% c("GMP", "CD56+ pre-neu.", "MDP", "CD25+ LMPP", "CD10+ LMPP"), commitment:="immature"]


sox <- AddMetaData(sox, meta$commitment, 'commitment')
sox <- AddMetaData(sox, meta$cell.type, 'cell.type')
DimPlot(sox, group.by = 'commitment')

saveRDS(sox, file = "sctapestri/sctapestri_nov29.Rds")



########
sox <- readRDS("sctapestri/sctapestri_nov29.Rds")

sox <- cluster.protein(sox, res=.4)
DimPlot(sox)

FeaturePlot(sox, features = "CD138", min.cutoff = "q5", max.cutoff = "q95" )

m <- FindAllMarkers(sox)
markers <- as.data.table(m)
markers[cluster == 6][order(-avg_log2FC)]

umap <- as.data.table(sox@reductions$umap@cell.embeddings, keep.rownames = "barcode")
meta <- as.data.table(sox@meta.data, keep.rownames = "barcode")
DT <- merge(meta, umap, by="barcode")
DT[patient == "P05" & clone ==3, label:="NRAS"]


# DT[commitment == "myeloid", .N, by=label]
# 
# (128 + 336) / 10779 * 100


## by clone
#DT[label %in% c("NPM1", "pre-leukemic"), label:="Other"]
DT[, label:=factor(label, levels=c("WT", "pre-leukemic", "NPM1", "NRAS",  "FLT3-TKD", "FLT3-ITD", "FLT3-LOH"))]


## by stage
ggplot(DT, aes(x=umap_1, y=umap_2, color=sample_name)) + 
geom_point_rast(scale = .1, raster.dpi = 50) + 
scale_color_manual(values = c("#4955F3", '#82B300')) +
theme_void()
ggsave("sctapestri/img/UMAP_by_stage.pdf", width = 4, height =4)


## by commitment
ggplot(DT, aes(x=umap_1, y=umap_2, color=commitment)) + 
geom_point_rast(scale = .1, raster.dpi = 50) + 
scale_color_manual(values = c('#CCCCCC', '#00b6eb', '#f8766d')) +
theme_void()
ggsave("sctapestri/img/UMAP_by_commitment.pdf", width = 4, height =4)


ggplot(DT, aes(x=umap_1, y=umap_2,  color=label)) +
geom_point_rast(scale = .1, raster.dpi = 50) +
scale_color_manual(values = c('#4daf4a',  '#935503', '#888888', '#8c21ea',  '#377eb8', '#f2c406', '#e41a1c')) +
theme_void()
ggsave("sctapestri/img/UMAP_by_clone.pdf", width = 4, height = 4)



DTX <- DT[, .N/nrow(DT), by=label]

# Basic piechart
ggplot(DTX, aes(x="", y=V1, fill=label)) +
geom_bar(stat="identity", width=1, color="white") +
coord_flip() +
#coord_polar("y", start=0, direction = -1) +
scale_fill_manual(values = c('#4daf4a',  '#935503', '#888888', '#8c21ea',  '#377eb8', '#f2c406', '#e41a1c')) +
theme_void() # remove background, grid, numeric labels


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




DTX <- DT[, .N, by=.(commitment, sample_name, cell.type, label)]


for(lc in c("immature", "lymphoid", "myeloid"))
{

  ggplot(DTX[commitment == lc], aes(x=cell.type, y=N, fill=label)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#4daf4a',  '#935503', '#888888', '#8c21ea',  '#377eb8', '#f2c406', '#e41a1c'),  breaks = c("WT", "pre-leukemic", "NPM1", "NRAS",  "FLT3-TKD", "FLT3-ITD", "FLT3-LOH")) + 
  facet_wrap(~sample_name, nrow = 1, scales = "free_y") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  ggsave(sprintf("sctapestri/img/barplot_%s.pdf", lc), width = 6, height = 4)    
  
}



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
