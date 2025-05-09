##################################################################################
# 
# Wout Megchelenbrink
# May 08, 2025
#
# UMAP of the scRNAseq for paient 01 colored by stage and cluster
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scRNA_uncorrected.rds")
sox <- subset(so, patient == "NPM1-01")

# Normalization, dim reduction and clustering
sox <- NormalizeData(sox) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims=1:20) %>%
  FindNeighbors() %>%
  FindClusters(resolution = .3)

# Get cluster markers
# m <- FindAllMarkers(sox, logfc.threshold = .5, min.pct = .1, only.pos = T)
# markers <- as.data.table(m)
# markers[cluster == 1][order(pct.2 - pct.1)][1:30,]

# Selected cluster markers
genes <- c("CD38", "SELL", "TFRC", "PRTN3", "ELANE", "KIT", "CLEC11A",
           "TYMS", "UHRF1", "TOP2A", "MZB1", "CD69", "IL3RA", # Cycling
           "JCHAIN", "CD34", "IL2RA", "PTPRC", "CD99", "IGHM", "PIM1", "PDLIM1", "TNFRSF4",  # LMPP
           "CD3D", "IL32", # T-Cells
           "GP1BB", "GATA1", "GFI1B", "KLF1", # MEP/MKP 
           "HBD", "HBA1", "GYPA" # Eryth.
           )
setdiff(genes, rownames(sox))

# Scale data for marker genes
sox <- ScaleData(sox, features = genes)
mat <- sox@assays$RNA@scale.data[genes, ]


# Colors for the relative expression levels
col_fun <- colorRamp2(seq(-1, 2, length.out=20),  viridisLite::viridis(20))

# Annotation color for the stage
cols.stage <- c("#4955F3", '#82B300')
names(cols.stage) <- c("DX","RE")

# Get meta data
meta <- as.data.table(sox@meta.data, keep.rownames = "BC")

# Annotation color for the clusters
cols.cluster <- hue_pal()(length(unique(meta$seurat_clusters)))
names(cols.cluster) <- sort(unique(meta$seurat_clusters))

# Heatmap top annotation
ha <- HeatmapAnnotation(df=meta[, .(cluster=seurat_clusters, stage)],
                        which = "col",
                        col = list(cluster=cols.cluster, stage=cols.stage))

pdf(file = sprintf("scRNAseq/img/SupFig5c_P01_heatmap_markers.pdf"), width = 8, height = 7)
Heatmap(mat, 
        cluster_rows = T, 
        cluster_columns = F, 
        show_row_dend = F,
        show_column_dend = F, 
        show_column_names = F,
        row_names_side = "left",
        col = col_fun,
        column_split = meta$seurat_clusters,
        use_raster = T,
        raster_quality = 4,
        top_annotation = ha)
dev.off()