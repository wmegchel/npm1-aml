#############################################################################################
#
# Wout Megchelenbrink
# 
# May 05, 2025 
#
# Annotated heatmap of marker genes that discriminate the immature clusters 1, 2, 4 and 7
# Notice that these are cluster 0, 1, 3 and 6 in the Seurat object
#############################################################################################

# Clear workspace
rm(list=ls())

# Set RNG seed
set.seed(42)

# Load scRNA
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")
clusters <- c(0,1,3,6)
sox <- subset(so, seurat_clusters %in% clusters)

## Selection of cluster markers. Best visible when the relative abundance of cells expressing the gene (not FC) is 
# much higher in cluster A vs B. Also pval is considered.
# markers <- FindAllMarkers(sox)
# m <- as.data.table(markers)
# m[cluster == 0][order(pct.2 - pct.1)][1:30,]


# Selected gene list (many more can be selected)
genes <- c("CD52", 
           "CD74",
           "HLA-DPB1",
           "HLA-DPA1",
           "JAML",
           "JCHAIN",
           "CD36",
           "LTB",
           "TNFRSF4",
           "CD53",
           "FCER1G",
           "ANXA6",
           "FAM30A",
           "IGHM",
           "MEF2C",
           "CD79B",
           "PTPRC",
           "ANXA2",
           "LGALS1",
           "CRIP1",
           "CD9",
           "CTSG",
           "ELANE",
           "RNASE2",
           "AZU1",
           "PRTN3",
           "MPO",
           "ZBTB16",
           "CD38",
           "CD69",
           "ITM2A",
           "STMN1",
           "VPREB1",
           "IGLL1",
           "MZB1",
           "TRBC2",
           "PROM1")
          

# Scale data for the selected genes           
sox <- ScaleData(sox, features=genes)

# Subsample 25K cells, otherwise it will take forever to make the heatmap
sox.subsampled <- sox[, sample.int(ncol(sox), 25e3)]
mat <- sox.subsampled@assays$RNA$scale.data[genes, ]

meta <- as.data.table(sox.subsampled@meta.data)

# Colors for the scaled expression
col_fun <- colorRamp2(seq(-1, 2, length.out=10),  viridisLite::viridis(10))

# Colors for the clusters
cols.cluster <- hue_pal()(length(unique(meta$seurat_clusters)))
names(cols.cluster) <- sort(unique(meta$seurat_clusters))

# Colors for the patients
cols.patient <- hue_pal()(10)
names(cols.patient) <- meta[, sort(unique(donor))]

# Colors for the stage
cols.stage <- c("#4955F3", '#82B300')
names(cols.stage) <- c("DX", "RE")

# Heatmap top annotation
ha <- HeatmapAnnotation(df=meta[, .(cluster=seurat_clusters, stage, patient=donor)], 
                        which = "column", 
                        col = list('cluster'= cols.cluster, 'stage' = cols.stage, 'patient'= cols.patient))

# Create a PDF of the heatmap
pdf(file = "scRNAseq/img/SupFig2f_cluster_markers_immature_cells.pdf", width = 8, height = 8)
Heatmap(mat, 
        cluster_rows = T, 
        cluster_columns = F, 
        row_names_side = "left",
        show_row_dend = F,
        show_column_names = F, 
        show_column_dend = F,
        col = col_fun,
        top_annotation = ha,
        column_split = meta$seurat_clusters,
        use_raster = T,
        raster_quality = 4)
dev.off()
