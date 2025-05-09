
library(data.table)
library(HDF5Array)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(viridis)
library(scales)
library(Seurat)

rm(list=ls())
gc()

ht_opt$message <- FALSE
# Check the structure of the HDF5 file
# get.protein.metadata <- function(h5.file)
# {
#   
#   # Get labels
#   meta <- data.table(barcode=h5read(file=h5.file, name = "/assays/protein_read_counts/ra/barcode"),
#                         original_sample_name=h5read(file=h5.file, name = "/assays/protein_read_counts/ra/sample_name"))
# 
#   return (meta)
# }




get.dna.barcode.metadata <- function(h5.file, remove.ambiguous.clones=T, min.clone.pct=0)
{
  
  meta <- data.table(do.call(cbind, h5read(file=h5.file, name = "/assays/dna_variants/ra")))
  setnames(meta, "sample_name", "original_sample_name")
  meta[, clone:=str_sub(label, end=str_locate(label, "__")[,1]-1)]
  
  meta[, sample_name:="allo-SCT"]
  meta[str_sub(original_sample_name, start=-1) == "D", sample_name:="Dx"]
  meta[str_sub(original_sample_name, start=-1) == "R", sample_name:="Rl"]
  
  if(remove.ambiguous.clones)
    meta <- meta[clone != "Ambiguous"]
  
  stats <- meta[, .(.N, pct=.N/nrow(meta)*100), by=.(clone)][pct >= min.clone.pct]
  meta <- meta[clone %in% stats$clone]
  
  return(meta)
}



get.protein.counts <- function(h5.file, raw=F)
{

  # Column info
  cols <- data.table(do.call(cbind, h5read(file=h5.file, name = "/assays/protein_read_counts/ca")))
  cols[antibody_sequence == 'CTCGTTTCCGTATCG', id:="FCER1"]

  # Read normalized protein counts
  if(raw)
  {
    mtx = t(h5read(file=h5.file, name = "/assays/protein_read_counts/layers/read_counts"))

  } else
  {
    mtx = t(h5read(file=h5.file, name = "/assays/protein_read_counts/layers/normalized_counts"))
  }
  colnames(mtx) <- cols$id

  # Set barcodes as rownames
  rownames(mtx) <- h5read(file=h5.file, name = "/assays/protein_read_counts/ra/barcode")
  return (mtx)
}


get.dna.variants <- function(h5.file, chr="13", type="all")
{
  # DNA
  dna.ca <- data.table(do.call(cbind, h5read(file=h5.file, name = "/assays/dna_variants/ca")))
  dna.ra <- data.table(do.call(cbind, h5read(file=h5.file, name = "/assays/dna_variants/ra")))
  
  # Load the detected/selected somatic DNA variants
  mtx <- t(h5read(file=h5.file, name="/assays/dna_variants/layers/NGT")) # or NGT_FILTERED
  rownames(mtx) <- dna.ra$barcode
  colnames(mtx) <- dna.ca$id
  
  if(str_detect(h5.file, "P09D"))
  {
    snvs <-  c("chr13:28592641:T/A", "chr13:28592642:C/A", "chr1:91180155:G/A", "chr1:211843642:G/A", "chr1:228210573:G/A")
    if (type == "snv")
    {
      # retain the FLT3-ITD and FLT3-TKD as SNVs
      id <- dna.ca[!(CHROM %in% c("1", "13") & str_length(ALT) <= 8) | id %in% snvs, id]
      mtx <- mtx[, id]
    } else if (type == "cnv")
    {
      id <- dna.ca[CHROM == chr & str_length(ALT) <= 8 & !id %in% snvs, id]
      mtx <- mtx[, id]
    }
    
    
  }
  else
  {
  
    if (type == "snv")
    {
      # retain the FLT3-ITD and FLT3-TKD as SNVs
      id <- dna.ca[!(CHROM == "13" & str_length(ALT) <= 8) | id %in% c("chr13:28592641:T/A",  "chr13:28592642:C/A"), id]
      mtx <- mtx[, id]
    } else if (type == "cnv")
    {
      id <- dna.ca[CHROM == "13" & str_length(ALT) <= 8 & !id %in% c("chr13:28592641:T/A",  "chr13:28592642:C/A"), id]
      mtx <- mtx[, id]
    }
  }  
  
  return (mtx)
}


get.snv.heatmap <- function(so)
{
  
  # Colors for clones
  if (length(unique(so$clone)) <= 3)
  {
    clone.colors <- c("#666666", "#E41A1C", "#377EB8")
  } else
  {
    clone.colors <- c("#666666", RColorBrewer::brewer.pal(length(unique(so$clone))-1, 'Set1'))
  }
  names(clone.colors) <- sort(unique(so$clone))
  
  # Sample colors
  sample.colors <- c("#4955F3", "#82B300", "#B20044") 
  names(sample.colors) <- c("Dx", "Rl", "allo-SCT")
  
  # NGT colors (WT, hetero, homo, NA)
  ngt_colors <- structure(c('#3C4E73','#78A3BC', '#D7ECED', '#000000'), names = c('0', '1', '2', '3'))  

  # order cells by clone and protein cluster  
  meta <- as.data.table(so@meta.data, keep.rownames = "bc")
  setorder(meta, clone, seurat_clusters)
  
  
  # Top annotation
  top.annot <- HeatmapAnnotation(clone=meta$clone, sample=meta$sample_name, col=(list(clone=clone.colors, sample=sample.colors)), which= "column")
  
  # Get matrix and make heatmap
  mtx <- data.matrix(GetAssayData(so, assay = "SNV"))
  
  # Discard variants from other patients
  mtx <- mtx[rowSums(mtx) > 0, meta$bc]

  hm  <- Heatmap(mtx, name = "GT", col=ngt_colors, column_split = meta$clone, 
                      cluster_rows = F, row_names_side = "left",
                      cluster_columns = F, top_annotation = top.annot, column_names_max_height = unit(50, 'mm'),
                      show_row_dend = F, use_raster = T, 
                      show_column_dend = F, show_column_names = F, show_row_names = T, 
                      cluster_column_slices = F, cluster_row_slices = F)
  
  return (hm)
}




get.cnv.heatmap <- function(so)
{
  ngt_colors <- structure(c('#3C4E73','#78A3BC', '#D7ECED', '#000000'), names = c('0', '1', '2', '3'))  
  
  # order cells by clone and protein cluster  
  meta <- as.data.table(so@meta.data, keep.rownames = "bc")
  setorder(meta, clone, seurat_clusters)
  
  # Get matrix and make heatmap
  mtx <- data.matrix(GetAssayData(so, assay = "CNV"))
  
  # Discard variants from other patients
  mtx <- mtx[rowSums(mtx) > 0, meta$bc]
  
  
  # CNV heatmap
  hm <- Heatmap(mtx, name = "GT", col=ngt_colors, column_split = meta$clone, 
                cluster_rows = F, column_names_side = "bottom",
                cluster_columns = F, use_raster = T)

  return (hm)
}




cluster.protein <- function(so, res=.4, selected.proteins=rownames(so), n.dims=15)
{
  VariableFeatures(so) <- selected.proteins
  so <- ScaleData(so)
  so <- RunPCA(so, npcs = n.dims)
  so <- RunUMAP(so, dims=1:n.dims)
  so <- FindNeighbors(so, dims = 1:n.dims)
  so <- FindClusters(so, resolution = res)
  
  return(so)
}






get.protein.heatmap <- function(so, selected.proteins=rownames(so))
{

  cat(sprintf("Creating heatmap with %s ADTs", length(selected.proteins)))

  # Row annot
  sample.colors <- c("#4955F3", "#82B300", "#B20044") 
  names(sample.colors) <- c("Dx", "Rl", "allo-SCT")
  
  n_clusters <- length(unique(so$seurat_clusters))
  cluster.colors <- hue_pal()(n_clusters)
  names(cluster.colors) <- sort(unique(so$seurat_clusters))

  # order cells by clone and protein cluster  
  meta <- as.data.table(so@meta.data, keep.rownames = "bc")
  setorder(meta, clone, seurat_clusters)
  
  
  # Top annotation
  
  anno.df <- data.table(stage=meta$sample_name, cluster=meta$seurat_clusters)
  annot <- rowAnnotation(df=anno.df, col=list(cluster=cluster.colors, stage=sample.colors))
  
  # Color scale for the normalized counts
  normalized_counts.colors <- colorRamp2(seq(-6, 6, length.out = 20), viridis::viridis(n=20))

  # Get the data
  mtx <- GetAssayData(so)
  mtx <- mtx[selected.proteins, meta$bc]
  
  # Top annotation
  top.annot <- HeatmapAnnotation(sample=meta$sample_name, cluster=meta$seurat_clusters, col=(list(sample=sample.colors, cluster=cluster.colors)), which= "column")
  
  # Heatmap
  hm <- Heatmap(mtx, top_annotation = top.annot, col=normalized_counts.colors,
                use_raster=T, show_row_dend=F, show_column_dend=F,
                cluster_columns=F, cluster_rows=T, show_row_names=T, row_names_side = "left",
                show_column_names=F)

  return (hm)
}





