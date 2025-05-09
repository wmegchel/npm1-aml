######################################################################################
#
# Wout Megchelenbrink
# March 22, 2025
# Gene ontology (GO) enrichment of immature cell types
#
# Modality: Single cell RNA-seq. 
######################################################################################
library(scales)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

# Clear workspace
rm(list=ls())

make.go.enrichment.plot <- function(sox, gencode.annot, cluster.id)
{

  # find cluster marker genes
  dm <- FindMarkers(sox, min.pct = .1, min.diff.pct = 0, ident.1 = cluster.id,  only.pos = T, logfc.threshold = 0.5)
  m <- as.data.table(dm, keep.rownames = "gene")[pct.1 >= .1]
  
  DT.genes <- merge(gencode.annot, m, by.x="gene_name", by.y="gene")

  # GO enrichment using cluster profiler
  ego <- enrichGO(gene          = DT.genes$ensembl_gene,
                  universe      = unique(gencode.annot$ensembl_gene),
                  OrgDb         = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  minGSSize = 50,
                  maxGSSize = 250,
                  readable = TRUE)
  
  # Convert to data.table
  DT <- as.data.table(ego)
  
  # Barplot of non-duplicated terms
  DT[, Description:=factor(Description, levels=rev(Description))]
  ggplot(DT[1:10,], aes(x=Description, y=-log10(p.adjust))) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  scale_y_continuous(limits = c(0, 10), breaks = seq(0,10, by=2.5)) +
  theme_minimal()

  # save PDF
  ggsave(sprintf("scRNAseq/img/Fig2g_GOenrichment_cluster_%d_top_10_terms.pdf", cluster.id+1), width = 8, height = 3)
}

# Read seurat object for stem- and progenitor cell clusters
sox <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")

# Gencode reference annotation
gencode.annot <- fread("processed_data/gencode.v38.genes.tsv")

# make the pdf barplots for each cluster
clusters <- c(0, 3, 6)
lapply(clusters, function(cluster.id) make.go.enrichment.plot(sox, gencode.annot, cluster.id))
