######################################################################################
#
# Wout Megchelenbrink
# March 22, 2025
# Plot module scores
#
# Modules based on https://www.biorxiv.org/content/10.1101/2024.05.07.592983v1.full.pdf, 
# Figure 1E: 
# LMPP: Genes expressed in LMPP1_CD45RA+, ProgB_CD10- or ProgB CD10+
# GMP: Genes expressed in LMPP3_CLL1+ (inc. CLEC12A/CLL1)
######################################################################################

# Clear workspace
rm(list=ls())

# Read the RPCA reduced progenitor data
so <- readRDS("processed_data/scRNA_stem_and_progenitors_RPCA_integrated.rds")

# Add module scores
so <- AddModuleScore(so, 
                     features = list( c("SERPINB1", "IGLL1", "CTSG", "MPO","AZU1", "ELANE", "CALR", "PRTN3", "IGLC2", "CLEC12A"),
                                      c("SPINK2", "FAM30A", "ADA", "DNTT", "IGHM", "LTB", "CD69", "CD79A", "MS4A1", "JCHAIN", "CD37", "IGKC")), name = "module")

# UMAP with GMP module score
FeaturePlot(so, features = "module1",
            raster = T, raster.dpi = c(320,320), pt.size = 0, 
            min.cutoff = -0.5, max.cutoff = 1) + 
scale_color_viridis_c() +
theme_void()
ggsave("scRNAseq/img/Fig2h_module_score_gmp.pdf", width=5, height = 4)

# UMAP with LMPP module score
FeaturePlot(so, features = "module2",
            raster = T, raster.dpi = c(320,320), pt.size = 0, 
            min.cutoff = -0.5, max.cutoff = 1) + 
scale_color_viridis_c() +
theme_void()
ggsave("scRNAseq/img/Fig2h_module_score_lmpp.pdf", width=5, height = 4)

