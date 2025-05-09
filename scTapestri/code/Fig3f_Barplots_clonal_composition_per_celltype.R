##################################################################################
# 
# Wout Megchelenbrink
# May 02, 2025
#
# Barplots with cell counts per main clone within each cell.type
##################################################################################

# Clear workspace
rm(list = ls())

# Load data
so <- readRDS("processed_data/scTapestri_filtered_annotated_integrated.rds")

meta <- as.data.table(so@meta.data, keep.rownames = "BC")
DT <- meta[, .N, by=.(commitment, stage, cell.type, clone.label)]


meta[, .N/nrow(meta) * 100, by=clone.label]


for(lc in c("immature", "lymphoid", "myeloid"))
{
  ggplot(DT[commitment == lc], aes(x=cell.type, y=N, fill=clone.label)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = c('#4daf4a',  '#935503', '#888888', '#8c21ea',  '#377eb8', '#f2c406', '#e41a1c'),
                    breaks = c("WT", "pre-leukemic", "NPM1", "NRAS",  "FLT3-TKD", "FLT3-ITD", "FLT3-LOH")) + 
  facet_wrap(~stage, nrow = 1, scales = "free_y") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  ggsave(sprintf("scTapestri//img/Figure_3f_barplot_%s.pdf", lc), width = 6, height = 4)    
  
}

