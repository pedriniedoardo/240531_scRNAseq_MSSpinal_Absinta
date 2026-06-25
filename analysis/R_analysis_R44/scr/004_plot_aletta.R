# AIM ---------------------------------------------------------------------
# queries from aletta

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(ggalluvial)

# read in the data --------------------------------------------------------
# read in the object
sobj_update <- readRDS("../../out/object/analysis_R44/001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds")

DimPlot(sobj_update,label = T)

# analysis ----------------------------------------------------------------
# aletta asked to provide a dotplot with a panel of maker genes per cell type
Idents(sobj_update) <- "cell_id"

shortlist_features_list_long <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC","AIF1","HLA-DRA","TYROBP"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
  NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
  ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
  PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(sobj_update,
                       features = shortlist_features_list_long,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "cell_id") +
  RotatedAxis() +
  labs(title = "cell_id")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/plot/analysis_R44/004_DotplotLong_cellid.pdf",width = 30,height = 6)

# aletta suggested using another panel of marker genes

shortlist_features_list_long2 <- list(
  Active_MG = c("HLA-DRA", "CD74", "TREM2"),
  Astro = c("AQP4", "SLC1A2", "ALDH1L1", "GFAP", "S100B"),
  B_cells = c("MS4A1", "CD79A", "CD79B", "CD38", "IGHG1"),
  Endo = c("CLDN5", "PECAM1", "VWF", "FLT1", "KDR"),
  Ependymam = c("FOXJ1", "CFAP299", "DNAH9", "RSPH1"),
  Macro = c("LYVE1", "CD163", "MRC1"),
  MG = c("P2RY12", "TMEM119", "CX3CR1", "ADGRG1", "CSF1", "TYROBP", "AIF1", "C1QA", "C1QB"),
  Neu = c("SNAP25", "RBFOX3", "SYT1", "STX1A", "NEFL"),
  Oligo = c("PLP1", "MOG", "MBP", "MAG", "MOBP"),
  OPC = c("PDGFRA", "CSPG4", "PTPRZ1", "OLIG1", "SOX10"),
  Peri = c("PDGFRB", "RGS5", "ACTA2"),
  Schwann = c("MPZ", "PMP22", "PRX"),
  Stromal = c("COL1A1", "DCN", "LAMA2"),
  T_cells = c("CD3D", "CD3E", "CD2", "CD8A", "SKAP1", "TRBC1")
  )

test_long02 <- DotPlot(sobj_update,
                       features = shortlist_features_list_long2,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "cell_id") +
  RotatedAxis() +
  labs(title = "cell_id")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long02,"../../out/plot/analysis_R44/004_DotplotLong_cellid_panelAletta.pdf",width = 30,height = 6)

# alluvian plot -----------------------------------------------------------
# Aletta suggested including an alluvian plot to show how the proportion of the different cell types shift across locations

# pull the metadata from the object
df_meta <- sobj_update@meta.data

# generate the plot
df_summary <- df_meta %>%
  group_by(sample_id,location,diagnosis_short,pathological_stage_short,cell_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(sample_id) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  mutate(pathology_location = paste0(pathological_stage_short,"_",location)) %>%
  mutate(pathology_location = factor(pathology_location,levels = c("CTRL_cervical","CTRL_thoracic","CTRL_lumbar","IN_cervical","IN_thoracic","IN_lumbar","ACT_cervical","ACT_thoracic","ACT_lumbar")))
# save the table
df_summary %>%
  write_tsv("../../out/table/analysis_R44/004_table_summary_aletta.tsv")

# attempt generic plot per sample
# df_summary %>%
#   ggplot(aes(x = sample_id, y = prop, fill = cell_id)) + 
#   geom_col(position = "stack") +
#   facet_wrap(~ pathology_location, scales = "free_x",nrow=1) +  # Drops samples not in that location
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
#   labs(x = "Sample ID", y = "Proportion", fill = "Cell ID")

p03 <- df_summary %>%
  ggplot(aes(x = sample_id, y = prop, fill = cell_id)) +
  geom_col() + # default is position = "stack"
  facet_grid(. ~ pathology_location, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.text = element_text(angle = 90),
        panel.grid = element_blank())
ggsave(plot=p03,"../../out/plot/analysis_R44/004_porpCellID_cellid_panelAletta.pdf",width = 15,height = 6)

# generate the plot simplify the aaggreagtionk
df_summary2 <- df_meta %>%
  group_by(location,diagnosis_short,pathological_stage_short,cell_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(location,diagnosis_short,pathological_stage_short) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot) %>%
  mutate(pathology_location = paste0(pathological_stage_short,"_",location)) %>%
  mutate(pathology_location = factor(pathology_location,levels = c("CTRL_cervical","CTRL_thoracic","CTRL_lumbar","IN_cervical","IN_thoracic","IN_lumbar","ACT_cervical","ACT_thoracic","ACT_lumbar")))
# save the table
df_summary2 %>%
  write_tsv("../../out/table/analysis_R44/004_table_summary_aletta2.tsv")

# attempt generic plot per sample
# df_summary %>%
#   ggplot(aes(x = sample_id, y = prop, fill = cell_id)) + 
#   geom_col(position = "stack") +
#   facet_wrap(~ pathology_location, scales = "free_x",nrow=1) +  # Drops samples not in that location
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
#   labs(x = "Sample ID", y = "Proportion", fill = "Cell ID")

p04 <- df_summary2 %>%
  ggplot(aes(x = pathological_stage_short, y = prop, fill = cell_id)) +
  geom_col() + # default is position = "stack"
  facet_grid(. ~ location, scales = "free_x", space = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        strip.text = element_text(angle = 90),
        panel.grid = element_blank())
ggsave(plot=p04,"../../out/plot/analysis_R44/004_porpCellID_cellid_panelAletta2.pdf",width = 10,height = 6)
