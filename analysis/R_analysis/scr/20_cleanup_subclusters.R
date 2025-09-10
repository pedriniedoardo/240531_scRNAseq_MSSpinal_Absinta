# AIM ---------------------------------------------------------------------
# after looking at the subclustering, martina recommended to remove the following clusters from teh subcluster analysis:
# AST res0.1 removal of clusters 1 (oligo) and 4 (imm)
# IMM res0.1 removal of clusters 3 (oligo) and 4 (astro)
# LYM res 0.1 removal cluster 8 (stromal), 4 (astro), 5 (imm), 9 (oligo), 3 (oligo)
# OLIGO res0.1 removal cluster 2 (mix a bit of everything)
# OPC res 0.1 removal cluster 5 (imm),7 (stromal), 6 (lymp),3 (astro)
# VAS res0.1 removal cluster 3 (oligo), 6 (imm), 7(astro), 8 (oligo)
# NEU res 0.1 only 0 (oligo) and 7 (imm), since they are few neurons
# OTHER res0.1 removal cluster 8 (imm) and 4 (oligo)

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)

# input -------------------------------------------------------------------
# source file
# sobj <- readRDS("../../out/object/manualClean/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType_manualAnnotation.rds")

# vasc subcluster
# VAS res0.1 removal cluster 3 (oligo), 6 (imm), 7(astro), 8 (oligo)
sobj_VAS <- readRDS("../../out/object/17_VAS_subcluster_HarmonySample.rds")
DimPlot(sobj_VAS,group.by = "RNA_snn_res.0.1",label = T)

# neu subcluster
# NEU res 0.1 only 0 (oligo) and 7 (imm), since they are few neurons
sobj_NEU <- readRDS("../../out/object/17_NEU_subcluster_HarmonySample.rds")
DimPlot(sobj_NEU,group.by = "RNA_snn_res.0.1",label = T)

# other subcluster
# OTHER res0.1 removal cluster 8 (imm) and 4 (oligo)
sobj_OTHER <- readRDS("../../out/object/17_OTHER_subcluster_HarmonySample.rds")
DimPlot(sobj_OTHER,group.by = "RNA_snn_res.0.1",label = T)

# AST res0.1 removal of clusters 1 (oligo) and 4 (imm)
sobj_AST <- readRDS("../../out/object/17_AST_subcluster_HarmonySample.rds")
DimPlot(sobj_AST,group.by = "RNA_snn_res.0.1",label = T)

# IMM res0.1 removal of clusters 3 (oligo) and 4 (astro)
sobj_IMM <- readRDS("../../out/object/17_IMMUNE_subcluster_HarmonySample.rds")
DimPlot(sobj_IMM,group.by = "RNA_snn_res.0.1",label = T)

# LYM res 0.1 removal cluster 8 (stromal), 4 (astro), 5 (imm), 9 (oligo), 3 (oligo)
sobj_LYM <- readRDS("../../out/object/17_LYM_subcluster_HarmonySample.rds")
DimPlot(sobj_LYM,group.by = "RNA_snn_res.0.1",label = T)

# OLIGO res0.1 removal cluster 2 (mix a bit of everything)
sobj_OLIGO <- readRDS("../../out/object/17_OLIGO_subcluster_HarmonySample.rds")
DimPlot(sobj_OLIGO,group.by = "RNA_snn_res.0.1",label = T)

# OPC res 0.1 removal cluster 5 (imm),7 (stromal), 6 (lymp),3 (astro)
sobj_OPC <- readRDS("../../out/object/17_OPC_subcluster_HarmonySample.rds")
DimPlot(sobj_OPC,group.by = "RNA_snn_res.0.1",label = T)

# panel marker genes
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


# VAS ---------------------------------------------------------------------
# VAS res0.1 removal cluster 3 (oligo), 6 (imm), 7(astro), 8 (oligo)
# flag the barcodes
sobj_VAS$br_keep <- !sobj_VAS@meta.data$RNA_snn_res.0.1 %in% c(3,6,7,8)
sobj_VAS_clean <- subset(sobj_VAS,subset = br_keep == T)

test_long_VAS <- DotPlot(sobj_VAS,
                         features = shortlist_features_list_long,
                         dot.scale = 8,
                         cluster.idents = T,
                         group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_VAS_clean <- DotPlot(sobj_VAS_clean,
                               features = shortlist_features_list_long,
                               dot.scale = 8,
                               cluster.idents = T,
                               group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_VAS/test_long_VAS_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_VAS.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_VAS_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_VAS_clean.rds")

# NEU ---------------------------------------------------------------------
# NEU res 0.1 only 0 (oligo) and 7 (imm), since they are few neurons
# flag the barcodes
sobj_NEU$br_keep <- !sobj_NEU@meta.data$RNA_snn_res.0.1 %in% c(0,7)
sobj_NEU_clean <- subset(sobj_NEU,subset = br_keep == T)

test_long_NEU <- DotPlot(sobj_NEU,
                         features = shortlist_features_list_long,
                         dot.scale = 8,
                         cluster.idents = T,
                         group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_NEU_clean <- DotPlot(sobj_NEU_clean,
                               features = shortlist_features_list_long,
                               dot.scale = 8,
                               cluster.idents = T,
                               group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_NEU/test_long_NEU_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_NEU.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_NEU_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_NEU_clean.rds")

# OTHER -------------------------------------------------------------------
# OTHER res0.1 removal cluster 8 (imm) and 4 (oligo)
# flag the barcodes
sobj_OTHER$br_keep <- !sobj_OTHER@meta.data$RNA_snn_res.0.1 %in% c(4,8)
sobj_OTHER_clean <- subset(sobj_OTHER,subset = br_keep == T)

test_long_OTHER <- DotPlot(sobj_OTHER,
                         features = shortlist_features_list_long,
                         dot.scale = 8,
                         cluster.idents = T,
                         group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_OTHER_clean <- DotPlot(sobj_OTHER_clean,
                               features = shortlist_features_list_long,
                               dot.scale = 8,
                               cluster.idents = T,
                               group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_OTHER/test_long_OTHER_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_OTHER.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_OTHER_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_OTHER_clean.rds")

# AST ---------------------------------------------------------------------
# AST res0.1 removal of clusters 1 (oligo) and 4 (imm)
# flag the barcodes
sobj_AST$br_keep <- !sobj_AST@meta.data$RNA_snn_res.0.1 %in% c(1,4)
sobj_AST_clean <- subset(sobj_AST,subset = br_keep == T)

test_long_AST <- DotPlot(sobj_AST,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_AST_clean <- DotPlot(sobj_AST_clean,
                                 features = shortlist_features_list_long,
                                 dot.scale = 8,
                                 cluster.idents = T,
                                 group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_AST/test_long_AST_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_AST.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_AST_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_AST_clean.rds")


# IMM ---------------------------------------------------------------------
# IMM res0.1 removal of clusters 3 (oligo) and 4 (astro)
# flag the barcodes
sobj_IMM$br_keep <- !sobj_IMM@meta.data$RNA_snn_res.0.1 %in% c(3,4)
sobj_IMM_clean <- subset(sobj_IMM,subset = br_keep == T)

test_long_IMM <- DotPlot(sobj_IMM,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_IMM_clean <- DotPlot(sobj_IMM_clean,
                                 features = shortlist_features_list_long,
                                 dot.scale = 8,
                                 cluster.idents = T,
                                 group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_IMM/test_long_IMM_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_IMM.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_IMM_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_IMM_clean.rds")

# LYM ---------------------------------------------------------------------
# LYM res 0.1 removal cluster 8 (stromal), 4 (astro), 5 (imm), 9 (oligo), 3 (oligo)
# flag the barcodes
sobj_LYM$br_keep <- !sobj_LYM@meta.data$RNA_snn_res.0.1 %in% c(3,4,5,8,9)
sobj_LYM_clean <- subset(sobj_LYM,subset = br_keep == T)

test_long_LYM <- DotPlot(sobj_LYM,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_LYM_clean <- DotPlot(sobj_LYM_clean,
                                 features = shortlist_features_list_long,
                                 dot.scale = 8,
                                 cluster.idents = T,
                                 group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_LYM/test_long_LYM_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_LYM.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_LYM_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_LYM_clean.rds")

# OLIGO -------------------------------------------------------------------
# OLIGO res0.1 removal cluster 2 (mix a bit of everything)
# flag the barcodes
sobj_OLIGO$br_keep <- !sobj_OLIGO@meta.data$RNA_snn_res.0.1 %in% c(2)
sobj_OLIGO_clean <- subset(sobj_OLIGO,subset = br_keep == T)

test_long_OLIGO <- DotPlot(sobj_OLIGO,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_OLIGO_clean <- DotPlot(sobj_OLIGO_clean,
                                 features = shortlist_features_list_long,
                                 dot.scale = 8,
                                 cluster.idents = T,
                                 group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_OLIGO/test_long_OLIGO_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_OLIGO.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_OLIGO_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_OLIGO_clean.rds")

# OPC ---------------------------------------------------------------------
# OPC res 0.1 removal cluster 5 (imm),7 (stromal), 6 (lymp),3 (astro)
# flag the barcodes
sobj_OPC$br_keep <- !sobj_OPC@meta.data$RNA_snn_res.0.1 %in% c(3,5,6,7)
sobj_OPC_clean <- subset(sobj_OPC,subset = br_keep == T)

test_long_OPC <- DotPlot(sobj_OPC,
                           features = shortlist_features_list_long,
                           dot.scale = 8,
                           cluster.idents = T,
                           group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 full")+
  theme(strip.text = element_text(angle = 90))

test_long_OPC_clean <- DotPlot(sobj_OPC_clean,
                                 features = shortlist_features_list_long,
                                 dot.scale = 8,
                                 cluster.idents = T,
                                 group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1 cleanup")+
  theme(strip.text = element_text(angle = 90))

ggsave(plot = (test_long_OPC/test_long_OPC_clean),paste0("../../out/plot/subclusters/20_dotplot_cleanup_OPC.pdf"),
       width = 30,
       height = 10)

# save the metadata of the cleanup
saveRDS(sobj_OPC_clean@meta.data,file = "../../out/object/subclusters/20_meta_sobj_OPC_clean.rds")
