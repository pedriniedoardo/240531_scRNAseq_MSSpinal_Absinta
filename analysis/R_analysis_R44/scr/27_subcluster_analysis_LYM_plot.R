# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration of the subclusters

# libraries ---------------------------------------------------------------
library(harmony)
library(Seurat)
library(dplyr)
library(cowplot)
library(tidyverse)
library(ggrepel)
library(scales)
library(RColorBrewer)
library(SeuratWrappers)
library(dittoSeq)
library(clustree)
library(pals)
library(patchwork)
library(ComplexHeatmap)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds")
# data.combined2 <- readRDS("../../out/object/100_LYM_subcluster_HarmonyRun.rds")

(DimPlot(data.combined,group.by = "sample_id") + ggtitle("Harmony Sample"))

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_tree_LYM_subcluster.pdf",width = 10,height = 10)

# general UMAP with former clustering

# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
  sort()

list_plot <- lapply(id_resolution,function(x){
  plot <- DimPlot(data.combined,
                  reduction = "umap",
                  group.by = x,
                  label = T,
                  raster = F)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_resolutions_LYM_subcluster.pdf",width = 25,height = 15)

# 
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

# make a shortlist of the markers
shortlist_features_list_short <- list(
  IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
  B_CELLS = c("IGHG1", "CD38"),
  T_CELLS =  c("SKAP1", "CD8A", "CD2"),
  OLIGO = c("MOG","MBP","MAG"),
  OPC = c("NLGN4X","OLIG1","OLIG2"),
  ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
  NEURONS = c("CUX2","SYP", "NEFL","SYT1"),
  VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
  SCHWANN = c("PMP22","MPZ","PRX"),
  EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
  STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
)

# score the signatures, both long and short
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_short,
                                        name = "_score_short")

data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_long,
                                        name = "_score_long")

df_rename_short <- data.frame(names = data.combined@meta.data %>%
                                colnames() %>%
                                str_subset("_score_short"),
                              rename = paste0("scoreShort_",names(shortlist_features_list_short)))

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score_long"),
                             rename = paste0("scoreLong_",names(shortlist_features_list_long)))

lookup_short <- df_rename_short$names
names(lookup_short) <- df_rename_short$rename

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_short))
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_short <- lapply(df_rename_short$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = F,pt.size = 0.2) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_short)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_short_LYM_subcluster.pdf",width = 22,height = 12)

list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = F,pt.size = 0.2) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_long_LYM_subcluster.pdf",width = 22,height = 12)

# res 0.1 -----------------------------------------------------------------
# select a specific subset of resolutions
# general UMAP with new clustering
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",label = T,raster = F)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_res0.1_LYM_subcluster.pdf",width = 4,height = 3)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",split.by = "diagnosis_short",label = T,raster = F,ncol=4)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_splitDisease_res0.1_LYM_subcluster.pdf",width = 8,height = 4)
# DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(data.combined,
                       features = shortlist_features_list_long,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/plot/analysis_R44/27_DotplotLong_res0.1_LYM_subcluster.pdf",width = 30,height = 6)

# same as above but as violin plot
list_plot <- lapply(df_rename_long$rename, function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
  return(test)
})

# make it a dataframe
# x <- list_plot[[1]]
df_violin <- lapply(list_plot,function(x){ 
  df <- x[[1]]$data 
  
  # extract the name of the gene 
  feature <- colnames(df)[1] 
  
  df %>% 
    mutate(feature = feature) %>% 
    setNames(c("value","ident","feature")) 
}) %>% 
  bind_rows()

head(df_violin) 

# plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 400,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(y = ident, x = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,nrow = 1,scales = "free") + 
#   theme_bw() + 
#   geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/plot/analysis_R44/121_ViolinCluster_res0.2.pdf",width = 21,height = 7)

# marker per cluster ------------------------------------------------------
DefaultAssay(data.combined) <- "RNA"
Idents(data.combined) <- "RNA_snn_res.0.1"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.1_LYM_subcluster.tsv")

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.1_LYM_subcluster_top100.tsv")

sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.1_LYM_subcluster_top100_noRIBOandMT.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.1")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/plot/analysis_R44/27_Ditto_HarmonySample_res0.1_LYM_subcluster.pdf",width = 10,height = 5)

# # plot the proportions
# df_summary <- data.combined@meta.data %>%
#   group_by(sample_id,diagnosis_short,location,RNA_snn_res.0.1) %>%
#   summarise(n = n()) %>%
#   ungroup() %>%
#   group_by(sample_id,diagnosis_short,location) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot)
# 
# df_summary %>%
#   ggplot(aes(x=diagnosis_short,y=prop)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2),shape=1)+
#   theme_bw() +
#   theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
#   facet_wrap(~RNA_snn_res.0.1,scales = "free")+
#   scale_y_sqrt()
# ggsave("../../out/plot/analysis_R44/27_plot_clusterProp_res0.1_LYM_subcluster.pdf",height = 8,width = 12)

# try the same with propeller on the same data
# renv::install("phipsonlab/speckle")
# renv::install("statmod")
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)

# run the test ------------------------------------------------------------
meta_test <- data.combined@meta.data
# filter(treat != "CSF.MS_RAPA")

table(meta_test$sample_id,meta_test$diagnosis_short)
table(meta_test$diagnosis_short,meta_test$sample_id)
out_diagnosis <- propeller(clusters = meta_test$RNA_snn_res.0.1,
                           sample = paste0(meta_test$sample_id),
                           group = meta_test$diagnosis_short)

out_diagnosis %>%
  rownames_to_column("RNA_snn_res.0.1") %>%
  write_tsv("../../out/table/analysis_R44/27_propeller_res0.1_LYM_subcluster.tsv")

# plotting diagnosis ------------------------------------------------------
# default plot
speckle::plotCellTypeProps(x = data.combined,
                           clusters = data.combined$RNA_snn_res.0.1,
                           sample = data.combined$diagnosis_short)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/27_plot_propeller_res0.1_LYM_subcluster.pdf",height = 5,width = 5)

# custom plot
df_summary_diagnosis <- meta_test %>% 
  # mutate(group_id = dataset) %>%
  group_by(RNA_snn_res.0.1,
           sample_id,
           diagnosis_short) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_diagnosis %>%
  write_tsv("../../out/table/analysis_R44/27_df_summary_diagnosis_res0.1_LYM_subcluster.tsv")

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=diagnosis_short,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/27_propeller_plot01_res0.1_LYM_subcluster.pdf",width = 9,height = 6)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis_short),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis_short),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
# ggsave("../../out/plot/analysis_R44/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)


# res 0.4 -----------------------------------------------------------------
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.4",label = T,raster = F)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_res0.4_LYM_subcluster.pdf",width = 4,height = 3)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.4",split.by = "diagnosis_short",label = T,raster = F,ncol=4)
ggsave("../../out/plot/analysis_R44/27_UMAPCluster_splitDisease_res0.4_LYM_subcluster.pdf",width = 8,height = 4)
# DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(data.combined,
                       features = shortlist_features_list_long,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "RNA_snn_res.0.4") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.4")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/plot/analysis_R44/27_DotplotLong_res0.4_LYM_subcluster.pdf",width = 30,height = 6)

# same as above but as violin plot
list_plot <- lapply(df_rename_long$rename, function(x){ 
  test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.4",raster = T)
  return(test)
})

# make it a dataframe
# x <- list_plot[[1]]
df_violin <- lapply(list_plot,function(x){ 
  df <- x[[1]]$data 
  
  # extract the name of the gene 
  feature <- colnames(df)[1] 
  
  df %>% 
    mutate(feature = feature) %>% 
    setNames(c("value","ident","feature")) 
}) %>% 
  bind_rows()

head(df_violin) 

# plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 400,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(y = ident, x = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,nrow = 1,scales = "free") + 
#   theme_bw() + 
#   geom_vline(data = df_plot_violin_summary,aes(xintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/plot/analysis_R44/121_ViolinCluster_res0.2.pdf",width = 21,height = 7)

# marker per cluster ------------------------------------------------------
DefaultAssay(data.combined) <- "RNA"
Idents(data.combined) <- "RNA_snn_res.0.4"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.4_LYM_subcluster.tsv")

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.4_LYM_subcluster_top100.tsv")

sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv("../../out/table/analysis_R44/27_FindAllMarkers_HarmonySample_res0.4_LYM_subcluster_top100_noRIBOandMT.tsv")

# try plotting the top markers
top_specific_markers <- sobj_total_h.markers %>%
  # filter ribosomal and mt genes
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  group_by(cluster) %>%
  top_n(3, avg_log2FC)

# And generate e.g. a dotplot:
dittoSeq::dittoDotPlot(data.combined,
                       vars = unique(top_specific_markers$gene), 
                       group.by = "RNA_snn_res.0.4")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
ggsave("../../out/plot/analysis_R44/27_Ditto_HarmonySample_res0.4_LYM_subcluster.pdf",width = 15,height = 6)

# # plot the proportions
# df_summary <- data.combined@meta.data %>%
#   group_by(sample_id,diagnosis_short,location,RNA_snn_res.0.4) %>%
#   summarise(n = n()) %>%
#   ungroup() %>%
#   group_by(sample_id,diagnosis_short,location) %>%
#   mutate(tot = sum(n)) %>%
#   mutate(prop = n/tot)
# 
# df_summary %>%
#   ggplot(aes(x=diagnosis_short,y=prop)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(position = position_jitter(width = 0.2),shape=1)+
#   theme_bw() +
#   theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
#   facet_wrap(~RNA_snn_res.0.4,scales = "free")+
#   scale_y_sqrt()
# ggsave("../../out/plot/analysis_R44/27_plot_clusterProp_res0.1_LYM_subcluster.pdf",height = 8,width = 12)

# try the same with propeller on the same data
# renv::install("phipsonlab/speckle")
# renv::install("statmod")
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)

# run the test ------------------------------------------------------------
meta_test <- data.combined@meta.data
# filter(treat != "CSF.MS_RAPA")

table(meta_test$sample_id,meta_test$diagnosis_short)
table(meta_test$diagnosis_short,meta_test$sample_id)
out_diagnosis <- propeller(clusters = meta_test$RNA_snn_res.0.4,
                           sample = paste0(meta_test$sample_id),
                           group = meta_test$diagnosis_short)

out_diagnosis %>%
  rownames_to_column("RNA_snn_res.0.4") %>%
  write_tsv("../../out/table/analysis_R44/27_propeller_res0.4_LYM_subcluster.tsv")

# plotting diagnosis ------------------------------------------------------
# default plot
speckle::plotCellTypeProps(x = data.combined,
                           clusters = data.combined$RNA_snn_res.0.4,
                           sample = data.combined$diagnosis_short)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/27_plot_propeller_res0.4_LYM_subcluster.pdf",height = 5,width = 5)

# custom plot
df_summary_diagnosis <- meta_test %>% 
  # mutate(group_id = dataset) %>%
  group_by(RNA_snn_res.0.4,
           sample_id,
           diagnosis_short) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(sample_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_diagnosis %>%
  write_tsv("../../out/table/analysis_R44/27_df_summary_diagnosis_res0.4_LYM_subcluster.tsv")

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=diagnosis_short,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.4,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/plot/analysis_R44/27_propeller_plot01_res0.4_LYM_subcluster.pdf",width = 9,height = 6)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.4,y=prop,color=diagnosis_short),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.4,y=prop,color=diagnosis_short),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
# ggsave("../../out/plot/analysis_R44/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)
