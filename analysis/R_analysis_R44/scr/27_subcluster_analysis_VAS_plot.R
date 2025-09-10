# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration of the subclusters of VAS cells

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
data.combined <- readRDS("../../out/object/130_VAS_subcluster_HarmonySample.rds")
# data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")

(DimPlot(data.combined,group.by = "dataset") + ggtitle("Harmony Sample"))
(DimPlot(data.combined,group.by = "dataset",split.by = "origin") + ggtitle("Harmony Sample"))

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../../out/image/130_UMAPCluster_tree_VAS_subcluster.pdf",width = 10,height = 10)

# general UMAP with new clustering
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",label = T,raster = F)
ggsave("../../out/image/130_UMAPCluster_res0.1_VAS_subcluster.pdf",width = 4,height = 3)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.4",label = T,raster = F)
ggsave("../../out/image/130_UMAPCluster_res0.4_VAS_subcluster.pdf",width = 4,height = 3)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",split.by = "origin",label = T,raster = F,ncol=4)
ggsave("../../out/image/130_UMAPCluster_splitTreat_res0.1_VAS_subcluster.pdf",width = 6,height = 3)
# DimPlot(data.combined, reduction = "umap", split.by = "origin",label = T,raster = T,ncol=5)
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
ggsave("../../out/image/130_UMAPCluster_resolutions_VAS_subcluster.pdf",width = 25,height = 15)

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
ggsave(plot=test_long01,"../../out/image/130_DotplotLong_res0.1_VAS_subcluster.pdf",width = 30,height = 5)

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
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_short)
ggsave("../../out/image/130_UMAPCluster_short_VAS_subcluster.pdf",width = 22,height = 12)

list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
ggsave("../../out/image/130_UMAPCluster_long_VAS_subcluster.pdf",width = 22,height = 12)

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

# -------------------------------------------------------------------------
# try to plot the score also by heatmap
# plot the average score per signature per cluster as a heatmap
df_signatureAbsinta <- data.combined@meta.data %>%
  dplyr::select(RNA_snn_res.0.1,contains("scoreLong_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-RNA_snn_res.0.1) %>%
  group_by(RNA_snn_res.0.1,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",RNA_snn_res.0.1,"_res0.1")) %>%
  ungroup()

mat_signatureAbsinta_avg <- df_signatureAbsinta %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

hm01 <- Heatmap(mat_signatureAbsinta_avg,name = "scaled_avg_exp")

pdf("../../out/image/130_heatmap_signatureAbsinta_res0.1_VASsubset.pdf",width = 9,height = 4)
draw(hm01,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

mat_signatureAbsinta_avg2 <- df_signatureAbsinta %>%
  group_by(cluster_id) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = signature,values_from = scaled_avg_score) %>%
  column_to_rownames("cluster_id")

hm02 <- Heatmap(mat_signatureAbsinta_avg2,name = "scaled_avg_exp")

pdf("../../out/image/130_heatmap_signatureAbsinta_res0.1_VASsubset2.pdf",width = 8,height = 4)
draw(hm02,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# -------------------------------------------------------------------------
# load more signatures for this dataset
sig01 <- read_csv("../../data/signatures/Yang2022/Yang2022_marker_general.csv") %>% 
  # filter(avg_logFC < 0)
  split(.$`Cell Type`) %>% map(function(x){x %>% pull(Gene)})
sig01_02 <- setNames(object = sig01,paste0("gen_",names(sig01) %>% str_remove_all(",|\\s|\\/")))

sig02 <- read_csv("../../data/signatures/Yang2022/Yang2022_marker_endo_subset.csv") %>%
  # filter(avg_logFC < 0)
  split(.$`Cell subtype`) %>% map(function(x){x %>% pull(Gene)})
sig02_02 <- setNames(object = sig02,paste0("gen_",names(sig02) %>% str_remove_all(",|\\s|\\/")))

sig03 <- read_csv("../../data/signatures/Yang2022/Yang2022_marker_mural_subset.csv") %>%
  # filter(avg_logFC < 0)
  split(.$`Cell subtype`) %>% map(function(x){x %>% pull(Gene)})
sig03_02 <- setNames(object = sig03,paste0("mur_",names(sig03) %>% str_remove_all(",|\\s|\\/")))

shortlist_feature_Yang2022 <- c(sig01_02,sig02_02,sig03_02)

# score the signatures, both long and short
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_feature_Yang2022,
                                        name = "_score_Yang2022")

df_rename_Yang2022 <- data.frame(names = data.combined@meta.data %>%
                                   colnames() %>%
                                   str_subset("_score_Yang2022"),
                                 rename = paste0("scoreYang2022_",names(shortlist_feature_Yang2022) %>% str_remove_all(",|\\s|\\/")))

lookup_Yang2022 <- df_rename_Yang2022$names
names(lookup_Yang2022) <- df_rename_Yang2022$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_Yang2022))

# plot the scores from AddModuleScore
list_plot_02_Yang2022 <- lapply(df_rename_Yang2022$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_Yang2022)
ggsave("../../out/image/130_UMAPCluster_Yang2022_VAS_subcluster.pdf",width = 25,height = 25)

# try to plot the score also by heatmap
# plot the average score per signature per cluster as a heatmap
df_signatureYang2022 <- data.combined@meta.data %>%
  dplyr::select(RNA_snn_res.0.1,contains("scoreYang2022_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-RNA_snn_res.0.1) %>%
  group_by(RNA_snn_res.0.1,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",RNA_snn_res.0.1,"_res0.1")) %>%
  ungroup()

mat_signatureYang2022_avg <- df_signatureYang2022 %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

hm01_02 <- Heatmap(mat_signatureYang2022_avg,name = "scaled_avg_exp")

pdf("../../out/image/130_heatmap_signatureYang2022_res0.1_VASsubset.pdf",width = 9,height = 5)
draw(hm01_02,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# attemtp at higher resolution
df_signatureYang2022_higres <- data.combined@meta.data %>%
  dplyr::select(RNA_snn_res.0.4,contains("scoreYang2022_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-RNA_snn_res.0.4) %>%
  group_by(RNA_snn_res.0.4,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",RNA_snn_res.0.4,"_res0.4")) %>%
  ungroup()

mat_signatureYang2022_avg_higres <- df_signatureYang2022_higres %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

mat_signatureYang2022_avg_higres %>%
  rownames_to_column("signature") %>%
  write_tsv("../../out/table/130_tableScaled_signatureYang2022_res0.4_VASsubset_higres.tsv")

df_signatureYang2022_higres %>%
  write_tsv("../../out/table/130_tableUnScaled_signatureYang2022_res0.4_VASsubset_higres.tsv")

hm01_02_highres <- Heatmap(mat_signatureYang2022_avg_higres,name = "scaled_avg_exp")

pdf("../../out/image/130_heatmap_signatureYang2022_res0.4_VASsubset_higres.pdf",width = 10,height = 5)
draw(hm01_02_highres,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# -------------------------------------------------------------------------


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
# ggsave("../../out/image/121_ViolinCluster_res0.2.pdf",width = 21,height = 7)

# marker per cluster ------------------------------------------------------
DefaultAssay(data.combined) <- "RNA"
Idents(data.combined) <- "RNA_snn_res.0.1"

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/130_FindAllMarkers_HarmonySample_res0.1_VAS_subcluster.tsv")

# pick the top 100 markers per cluster
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/130_FindAllMarkers_HarmonySample_res0.1_VAS_subcluster_top100.tsv")

sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
  filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
  filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
  write_tsv("../../out/table/130_FindAllMarkers_HarmonySample_res0.1_VAS_subcluster_top100_noRIBOandMT.tsv")

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
ggsave("../../out/image/130_Ditto_HarmonySample_res0.1_VAS_subcluster.pdf",width = 15,height = 5)

# plot the proportions
df_summary <- data.combined@meta.data %>%
  group_by(dataset,pathology_class,origin,RNA_snn_res.0.1) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(dataset,pathology_class,origin) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

df_summary %>%
  ggplot(aes(x=pathology_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.2),shape=1)+
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 45)) +
  facet_wrap(~RNA_snn_res.0.1,ncol = 1,scales = "free")+
  scale_y_sqrt()
ggsave("../../out/image/130_plot_clusterProp_res0.1_VAS_subcluster.pdf",height = 20,width = 4)

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

table(meta_test$dataset,meta_test$orig.ident)
table(meta_test$pathology_class,meta_test$dataset)
out_diagnosis <- propeller(clusters = meta_test$RNA_snn_res.0.1,
                           sample = paste0(meta_test$dataset),
                           group = meta_test$pathology_class)

out_diagnosis %>%
  rownames_to_column("RNA_snn_res.0.1") %>%
  write_tsv("../../out/table/130_propeller_res0.1_VAS_subcluster.tsv")

# plotting diagnosis ------------------------------------------------------
# default plot
speckle::plotCellTypeProps(x = data.combined,
                           clusters = data.combined$RNA_snn_res.0.1,
                           sample = data.combined$pathology_class)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/130_plot_propeller_res0.1_VAS_subcluster.pdf",height = 5,width = 5)

# custom plot
df_summary_diagnosis <- meta_test %>% 
  mutate(group_id = dataset) %>%
  group_by(RNA_snn_res.0.1,
           group_id,
           pathology_class) %>% 
  summarise(n=n()) %>% 
  ungroup() %>% 
  group_by(group_id) %>% 
  mutate(tot = sum(n),
         prop = n/tot)

df_summary_diagnosis %>%
  write_tsv("../../out/table/130_df_summary_diagnosis_res0.1_VAS_subcluster.tsv")

# plot 01
df_summary_diagnosis %>%
  ggplot(aes(x=pathology_class,y=prop))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
  facet_wrap(~RNA_snn_res.0.1,scales = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
ggsave("../../out/image/130_propeller_plot01_VAS_subcluster.pdf",width = 10,height = 10)

# plot 02
df_summary_diagnosis %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=pathology_class),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=pathology_class),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
# ggsave("../../out/image/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)
