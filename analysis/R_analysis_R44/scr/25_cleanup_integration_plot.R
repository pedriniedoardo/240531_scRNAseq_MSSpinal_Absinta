# AIM ---------------------------------------------------------------------
# the aim is to plot the data after martina's cleanup and integration.

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
library(presto)
library(glmGamPoi)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/24_sobj_integrated_cleanup.rds")
Idents(data.combined) <- "cell_id_old"

# show the UMAP with the former anntoation
DimPlot(data.combined,label = T)

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../../out/plot/analysis_R44/25_UMAPCluster_tree.pdf",width = 10,height = 10)

# plot the UMAP with all the resolutions runs
id_resolution <- str_subset(colnames(data.combined@meta.data),pattern = "RNA_snn_res") %>%
  sort()

list_plot <- lapply(id_resolution,function(x){
  plot <- DimPlot(data.combined,
                  reduction = "umap",
                  group.by = x,
                  label = T,
                  raster = T)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../../out/plot/analysis_R44/25_UMAPCluster_resolutions.pdf",width = 25,height = 15)

# save the panel martina uses for the annotation
# shortlist_features_list_long <- list(
#   IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
#   B_CELLS = c("IGHG1", "CD38"),
#   T_CELLS =  c("SKAP1", "CD8A", "CD2"),
#   OLIGOLINEAGE = c("PLP1","MOG","PPP1R16B","TNS3","HMGB1","CD81","B2M","C1QL1","HLA-A","HLA-C","NLGN4X","OLIG1","OLIG2","CSPG4", "PDGFRA", "SOX6", "SOX10", "MBP","MAG"),
#   ASTRO = c("AQP4", "GFAP", "CD44", "AQP1", "VIM","APOE", "VCAN", "STAT3", "ABCA1", "TNC", "SDC4","SLC1A2","S100B"),
#   NEURONS = c("GAD2", "PVALB", "SV2C", "VIP", "TLE4", "CUX2", "THY1", "SLC17A7", "NRGN", "SATB2", "RORB", "SST", "STX1A", "STX1B", "SYP", "TH", "NEFL","SYT1"),
#   NPC = c("NES", "PAX6", "SOX1"),
#   ENDO = c("VWF","CDH5","TEK","PECAM1","FLT1","KDR","NOS3","MCAM","MMRN1","CLDN5","BMX","ANGPT2","GJA4","TIE1","ROBO4","ECSCR"),
#   PERICYTE = c("PDGFRB","DES","ACTA2","ANPEP","RGS5","ABCC9","KCNJ8","CD248","DLK1","NT5E","ANGPT1"),
#   SCHWANN = c("PMP22","MPZ","PRX"),
#   EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
#   STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
# )

# martina also suggested to reomve the NPC
shortlist_features_list_long2 <- list(
  IMMUNE = c("LYVE1","CD163","MRC1","LINGO1","HSPA1A","MOBP","CD22","CD83","HIF1A","VEGFA","SOD1","TREM2","CX3CR1","P2RY12","C3","CSF1R", "CD74", "RUNX1","C1QB","PTPRC"),
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

# # make a shortlist of the markers
# shortlist_features_list_short <- list(
#   IMMUNE = c("CX3CR1","P2RY12","C3","CSF1R", "CD74","C1QB"),
#   B_CELLS = c("IGHG1", "CD38"),
#   T_CELLS =  c("SKAP1", "CD8A", "CD2"),
#   OLIGO = c("MOG","MBP","MAG"),
#   OPC = c("NLGN4X","OLIG1","OLIG2"),
#   ASTRO = c("AQP4", "GFAP", "CD44", "AQP1"),
#   NEURONS = c("CUX2","SYP", "NEFL","SYT1"),
#   VAS = c("VWF","FLT1","CLDN5","PDGFRB"),
#   SCHWANN = c("PMP22","MPZ","PRX"),
#   EPENDYMA = c("CFAP299","DNAH7","DNAH9"),
#   STROMAL = c("LAMA2","RBMS3","CEMIP","GPC6")
# )

# plot the shortlisted feature per cluster
# notice that this is done only on the subset of the young (control) cells
test_long01 <- DotPlot(data.combined,
                       features = shortlist_features_list_long2,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "RNA_snn_res.0.1") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.1")+
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long01,"../../out/plot/analysis_R44/25_DotplotLong_sobj_integrated_cleanup_res0.1.pdf",width = 30,height = 6)

test_long04 <- DotPlot(data.combined,
                       features = shortlist_features_list_long2,
                       dot.scale = 8,
                       cluster.idents = T,
                       group.by = "RNA_snn_res.0.4") +
  RotatedAxis() +
  labs(title = "RNA_snn_res.0.4") +
  theme(strip.text = element_text(angle = 90))
ggsave(plot=test_long04,"../../out/plot/analysis_R44/25_DotplotLong_sobj_integrated_cleanup_res0.4.pdf",width = 30,height = 6)

# score the signatures, both long and short
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_long2,
                                        name = "_score_long")

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score_long"),
                             rename = paste0("scoreLong_",names(shortlist_features_list_long2)))

lookup_long <- df_rename_long$names
names(lookup_long) <- df_rename_long$rename

# rename the columns
data.combined@meta.data <- dplyr::rename(data.combined@meta.data,all_of(lookup_long))

# plot the scores from AddModuleScore
list_plot_02_long <- lapply(df_rename_long$rename,function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_02_long)
ggsave("../../out/plot/analysis_R44/25_UMAPCluster_ExpertAnnotaiton_long.pdf",width = 22,height = 15)

# -------------------------------------------------------------------------
# # same as above but as violin plot
# list_plot <- lapply(df_rename_long$rename, function(x){ 
#   test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
#   return(test)
# })
# 
# # make it a dataframe
# # x <- list_plot[[1]]
# df_violin <- lapply(list_plot,function(x){ 
#   df <- x[[1]]$data 
#   
#   # extract the name of the gene 
#   feature <- colnames(df)[1] 
#   
#   df %>% 
#     mutate(feature = feature) %>% 
#     setNames(c("value","ident","feature")) 
# }) %>% 
#   bind_rows()
# 
# head(df_violin) 
# 
# # how many cells per ident
# df_violin %>%
#   group_by(ident,feature) %>%
#   summarise(n = n()) %>%
#   filter(feature %in% c("score_ASTRO"))
# 
# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 800,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_summary <- df_plot_violin %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin %>%
#   ggplot(aes(x = ident, y = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,ncol = 1,scales = "free") + 
#   theme_bw() + 
#   geom_hline(data = df_plot_violin_summary,aes(yintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/plot/analysis_R44/25_ViolinCluster_ExpertAnnotaiton_res0.2.pdf",width = 7,height = 21)
# -------------------------------------------------------------------------

# use a similar approach to score potential technical clusters
# plot the scores from AddModuleScore
list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"),function(x){
  plot <- FeaturePlot(data.combined,features = x,order = T,
                      reduction = "umap",
                      raster = T) + scale_color_viridis_c(option = "turbo")
  return(plot)
})

wrap_plots(list_plot_technical)
ggsave("../../out/plot/analysis_R44/25_UMAPCluster_technical.pdf",width = 10,height = 8)

# -------------------------------------------------------------------------
# list_plot_technical <- lapply(c("nFeature_RNA","percent.mt","percent.ribo","percent.globin"), function(x){ 
#   test <- VlnPlot(object = data.combined,features = x, group.by = "RNA_snn_res.0.1",raster = T)
#   return(test)
# })
# 
# # make it a dataframe
# # x <- list_plot[[1]]
# df_violin_technical <- lapply(list_plot_technical,function(x){ 
#   df <- x[[1]]$data 
#   
#   # extract the name of the gene 
#   feature <- colnames(df)[1] 
#   
#   df %>% 
#     mutate(feature = feature) %>% 
#     setNames(c("value","ident","feature")) 
# }) %>% 
#   bind_rows()
# 
# head(df_violin_technical) 
# 
# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin_technical <- df_violin_technical %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 1000,replace = F) %>%
#   ungroup()
# 
# df_plot_violin_technical_summary <- df_plot_violin_technical %>%
#   group_by(feature) %>%
#   summarise(med_score = median(value))
# 
# df_plot_violin_technical %>%
#   ggplot(aes(x = ident, y = value)) + 
#   geom_violin(scale = "width")+ 
#   #geom_boxplot(outlier.shape = NA,position = position_dodge(width=0.9),width=0.05) + 
#   geom_point(position = position_jitter(width = 0.2),alpha = 0.05,size = 0.5) + 
#   facet_wrap(~feature,ncol = 1,scales = "free") + 
#   theme_bw() + 
#   geom_hline(data = df_plot_violin_technical_summary,aes(yintercept = med_score),linetype="dashed",col="red") +
#   theme(strip.background = element_blank(),
#         axis.text.x = element_text(hjust = 1,angle = 45))
# ggsave("../../out/plot/analysis_R44/25_ViolinCluster_technical_res0.1.pdf",width = 7,height = 12)
# -------------------------------------------------------------------------

# pick one resolution onward 0.1

# main umap
plot03 <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,raster = T)
ggsave(plot=plot03,"../../out/plot/analysis_R44/25_UMAPCluster_sobj_integrated_cleanup_res0.1.pdf",width = 6,height = 5)

plot03a <- DimPlot(data.combined, reduction = "umap", group.by = "RNA_snn_res.0.1",label = T,raster = T,split.by = "original_sample_name",ncol = 5)
ggsave(plot=plot03a,"../../out/plot/analysis_R44/25_UMAPClusterSplit_sobj_integrated_cleanup_res0.1_02.pdf",width = 13,height = 12)

# main umap cell cycle
plot04 <- DimPlot(data.combined, reduction = "umap", group.by = "Phase",raster = T,order = T)
ggsave(plot=plot04,"../../out/plot/analysis_R44/25_UMAPPhase_sobj_integrated_cleanup.pdf",width = 5,height = 4)

# split by sample
Idents(data.combined) <- "RNA_snn_res.0.1"
df_meta <- data.combined@meta.data %>%
  rownames_to_column("rowname")
df_UMAP <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column("rowname")

# plot the proporition for the phase per cluster
df_summary_phase <- df_meta %>%
  group_by(RNA_snn_res.0.1,Phase) %>%
  summarise(n = n()) %>%
  group_by(RNA_snn_res.0.1) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:8))) %>%
  ggplot() +
  geom_col(aes(x=RNA_snn_res.0.1,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/analysis_R44/25_BarplotPhase_summary_sobj_integrated_cleanup_res0.1.pdf",width = 7,height = 6)

# split by donor
df_summary_phase2 <- df_meta %>%
  group_by(original_sample_name,RNA_snn_res.0.1,Phase) %>%
  summarise(n = n()) %>%
  group_by(original_sample_name,RNA_snn_res.0.1) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)

df_summary_phase2 %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:8))) %>%
  ggplot() +
  geom_col(aes(x=original_sample_name,y=prop,fill=Phase))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  facet_wrap(~RNA_snn_res.0.1,ncol = 1)+theme(strip.background = element_blank())
ggsave("../../out/plot/analysis_R44/25_BarplotPhase2_summary_sobj_integrated_cleanup_res0.1.pdf",width = 7,height = 20)

# proportion of cell per cluster
df_summary <- df_meta %>%
  group_by(original_sample_name,RNA_snn_res.0.1) %>%
  summarise(n = n()) %>%
  group_by(original_sample_name) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary,"../../out/table/analysis_R44/25_summary_sobj_integrated_cleanup_res0.1.tsv")

color_id <- alphabet(length(unique(df_summary$RNA_snn_res.0.1)))
# check the colors
show_col(color_id)

df_summary %>%
  ggplot() +
  geom_col(aes(x=original_sample_name,y=prop,fill=RNA_snn_res.0.1))+theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  scale_fill_manual(values = unname(color_id))
ggsave("../../out/plot/analysis_R44/25_Barplot_summary_sobj_integrated_cleanup_res0.1.pdf",width = 7,height = 6)

# following Martina's reccomandation I added the info for diagnosis or location
df_summary_02 <- df_meta %>%
  group_by(original_sample_name,RNA_snn_res.0.1,location,diagnosis) %>%
  summarise(n = n()) %>%
  group_by(original_sample_name) %>%
  mutate(tot=sum(n)) %>%
  ungroup() %>%
  mutate(prop=n/tot)
write_tsv(df_summary_02,"../../out/table/analysis_R44/25_summary02_sobj_integrated_cleanup_res0.1.tsv")

df_summary_02 %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:8))) %>%
  mutate(location = fct_relevel(location,c("cervical","thoracic","lumbar"))) %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=location),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=location),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/25_BoxplotLocation_summary_sobj_integrated_cleanup_res0.1.pdf",width = 8,height = 5)

df_summary_02 %>%
  mutate(RNA_snn_res.0.1 = fct_relevel(RNA_snn_res.0.1,as.character(0:8))) %>%
  ggplot() +
  geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis),outlier.shape = NA) +
  geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=diagnosis),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
  theme_cowplot()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  scale_y_sqrt()
ggsave("../../out/plot/analysis_R44/25_BoxplotDiagnosis_summary_sobj_integrated_cleanup_res0.1.pdf",width = 8,height = 5)

# render the same plot as an heatmap
sample_prop_wide <- df_summary_02 %>%
  # scale by rows
  group_by(RNA_snn_res.0.1) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(original_sample_name,RNA_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = original_sample_name,values_from = zscore,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.1")

rowSums(sample_prop_wide)
colSums(sample_prop_wide)

# plot the data as heatmap
meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
  left_join(df_meta %>%
              group_by(original_sample_name,sex,diagnosis,location) %>%
              summarise(),by=c("colname" = "original_sample_name"))

column_meta_sample_prop <- HeatmapAnnotation(gender = meta_sample_prop$sex,
                                             location = meta_sample_prop$location,
                                             diagnosis = meta_sample_prop$diagnosis,
                                             col = list(gender = c("m" = "blue",
                                                                   "f" = "pink"),
                                                        location = c("lumbar" = "gray90",
                                                                     "thoracic" = "gray50",
                                                                     "cervical" = "black"),
                                                        diagnosis = c("Non-demented control" = "green",
                                                                      "Multiple sclerosis" = "red")))

ht2_shr_MG2 <- Heatmap(sample_prop_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                       name = "zscore \nprop_cell_type \nscale cluster",
                       column_title = "sample",
                       # col = viridis::viridis(option = "turbo",n = 10),
                       
                       # row_names_gp = gpar(fontsize = 3),
                       top_annotation = column_meta_sample_prop,show_row_names = T
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
)
pdf("../../out/plot/analysis_R44/25_HeatmapCluster_summary_sobj_integrated_cleanup_res0.1.pdf",width = 10,height = 6)
draw(ht2_shr_MG2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# render the same plot as an heatmap
sample_prop_wide2 <- df_summary_02 %>%
  # scale by rows
  group_by(original_sample_name) %>%
  mutate(zscore = (prop-mean(prop))/sd(prop)) %>%
  # make it long
  dplyr::select(original_sample_name,RNA_snn_res.0.1,zscore) %>%
  pivot_wider(names_from = original_sample_name,values_from = zscore,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.1")

rowSums(sample_prop_wide2)
colSums(sample_prop_wide2)

# plot the data as heatmap
meta_sample_prop2 <- data.frame(colname = colnames(sample_prop_wide2)) %>%
  left_join(df_meta %>%
              group_by(original_sample_name,sex,diagnosis,location) %>%
              summarise(),by=c("colname" = "original_sample_name"))

column_meta_sample_prop2 <- HeatmapAnnotation(gender = meta_sample_prop2$sex,
                                              location = meta_sample_prop2$location,
                                              diagnosis = meta_sample_prop2$diagnosis,
                                              col = list(gender = c("m" = "blue",
                                                                    "f" = "pink"),
                                                         location = c("lumbar" = "gray90",
                                                                      "thoracic" = "gray50",
                                                                      "cervical" = "black"),
                                                         diagnosis = c("Non-demented control" = "green",
                                                                       "Multiple sclerosis" = "red")))

ht2_shr_MG22 <- Heatmap(sample_prop_wide2, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                        name = "zscore \nprop_cell_type \nscale sample",
                        column_title = "sample",
                        # col = viridis::viridis(option = "turbo",n = 10),
                        
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_meta_sample_prop2,show_row_names = T
                        # cluster_rows = F,
                        # right_annotation = row_ha,
                        # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                        
)
pdf("../../out/plot/analysis_R44/25_HeatmapSample_summary_sobj_integrated_cleanup_res0.1.pdf",width = 10,height = 6)
draw(ht2_shr_MG22,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# render the same plot as an heatmap
sample_prop_wide3 <- df_summary_02 %>%
  # make it long
  dplyr::select(original_sample_name,RNA_snn_res.0.1,prop) %>%
  pivot_wider(names_from = original_sample_name,values_from = prop,values_fill = 0) %>%
  column_to_rownames("RNA_snn_res.0.1")

rowSums(sample_prop_wide3)
colSums(sample_prop_wide3)

# plot the data as heatmap
meta_sample_prop3 <- data.frame(colname = colnames(sample_prop_wide3)) %>%
  left_join(df_meta %>%
              group_by(original_sample_name,sex,diagnosis,location) %>%
              summarise(),by=c("colname" = "original_sample_name"))

column_meta_sample_prop3 <- HeatmapAnnotation(gender = meta_sample_prop3$sex,
                                              location = meta_sample_prop3$location,
                                              diagnosis = meta_sample_prop3$diagnosis,
                                              col = list(gender = c("m" = "blue",
                                                                    "f" = "pink"),
                                                         location = c("lumbar" = "gray90",
                                                                      "thoracic" = "gray50",
                                                                      "cervical" = "black"),
                                                         diagnosis = c("Non-demented control" = "green",
                                                                       "Multiple sclerosis" = "red")))

ht2_shr_MG23 <- Heatmap(sample_prop_wide3, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                        name = "prop_cell_type \nscale sample",
                        column_title = "sample",
                        col = viridis::viridis(option = "turbo",n = 10),
                        
                        # row_names_gp = gpar(fontsize = 3),
                        top_annotation = column_meta_sample_prop3,show_row_names = T
                        # cluster_rows = F,
                        # right_annotation = row_ha,
                        # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                        
)
pdf("../../out/plot/analysis_R44/25_HeatmapSampleProp_summary_sobj_integrated_cleanup_res0.1.pdf",width = 10,height = 6)
draw(ht2_shr_MG23,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 2), "mm"))
dev.off()

# # define a convenient palette of colors
# show_col(hue_pal()(9))
# RColorBrewer::display.brewer.all()
# col_pal <- RColorBrewer::brewer.pal(name = "Paired",n = 9)
# # col_pal <- c("#E6E6E6","#ffff09","#c1ce08","#446d05","#053c03","#4D4D4D","#06a8ce","#033b6d","#ff0ed7","#9a0404")
# show_col(col_pal)

# Identify conserved cell type markers ------------------------------------
# data
DefaultAssay(data.combined) <- "RNA"
# notice that in this case the data object of the RNA slot is already filled with the normalzied data, therefore in this case (following the Normalize workfloe for the integration) there is no need to run the NormalizeData on the RNA slot of the integrated object
# sobj_total_h@assays$RNA@data[20:50,1:10]
# dim(sobj_total_h@assays$RNA@data)
# 
# # scale the data see the note on evernote on why this can be done also at this point. this is needed because the scale.data is empty
# sobj_total_h@assays$RNA@scale.data
# all.genes <- rownames(sobj_total_h)
# sobj_total_h <- ScaleData(sobj_total_h,vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
# # confirm now the slot is loaded
# sobj_total_h@assays$RNA@scale.data[1:10,1:10]

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
Idents(data.combined) <- "RNA_snn_res.0.1"
# sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj_total_h.markers <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers %>%
  write_tsv("../../out/table/analysis_R44/25_FindAllMarkers_sobj_integrated_cleanup_res0.1.tsv")

# save the top 100
sobj_total_h.markers %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/analysis_R44/25_FindAllMarkers_sobj_integrated_cleanup_res0.1_top100.tsv")

# sobj_total_h.markers <- read_tsv("../../out/table/FindAllMarkers_harmonySkipIntegration_AllSoupX_01000_06000_15.tsv")


Idents(data.combined) <- "RNA_snn_res.0.4"
# sobj_total_h.markers2 <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sobj_total_h.markers2 <- FindAllMarkers(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# save the table of all markers
sobj_total_h.markers2 %>%
  write_tsv("../../out/table/analysis_R44/25_FindAllMarkers_sobj_integrated_cleanup_res0.4.tsv")

# save the top 100
sobj_total_h.markers2 %>%
  group_by(cluster) %>%
  dplyr::slice(1:100) %>%
  write_tsv("../../out/table/analysis_R44/25_FindAllMarkers_sobj_integrated_cleanup_res0.4_top100.tsv")

# -------------------------------------------------------------------------
# # try plotting the top markers
# top_specific_markers <- sobj_total_h.markers %>%
#   # filter ribosomal and mt genes
#   filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
#   filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
#   filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
#   group_by(cluster) %>%
#   top_n(5, avg_log2FC)
# 
# # And generate e.g. a dotplot:
# dittoSeq::dittoDotPlot(data.combined,
#                        vars = unique(top_specific_markers$gene), 
#                        group.by = "RNA_snn_res.0.2")+scale_color_viridis_c(option = "turbo",name="relative \nexpression")
# ggsave("../../out/plot/analysis_R44/manualClean/TopMarkersDitto_sobj_integrated_cleanup_res0.2.pdf",width = 15,height = 6)

# save the object ---------------------------------------------------------
# save the object with the full annotation
# saveRDS(data.combined,"../../out/object/manualClean/25_data.combined_sobj_integrated_cleanup_Annotation.rds")
