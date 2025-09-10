# AIM ---------------------------------------------------------------------
# the aim is to plot the data after integration. in particular to plot the sognature score for Martina's cell assignament.

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
library(circlize)

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/27_IMMUNE_subcluster_HarmonySample.rds")

(DimPlot(data.combined,group.by = "sample_id") + ggtitle("Harmony Sample"))

# load the signatures defined by Martina for the immune cell types
# use the top 100 DEGs she has defined  (the file was provided by martina)
df_LUT <- data.frame(type = c("homeostatic","MIMs foamy","MIMs iron","DC_2","DC_6","stressed MG"),
                     cluster = c(0,1,8,2,6,3))

list_sig <- read_csv("../../data/signatures/Signature_IMM_subset_Absinta2019.csv") %>%
  left_join(df_LUT,by = "cluster") %>%
  mutate(type = case_when(is.na(type)~paste0("clu_",cluster),
                          T~type)) %>%
  split(.$type) %>%
  lapply(function(x){
    x %>%
      pull(gene) %>%
      unique()
  })

# plots -------------------------------------------------------------------
# score the signatures for the custom signature
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = list_sig,
                                        name = "_score")

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score"),
                             rename = paste0("score_",names(list_sig)))

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
ggsave("../../out/plot/analysis_R44/27_UMAP_SigAbsinta2019_IMMUNE_subset_highres.pdf",width = 22,height = 15)

# res 0.1 -----------------------------------------------------------------
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

# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 150,replace = F) %>%
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
# ggsave("../../out/image/122_Violin_senescence_res0.2.pdf",width = 30,height = 4)

# plot the average score per signature per cluster as a heatmap
df_custom <- data.combined@meta.data %>%
  dplyr::select(RNA_snn_res.0.1,contains("score_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-RNA_snn_res.0.1) %>%
  group_by(RNA_snn_res.0.1,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",RNA_snn_res.0.1,"_res0.1")) %>%
  ungroup()

mat_custom_avg <- df_custom %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

# hm01 <- Heatmap(mat_custom_avg,name = "scaled_avg_exp",col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")))
hm01 <- Heatmap(mat_custom_avg,name = "scaled_avg_exp")

mat_custom_avg %>%
  rownames_to_column("signatures") %>%
  write_tsv("../../out/table/analysis_R44/27_tableScaled_SigAbisnta2019_res0.1_IMMUNE_subset_highres.tsv")

df_custom %>%
  write_tsv("../../out/table/analysis_R44/27_tableUnScaled_SigAbisnta2019_res0.1_IMMUNE_subset_highres.tsv")

pdf("../../out/plot/analysis_R44/27_heatmap_SigAbisnta2019_res0.1_IMMUNE_subset_highres.pdf",width = 8,height = 5)
draw(hm01,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# try also the scaling in the other direction, cluster-wise
mat_custom_avg2 <- df_custom %>%
  group_by(cluster_id) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = signature,values_from = scaled_avg_score) %>%
  column_to_rownames("cluster_id")

hm02 <- Heatmap(mat_custom_avg2,name = "scaled_avg_exp")

pdf("../../out/plot/analysis_R44/27_heatmap_SigAbisnta2019_res0.1_IMMUNE_subset_highres2.pdf",width = 9,height = 4)
draw(hm02,heatmap_legend_side = "left",padding = unit(c(4, 2, 2, 80), "mm"))
dev.off()

# res 0.4 -----------------------------------------------------------------
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

# # plot at maximum 1000 cells per group
# set.seed(123)
# df_plot_violin <- df_violin %>% 
#   group_by(ident,feature) %>%
#   sample_n(size = 150,replace = F) %>%
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
# ggsave("../../out/image/122_Violin_senescence_res0.2.pdf",width = 30,height = 4)

# plot the average score per signature per cluster as a heatmap
df_custom <- data.combined@meta.data %>%
  dplyr::select(RNA_snn_res.0.4,contains("score_")) %>%
  pivot_longer(names_to = "signature",values_to = "score",-RNA_snn_res.0.4) %>%
  group_by(RNA_snn_res.0.4,signature) %>%
  summarise(avg_score = mean(score),
            med_score = median(score)) %>%
  mutate(cluster_id = paste0("clu_",RNA_snn_res.0.4,"_res0.4")) %>%
  ungroup()

mat_custom_avg <- df_custom %>%
  group_by(signature) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
  column_to_rownames("signature")

# hm01 <- Heatmap(mat_custom_avg,name = "scaled_avg_exp",col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")))
hm01 <- Heatmap(mat_custom_avg,name = "scaled_avg_exp")

mat_custom_avg %>%
  rownames_to_column("signatures") %>%
  write_tsv("../../out/table/analysis_R44/27_tableScaled_SigAbisnta2019_res0.4_IMMUNE_subset_highres.tsv")

df_custom %>%
  write_tsv("../../out/table/analysis_R44/27_tableUnScaled_SigAbisnta2019_res0.4_IMMUNE_subset_highres.tsv")

pdf("../../out/plot/analysis_R44/27_heatmap_SigAbisnta2019_res0.4_IMMUNE_subset_highres.pdf",width = 10,height = 5)
draw(hm01,heatmap_legend_side = "left",padding = unit(c(2, 2, 2, 80), "mm"))
dev.off()

# try also the scaling in the other direction, cluster-wise
mat_custom_avg2 <- df_custom %>%
  group_by(cluster_id) %>%
  mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
  ungroup() %>%
  select(signature,cluster_id,scaled_avg_score) %>%
  pivot_wider(names_from = signature,values_from = scaled_avg_score) %>%
  column_to_rownames("cluster_id")

hm02 <- Heatmap(mat_custom_avg2,name = "scaled_avg_exp")

pdf("../../out/plot/analysis_R44/27_heatmap_SigAbisnta2019_res0.4_IMMUNE_subset_highres2.pdf",width = 9,height = 5)
draw(hm02,heatmap_legend_side = "left",padding = unit(c(4, 2, 2, 80), "mm"))
dev.off()
