# AIM ---------------------------------------------------------------------
# the aim is to plot the data of the subclusters

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
data.combined <- readRDS("../../out/object/17_NEU_subcluster_HarmonySample.rds")
# data.combined2 <- readRDS("../../out/object/100_MG_subcluster_HarmonyRun.rds")

(DimPlot(data.combined,group.by = "sample_id") + ggtitle("Harmony Sample"))
(DimPlot(data.combined,group.by = "dataset",split.by = "location") + ggtitle("Harmony Sample"))

# plots -------------------------------------------------------------------
# plot the tree of the cluster dependencies. this will justify the choice of the resolution, not too granula for the moment.
# library(clustree)
clustree::clustree(data.combined@meta.data[,grep("RNA_snn_res", colnames(data.combined@meta.data))],
                   prefix = "RNA_snn_res.")
ggsave("../../out/plot/subclusters/18_UMAPCluster_tree_NEU_subcluster.pdf",width = 10,height = 10)

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
                  raster = T)
  return(plot)
})

wrap_plots(list_plot)
ggsave("../../out/plot/subclusters/18_UMAPCluster_resolutions_NEU_subcluster.pdf",width = 25,height = 15)

# general UMAP with new clustering
DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",label = T,raster = T)
ggsave("../../out/plot/subclusters/18_UMAPCluster_res0.1_NEU_subcluster.pdf",width = 6,height = 5)

DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.4",label = T,raster = T)
ggsave("../../out/plot/subclusters/18_UMAPCluster_res0.4_NEU_subcluster.pdf",width = 6,height = 5)

# DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",split.by = "location",label = T,raster = T,ncol=4)
# DimPlot(data.combined, reduction = "umap",group.by = "RNA_snn_res.0.1",split.by = "pathological_stage",label = T,raster = T,ncol=4)
# ggsave("../../out/plot/subclusters/18_UMAPCluster_splitTreat_res0.1_NEU_subcluster.pdf",width = 6,height = 3)

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

# score the signatures
data.combined <- Seurat::AddModuleScore(data.combined,
                                        features = shortlist_features_list_long,
                                        name = "_score_long")

df_rename_long <- data.frame(names = data.combined@meta.data %>%
                               colnames() %>%
                               str_subset("_score_long"),
                             rename = paste0("scoreLong_",names(shortlist_features_list_long)))

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
ggsave("../../out/plot/subclusters/18_UMAPCluster_long_NEU_subcluster.pdf",width = 22,height = 12)

# -------------------------------------------------------------------------
# define the set of resolutions to test
res_ids <- c("res.0.1","res.0.4")
# rs <- "res.0.1"
# rs <- "res.0.4"

lapply(res_ids,function(rs){
  # plot the dotplot
  # plot the shortlisted feature per cluster
  # notice that this is done only on the subset of the young (control) cells
  test_long01 <- DotPlot(data.combined,
                         features = shortlist_features_list_long,
                         dot.scale = 8,
                         cluster.idents = T,
                         group.by = paste0("RNA_snn_",rs)) +
    RotatedAxis() +
    labs(title = paste0("RNA_snn_",rs))+
    theme(strip.text = element_text(angle = 90))
  ggsave(plot=test_long01,paste0("../../out/plot/subclusters/18_DotplotLong_",rs,"_NEU_subcluster.pdf"),
         width = 30,
         height = length(table(data.combined[[paste0("RNA_snn_",rs)]]))*0.2+4)
  
  # try to plot the score also by heatmap
  # plot the average score per signature per cluster as a heatmap
  df_signatureAbsinta <- data.combined@meta.data %>%
    dplyr::select(paste0("RNA_snn_",rs),contains("scoreLong_")) %>%
    pivot_longer(names_to = "signature",values_to = "score",-paste0("RNA_snn_",rs)) %>%
    
    group_by(!!sym(paste0("RNA_snn_",rs)),signature) %>%
    summarise(avg_score = mean(score),
              med_score = median(score)) %>%
    dplyr::rename(cluster = paste0("RNA_snn_",rs)) %>%
    mutate(cluster_id = paste0("clu_",cluster,"_",rs)) %>%
    ungroup()
  
  mat_signatureAbsinta_avg <- df_signatureAbsinta %>%
    group_by(signature) %>%
    mutate(scaled_avg_score = (avg_score - mean(avg_score))/sd(avg_score)) %>%
    ungroup() %>%
    select(signature,cluster_id,scaled_avg_score) %>%
    pivot_wider(names_from = cluster_id,values_from = scaled_avg_score) %>%
    column_to_rownames("signature")
  
  hm01 <- Heatmap(mat_signatureAbsinta_avg,name = "scaled_avg_exp")
  
  pdf(paste0("../../out/plot/subclusters/18_heatmap_signatureLong_",rs,"_NEUsubset.pdf"),
      width = dim(mat_signatureAbsinta_avg)[2]*0.5+6,
      height = 4)
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
  
  pdf(paste0("../../out/plot/subclusters/18_heatmap_signatureLong_",rs,"_NEUsubset2.pdf"),
      width = 6,
      height = dim(mat_signatureAbsinta_avg2)[1]*0.5+5)
  draw(hm02,heatmap_legend_side = "left",padding = unit(c(80, 2, 2, 2), "mm"))
  dev.off()
  
  # DE for markers
  DefaultAssay(data.combined) <- "RNA"
  Idents(data.combined) <- paste0("RNA_snn_",rs)
  
  # find markers for every cluster compared to all remaining cells, report only the positive
  # ones
  sobj_total_h.markers <- RunPrestoAll(data.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
    mutate(res = paste0("RNA_snn_",rs))
  
  # save the table of all markers
  sobj_total_h.markers %>%
    write_tsv(paste0("../../out/table/subclusters/18_FindAllMarkers_HarmonySample_",rs,"_NEU_subcluster.tsv"))
  
  # pick the top 100 markers per cluster
  sobj_total_h.markers %>%
    group_by(cluster) %>%
    dplyr::slice(1:100) %>%
    write_tsv(paste0("../../out/table/subclusters/18_FindAllMarkers_HarmonySample_",rs,"_NEU_subcluster_top100.tsv"))
  
  sobj_total_h.markers %>%
    group_by(cluster) %>%
    dplyr::slice(1:100) %>%
    filter(str_detect(gene,pattern = "^MT-",negate=T)) %>%
    filter(str_detect(gene,pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA",negate=T)) %>%
    filter(str_detect(gene,pattern = "^HB[^(P)]",negate=T)) %>%
    write_tsv(paste0("../../out/table/subclusters/18_FindAllMarkers_HarmonySample_",rs,"_NEU_subcluster_top100_noRIBOandMT.tsv"))
  
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
                         group.by = paste0("RNA_snn_",rs))+
    scale_color_viridis_c(option = "turbo",name="relative \nexpression")
  ggsave(paste0("../../out/plot/subclusters/18_Ditto_HarmonySample_",rs,"_NEUsubset.pdf"),
         width = dim(top_specific_markers)[1]*0.5+1,
         height = dim(top_specific_markers)[2]*0.4+1)
  
  # plot the proportions
  df_summary <- data.combined@meta.data %>%
    group_by(sample_id,pathological_stage,location,!!sym(paste0("RNA_snn_",rs))) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(sample_id,pathological_stage,location) %>%
    mutate(tot = sum(n)) %>%
    mutate(prop = n/tot)
  
  write_tsv(df_summary,paste0("../../out/table/subclusters/18_summary_prop_",rs,"_NEU_subcluster.tsv"))
  
  df_summary %>%
    dplyr::rename(cluster = paste0("RNA_snn_",rs)) %>%
    ggplot(aes(x=pathological_stage,y=prop,color=location))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(position = position_jitterdodge(dodge.width = 0.8,jitter.width = 0.2),shape=1)+
    theme_bw() +
    theme(strip.background = element_blank(), axis.text.x = element_text(hjust = 1,angle = 90)) +
    facet_wrap(~cluster,scales = "free")+
    scale_y_sqrt()
  ggsave(paste0("../../out/plot/subclusters/18_plot_clusterProp_",rs,"_NEUsubset1.pdf"),
         width = length(table(df_summary[[paste0("RNA_snn_",rs)]]))*1.5+1,
         height = length(table(df_summary[[paste0("RNA_snn_",rs)]]))*1.5)
})


# # try the same with propeller on the same data
# # renv::install("phipsonlab/speckle")
# # renv::install("statmod")
# library(speckle)
# library(limma)
# library(statmod)
# library(cowplot)
# library(ggrepel)
# library(finalfit)
# 
# # run the test ------------------------------------------------------------
# meta_test <- data.combined@meta.data
# # filter(treat != "CSF.MS_RAPA")
# 
# table(meta_test$dataset,meta_test$sample_id)
# table(meta_test$pathology_class,meta_test$pathological_stage)
# out_diagnosis <- propeller(clusters = meta_test$RNA_snn_res.0.1,
#                            sample = paste0(meta_test$sample_id),
#                            group = meta_test$pathological_stage)
# 
# # out_diagnosis %>%
# #   rownames_to_column("RNA_snn_res.0.1") %>%
# #   write_tsv("../../out/table/subclusters/131_propeller_res0.1_AST_subcluster.tsv")
# 
# # plotting diagnosis ------------------------------------------------------
# # default plot
# speckle::plotCellTypeProps(x = data.combined,
#                            clusters = data.combined$RNA_snn_res.0.1,
#                            sample = data.combined$pathology_class)+theme_minimal()+theme(panel.grid = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# # ggsave("../../out/image/131_plot_propeller_res0.1_AST_subcluster.pdf",height = 5,width = 5)
# 
# # custom plot
# df_summary_diagnosis <- meta_test %>% 
#   mutate(group_id = dataset) %>%
#   group_by(RNA_snn_res.0.1,
#            group_id,
#            pathology_class) %>% 
#   summarise(n=n()) %>% 
#   ungroup() %>% 
#   group_by(group_id) %>% 
#   mutate(tot = sum(n),
#          prop = n/tot)
# 
# # df_summary_diagnosis %>%
# #   write_tsv("../../out/table/subclusters/131_df_summary_diagnosis_res0.1_AST_subcluster.tsv")
# 
# # plot 01
# df_summary_diagnosis %>%
#   ggplot(aes(x=pathology_class,y=prop))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_point(position = position_jitter(width = 0.1),shape=1,alpha =0.7)+
#   facet_wrap(~RNA_snn_res.0.1,scales = "free")+
#   theme_bw()+
#   theme(strip.background = element_blank(),axis.text.x = element_text(hjust = 1,angle = 45))
# # ggsave("../../out/image/131_propeller_plot01_AST_subcluster.pdf",width = 10,height = 10)
# 
# # plot 02
# df_summary_diagnosis %>%
#   ggplot() +
#   geom_boxplot(aes(x=RNA_snn_res.0.1,y=prop,color=pathology_class),outlier.shape = NA) +
#   geom_point(aes(x=RNA_snn_res.0.1,y=prop,color=pathology_class),position = position_jitterdodge(jitter.width = 0.1,dodge.width = 0.8),alpha=0.7) +
#   theme_cowplot()+
#   theme(axis.text.x = element_text(hjust = 1,angle = 90))+
#   scale_y_sqrt()
# # ggsave("../../out/image/manualClean/propeller_plot02_diagnosis_cellid.pdf",width = 8,height = 5)
