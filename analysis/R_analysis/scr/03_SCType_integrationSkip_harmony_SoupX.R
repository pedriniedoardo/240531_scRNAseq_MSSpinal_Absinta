# AIM ---------------------------------------------------------------------
# attempat the automatic annotation using sctype algorithm

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggraph)
library(igraph)
library(data.tree)
library(HGNChelper)
library(openxlsx)
library(cowplot)
library(ComplexHeatmap)

# read in the functions ---------------------------------------------------
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R",destfile = "scr/R/gene_sets_prepare.R")
source("scr/SCType/gene_sets_prepare.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R",destfile = "scr/R/sctype_score_.R")
source("scr/SCType/sctype_score_.R")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R",destfile = "scr/R/auto_detect_tissue_type.R")
source("scr/SCType/auto_detect_tissue_type.R")

# download the database file
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",destfile = "../../data/ScTypeDB_full.xlsx")
# download.file("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",destfile = "../../data/ScTypeDB_short.xlsx")

# read in the final object ------------------------------------------------
scobj <- readRDS("../../out/object/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_Annotation.rds")

# confirm the identity of the object
DimPlot(scobj,label = T,raster = T)

# DimPlot(sobj_test_cluster,label = T)

# confirm the object has a scaled matrix in it. this is the one that will bu used for the scoring of the cell types
# data.combined@assays$RNA@scale.data[1:10,1:10]
# dim(data.combined@assays$RNA@scale.data)
scobj@assays$RNA@scale.data[1:10,1:10]
dim(scobj@assays$RNA@scale.data)

# read in the database file
db_ <- "../../data/ScTypeDB_full.xlsx"
tissue <- "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)
str(gs_list)

# Finally, let's assign cell types to each cluster:
# get cell-type by cell matrix
# -------------------------------------------------------------------------
# this is the official version using the scaled.data. but if not all the features have been scaled, some genes might be missing
# es.max <- sctype_score(scRNAseqData = scobj[["RNA"]]@scale.data, scaled = TRUE, 
#                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# -------------------------------------------------------------------------
# run an ad hoc scaling to include the genes for the cell type annotation
scobj_test <- scobj %>%
  # I can scale the missing features afterwards now focus on the highly variable one for speed purposes
  ScaleData(vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"), verbose = T,features = unique(unlist(gs_list))) %>% 
  identity()

dim(scobj_test@assays$RNA@scale.data)

es.max2 <- sctype_score(scRNAseqData = scobj_test[["RNA"]]@scale.data, scaled = TRUE, 
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
# es.max[1:11,1:10]
es.max2[1:24,1:10]

# wrangling ---------------------------------------------------------------
# after running the annotation in theory all the cells can have a specific annotation
df_score <- es.max2 %>% 
  data.frame() %>% 
  rownames_to_column("annotation") %>% 
  pivot_longer(names_to = "rowname",values_to = "score",-annotation) %>% 
  mutate(rowname = str_replace(rowname,pattern = "\\.",replacement = "-"))

df_score

# add the new annotation to the original object and show the umap
df_meta <- scobj@meta.data %>% 
  rownames_to_column("rowname")

head(df_meta)

# for each cluster in the dataset sum the scores per cell for each cluster
# cl <- "0"
df_score_cluster <- lapply(unique(df_meta$RNA_snn_res.0.1),function(cl){
  # pull the cells belonging to the clster
  id <- df_meta %>% 
    filter(RNA_snn_res.0.1 %in% cl) %>% 
    pull(rowname)
  
  # sum the score from eah cell for each category
  df_tot <- df_score %>% 
    filter(rowname %in% id) %>% 
    # add the cluster information
    mutate(cluster = cl) %>% 
    # add the total number of cells
    # mutate(n_cell = length(id)) %>% 
    group_by(annotation,cluster) %>% 
    summarise(sum_score = sum(score),
              n_cell = n())
  
  return(df_tot)
  
}) %>% 
  bind_rows() %>%
  ungroup()


# plot the results --------------------------------------------------------
# define the summary per annotation
df_score_cluster_summary <- df_score_cluster %>%
  group_by(annotation) %>%
  summarise(med_score = median(sum_score),
            max_score = max(sum_score)) %>%
  ungroup() %>%
  arrange(-max_score) %>%
  mutate(annotation = fct_reorder(annotation,max_score,.desc = T))

# plot the score per annotation
df_score_cluster %>%
  left_join(df_score_cluster_summary,by = "annotation") %>%
  mutate(annotation = fct_reorder(annotation,max_score,.desc = T)) %>%
  mutate(cluster = fct_relevel(cluster,as.character(0:12))) %>%
  ggplot(aes(x=cluster,y=sum_score)) +
  geom_point() +
  facet_wrap(~annotation,ncol = 4) +
  geom_hline(data = df_score_cluster_summary,aes(yintercept = med_score),col="red",linetype = "dashed") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/Dotplot_SCTypeClusterAnnotation_harmonySkipIntegration_AllSoupX_00500_07000_05.pdf",height = 12,width = 12)

# render the same plot as an heatmap
sample_score_cluster_wide_01 <- df_score_cluster %>%
  # normalize per cluster
  group_by(annotation) %>%
  mutate(norm = (sum_score-mean(sum_score))/sd(sum_score)) %>%
  # make it wide
  dplyr::select(annotation,cluster,norm) %>%
  pivot_wider(names_from = annotation,values_from = norm) %>%
  column_to_rownames("cluster")

# # plot the data as heatmap
# meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
#   left_join(df_meta %>%
#               group_by(original_sample_name,sex) %>%
#               summarise(),by=c("colname" = "original_sample_name"))
# 
# column_meta_sample_prop <- HeatmapAnnotation(gender = meta_sample_prop$sex,
#                                              col = list(location = c("m" = "blue",
#                                                                      "f" = "pink")))

ht2_shr_MG2_01 <- Heatmap(sample_score_cluster_wide_01, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                       name = "SCType score \nannotation scaled",
                       column_title = "annotation",
                       # col = viridis::viridis(option = "turbo",n = 10),
                       
                       # row_names_gp = gpar(fontsize = 3),
                       # top_annotation = column_meta_sample_prop,
                       show_row_names = T
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
)
pdf("../../out/plot/Heatmap_SCTypeClusterAnnotation_harmonySkipIntegration_AllSoupX_00500_07000_05.pdf",width = 10,height = 8)
draw(ht2_shr_MG2_01,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2, 30), "mm"))
dev.off()

# define the summary per annotation
df_score_cluster_summary2 <- df_score_cluster %>%
  group_by(cluster) %>%
  summarise(med_score = median(sum_score),
            max_score = max(sum_score)) %>%
  ungroup() %>%
  arrange(-max_score) %>%
  mutate(cluster = fct_relevel(cluster,as.character(0:12)))

# plot the score per annotation
df_score_cluster %>%
  left_join(df_score_cluster_summary2,by = "cluster") %>%
  # mutate(annotation = fct_reorder(annotation,max_score,.desc = T)) %>%
  mutate(cluster = fct_relevel(cluster,as.character(0:12))) %>%
  ggplot(aes(x=annotation,y=sum_score)) +
  geom_point() +
  facet_wrap(~cluster,ncol = 4) +
  geom_hline(data = df_score_cluster_summary2,aes(yintercept = med_score),col="red",linetype = "dashed") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1,angle = 90))
ggsave("../../out/plot/Dotplot_SCTypeAnnotationCluster_harmonySkipIntegration_AllSoupX_00500_07000_05.pdf",height = 8,width = 15)

# render the same plot as an heatmap
sample_score_cluster_wide <- df_score_cluster %>%
  # normalize per cluster
  group_by(cluster) %>%
  mutate(norm = (sum_score-mean(sum_score))/sd(sum_score)) %>%
  # make it wide
  dplyr::select(annotation,cluster,norm) %>%
  pivot_wider(names_from = cluster,values_from = norm) %>%
  column_to_rownames("annotation")

# # plot the data as heatmap
# meta_sample_prop <- data.frame(colname = colnames(sample_prop_wide)) %>%
#   left_join(df_meta %>%
#               group_by(original_sample_name,sex) %>%
#               summarise(),by=c("colname" = "original_sample_name"))
# 
# column_meta_sample_prop <- HeatmapAnnotation(gender = meta_sample_prop$sex,
#                                              col = list(location = c("m" = "blue",
#                                                                      "f" = "pink")))

ht2_shr_MG2 <- Heatmap(sample_score_cluster_wide, show_column_names = T,raster_by_magick = T,show_row_dend = F, use_raster = T,
                       name = "SCType score \ncluster scaled",
                       column_title = "cluster",
                       # col = viridis::viridis(option = "turbo",n = 10),
                       
                       # row_names_gp = gpar(fontsize = 3),
                       # top_annotation = column_meta_sample_prop,
                       show_row_names = T
                       # cluster_rows = F,
                       # right_annotation = row_ha,
                       # row_split = rep(c(1,2,3,4),c(2,3,4,7))
                       
)
pdf("../../out/plot/Heatmap_SCTypeAnnotationCluster_harmonySkipIntegration_AllSoupX_00500_07000_05.pdf",width = 10,height = 6)
draw(ht2_shr_MG2,heatmap_legend_side = "left",annotation_legend_side = "left",padding = unit(c(30,2,2,30), "mm"))
dev.off()


# fill the annotation -----------------------------------------------------
# pull the top score per cluster
df_score_cluster_top <- df_score_cluster %>%
  group_by(cluster) %>%
  top_n(n = 1, wt = sum_score) %>% 
  # set low-confident (low ScType score) clusters to "unknown"
  mutate(annotation_confident = case_when(sum_score < n_cell/4 ~"Unknown",
                                          T~annotation)) %>% 
  ungroup()

df_score_cluster_top

# add the annotation to the original dataset.
# add another costum annotation for the origin metadata
df_meta_full <- left_join(df_meta,df_score_cluster_top,c("RNA_snn_res.0.1"="cluster"))

# df_meta_full <- left_join(df_meta,df_score_cluster_top,c("seurat_clusters"="cluster")) %>% 
#   mutate(orig_alt = case_when(orig.ident %in% c("s31","s27","s26","s25")~"wm_new",
#                               T ~ origin))
head(df_meta_full)

# add it back to the object
scobj@meta.data <- df_meta_full %>% 
  column_to_rownames("rowname")

head(scobj@meta.data)

# plotting ----------------------------------------------------------------
DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation',raster = T)
ggsave("../../out/plot/UMAPSCTypeRNN01_harmonySkipIntegration_AllSoupX_00500_07000_05.pdf",width = 8,height = 5)

# for francesca plot also the split between cortex and wm
# DimPlot(scobj, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'annotation',split.by = "origin",raster = T)
# ggsave("../../out/image/ManualClean/UMAP_data.combined_WM_CX_harmonySkipIntegAllSoupX_SCType2.pdf",width = 15,height = 7)

# save the annotated object -----------------------------------------------
saveRDS(scobj,"../../out/object/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.rds")
