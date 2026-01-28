# AIM ---------------------------------------------------------------------
# run the label transfer from the martina's manual annotation of the ASTRO subcluster brain dataset, onto the subcluster of the ASTRO of the spinal dataset

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
# library(SeuratData)
library(ggridges)
library(ComplexHeatmap)
library(SeuratWrappers)
library(cowplot)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# read in the data --------------------------------------------------------
# read in a sample object to use as reference
ref_WM <- readRDS("/beegfs/scratch/ric.cosr/pedrini.edoardo/project_edoardo/220501_scRNAseq_MSbrain_Absinta/out/object/132_AST_subcluster_HarmonySample_martina.rds")

# check the reference
DimPlot(ref_WM,label = T,group.by = "cellid")
DimPlot(ref_WM,label = T,group.by = "cellid",split.by = "origin")

# read in the query dataset this is the clean3 version of the dataset
query <- readRDS("../../out/object/analysis_R44/27_ASTRO_subcluster_HarmonySample.rds")

# check the query
DimPlot(query,label = T,raster = T,group.by = "RNA_snn_res.0.4")

# if needed fix the meta to add the general cell annotation
query$seurat_clusters_fix1 <- paste0("clu_",query$RNA_snn_res.0.4)
query@meta.data

# run the label transfer: ref over query ----------------------------------

# use the RNA slot to avoid issue in the FindTransferAnchors call
DefaultAssay(ref_WM)<-"RNA"

dim(ref_WM)
ref_WM@meta.data

# in the vignetted the query isn't the RNA slot
# DefaultAssay(query) <- "integrated"
# use the RNA slot in general to make sure all the features are accounted
DefaultAssay(query) <- "RNA"
dim(query)

# find the anchors
# FindTransferAnchors is highly sensitive to how the data was normalized. If your reference was normalized with LogNormalize and your query with SCTransform (or vice versa), the anchors will be poor.
# Ensure both objects have been processed with the same pipeline. If using the "RNA" assay, both should have gone through NormalizeData and FindVariableFeatures.
test.anchors_0 <- FindTransferAnchors(reference = ref_WM,
                                      # k.filter = NA,
                                      query = query,
                                      dims = 1:30,
                                      reference.reduction = "pca")

# predict the labels from the reference onto the query
predictions_0 <- TransferData(anchorset = test.anchors_0,
                              refdata = ref_WM$cellid,
                              dims = 1:30)

# notice that the dimension of the prediction is the same as the query object (in terms of number of cells)
dim(predictions_0)
dim(query)

# wrangling ---------------------------------------------------------------
# add the coordinates to the metadata
query_transfer <- AddMetaData(query, metadata = query@reductions$umap@cell.embeddings)

# add the predictions ot the query dataset
query_transfer <- AddMetaData(query_transfer, metadata = predictions_0)

# add the robust score metric to the object. use the 0.75 cut off as referenced in Azimuth https://azimuth.hubmapconsortium.org/
query_transfer$robust_score <- query_transfer@meta.data %>%
  rownames_to_column(var = "barcodes") %>%
  # if above the threshold, confirm the prediction, otherwise call uncertain
  mutate(robust_score = case_when(prediction.score.max>0.75~predicted.id,
                                  T ~ "uncertain")) %>%
  pull(robust_score)

# save the object
saveRDS(query_transfer,"../../out/object/28_query_ASTRO_transfer_RefOverQuery_res0.4.rds")

# EDA of the tranferring --------------------------------------------------

# plot the predicition score =============================================

# show the distribution of the max scores
predictions_0 %>%
  ggplot(aes(x=prediction.score.max))+geom_histogram()+theme_bw()
# ggsave("../../out/image/129_histo_prediction_label_transfer_refWMImmune_queryIMMLYMsubset.pdf",width = 5,height = 4)

# check the table of the robust score
table(query_transfer$robust_score)

# default plot: plot the raw prediction based on the reference dataset
DimPlot(query_transfer,label = T,group.by = "predicted.id")

# default plot: plot the robust prediction
DimPlot(query_transfer,label = T,group.by = "robust_score")

# average the position of the clusters
data_transfer_avg_raw <- query_transfer@meta.data %>% group_by(predicted.id) %>% select(umap_1, umap_2) %>% summarize_all(mean)

# custom plot: plot the raw prediction based on the reference dataset
ggplot(label= TRUE)+
  geom_point(data = query_transfer@meta.data,aes(x = umap_1,y = umap_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_raw,aes(x = umap_1,y = umap_2,label = predicted.id),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_raw.pdf",width = 5,height = 4)

# plot the robust prediction based on the reference dataset
# split the data for convinience
data_transfer_unc <- query_transfer@meta.data %>% 
  filter(robust_score == "uncertain")
#
data_transfer_defined <- query_transfer@meta.data %>% 
  filter(robust_score != "uncertain")

data_transfer_avg_robust <- data_transfer_defined %>% group_by(predicted.id) %>% select(umap_1, umap_2) %>% summarize_all(mean)

# custom plot: plot the robust prediction based on the reference dataset
ggplot(label= TRUE)+
  # plot the uncertain only in the background as gray
  geom_point(data = data_transfer_unc,aes(x = umap_1,y = umap_2),size=0.1,alpha=0.5,col="gray") +
  geom_point(data = data_transfer_defined,aes(x = umap_1,y = umap_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_robust,aes(x = umap_1,y = umap_2,label = predicted.id),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
# ggsave("../../out/image/129_UMAP_label_transfer_refWMImmune_queryIMMLYMsubset_robust.pdf",width = 5,height = 4)

# plot the predicted score per per cell type =============================

# since I am placing my cells on a reference it makes sense to me to show the score of how well each of my cells i mapped on the reference
id_factor <- query_transfer@meta.data %>% 
  group_by(seurat_clusters_fix1) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(seurat_clusters_fix1)

# plot the distrobution of the prediction score per cluster to show clusters with low average prediction score
query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
# ggsave("../out/plot/129_histo_label_transfer_refWMImmune_queryIMMLYMsubset_res0.9.pdf",width = 3,height = 2)

# do the same as above but split by condition
query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)+facet_wrap(~diagnosis_short)+theme(strip.background = element_blank())
# ggsave("../out/plot/129_histo_label_transfer_refWMImmune_queryIMMLYMsubset_split_res0.9.pdf",width = 9,height = 6)

# calculate the proportion os cells assigned to each annotation topped by the total number of cells in the cluster per condition
prop_table_tot_test <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,diagnosis_short,predicted.id,robust_score) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(seurat_clusters_fix1,diagnosis_short) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>% 
  ungroup() %>% 
  mutate(robust = case_when(robust_score!="uncertain"~"robust",
                            T~robust_score))

prop_table_tot_test %>% 
  filter(robust == "robust")

# plot the proprtion of cells split by annotation
prop_table_tot_test %>% 
  ggplot(aes(x=diagnosis_short,y=prop,col=seurat_clusters_fix1))+
  geom_point(position = position_jitter(width = 0.2))+
  facet_grid(robust~predicted.id)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),strip.background = element_blank())
# ggsave("../out/plot/label_transfer_refWMImmune_queryBrainsphere_prop_table_tot_test.pdf",width = 10,height = 6)

# same as above, but remove the scales per panel
prop_table_tot_test %>% 
  ggplot(aes(x=diagnosis_short,y=prop,col=seurat_clusters_fix1))+
  geom_point(position = position_jitter(width = 0.2))+
  facet_grid(predicted.id~robust,scales="free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),strip.background = element_blank())

# treat the uncertain as another annotation
prop_table_tot_test2 <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,diagnosis_short,robust_score) %>%
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(seurat_clusters_fix1,diagnosis_short) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>% 
  ungroup()
# mutate(robust = case_when(robust_score!="uncertain"~"robust",
#                           T~robust_score))

prop_table_tot_test2 %>% 
  ggplot(aes(x=diagnosis_short,y=prop,col=seurat_clusters_fix1))+
  geom_point(position = position_jitter(width = 0.2))+
  facet_grid(~robust_score)+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 45),strip.background = element_blank())
# ggsave("../../out/image/label_transfer_refWMImmune_queryBrainsphere_prop_table_tot_test2.pdf",width = 7,height = 3)

# heatmap for the assignament =============================================

# plot the proportion data as an heatmap regardless of the disease annotation

# regardless of the rubust score
prop_table_tot <- query_transfer@meta.data %>%
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

# only the rubust score
prop_table_robust <- query_transfer@meta.data %>%
  filter(robust_score != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

# raw assignament
# pdf("../out/plot/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_raw_res0.9.pdf",height = 4,width = 5)
hm1 <- Heatmap(prop_table_tot,
        name = "prop", 
        column_title = "predicted ASTRO id raw",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()

# only on robust cell assignament
# pdf("../out/plot/129_heatmap_label_transfer_refWMImmune_queryIMMLYMsubset_robust_res0.9.pdf",height = 4,width = 5)
hm2 <- Heatmap(prop_table_robust,
        name = "prop", 
        column_title = "predicted ASTRO id robust",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(11))
# dev.off()

hm1 + hm2

# use the Jaccard score ==================================================
# similarly to the proportion of cells per annotation. I could also measure the cross cluster similarity per cell (how similar are the cluster from the query comapred to the annotation derived from the reference), using the Jaccard score.

# define the jaccard score function
jaccard <- function(a, b) {
  intersection <- length(intersect(a, b))
  union <- length(a) + length(b) - intersection
  return (intersection/union)
}

# test
a <- c('potato', 'tomotto', 'chips', 'baloon')
b <- c('car', 'chips', 'bird', 'salt')

jaccard(a, b)

# build the dataset for the correlatino plot cross the two annotation per cells
df_crossing <- crossing(id_ref = unique(query_transfer@meta.data$predicted.id),
                        id_query = unique(query_transfer@meta.data$seurat_clusters_fix1))

# df_crossing

# build the scatter plot
df_jaccard_score <- pmap(list(id_ref = df_crossing$id_ref,
                              id_query = df_crossing$id_query), function(id_ref,id_query){
                                
                                # calculate the jaccard score
                                a <- query_transfer@meta.data %>%
                                  # rownames_to_column("barcode") %>%
                                  filter(predicted.id == id_ref) %>% pull(barcode)
                                
                                b <- query_transfer@meta.data %>%
                                  # rownames_to_column("barcode") %>%
                                  filter(seurat_clusters_fix1 == id_query) %>% pull(barcode)
                                
                                jaccard_score <- jaccard(a,b)
                                
                                # build a data.frame
                                df <- data.frame("id_ref" = id_ref,
                                                 "id_query" = id_query,
                                                 "jaccard_score" = jaccard_score)
                                return(df)
                              }) %>%
  bind_rows()

# check the table
head(df_jaccard_score)

# shape it as a matrix
mat_jaccard_score <- df_jaccard_score %>%
  pivot_wider(names_from = id_ref,values_from = jaccard_score) %>%
  column_to_rownames("id_query")

mat_jaccard_score

# plot the matrix
ht_02 <- Heatmap(mat_jaccard_score,
                 name = "Jaccard\nscore",
                 # col = colorRamp2(c(-1, 0, 1), colors = c("blue", "white", "red")),
                 col = viridis::viridis(option = "turbo",n = 20),
                 row_names_side = "right",
                 # row_names_gp = gpar(fontsize = 8),
                 column_names_side = "bottom",
                 # column_names_gp = gpar(fontsize = 8),
                 row_dend_reorder = FALSE,
                 column_dend_reorder = FALSE,
                 row_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 column_title_gp = gpar(fontsize = 10, fontface = "bold"),
                 show_column_names = T,
                 show_row_names = T)

# pdf("../out/plot/129_heatmap_jaccard_res0.9.pdf",height = 4,width = 5)
ht_02
# dev.off()
