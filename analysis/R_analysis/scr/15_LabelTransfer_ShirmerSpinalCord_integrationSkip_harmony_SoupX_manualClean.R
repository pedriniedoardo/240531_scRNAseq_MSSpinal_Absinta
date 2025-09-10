# AIM ---------------------------------------------------------------------
# run a label transfer to help identify where the cells of an unannotated object (query) are falling with respect with another reference object (reference)

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratData)
library(ggridges)
library(ComplexHeatmap)
library(cowplot)

# read in the reference data ----------------------------------------------
# read in the reference dataset from a spinal cord dataset https://cells.ucsc.edu/?ds=ms-cross-regional+sc
ref <- readRDS("../../data/Shirmer_spinal_2021_seurat.rds")
dim(ref)

DimPlot(ref,group.by = "cell_subtypes",label = T)
ggsave("../../out/plot/manualClean/UMAP_refSCShirmer.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_refSCShirmer.png",width = 10,height = 8,bg="white")

# define the reference Assay and identity
# DefaultAssay(ref) <- "RNA"
Idents(ref) <- "cell_subtypes"

# wrangling ---------------------------------------------------------------
# data.combined.markers <- FindAllMarkers(ref, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# # save the table of all markers
# data.combined.markers %>%
#   write_tsv("out/table/FindAllMarkers_ref.tsv")
# 
# # top 100 per cluster
# data.combined.markers %>%
#   group_by(cluster) %>%
#   mutate(rank = rank(order(p_val_adj, -abs(avg_log2FC)), ties.method='first')) %>%
#   arrange(cluster,rank) %>%
#   filter(rank < 101) %>%
#   write_tsv("out/table/FindAllMarkers_ref_top100.tsv")


# read in the query dataset -----------------------------------------------
query <- readRDS("../../out/object/manualClean/data.combined_harmonySkipIntegration_AllSoupX_00500_07000_05_AnnotationSCType.rds")

# fix the meta to add the general cell annotation
query$seurat_clusters_fix1 <- paste0("clu_",query$RNA_snn_res.0.2)
DimPlot(query,group.by = "seurat_clusters_fix1",label = T,raster = T)
ggsave("../../out/plot/manualClean/UMAP_querySC_res0.2.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_querySC_res0.2.png",width = 10,height = 8,bg="white")

# label transfer ----------------------------------------------------------
# define the reference Assay and identity
DefaultAssay(ref) <- "RNA"
Idents(ref) <- "cell_subtypes"

# in the vignetted the query is the RNA slot
DefaultAssay(query) <- "RNA"
Idents(ref) <- "seurat_clusters_fix1"

# the object is missing the HVG. add them to the slot
ref <- ref %>%
  FindVariableFeatures(selection.method = "vst")

# find tranfer anchor
test.anchors_0 <- FindTransferAnchors(reference = ref,
                                      query = query,
                                      dims = 1:30,
                                      reference.reduction = "pca")

# tranfer the lables
predictions_0 <- TransferData(anchorset = test.anchors_0,
                              # specify the label from the reference to transfer
                              refdata = ref$cell_subtypes,
                              dims = 1:30)

# add the predictions ot the query dataset
query_transfer <- AddMetaData(query, metadata = predictions_0)

# wrangling ---------------------------------------------------------------
# add the meta to the coordinates
data_transfer <- left_join(
  # get the coordinates of the query dataset
  query_transfer@reductions$umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column(var = "rowname"),
  
  # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
  query_transfer@meta.data %>%
    rownames_to_column(var = "rowname") %>%
    # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
    mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                    T ~ "uncertain"))
  ,"rowname")

# table of the assigned labels per cluster
table(data_transfer$robust_score,data_transfer$seurat_clusters_fix1)

# plotting ----------------------------------------------------------------
# show the distribution of the max scores
predictions_0 %>%
  ggplot(aes(x=prediction.score.max))+geom_histogram()+theme_bw()
ggsave("../../out/plot/manualClean/histo_prediction_label_transfer_refSCShirmer_querySC_res0.2.pdf",width = 5,height = 4)

dim(predictions_0)

# plot on the original query dataset
# DimPlot(query, reduction = "umap", group.by = "predicted.id")

# see how is different from the integrated one
# DimPlot(ref_WM, reduction = "umap", group.by = "clusterCellType")

# average the position of the clusters
data_transfer_avg_raw <- data_transfer %>% group_by(predicted.id) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
ggplot(label= TRUE)+
  geom_point(data = data_transfer,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_raw,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_raw_res0.2.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_raw_res0.2.png",width = 10,height = 8,bg="white")

# plot the raw prediction based on the reference dataset robust
data_transfer_unc <- data_transfer %>% 
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score == "uncertain")
#
data_transfer_defined <- data_transfer %>% 
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score != "uncertain")

data_transfer_avg_robust <- data_transfer_defined %>% group_by(predicted.id) %>% dplyr::select(UMAP_1, UMAP_2) %>% summarize_all(mean)

ggplot(label= TRUE)+
  geom_point(data = data_transfer_unc,aes(x = UMAP_1,y = UMAP_2),size=0.1,alpha=0.1,col="gray") +
  geom_point(data = data_transfer_defined,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_robust,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2.png",width = 10,height = 8,bg="white")

ggplot(label= TRUE)+
  geom_point(data = data_transfer_unc,aes(x = UMAP_1,y = UMAP_2),size=0.1,alpha=0.1,col="gray") +
  geom_point(data = data_transfer_defined,aes(x = UMAP_1,y = UMAP_2, col = predicted.id),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_transfer_avg_robust,aes(x = UMAP_1,y = UMAP_2,label = predicted.id),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme(legend.position = "none")
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_noLab.pdf",width = 8,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_noLab.png",width = 8,height = 8,bg="white")

# plot the predicted score per per cell type
# since I am placing my cells on a reference it makes sense to me to show the score of how well each of my cells i mappe on the reference
id_factor <- query_transfer@meta.data %>% 
  group_by(seurat_clusters_fix1) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(seurat_clusters_fix1)

# plot the mapping scores for the original label ID
query_transfer@meta.data %>%
  mutate(seurat_clusters_fix1= factor(seurat_clusters_fix1,levels = id_factor)) %>% 
  ggplot(aes(x=prediction.score.max,y=seurat_clusters_fix1)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
ggsave("../../out/plot/manualClean/histo_label_transfer_refSCShirmer_querySC_cluster_res0.2.pdf",width = 3,height = 4)

# plot the mapping scores for the ref labels
id_factor2 <- query_transfer@meta.data %>% 
  group_by(predicted.id) %>% 
  summarise(med = median(prediction.score.max)) %>% 
  arrange(med) %>% 
  pull(predicted.id)

query_transfer@meta.data %>%
  mutate(predicted.id= factor(predicted.id,levels = id_factor2)) %>% 
  ggplot(aes(x=prediction.score.max,y=predicted.id)) +
  geom_density_ridges(scale = 0.9) +
  theme_bw() +
  geom_vline(xintercept = 0.70,linetype="dashed",col="red",alpha=0.5)
ggsave("../../out/plot/manualClean/histo_label_transfer_refSCShirmer_querySC_predictedID_res0.2.pdf",width = 3,height = 4)

# render the assignaments as an heatmap
prop_table_tot <- query_transfer@meta.data %>%
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  # filter(robust_score_subclass != "uncertain") %>%
  # filter(robust_score_subclass != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>% 
  dplyr::select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

colSums(prop_table_tot)
rowSums(prop_table_tot)

prop_table_robust <- query_transfer@meta.data %>%
  mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                  T ~ "uncertain")) %>% 
  filter(robust_score != "uncertain") %>%
  group_by(seurat_clusters_fix1,predicted.id) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  dplyr::select(seurat_clusters_fix1,predicted.id,prop) %>%
  pivot_wider(names_from = seurat_clusters_fix1,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.id") %>%
  as.matrix()

colSums(prop_table_robust)
rowSums(prop_table_robust)

pdf("../../out/plot/manualClean/heatmap_label_transfer_refSCShirmer_querySC_raw_res0.2.pdf",height = 6,width = 6)
Heatmap(prop_table_tot,
        name = "prop \ncolumn scaled", 
        column_title = "predicted id raw",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

pdf("../../out/plot/manualClean/heatmap_label_transfer_refSCShirmer_querySC_robust_res0.2.pdf",height = 6,width = 6)
Heatmap(prop_table_robust,
        name = "prop \ncolumn scaled", 
        column_title = "predicted id robust",
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

# Unimodal UMAP Projection ------------------------------------------------
# In Seurat v4, we also enable projection of a query onto the reference UMAP structure. This can be achieved by computing the reference UMAP model and then calling MapQuery() instead of TransferData().
# this one run on the subset will overwrite the model
# if I don't run it the result for the mapping looks weird
ref2 <- RunUMAP(ref,
                dims = 1:30,
                reduction = "pca",
                return.model = TRUE)

query_transfer <- MapQuery(anchorset = test.anchors_0,
                           reference = ref2,
                           query = query_transfer,
                           refdata = list(celltype = "cell_subtypes"),
                           reference.reduction = "pca", reduction.model = "umap")

# We can now visualize the query cells alongside our reference.
# add the meta to the coordinates
data2_transfer <- left_join(
  # get the coordinates of the query dataset
  query_transfer@reductions$ref.umap@cell.embeddings %>%
    data.frame() %>%
    rownames_to_column(var = "barcodes"),
  
  # get the metadata from the query dataset, after adding the prediciton information. use a threshold for the score
  query_transfer@meta.data %>%
    rownames_to_column(var = "barcodes") %>%
    # this is based on martina's table, I need to ask her how she provided the accuracy of the imputation
    mutate(robust_score = case_when(prediction.score.max>0.70~predicted.id,
                                    T ~ "uncertain"))
  ,"barcodes")

# average the position of the clusters
data2_transfer_avg_raw <- data2_transfer %>% group_by(seurat_clusters_fix1) %>% dplyr::select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# plot the raw prediction based on the reference dataset
# save the reference UMAPS
UMAP_ref_WM <- DimPlot(ref2,group.by = "cell_subtypes",label = T)
UMAP_ref_WM$data %>%
  ggplot(label= TRUE)+
  UMAP_ref_WM$layers[[1]]+
  UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/plot/manualClean/UMAP_refSCShirmerRunUMAP.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_refSCShirmerRunUMAP.png",width = 10,height = 8,bg="white")

# save the image with no label and no legend
UMAP_ref_WM$data %>%
  ggplot(label= F)+
  UMAP_ref_WM$layers[[1]]+
  # UMAP_ref_WM$layers[[2]]+
  theme_cowplot()+
  theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave("../../out/plot/manualClean/UMAP_refSCShirmerRunUMAP_noLab.pdf",width = 8,height = 8)
ggsave("../../out/plot/manualClean/UMAP_refSCShirmerRunUMAP_noLab.png",width = 8,height = 8,bg="white")

ggplot(label= TRUE) +
  geom_point(data = ref2@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = data2_transfer,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_raw_res0.2_inverse.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_raw_res0.2_inverse.png",width = 10,height = 8,bg="white")


# save the image with no label and no legend
ggplot(label= F)+
  geom_point(data = ref2@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  geom_point(data = data2_transfer,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg_raw,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
# ggsave("../../out/image/UMAP_label_transfer_refWMOLIGO_queryWMCX_raw_inverse_no_lab.pdf",width = 5,height = 4)

data2_transfer %>% 
  filter(robust_score != "uncertain") %>% 
  group_by(seurat_clusters_fix1,robust_score) %>% 
  summarise(n = n()) %>% 
  print(n=70)

# divide the dataset into uncertain and not
data2_transfer_unc <- data2_transfer %>%
  filter(robust_score == "uncertain")
#
data2_transfer_defined <- data2_transfer %>%
  filter(robust_score != "uncertain")
# # focus only on cluster with more than 20 cells with robust prediction
# group_by(seurat_clusters_fix1,robust_score) %>% 
# mutate(tot_n = n()) %>% 
# filter(tot_n>20)

# # focus only on cluster with more than 20 cells with robust prediction
# data2_transfer_defined  %>% 
#   group_by(seurat_clusters_fix1,robust_score) %>% 
#   summarise(tot_n = n()) %>% 
#   filter(seurat_clusters=="cluster_10")

# average the position of the clusters
data2_transfer_avg <- data2_transfer_defined %>% group_by(seurat_clusters_fix1) %>% dplyr::select(refUMAP_1, refUMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot(label= TRUE) +
  geom_point(data = ref2@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_inverse.pdf",width = 10,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_inverse.png",width = 10,height = 8,bg="white")

# print also without the labels
ggplot(label= F) +
  geom_point(data = ref2@reductions$umap@cell.embeddings %>%
               data.frame(),aes(x=UMAP_1,y=UMAP_2),col="gray",alpha=0.3,size=0.1) +
  # geom_point(data = data2_transfer_unc,aes(x = refUMAP_1,y = refUMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data2_transfer_defined,aes(x = refUMAP_1,y = refUMAP_2, col = seurat_clusters_fix1),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data2_transfer_avg,aes(x = refUMAP_1,y = refUMAP_2,label = seurat_clusters_fix1),col="black")+
  theme_cowplot()+
  guides(colour = guide_legend(override.aes = list(size=5),ncol = 1))+
  theme(legend.position = "none")
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_inverse_noLab.pdf",width = 8,height = 8)
ggsave("../../out/plot/manualClean/UMAP_label_transfer_refSCShirmer_querySC_robust_res0.2_inverse_noLab.png",width = 8,height = 8,bg="white")

