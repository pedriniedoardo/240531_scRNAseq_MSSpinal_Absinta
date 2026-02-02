# AIM ---------------------------------------------------------------------
# sample processing of Azimuth processing for reference labelling
# this comes after running the sample_Azimuth_01_preprocessing.R script
# this part run the actual label tranfer
# this is for the LYM subset dataset

# libraries ---------------------------------------------------------------
library(Seurat)
library(Azimuth)
library(SeuratData)
library(patchwork)
library(tidyverse)
library(ComplexHeatmap)
library(cowplot)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v5")
options(future.globals.maxSize = 1000 * 1024^2)

# run azimuth -------------------------------------------------------------
# read in the query
query <- LoadFileInput("../../out/object/analysis_R45/28_LYM_subcluster_HarmonySample_diet.rds")
head(query@meta.data)

# The RunAzimuth function can take a Seurat object as input
# pass the folder with the reference data
# Name of reference to map to or a path to a directory containing ref.Rds and idx.annoy
query <- RunAzimuth(query, reference = "renv/library/linux-rocky-9.5/R-4.5/x86_64-conda-linux-gnu/pbmcref.SeuratData/azimuth/")

# can also be passed
# query <- RunAzimuth(query, reference = "humancortexref")
saveRDS(query,"../../out/object/analysis_R45/28_Azimuth_LYM_subcluster_HarmonySample.rds")

# the metadata have all the info of interest
head(query@meta.data)

# specify the tags for the processing -------------------------------------
# define the annotation covariate
class <- "celltype.l2"
class_label <- paste0("predicted.",class)

# Construct the exact column names dynamically
col_class <- paste0("predicted.", class)
col_class_score <- paste0("predicted.", class, ".score")

# define the tag labels of the file
sample_id <- "LYM_subcluster_HarmonySample"
reference_id <- "HumanPBMCRef"

# define the cluster covariate
cluster_query <- "RNA_snn_res.0.4"

# condition splitting covariate
var_split <- "diagnosis_short"

# wrangling ---------------------------------------------------------------
DimPlot(query, group.by = class_label, label = TRUE, label.size = 3,raster = T) + NoLegend()
DimPlot(query, group.by = class_label, label = TRUE, label.size = 3,reduction = "ref.umap",raster = T) + NoLegend()
DimPlot(query, group.by = class_label, label = TRUE, label.size = 3,split.by = "orig.ident",raster = T) + NoLegend()

# save the UMAP coordinates of the new reference
query@reductions$ref.umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode") %>% 
  write_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_CoordUMAP_Azimuth_",reference_id,"_",class,".tsv"))

# save the new annotation from azimuth
query@meta.data %>% 
  data.frame() %>% 
  # rownames_to_column("barcode") %>% 
  write_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_meta_Azimuth_",reference_id,"_",class,".tsv"))

# custome plots threshold 0.75 --------------------------------------------
# costume visualization
UMAP_query <- read_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_CoordUMAP_Azimuth_",reference_id,"_",class,".tsv"))
meta_query <- read_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_meta_Azimuth_",reference_id,"_",class,".tsv")) %>%
  mutate(
    robust_class = case_when(.data[[col_class_score]] > 0.75 & mapping.score > 0.75 ~ .data[[col_class]],
                             TRUE ~ "uncertain"),
    predicted.class = .data[[col_class]],
    predicted.class.score = .data[[col_class_score]],
    cluster.query = .data[[cluster_query]],
    var_split = .data[[var_split]])

# define the levels of the cluster variable
level_class <- meta_query %>%
  group_by(predicted.class) %>%
  summarise(med = median(predicted.class.score)) %>%
  mutate(predicted.class = fct_reorder(predicted.class, med,.desc = T)) %>%
  pull(predicted.class) %>%
  levels()

# for each assignment what is the distribution of the scores
meta_query %>%
  mutate(predicted.class = factor(predicted.class, levels = level_class)) %>%
  ggplot(aes(x=predicted.class,y=predicted.class.score))+
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.2),alpha=0.01,size=0.1)+
  theme_bw() + 
  theme(axis.text.x = element_text(hjust = 1,angle = 45)) + 
  geom_hline(yintercept = 0.75,col="red")
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_score_075.pdf"),height = 4,width = 4)

meta_query %>%
  mutate(predicted.class = factor(predicted.class, levels = level_class)) %>%
  ggplot(aes(y=predicted.class,x=predicted.class.score))+
  ggridges::geom_density_ridges()+
  theme_bw() + 
  geom_vline(xintercept = 0.75,col="red",linetype="dashed")
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_ridges_score_075.pdf"),height = 5,width = 4)

# identify the most likely assignment for each seurat cluster
# first using all the subcluster annotation, not filtering for threshold of scores
prop_table_subclass <- meta_query %>%
  group_by(cluster.query,predicted.class) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(cluster.query,predicted.class,prop) %>%
  pivot_wider(names_from = cluster.query,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.class") %>% 
  as.matrix()

pdf(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_",cluster_query,"_Heatmap_score_075.pdf"),height = 4,width = 5)
Heatmap(prop_table_subclass,
        name = "prop", 
        column_title = paste(cluster_query,class,"score all",sep = " "),
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_class != "uncertain") %>%
  group_by(cluster.query,predicted.class) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot)

# first using all the subcluster annotation, filtering for threshold of scores
prop_table_subclass_filter <- meta_query %>%
  # filter(robust_score_subclass != "uncertain") %>%
  filter(robust_class != "uncertain") %>%
  group_by(cluster.query,predicted.class) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(cluster.query,predicted.class,prop) %>%
  pivot_wider(names_from = cluster.query,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.class") %>%
  as.matrix()

pdf(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_",cluster_query,"_HeatmapFilter_score_075.pdf"),height = 4,width = 5)
Heatmap(prop_table_subclass_filter,
        name = "prop", 
        column_title = paste(cluster_query,class,"score robust",sep = " "),
        cluster_rows = T,cluster_columns = T,col = viridis::turbo(10))
dev.off()

meta_query %>%
  filter(robust_class != "uncertain") %>%
  group_by(cluster.query,predicted.class) %>%
  summarise(n = n()) %>%
  mutate(tot = sum(n)) %>%
  mutate(prop = n/tot) %>%
  select(cluster.query,predicted.class,prop) %>%
  pivot_wider(names_from = cluster.query,values_from = prop,values_fill=0) %>%
  column_to_rownames("predicted.class")

# add the meta to the coordinates
data_query <- left_join(UMAP_query,meta_query,"barcode")
# add the groupign suggested by Martina
# data_query$group_martina <- case_when(data_query$orig.ident %in% c("pool_MOCK_CTRL")~"MOCK",
#                                       data_query$orig.ident %in% c("CTRL8_WT_CSF","CTRL8_WT_CTRL")~"CTRL",
#                                       T~"SOX10")

# divide the dataset into uncertain and not
data_query_unc <- data_query %>%
  filter(robust_class == "uncertain")
#
data_query_defined <- data_query %>%
  filter(robust_class != "uncertain")

# average the position of the clusters
data_query_avg <- data_query_defined %>% group_by(robust_class) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)
data_query_avg2 <- data_query_defined %>% group_by(cluster.query) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot()+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_class),col="black")+
  theme_cowplot()
# ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_UMAP_075.pdf",height = 5,width = 3)

# build the plot using both info
ggplot()+
  geom_point(data = data_query_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_query_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_query_avg,aes(x = UMAP_1,y = UMAP_2,label = robust_class),col="black")+
  theme_cowplot()+
  facet_wrap(~var_split,ncol=4)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthHumanCortexRef_subclass_split_UMAP_075.pdf",height = 15,width = 20)

# costum plotting ---------------------------------------------------------
# costume visualization
# umap coordinates of the reference onbject
UMAP_ref <- read_tsv("../../out/table/analysis_R45/azimuth_pbmcref_CoordUMAP.tsv")

# umap coordinates of the query object
UMAP_ref2 <- read_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_CoordUMAP_Azimuth_",reference_id,"_",class,".tsv"))

# original metadata of the reference object
meta_ref <- read_tsv("../../out/table/analysis_R45/azimuth_pbmcref_Metadata.tsv") %>%
  mutate(class = .data[[class]])
  
meta_query <- read_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_meta_Azimuth_",reference_id,"_",class,".tsv")) %>%
  mutate(
    robust_class = case_when(.data[[col_class_score]] > 0.75 & mapping.score > 0.75 ~ .data[[col_class]],
                             TRUE ~ "uncertain"),
    predicted.class = .data[[col_class]],
    predicted.class.score = .data[[col_class_score]],
    cluster.query = .data[[cluster_query]],
    var_split = .data[[var_split]])

# add the meta to the coordinates
data_ref <- left_join(UMAP_ref,meta_ref,"barcode")
data_ref2 <- left_join(UMAP_ref2,meta_query,"barcode") %>%
  mutate(cluster.query=factor(cluster.query))

# average the position of the clusters
data_ref_avg <- data_ref %>% group_by(class) %>% select(refumap_1,refumap_2) %>% summarize_all(mean)

# average the position of the clusters
data_ref_avg2 <- data_ref2 %>%
  filter(robust_class != "uncertain") %>%
  group_by(cluster.query) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# plot the reference
ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2,col=class),size=0.3,alpha=1) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg,aes(x = refumap_1,y = refumap_2,label = class),col="black",max.overlaps = 12)+theme_cowplot()
ggsave(paste0("../../out/plot/analysis_R45/",reference_id,"_",class,"_label.pdf"),width = 8,height = 6)

ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2,col=class),size=0.3,alpha=1) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=2)) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = data_ref_avg,aes(x = refumap_1,y = refumap_2,label = class),col="black",max.overlaps = 12)+
  theme_cowplot()
ggsave(paste0("../../out/plot/analysis_R45/",reference_id,"_",class,"_Nolabel.pdf"),width = 8,height = 6)

# build the plot using both info
ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_class != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = cluster.query),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = cluster.query),col="black",max.overlaps = 12)+
  theme_cowplot()
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAP_Cluster_075.pdf"),height = 6,width = 7)
# ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075.pdf",width = 9,height = 6)

ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_class != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = cluster.query),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2,aes(x = UMAP_1,y = UMAP_2,label = cluster.query),col="black",max.overlaps = 15)+
  theme_cowplot()+
  facet_wrap(~var_split,ncol=4)+theme(strip.background = element_blank())
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAP_Cluster_slit_075.pdf"),height = 6,width = 15)
# ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_split.pdf",width = 20,height = 15)

# average the position of the clusters
data_ref_avg2subclass <- data_ref2 %>%
  filter(robust_class != "uncertain") %>%
  group_by(predicted.class) %>% select(UMAP_1, UMAP_2) %>% summarize_all(mean)

# build the plot using both info
ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_class != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.class),col="black")+
  theme_cowplot()
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAP_class_075.pdf"),height = 6,width = 7)
# ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass.pdf",width = 8,height = 6)

ggplot()+
  geom_point(data = data_ref,aes(x = refumap_1,y = refumap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = data_ref2 %>%
               filter(robust_class != "uncertain"),aes(x = UMAP_1,y = UMAP_2, col = predicted.class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5),ncol=1)) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = data_ref_avg2subclass,aes(x = UMAP_1,y = UMAP_2,label = predicted.class),col="black")+
  theme_cowplot()+
  facet_wrap(~var_split,ncol=4)+
  theme(strip.background = element_blank())
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAP_class_split_075.pdf"),height = 6,width = 14)
# ggsave("../out/plot/human_motorcortexv1.0.0_subclass_data.combined_annotated_norm_fix_resolution_UMAP_075_subclass_split.pdf",width = 20,height = 15)

# inverse transfer, see the ref on the query ------------------------------
# plot the annotation over the original data
# read in the original dataset
# "../out/test_introns/object/data.combined_harmonySkipIntegration_AllSoupX_01000_06000_15_V3_diet.rds"
sobj <- readRDS("../../out/object/analysis_R44/27_LYM_subcluster_HarmonySample.rds")

# load the meta with annotations
meta_query <- read_tsv(paste0("../../out/table/analysis_R45/28_",sample_id,"_meta_Azimuth_",reference_id,"_",class,".tsv")) %>%
  mutate(
    robust_class = case_when(.data[[col_class_score]] > 0.75 & mapping.score > 0.75 ~ .data[[col_class]],
                             TRUE ~ "uncertain"),
    predicted.class = .data[[col_class]],
    predicted.class.score = .data[[col_class_score]],
    cluster.query = .data[[cluster_query]],
    var_split = .data[[var_split]])

test_UMAP_query <- sobj@reductions$umap@cell.embeddings %>% 
  data.frame() %>% 
  rownames_to_column("barcode")

test_data_query <- left_join(test_UMAP_query,meta_query,"barcode")
# add the groupign suggested by Martina
# data_query$group_martina <- case_when(data_query$orig.ident %in% c("pool_MOCK_CTRL")~"MOCK",
#                                       data_query$orig.ident %in% c("CTRL8_WT_CSF","CTRL8_WT_CTRL")~"CTRL",
#                                       T~"SOX10")

# divide the dataset into uncertain and not
test_data_query_unc <- test_data_query %>%
  filter(robust_class == "uncertain")

#
test_data_query_defined <- test_data_query %>%
  filter(robust_class != "uncertain")

# average the position of the clusters
test_data_query_avg <- test_data_query_defined %>% group_by(robust_class) %>% select(umap_1, umap_2) %>% summarize_all(mean)
test_data_query_avg2 <- test_data_query_defined %>% group_by(cluster.query) %>% select(umap_1, umap_2) %>% summarize_all(mean)

# build the plot using both info
ggplot()+
  geom_point(data = test_data_query_unc,aes(x = umap_1,y = umap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = test_data_query_defined,aes(x = umap_1,y = umap_2, col = robust_class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = test_data_query_avg,aes(x = umap_1,y = umap_2,label = robust_class),col="black")+
  theme_cowplot()
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAPinverse_class_075.pdf"),height = 6,width = 7)
# ggsave("../out/plot/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthBS_subclass_UMAP_075.pdf",height = 3,width = 5)

# build the plot using both info
ggplot()+
  geom_point(data = test_data_query_unc,aes(x = umap_1,y = umap_2),size=0.3,alpha=0.1,col="gray") +
  geom_point(data = test_data_query_defined,aes(x = umap_1,y = umap_2, col = robust_class),size=0.3,alpha=0.8) +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = test_data_query_avg,aes(x = umap_1,y = umap_2,label = robust_class),col="black")+
  theme_cowplot() +
  facet_wrap(~var_split,ncol=4)+
  theme(strip.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA))
ggsave(paste0("../../out/plot/analysis_R45/28_",sample_id,"_",class,"_",reference_id,"_UMAPinverse_class_split_075.pdf"),height = 6,width = 14)
#ggsave("../../out/image/harmonySkipIntegration_AllSoupX_01000_06000_15_CoordUMAP_AzimuthBS_subclass_split_UMAP_075.pdf",height = 15,width = 20)