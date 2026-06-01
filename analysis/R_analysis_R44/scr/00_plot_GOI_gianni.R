# AIM ---------------------------------------------------------------------
# explore the expression of some GOIs in the final dataset
# these are the genes asked by Gianni

# libraries ---------------------------------------------------------------
library(Seurat)
library(SeuratWrappers)
library(tidyverse)
library(scales)
library(ComplexHeatmap)
library(lemon)
library(finalfit)
library(cowplot)
library(patchwork)
library(Nebulosa)

# read in the dataset -----------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")

# specify the covariate of the annotation
cov_anno <- "cell_id"
cov_pathology <- "pathological_stage"
cov_sample <- "sample_id"

# check the object version
class(data.combined@assays$RNA)

# check the dim reduction of the object of interest
DimPlot(data.combined,label = T,raster = T,group.by = cov_anno)
ggsave(paste0("../../out/plot/analysis_R44/00_UMAPSeurat_",cov_anno,".pdf"),width = 6,height = 5)

# define the gene of interest GOI
# GOI <- c("Irf7","Ddx58")
GOI <- c("LHFPL3")

table(data.combined@meta.data$expertAnno.l1)

# add a broader cell grouping
# data.combined@meta.data$cell_type2 <- data.combined@meta.data |> 
#   mutate(cell_type2 = case_when(# according to label transfer cluster 12 and 6 are Bipolar cells
#     cell_type %in% c("BP_0","BP_10","BP_11","BP_13","BP_16","BP_17","BP_5","BP_8","12","6")~"BP",
#     T~cell_type)) |> 
#   pull(cell_type2)


# run DE analysis at the annotation level ---------------------------------
# # run the DE analysis on the defined annotation
# DefaultAssay(data.combined)<-"RNA"
# Idents(data.combined) <- cov_anno
# 
# # define the filters
# filt.pct <- 0.25
# filt.logfc <- 0.25
# 
# # calculate the markers
# data.combined.markers <- FindAllMarkers(data.combined,
#                                         only.pos = TRUE,
#                                         min.pct = filt.pct,
#                                         logfc.threshold = filt.logfc)
# 
# # save the table of all markers
# data.combined.markers %>%
#   write_tsv(paste0("../../out/table/analysis_R44/00_FindAllMarkers_",cov_anno,"_pct",filt.pct,"_logfc",filt.logfc,".tsv"))
# 
# # top 100 per cluster
# data.combined.markers %>%
#   group_by(cluster) %>%
#   mutate(rank = rank(order(p_val_adj, -abs(avg_log2FC)), ties.method='first')) %>%
#   arrange(cluster,rank) %>%
#   filter(rank < 101) %>%
#   write_tsv(paste0("../../out/table/analysis_R44/00_FindAllMarkers_",cov_anno,"_pct",filt.pct,"_logfc",filt.logfc,"_top100.tsv"))

# read in the output of the DE at the level of the cell anntation
data.combined.markers <- read_tsv("../../out/table/analysis_R44/00_FindAllMarkers_cell_id_pct0.25_logfc0.25.tsv")

# load the dataset and check if the gene is expressed in the top markers per cluster
data.combined.markers %>%
  filter(gene %in% GOI)

# generate the table for the plots ----------------------------------------
# get the metadata from the other object
meta <- data.combined@meta.data %>%
  rownames_to_column(var = "barcodes")

# extrac the expression value
df_exp <- FetchData(data.combined, vars = GOI,layer = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "exp",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(exp_min_max = ((exp - min(exp))/(max(exp)-min(exp))),
         exp_cat = case_when(exp > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(exp_fix = exp + rnorm(nrow(.))/100000)

# get the coordinates
UMAP1_df <- data.combined@reductions$umap@cell.embeddings %>%
  data.frame() %>%
  rownames_to_column(var = "barcodes")

# generate the dataset for mapping the data in the umamp
dim(UMAP1_df)
dim(df_exp)
dim(meta)

# put everithing in a single dataset
df_tot <- reduce(list(meta,UMAP1_df,df_exp),left_join, by="barcodes")
df_tot_avg <- df_tot %>%
  group_by(.data[[cov_anno]]) %>%
  dplyr::select(umap_1, umap_2) %>%
  summarize_all(mean)

dim(df_tot)

head(data.combined@meta.data)

# plot the average expression per sample use the variable cell tyep per donor as grouping
# data.combined$group <- paste0(data.combined$orig.ident,".",data.combined$cell_type2)
data.combined$group <- paste0(data.combined@meta.data[[cov_pathology]],"|",
                              data.combined@meta.data[[cov_anno]],"|",
                              data.combined@meta.data[[cov_sample]])
# data.combined$group2 <- paste0(data.combined$orig.ident,".",data.combined$treat,".",data.combined$cell_type2)
Idents(data.combined) <- "group"
DefaultAssay(data.combined) <- "RNA"

average_GOI <- AverageExpression(data.combined,features = GOI,group.by = c("group"))

# plot general UMAP -------------------------------------------------------
# build the plot using both info
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = umap_1,y = umap_2,col=.data[[cov_anno]]),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  ggrepel::geom_text_repel(data = df_tot_avg,aes(x = umap_1,y = umap_2,label = .data[[cov_anno]]),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident.pdf",width = 7,height = 5)

# no lab
ggplot(label= TRUE) +
  # geom_point(data = data2,aes(x = UMAP_1,y = UMAP_2,col=cell_type),size=0.3,alpha=0.1) +
  geom_point(data = df_tot,aes(x = umap_1,y = umap_2,col=.data[[cov_anno]]),size=0.3) +
  # geom_point(data = data2_unc,aes(x = UMAP_1,y = UMAP_2),size=0.3,alpha=0.1,col="gray") +
  # geom_point(data = data2_defined,aes(x = UMAP_1,y = UMAP_2, col = robust_score),size=0.3,alpha=0.8) +
  # labs(color= "Clusters") +
  # ggrepel::geom_text_repel(data = df_tot_avg,aes(x = UMAP_1,y = UMAP_2,label = cell_type2),col="black")+
  guides(colour = guide_legend(override.aes = list(size=5)))+
  theme_bw()
# facet_wrap(~infection)
# ggsave("../../out/image/ManualClean/UMAP_38_annotationConfident_noLab.pdf",width = 7,height = 5)

# expression distribution -------------------------------------------------
# crop the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_grid(~gene)+theme_bw()+scale_x_log10()+geom_vline(xintercept = 1,col="red",linetype="dotted")

# keep the 0 expressing cells
df_exp %>%
  ggplot(aes(x=exp))+geom_histogram()+facet_wrap(~gene)+theme_bw()+
  # scale_x_log10()+
  geom_vline(xintercept = 2.5,col="red",linetype="dotted")

# library(scales)
# show_col(c("#4662D7FF","#FABA39FF","#7A0403FF"))

# plotting expression -----------------------------------------------------
# by counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~.data[[cov_pathology]]) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_count.pdf",width = 13,height = 12)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,split.by = cov_pathology,raster = T,order = T,ncol = 3)
# ggsave("../../out/image/06_UMAPSeurat_annotationConfident_DNMT3A_brain_count.pdf",width = 25,height = 3)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(~gene) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_count2.pdf",width = 6,height = 5)

# do the same using Seurat
FeaturePlot(data.combined,features = GOI,raster = T,order = T)
ggsave("../../out/plot/analysis_R44/00_UMAPSeurat_GOI_Gainni.pdf",width = 6,height = 5)

# try nebulosa option
plot_density(data.combined, GOI,reduction = "umap") + scale_color_viridis_c(option = "turbo")
ggsave("../../out/plot/analysis_R44/00_UMAPSeurat_GOI_Gainni_nebulosa.pdf",width = 6,height = 5)
# ggsave("../../out/image/00_UMAPSeurat_annotationConfident_Valentina_brain_countNebulosa.pdf",width = 6,height = 5)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp)) + geom_point(alpha = 0.5,size = 0.05) +
  facet_wrap(gene~.data[[cov_pathology]]) +
  theme_cowplot() +
  # scale_color_gradient(low = "gray",high = "blue") +
  scale_color_viridis_c(option = "turbo") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_count_alt.pdf",width = 13,height = 12)

# by min max normalized counts
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp_min_max)) + geom_point(alpha = 0.5,size = 0.2) +
  facet_wrap(gene~.data[[cov_pathology]]) +
  theme_cowplot() +
  scale_color_gradient(low = "gray",high = "blue") +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_minmax.pdf",width = 13,height = 12)

# plot the category. being 0 or non zero per cell
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(gene~.data[[cov_pathology]])+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_proppos.pdf",width = 13,height = 12)

df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  mutate(exp_cat = factor(exp_cat,levels = c("neg","pos"))) %>%
  arrange(exp_cat) %>%
  ggplot(aes(x = umap_1, y = umap_2,col = exp_cat)) + geom_point(alpha = 0.5,size = 0.05) +
  # facet_wrap(gene~NMDA_time,nrow = 2) +
  # facet_rep_wrap(gene~treat,repeat.tick.labels = "all",nrow=3)+
  facet_wrap(~gene)+
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_cowplot() +
  scale_color_manual(values = c("gray","blue")) +
  # theme(strip.background = element_blank(),
  #       panel.border = element_rect(colour = "black", fill = NA))
  theme(strip.background = element_blank())
# ggsave("../../out/image/00_UMAPggplot_annotationConfident_Valentina_brain_proppos2.pdf",width = 6,height = 5)

# violin plot for GOI expression use macro categories
df_tot %>%
  # filter(NMDA_time%in%c("NMDA_00","NMDA_03","NMDA_06","NMDA_12","NMDA_24","NMDA_36")) %>%
  # this is the processing shown in the violinplot function
  # mutate(exp_fix = exp + rnorm(nrow(.))/100000) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=.data[[cov_pathology]],y=exp_fix)) + 
  geom_violin(scale = "width")+
  geom_point(position=position_jitter(width = 0.2),alpha=0.01) +
  facet_wrap(~.data[[cov_anno]]) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,angle = 90)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
# ggsave("../../out/image/00_violin_annotationConfident_Valentina_brain.pdf",width = 15,height = 10)

# try to depict the average expression there is roughly one sample per condition
pattern_pathology <- data.combined@meta.data[[cov_pathology]] %>% unique() %>% str_replace_all(pattern = "\\s|_",replacement = "\\.") %>% paste0(collapse = "|")
pattern_anno <- data.combined@meta.data[[cov_anno]] %>% unique() %>% str_replace_all(pattern = "\\s|_",replacement = "\\.") %>% paste0(collapse = "|")
pattern_sample <- data.combined@meta.data[[cov_sample]] %>% unique() %>% str_replace_all(pattern = "\\s|_",replacement = "\\.") %>% paste0(collapse = "|")

df_avg <- average_GOI$RNA %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(gene = GOI) %>%
  pivot_longer(names_to = "group",values_to = "avg_exp",-gene) %>%
  # filter(!str_detect(group,pattern="doublet|unassigned")) |> 
  mutate(path_anno = str_extract(group,pattern = pattern_pathology)) %>%
  mutate(donor = str_extract(group,pattern = pattern_sample)) %>%
  mutate(cell_anno = str_extract(group,pattern = pattern_anno))

# plot the average expresison by cell annotation
df_avg %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=cell_anno,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free")+
  scale_y_continuous(trans = "log1p")
ggsave("../../out/plot/analysis_R44/00_boxplot_avgExp_manualAnnotation_GOI_Gainni.pdf",width = 6,height = 6)

df_avg %>%
  filter(cell_anno == "ASTRO") %>%
  group_by(path_anno) %>%
  summarise(n = n())

# do the same as above but split by condition
df_avg %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=path_anno,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~gene,scales = "free") +
  scale_y_continuous(trans = "log1p")

# plot splitting by treat full
df_avg %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=path_anno,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(~cell_anno,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/06_dotplot_annotationConfident_DNMT3A_brain_expressionAvg_treatFull.pdf",width = 9,height = 9)

# plot splitting by treat
df_avg %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=path_anno,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(gene~cell_anno,scales = "free")
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_Valentina_brain_expressionAvg_treat.pdf",width = 25,height = 15)

# try to keep the same scale
# calculate the median per annotation
df_avg_summary <- df_avg %>%
  group_by(cell_anno) %>%
  summarise(med = median(avg_exp)) %>%
  ungroup() %>%
  mutate(cell_anno = fct_reorder(cell_anno,desc(med)))

df_avg %>%
  mutate(cell_anno = factor(cell_anno,levels = levels(df_avg_summary$cell_anno))) %>%
  # ggplot(aes(x=NMDA_time,y=count)) + 
  ggplot(aes(x=path_anno,y=avg_exp))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position = position_jitter(width = 0.1),alpha = 0.6)+
  geom_hline(data = df_avg_summary,aes(yintercept = med),col="red",linetype="dashed") +
  # geom_col()+
  # facet_wrap(~cell_type2,scales = "free")+
  theme_bw()+theme(axis.text.x = element_text(hjust = 1,angle = 90))+
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))+
  facet_wrap(gene~cell_anno)
# scale_fill_viridis_c(option = "plasma",name="log10 number \nof cells")
# ggsave("../../out/image/00_dotplot_annotationConfident_Valentina_brain_expressionAvg_treat_scale.pdf",width = 25,height = 15)
