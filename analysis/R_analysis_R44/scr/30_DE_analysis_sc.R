# AIM ---------------------------------------------------------------------
# run the DE analysis on all the levels of the main annotation for the comparion MS vs CTRL

# libraries ---------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(limma)
library(ggrepel)
library(presto)

# specify the version of Seurat Assay -------------------------------------
# set seurat compatible with seurat4 workflow
# options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
data.combined <- readRDS("../../out/object/analysis_R44/26_sobj_integrated_cleanup_manualAnnotation.rds")

# check the object version
class(data.combined@assays$RNA)

# compare the clusters per treatment --------------------------------------
# make suer the correct defalult dataset il loaded should be RNA
DefaultAssay(data.combined)

# check that in the RNA slot the data object is indeed loaded with normalized values
data.combined@assays$RNA$data[1:10,1:10]

# define the grouping variables for the comparison of stim
head(data.combined@meta.data)

# define the new grouping
table(data.combined@meta.data$cell_id,data.combined@meta.data$diagnosis_short)
table(data.combined@meta.data$sample_id,data.combined@meta.data$diagnosis_short)

data.combined$condition.annotation <- paste(data.combined$diagnosis_short,data.combined$cell_id, sep = "_")
head(data.combined@meta.data)

# update the idents of the object
Idents(data.combined) <- "condition.annotation"

# avg_logFC: log fold-chage of the average expression between the two groups. Positive values indicate that the gene is more highly expressed in the first group
# run the differential expression over the WT
# clusters_id <- as.character(sort(unique(data.combined$annotation_confident))) %>%
#   str_subset(pattern = "P09",negate = T)
clusters_id <- as.character(sort(unique(data.combined$cell_id)))

# NOLD vs control all
# run presto to make the comparison faster
list_MS_vs_CTRL <- lapply(clusters_id,function(x){
  id_1 <- paste0("MS_",x)
  id_2 <- paste0("CTRL_",x)
  # latest version of seurat might have implemented the presto methods by default, therefore the regular FindMarkers function might work well
  # response <- FindMarkers(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response <- RunPresto(data.combined, ident.1 = id_1, ident.2 = id_2, verbose = T,logfc.threshold = 0)
  response %>%
    rownames_to_column("gene") %>%
    mutate(id_1 = id_1,
           id_2 = id_2) %>%
    mutate(cluster = x)
})

list_MS_vs_CTRL %>%
  setNames(paste0(clusters_id,"_MS_vs_CTRL")) %>%
  bind_rows(.id = "annotation") %>%
  write_tsv("../../out/table/analysis_R44/30_response_MS_vs_CTRL_seurat_sc.tsv")

# plot the volcanos per cluster -------------------------------------------
folder <- "../../out/table/analysis_R44/"
file <- dir(folder) %>%
  str_subset(pattern = "30_response_MS_vs_CTRL_seurat_sc.tsv")

df_res <- lapply(file, function(x){
  test_plot <- read_tsv(paste0(folder,x))
}) %>%
  bind_rows()

# show the distribution of the pvalues
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "MS_vs_CTRL")) %>%
  ggplot(aes(x = p_val)) + geom_histogram()+theme_bw()+facet_wrap(~cluster)+theme(strip.background = element_blank())
ggsave("../../out/plot/analysis_R44/30_dist_p_value_MS_vs_CTRL_seurat_sc.pdf",width = 12,height = 9)

# render all of them as a volcano plot
test_significant <- df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "MS_vs_CTRL")) %>%
  mutate(threshold = case_when(abs(avg_log2FC) > 1 & p_val_adj<0.05~1,
                               T~0)) %>%
  filter(threshold == 1)

# library(ggrepel)
df_res %>%
  mutate(comparison = str_extract(annotation,pattern = "MS_vs_CTRL")) %>%
  # filter(symbol %in% setdiff(GOI_SLC,GOI)) %>%
  ggplot(aes(x = avg_log2FC,y = -log(p_val_adj))) +
  geom_point(alpha = 0.01) +
  geom_point(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj)),col="red",alpha = 0.5) +
  geom_vline(xintercept = c(-1,1),col="gray",linetype="dashed")+
  geom_hline(yintercept = (-log(0.05)),col="gray",linetype="dashed")+
  geom_text_repel(data = test_significant,aes(x = avg_log2FC,y = -log(p_val_adj),label = gene)) +
  facet_wrap(~cluster) +
  theme_bw()+theme(strip.background = element_blank())
ggsave("../../out/plot/analysis_R44/30_volcano_MS_vs_CTRL_seurat_sc.pdf",width = 16,height = 12)

# sample plotting specific GOI --------------------------------------------
GOI <- c("IGKC")
# sGOI <- c("IFITM3", "BST2", "MT2A", "PLSCR1", "BTG1", "MX2", "MX1", "FPR1", "CASP1")

# wrangling ---------------------------------------------------------------
# tidy up the metadata
data.combined$sample_id <- str_replace_all(data.combined$sample_id,pattern = "_","-")
# data.combined$cell_id <- str_replace_all(data.combined$cell_id,pattern = "\\s","-")

# EDA ---------------------------------------------------------------------
# check the global expression
VlnPlot(object = data.combined,features = GOI,group.by = "cell_id")
# ggsave("../../out/plot/VlnProcr_cluster.pdf",width = 6,height = 4)

# manual plotting ---------------------------------------------------------
# extract the expression data
df_exp <- FetchData(data.combined, vars = GOI,layer = "data") |> 
  rownames_to_column("barcodes") |> 
  pivot_longer(names_to = "gene",values_to = "count",-barcodes) |> 
  # try to min/max normalize the count varaible per gene in order to rescale the difference in terms of expression
  group_by(gene) %>%
  # threshold of positiveness is based on the distriubtion of the expression of the signal in tihs case
  mutate(norm_min_max = ((count - min(count))/(max(count)-min(count))),
         exp_cat = case_when(count > 0~"pos",
                             T~"neg")) %>%
  ungroup() %>%
  mutate(count_fix = count + rnorm(nrow(.))/100000)
# separate(barcodes,into = c("barcode","barcode_id"),sep = "-",remove = F)

head(df_exp)

exp_sc_wide <- df_exp %>%
  dplyr::select(barcodes,gene,count_fix) %>%
  pivot_wider(names_from = gene,values_from = count_fix)

# ggpairs(exp_sc_wide[,-1])

# average expression
# build the grouping variable
group_id_treat <- paste(data.combined@meta.data$sample_id,
                        data.combined@meta.data$cell_id,
                        data.combined@meta.data$diagnosis_short,
                        sep = "|")
# add it to the meta
data.combined$group_id_treat <- group_id_treat

# set the idents
Idents(data.combined) <- "group_id_treat"
# DefaultAssay(data.combined) <- "RNA"

df_average_treat <- AverageExpression(data.combined,features = GOI)$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# build the pattern to extract the metadata
pattern_cellid <- paste0(unique(data.combined@meta.data$cell_id),collapse = "|")
pattern_sample <- paste0(unique(data.combined@meta.data$sample_id),collapse = "|")
pattern_treat <- paste0(unique(data.combined@meta.data$diagnosis_short),collapse = "|")

df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = cellid,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(treat~gene)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/plot/analysis_R44/30_DE_sc_GOI01_treat.pdf",width = 8,height = 4)

# try to plot the panel facetting by cell_id in order to rescale the range
df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = treat,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(cellid~gene,scales="free")+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/plot/analysis_R44/30_DE_sc_GOI01_treat_01.pdf",width = 9,height = 8)

# try to color by treat rather than splitting
df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = cellid,y=avg_exp,col=treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw() +
  facet_wrap(~gene)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/plot/analysis_R44/30_DE_sc_GOI01_treat_02.pdf",width = 4,height = 4)

# test if it is significant
df_average_treat %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  lm(data = .,avg_exp~cellid+treat) %>%
  summary()

# -------------------------------------------------------------------------
# Aletta suggested to also explore a panel of known markers of B cells to see if they have a trend in MS
GOIs <- c("CD38","IGKC","IGHG1","MZB1","CD79A")

df_average_treat2 <- AverageExpression(data.combined,features = GOIs)$RNA %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(names_to = "group_id",values_to = "avg_exp",-gene)

# try to plot the panel facetting by cell_id in order to rescale the range
df_average_treat2 %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = treat,y=avg_exp)) + geom_boxplot(outlier.shape = NA) + geom_point(position = position_jitter(width = 0.2),shape = 1)+
  theme_bw() +
  facet_wrap(cellid~gene,scales="free",ncol=5)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/plot/analysis_R44/30_DE_sc_GOI02_treat_01.pdf",width = 10,height = 22)

# try to color by treat rather than splitting
df_average_treat2 %>%
  mutate(cellid = str_extract_all(group_id,pattern = pattern_cellid) %>% unlist()) %>%
  mutate(sample = str_extract_all(group_id,pattern = pattern_sample) %>% unlist()) %>%
  mutate(treat = str_extract_all(group_id,pattern = pattern_treat) %>% unlist()) %>%
  # group_by(RNA_snn_res.0.1) %>%
  # mutate(group_avg = mean(avg_exp)) %>%
  # ungroup() %>%
  # mutate(RNA_snn_res.0.1 = fct_reorder(RNA_snn_res.0.1,group_avg,.desc = T)) %>%
  ggplot(aes(x = cellid,y=avg_exp,col=treat)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8),shape = 1)+
  theme_bw() +
  facet_wrap(~gene,scales="free",nrow = 1)+
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 90,hjust = 1))
ggsave("../../out/plot/analysis_R44/30_DE_sc_GOI02_treat_02.pdf",width = 11,height = 3)
