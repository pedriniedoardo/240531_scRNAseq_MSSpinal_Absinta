# AIM ---------------------------------------------------------------------
# Apply the filter stated in the label of the file.
# This is either a testing of the filtering results, or the actual filtering choice to go forward

# libraries ---------------------------------------------------------------
library(scater)
library(Seurat)
library(tidyverse)
library(robustbase)
library(patchwork)
library(DoubletFinder)

# define the filtering parameters -----------------------------------------
featureLow_thr <- 500
featureHigh_thr <- 7000
mito_thr <- 5
label <- "00500_07000_05"

# read in the data --------------------------------------------------------
# location of all the raw matrices
id_sample <- dir("../../data/SoupX_default_cellranger7/")

# load the LUT of the dataset
LUT <- read_csv("../../data/LUT_sample.csv")

# load the LUT of the doublet rate estimate
df_doublet <- read_csv("../../data/dublets_rate_2023.csv")

# run the processing ------------------------------------------------------
# do the preprocessing over all the dataset and save the objects
# x <- "Sample_1"
# id_sample <- "Sample_1"
list_datasc <- lapply(id_sample,function(x){
  # to track the processing of the progress of the lapply
  print(x)
  data <- Read10X(data.dir = paste0("../../data/SoupX_default_cellranger7/",x))
  
  datasc <- CreateSeuratObject(counts = data, project = LUT %>%
                                 filter(sample_id == x) %>%
                                 pull(sample), min.cells = 20, min.features = 200)
  
  # datasc <- CreateSeuratObject(counts = data, min.cells = 20, min.features = 200)
  
  # datasc@meta.data
  datasc$percent.mt <- PercentageFeatureSet(datasc, pattern = "^MT-")
  datasc$percent.ribo <- PercentageFeatureSet(datasc, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  # add also the percentage of globin. in this dataset it is not meaningful as there is no blood
  datasc$percent.globin <- Seurat::PercentageFeatureSet(datasc,pattern = "^HB[^(P)]")
  
  # label the cells based on the mt reads content
  datasc$mt_bin <- datasc@meta.data %>%
    mutate(test = case_when(percent.mt < 1~"low",
                            percent.mt < 10~"mid",
                            T ~ "high")) %>%
    pull(test)
  
  # datasc$treat <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(treat)
  # 
  # datasc$clone <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(clone)
  
  # datasc$doxy <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(doxy)
  
  # datasc$exposure <- LUT %>%
  #   filter(sample_id == x) %>%
  #   pull(exposure)
  
  # datasc$ID <- LUT %>%
  #   filter(sample == x) %>%
  #   pull(ID)
  
  # add the filtering variable based on the fixed threshold
  datasc$test <- datasc@meta.data %>%
    # mutate(test = percent.mt < 20 & nFeature_RNA > 700 & nFeature_RNA < 9000) %>%
    # mutate(test = percent.mt > 1 & percent.mt < 10 & nFeature_RNA > 1000 & nFeature_RNA < 9000) %>%
    mutate(test = percent.mt < mito_thr & nFeature_RNA > featureLow_thr & nFeature_RNA < featureHigh_thr) %>% 
    pull(test)
  
  # add the filtering variable based on the
  stats <- cbind(log10(datasc@meta.data$nCount_RNA), log10(datasc@meta.data$nFeature_RNA),
                 datasc@meta.data$percent.mt)
  
  # library(robustbase)
  outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
  #library(scater)
  multi.outlier <- isOutlier(outlying, type = "higher")
  # summary(multi.outlier)
  
  datasc$not_outlier <- !as.vector(multi.outlier)
  
  return(datasc)
}) %>%
  setNames(id_sample)

# before filtering --------------------------------------------------------
# this step is not modifine the object and is already run tin the 00_explore_fixed_threshold
meta_total <- lapply(list_datasc, function(x){
  x@meta.data %>%
    rownames_to_column("barcode") %>%
    mutate(barcode = paste0(barcode,"|",orig.ident))
}) %>%
  bind_rows(.id = "dataset")

# count the cells that will pass the filter
meta_total %>%
  dplyr::count(orig.ident)

# filtering of the cells --------------------------------------------------
# perform the filtering based on the fixed threshold
list_datasc_fixed <- lapply(list_datasc,function(x){
  datasc_filter <- subset(x, subset = test == 1)
})

# Pre-process the data for DoubletFinder ----------------------------------
# use the Normalize function. for the doublet removla I need to run the full preprocessing steps
list_datasc_fixed_norm <- lapply(list_datasc_fixed, function(x){
  datasc_filter <- x
  datasc_filter <- NormalizeData(datasc_filter)
  datasc_filter <- FindVariableFeatures(datasc_filter, selection.method = "vst", nfeatures = 2000)
  datasc_filter <- ScaleData(datasc_filter) 
  datasc_filter <- RunPCA(datasc_filter) 
  datasc_filter <- FindNeighbors(datasc_filter,dims = 1:30) 
  datasc_filter <- FindClusters(datasc_filter) 
  datasc_filter <- RunUMAP(datasc_filter,dims = 1:30) 
  
  datasc_filter
})

# DoubletFinder -----------------------------------------------------------
# determine the number of cell recovered
# df_doublet <- read_csv("data/doublets.csv")
# notice that I have update the file with doublets rate
# df_doublet <- read_csv("../../data/dublets_rate_2023.csv")

# remove the sample 6 from the analysis
# list_datasc_fixed_norm_fix <- list_datasc_fixed_norm[!names(list_datasc_fixed_norm)%in%c("06_cr_61")]
list_datasc_fixed_norm_fix <- list_datasc_fixed_norm

# needed for barcode.csv --------------------------------------------------
# list_nExp <- lapply(list_datasc_fixed_norm_fix,function(x){
#   recover_cell <- round(dim(x)[2]/1000,digits = 0)*1000
#   recover_cell
#   ifelse(test = recover_cell > 10000,
#          yes = 0.076,
#          no = df_doublet[df_doublet$CellRecovered == recover_cell,][["MultipletRate"]]
#          )
# })
# -------------------------------------------------------------------------
# notice that in this case we should evaluate carefully how many doublets to assign to sample prepared with more than 30k cells. In this cases we are using the same rate applied to 30k
list_nExp <- lapply(list_datasc_fixed_norm_fix,function(x){
  recover_cell <- dim(x)[2]
  recover_cell
  # pick the rate with the lowest distance form the reference
  df_delta <- df_doublet %>% 
    mutate(abs_delta = abs(CellRecovered-recover_cell)) %>% 
    arrange(abs_delta)
  
  df_delta$MultipletRate[1]
})

# runt he simulation for doublet finder
list_datasc_fixed_norm_doublet <- pmap(list(list_datasc_fixed_norm_fix,list_nExp),function(x,y){
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep_v3(x, PCs = 1:30, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE) 
  bcmvn <- find.pK(sweep.stats) 
  # plot the results to justify the pK choice 
  # bcmvn %>% 
  #   ggplot(aes(pK,BCmetric,group=1)) + 
  #   geom_point()+ 
  #   geom_line() 
  # save the max BCmetrics 
  pK <- bcmvn %>% 
    filter(BCmetric == max(BCmetric)) %>% 
    pull(pK) %>% 
    as.character() %>% 
    as.numeric() 
  
  # homotipic doublet proportion estimate ----------------------------------- 
  annotation <- x$seurat_clusters 
  homotipic.prop <- modelHomotypic(annotation) 
  
  # 0.076 is the expected nbumber of doublets based on the number of cells recovered 
  nExp.poi <- round(y*nrow(x@meta.data)) 
  nExp.poi.adj <- round(nExp.poi*(1-homotipic.prop)) 
  
  # run the DoubletFinder --------------------------------------------------- 
  pbmc.seurat.filteres <- doubletFinder_v3(x,PCs = 1:30,pN = 0.25,pK = pK,nExp = nExp.poi.adj,reuse.pANN = F,sct = F)
  pbmc.seurat.filteres
})

# trimm the datasets and save the filtered objects ------------------------
colnames(list_datasc_fixed_norm_doublet$Sample_1@meta.data)
# "orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN_0.25_0.005_441","DF.classifications_0.25_0.005_441"
# "orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF"

# count the cells per object
lapply(list_datasc_fixed_norm_doublet, function(x){
  meta <- x@meta.data
  c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  x@meta.data <- meta
  head(x@meta.data)
  dim(x@meta.data)[1]
})

# save the individula filtered and normalized objects
list_norm_doublet <- pmap(list(names(list_datasc_fixed_norm_doublet),list_datasc_fixed_norm_doublet),function(x,y){
  # fix the meta and save the object
  meta <- y@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  y@meta.data <- meta
  # head(y@meta.data)
  return(y)
}) %>% 
  setNames(names(list_datasc_fixed_norm_doublet))

saveRDS(object = list_norm_doublet,file = paste0("../../out/object/list_datasc_fix_filter_norm_doublet_SoupX_",label,".rds"))

# save the taola meta with the doublet imputation
meta_total_doublet <- lapply(list_datasc_fixed_norm_doublet, function(x){
  meta <- x@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  x@meta.data <- meta
  # head(x@meta.data)
  x@meta.data
}) %>%
  bind_rows(.id = "dataset") %>%
  rownames_to_column("barcode")

# save the total meta
meta_total_doublet %>%
  write_tsv(paste0("../../out/table/meta_datasc_fix_filter_norm_doublet_SoupX_",label,".tsv"))

# subset the datasets after removal of the doublets befreo the integration
list_norm_singlets <- pmap(list(names(list_datasc_fixed_norm_doublet),list_datasc_fixed_norm_doublet),function(x,y){
  # fix the meta
  meta <- y@meta.data
  colnames(meta) <- c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.globin","mt_bin","test","not_outlier","RNA_snn_res.0.8","seurat_clusters","pANN","DF")
  y@meta.data <- meta
  # head(y@meta.data)
  # filter only the singlets
  datasc_filter <- subset(y, subset = DF == "Singlet")
  return(datasc_filter)
}) %>%
  setNames(names(list_datasc_fixed_norm_doublet))

saveRDS(object = list_norm_singlets,file = paste0("../../out/object/list_datasc_fix_filter_norm_doubletSinglet_SoupX_",label,".rds"))

# # use the SCTransform function
# list_datasc_fixed_SC <- lapply(list_datasc_fixed, function(x){
#   datasc_filter <- x
#   datasc_filterSC <- SCTransform(datasc_filter, vars.to.regress = c("percent.mt","nCount_RNA"), verbose = FALSE)
#   
#   datasc_filterSC
# })
# 
# # save the individula filtered and normzlized objects
# pmap(list(names(list_datasc_fixed_SC),list_datasc_fixed_SC),function(x,y){
#   saveRDS(object = y,file = paste0("out/datasc_fix_filter_SCnorm_",x,".rds"))
# })

# -------------------------------------------------------------------------
# test <- readRDS("out/datasc_fix_filter_norm_PN0242_0004.rds")
# test <- readRDS("out/datasc_fix_filter_SCnorm_PN0242_0004.rds")
