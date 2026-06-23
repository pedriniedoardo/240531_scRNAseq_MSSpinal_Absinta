# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions
# Aletta suggesetd to add also the forest plot for the estimates of the statistics
# I will be using scProportionTest as tool for generating the estimate

# libraries ---------------------------------------------------------------
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(scProportionTest)

# read in the data --------------------------------------------------------
scobj <- readRDS("../../out/object/analysis_R44/001_sobj_integrated_cleanup_manualAnnotation_updateDemyelination.rds")

# confirm the object and annotations
DimPlot(scobj,group.by="cell_id")

# wrangling ---------------------------------------------------------------
# the tool needs to convert the datafarme into its own format, but add the barcodes as cell_id.
# rename the cellid
scobj$cellid <- scobj$cell_id

meta_ref <- scobj@meta.data %>%
  rownames_to_column()

# generate the correct object for the test
meta_ref2 <- sc_utils(scobj)

# confirm the numbers from the tissue dataset
meta_ref %>% 
  group_by(sample_id,cellid,diagnosis) %>% 
  summarise(n = n())

# run the proportion test diagnosis ---------------------------------------
# check the covariates
meta_ref$sample_id %>% table()
meta_ref$diagnosis_short %>% table()

# the function run the test camparing
# arguments
# sc_utils_obj: sc_utils object
# cluster_identity: Column that has cluster names
# sample_1: First sample to compare (ie. control)
# sample_2: Sample to compare to first sample (ie. treatment)
# sample_identity: Column that has sample names
# n_permutations: Number of permutations
prop_test_results <- permutation_test(meta_ref2,
                                      cluster_identity = "cellid",
                                      sample_1 = "CTRL",
                                      sample_2 = "MS",
                                      sample_identity = "diagnosis_short")

# extract the estimates
prop_test_results@results

# plot the data
permutation_plot(prop_test_results)
