# AIM ---------------------------------------------------------------------
# run a statistical test on the proportion differences across clusters, between conditions
# Aletta suggesetd to add also the forest plot for the estimates of the statistics
# I will be using scProportionTest as tool for generating the estimate

# libraries ---------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(speckle)
library(limma)
library(statmod)
library(cowplot)
library(ggrepel)
library(finalfit)
library(scProportionTest)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(emmeans)

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

# refine test -------------------------------------------------------------
# let's try to implemente a better test for comparing proportions acocunting for the nesting of the single cells per sample
# make sure to include 0 in case of missing combinations fo cell tyepes per sample

df_summary <- meta_ref %>%
  group_by(sample_id,diagnosis_short,cell_id) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  # fill the missing combinations
  pivot_wider(names_from = cell_id,values_from = n,values_fill = 0) %>%
  pivot_longer(names_to = "cell_id",values_to = "n",-c(sample_id,diagnosis_short)) %>%
  group_by(sample_id) %>%
  mutate(tot = sum(n)) %>%
  ungroup() %>%
  mutate(prop = n/tot)

# for the beta family modelling tweak the absolute 0 and 1 into sometihng slightly below
# pull the minimumnb value
min_prop <- df_summary %>%
  filter(prop > 0) %>%
  summarise(min_prop = min(prop)) %>%
  pull(min_prop)

df_summary <- df_summary %>%
  mutate(prop_fix = case_when(prop == 0 ~ (0 + min_prop/100),
                          prop == 1 ~ (1 - min_prop/100),
                          T ~ prop))

# compare prop and prop_fix
df_summary %>%
  ggplot(aes(x = prop,y = prop_fix)) + geom_point()

# how many where changed
df_summary %>%
  filter(prop == 0 | prop == 1) %>%
  summarise(n = n())

# make the model using the information about the nesting of the sample
# Proportions are strictly bounded between 0 and 1 (or 0% and 100%). The Beta distribution is mathematically defined only on the open interval $(0, 1)$
# Because there are multiple cell types measured from the exact same physical sample, those data points are dependent on one another. To account for this sample-level grouping, use sample_id as a random intercept:

model <- glmmTMB(prop_fix ~ cell_id * diagnosis_short + (1 | sample_id),
                 data = df_summary,
                 family = beta_family(link = "logit"))

summary(model)

# make the stiamates base on the model
emm <- emmeans(model, ~ diagnosis_short | cell_id, type = "response")

# make the pairs split by the comparisons of interest
res <- pairs(emm,reverse = T)

# add the confidence interval
df_res <- left_join(
  res %>%
    # Adds lower.CL and upper.CL columns
    confint() %>%
    data.frame(),
  res %>%
    data.frame() %>%
    select(contrast,cell_id,null,z.ratio,p.value), by = c("contrast","cell_id")) %>%
  arrange(p.value)

# make the plot
df_res %>%
  mutate(cell_id = fct_reorder(cell_id,odds.ratio)) %>%
  ggplot(aes(y=cell_id,x=odds.ratio)) +
  geom_point() +
  geom_errorbar(aes(xmin = asymp.LCL,
                    xmax = asymp.UCL),width = 0.15, linewidth = 0.5) +
  geom_vline(col="red",linetype = "dashed",xintercept = 1) +
  theme_bw() +
  scale_x_log10()

# test an alternative model -----------------------------------------------
# usage of the betabinomial() family. cbind(successes, failures) is R's convention for specifying a binomial-style response when you have the raw counts. It's equivalent to modeling the proportion n/tot, but the model now knows the total and can properly weight the precision. A sample with 5000 total cells gives a much more certain proportion estimate than one with 50, and the betabinomial likelihood accounts for that automatically.
model2 <- glmmTMB(cbind(n, tot - n) ~ cell_id * diagnosis_short + (1 | sample_id),
                 data = df_summary,
                 family = betabinomial())

summary(model2)

# make the stiamates base on the model
emm2 <- emmeans(model2, ~ diagnosis_short | cell_id, type = "response")

# make the pairs split by the comparisons of interest
res2 <- pairs(emm2,reverse = T)

# add the confidence interval
df_res2 <- left_join(
  res2 %>%
    # Adds lower.CL and upper.CL columns
    confint() %>%
    data.frame(),
  res2 %>%
    data.frame() %>%
    select(contrast,cell_id,null,z.ratio,p.value), by = c("contrast","cell_id")) %>%
  arrange(p.value)

# make the plot
df_res2 %>%
  mutate(cell_id = fct_reorder(cell_id,odds.ratio)) %>%
  ggplot(aes(y=cell_id,x=odds.ratio)) +
  geom_point() +
  geom_errorbar(aes(xmin = asymp.LCL,
                    xmax = asymp.UCL),width = 0.15, linewidth = 0.5) +
  geom_vline(col="red",linetype = "dashed",xintercept = 1) +
  theme_bw() +
  scale_x_log10()

# compare with propeller estimate -----------------------------------------
# remember that this works on the full dataset

# group is the grouping id
res_propeller_disease <- propeller(clusters = meta_ref$cell_id,
                                   sample = meta_ref$sample_id,
                                   group = meta_ref$diagnosis_short)


# here is the implementation fot provide the confidence interval.

# -------------------------------------------------------------------------
# NOTICE: this processing below has been suggested by AI
# same logit-transformed proportions propeller uses internally
props <- getTransformedProps(clusters = meta_ref$cell_id,
                             sample  = meta_ref$sample_id,
                             transform = "logit")

# align sample order to group labels
sample_info <- data.frame(sample_id = colnames(props$TransformedProps)) %>%
  left_join(meta_ref %>% distinct(sample_id, diagnosis_short), by = "sample_id")

# fit limma (CTRL as reference)
group  <- factor(sample_info$diagnosis_short, levels = c("CTRL", "MS"))
design <- model.matrix(~ group)

fit <- lmFit(props$TransformedProps, design) %>% eBayes()

# logFC is on logit scale = log odds ratio MS vs CTRL
df_propeller_ci <- topTable(fit, coef = 2, n = Inf, confint = TRUE) %>%
  rownames_to_column("cell_id") %>%
  mutate(
    odds.ratio = exp(logFC),
    asymp.LCL  = exp(CI.L),
    asymp.UCL  = exp(CI.R)
  )

df_propeller_ci %>%
  mutate(cell_id = fct_reorder(cell_id, odds.ratio)) %>%
  ggplot(aes(y = cell_id, x = odds.ratio)) +
  geom_point() +
  geom_errorbar(aes(xmin = asymp.LCL, xmax = asymp.UCL), width = 0.15, linewidth = 0.5) +
  geom_vline(col = "red", linetype = "dashed", xintercept = 1) +
  theme_bw() +
  scale_x_log10()

# -------------------------------------------------------------------------


