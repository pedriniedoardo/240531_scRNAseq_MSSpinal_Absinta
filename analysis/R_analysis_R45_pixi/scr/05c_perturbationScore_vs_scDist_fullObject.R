# AIM ---------------------------------------------------------------------
# compare the two measures of the response to the condition, one point per group:
# 1. the perturbation score (n DEGs / n cells), from 05_DE_sc_fullObject.R (05_DE_sc_summary.tsv)
# 2. the scDist distance, from 05b_scDist_fullObject.R. By default the model adjusted for sex and location (05_scDist_allComparisons_adj.tsv); set tag_model <- "" to use the model with the condition only (05_scDist_allComparisons.tsv)
# run it after both scripts are done (no object is loaded, it only reads the two tables).
# Spearman correlation per comparison and scatter plot of the two measures per level.

# LIBRARIES ---------------------------------------------------------------
library(tidyverse)
library(ggrepel)
library(cowplot)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)

folder_table <- "../../out/table/analysis_R45_pixi/05_DE_sc/"
folder_plot <- "../../out/plot/analysis_R45_pixi/05_DE_sc/"

# suffix of the scDist model: "_adj" = cond + sex + location + (1|sample), "" = cond + (1|sample). Also appended to the names of the output files
tag_model <- "_adj"

# read in the data --------------------------------------------------------
file_summary <- paste0(folder_table, "05_DE_sc_summary.tsv")
file_scdist <- paste0(folder_table, "05_scDist_allComparisons", tag_model, ".tsv")

df_summary <- read_tsv(file_summary, col_types = cols(cluster = col_character(), .default = col_guess()))
df_scdist_all <- read_tsv(file_scdist, col_types = cols(cluster = col_character(), .default = col_guess()))

# perturbation score vs scDist --------------------------------------------
df_corr <- df_summary %>%
  select(level, contrast, comparison, cluster, score) %>%
  inner_join(df_scdist_all %>% select(level, contrast, comparison, cluster, dist = `Dist.`), by = c("level", "contrast", "comparison", "cluster"))

df_corr %>%
  write_tsv(paste0(folder_table, "05_perturbationScore_vs_scDist", tag_model, ".tsv"))

# Spearman correlation per comparison (needs at least 3 groups)
df_cor <- df_corr %>%
  group_by(level, contrast, comparison) %>%
  filter(n() >= 3) %>%
  summarise(rho = cor.test(score, dist, method = "spearman", exact = FALSE)$estimate,
            p = cor.test(score, dist, method = "spearman", exact = FALSE)$p.value,
            n_group = n(),
            .groups = "drop")

df_cor %>%
  write_tsv(paste0(folder_table, "05_perturbationScore_vs_scDist_spearman", tag_model, ".tsv"))

df_corr %>%
  split(f = ~level) %>%
  iwalk(function(df, level){
    df <- df %>%
      left_join(df_cor, by = c("level", "contrast", "comparison")) %>%
      mutate(panel = paste0(contrast, " ", comparison, "\nrho = ", round(rho, 2), ", p = ", signif(p, 2)))

    n_panel <- n_distinct(df$panel)
    n_col <- min(3, n_panel)
    n_row <- ceiling(n_panel / n_col)

    p_corr <- df %>%
      ggplot(aes(x = score, y = dist)) +
      geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "grey50") +
      geom_point(size = 2, color = "grey30") +
      geom_text_repel(aes(label = cluster), size = 3, max.overlaps = 5) +
      facet_wrap(~panel, ncol = n_col, scales = "free") +
      theme_cowplot() +
      theme(strip.background = element_blank()) +
      labs(title = paste0("perturbation score vs scDist - ", level),
           x = paste0("perturbation score (n DEGs / n cells)"), y = "scDist Dist.")

    ggsave(paste0(folder_plot, "05_scatter_perturbationScore_vs_scDist_", level, tag_model, ".pdf"),
           plot = p_corr, width = 7 * n_col, height = 6 * n_row, limitsize = FALSE)
  })

message("\nAll done.")
