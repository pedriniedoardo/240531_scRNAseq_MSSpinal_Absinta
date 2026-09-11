# AIM ---------------------------------------------------------------------
# EDA on the final subcluster-filtered full object (output of R44's 29_apply_filter_annotate_fullObject.R):
# 1) DC vs microglia vs BAM markers in IMMUNE clusters
# 2) ciliated genes in ASTRO clusters
# 3) lymphocyte subtype markers in LYM clusters
# 4) IFN-signature validity across OPC + OLIGO + IMMUNE clusters (OPC6/OLIGO9/IMMUNE10)
# 5) OPALIN/oligodendrocyte maturation panel (Jakel 2019, Castelo-Branco lab) in the full object,
#    OLIGO shown at subcluster resolution against all other cell types pooled as reference

# LIBRARIES -----------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)

# PARAMETERS ------------------------------------------------------------------
options(stringsAsFactors = FALSE)
options(Seurat.object.assay.version = "v5")

# read in the data --------------------------------------------------------
sobj <- readRDS("../../out/object/analysis_R44/29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.rds")
DefaultAssay(sobj) <- "RNA"

# save the full metadata with the current draft anntoation from Martina
sobj@meta.data %>% group_by(cell_id_subcluster,cell_id_subcluster2) %>% summarise(n = n()) %>%
  write_tsv("../../out/table/analysis_R45_pixi/01_meta_29_sobj_integrated_cleanup_manualAnnotation_subclusterFiltered.tsv")

out_plot_dir <- "../../out/plot/analysis_R45_pixi"

# helper A: DotPlot a named gene-panel list for the target populations only, at fine resolution (group_col) -- the original focal-only view --------
plot_marker_dotplot_focal <- function(obj,
                                      gene_list,
                                      populations,
                                      group_col = "cell_id_subcluster2",
                                      pop_col = "cell_id",
                                      title, filename, width, height) {
  genes_present <- lapply(gene_list, function(g) intersect(g, rownames(obj)))
  missing <- setdiff(unlist(gene_list), unlist(genes_present))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))

  sub <- subset(obj, cells = colnames(obj)[obj[[pop_col, drop = TRUE]] %in% populations])

  p <- DotPlot(sub, features = genes_present, group.by = group_col,
               dot.scale = 8, cluster.idents = TRUE) +
    RotatedAxis() +
    labs(title = title) +
    theme(strip.text = element_text(angle = 90))

  ggsave(plot = p, filename = filename, width = width, height = height)
  return(p)
}

# helper B: DotPlot a named gene-panel list for the target populations at fine resolution (group_col), plus every non-target cell_id pooled as "Ref_<cell_id>" below a dashed separator -- so the panel is judged against general background too, not just within the target subclusters (same focal-vs-reference idiom as task 5) --------
plot_marker_dotplot_vs_ref <- function(obj,
                                       gene_list,
                                       populations,
                                       group_col = "cell_id_subcluster2",
                                       pop_col = "cell_id",
                                       title, filename, width, height) {
  genes_present <- lapply(gene_list, function(g) intersect(g, rownames(obj)))
  missing <- setdiff(unlist(gene_list), unlist(genes_present))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))

  is_focal <- obj[[pop_col, drop = TRUE]] %in% populations
  obj$plot_id <- if_else(is_focal,
                          as.character(obj[[group_col, drop = TRUE]]),
                          paste0("Ref_", obj[[pop_col, drop = TRUE]]))

  focal_ids <- sort(unique(obj$plot_id[is_focal]))
  ref_ids   <- sort(unique(obj$plot_id[!is_focal]))
  id_levels <- c(focal_ids, ref_ids)
  obj$plot_id <- factor(obj$plot_id, levels = id_levels)
  Idents(obj) <- "plot_id"

  sep_y <- length(focal_ids) + 0.5

  p <- DotPlot(obj, features = genes_present, dot.scale = 8, cluster.idents = FALSE) +
    RotatedAxis() +
    geom_hline(yintercept = sep_y, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(title = title,
         caption = paste0(paste(populations, collapse = "/"),
                           " subclusters below the line | Ref_* = pooled other cell_id above")) +
    theme(strip.text = element_text(angle = 90),
          plot.caption = element_text(hjust = 0, colour = "grey40", size = 7))

  ggsave(plot = p, filename = filename, width = width, height = length(id_levels) * 0.35 + 3)
  return(p)
}

# 1. DC vs microglia vs BAM -- IMMUNE clusters -------------------------------
dc_panel <- list(
  cDC1                       = c("CLEC9A", "XCR1", "BATF3"),
  `cDC2/conventional_DC`     = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  `DC_activation/migratory`  = c("CD83", "COTL1", "CCR7"),
  homeostatic_microglia      = c("P2RY12", "TMEM119", "SALL1", "GPR34", "CX3CR1"),
  `BAM/perivascular_mac`     = c("MRC1", "LYVE1", "SELENOP", "CD163", "STAB1")
)

plot_marker_dotplot_focal(
  sobj, dc_panel, populations = "IMMUNE",
  title = "DC vs microglia vs BAM markers -- IMMUNE subclusters (focal only)",
  filename = file.path(out_plot_dir, "01_dotplot_DC_markers_IMMUNE_focalOnly.pdf"),
  width = 12, height = 6
)

plot_marker_dotplot_vs_ref(
  sobj, dc_panel, populations = "IMMUNE",
  title = "DC vs microglia vs BAM markers -- IMMUNE subclusters vs reference",
  filename = file.path(out_plot_dir, "01_dotplot_DC_markers_IMMUNE_vsRef.pdf"),
  width = 12, height = 6
)

# 2. Ciliated astrocyte genes -- ASTRO clusters ------------------------------
ciliated_panel <- list(
  ciliated_astrocyte = c("FRMPD2", "ADGB", "SPAG17", "FOXJ1", "DNAH11", "DNAH6", "CFAP54", "CFAP299")
)

plot_marker_dotplot_focal(
  sobj, ciliated_panel, populations = "ASTRO",
  title = "Ciliated astrocyte genes -- ASTRO subclusters (focal only)",
  filename = file.path(out_plot_dir, "01_dotplot_ciliatedGenes_ASTRO_focalOnly.pdf"),
  width = 8, height = 5
)

plot_marker_dotplot_vs_ref(
  sobj, ciliated_panel, populations = "ASTRO",
  title = "Ciliated astrocyte genes -- ASTRO subclusters vs reference",
  filename = file.path(out_plot_dir, "01_dotplot_ciliatedGenes_ASTRO_vsRef.pdf"),
  width = 8, height = 5
)

# 3. Lymphocyte subtype markers -- LYM clusters ------------------------------
# note: cell_id_subcluster(2) uses the pooled "LYM|..." prefix for these cells even though cell_id itself is split into "B CELLS"/"T CELLS"
lym_panel <- list(
  `MAIT/mucosal`  = c("TRAV1-2", "KLRB1", "SLC4A10", "DPP4"),
  pan_T           = c("TRBC1", "TRAC", "IL7R", "LTB"),
  `CD4/CD8`       = c("CD4", "CD8A", "CD8B"),
  trafficking     = c("CCR6", "CXCR6"),
  `NK/cytotoxic`  = c("KLRD1", "NCR1", "FCGR3A", "TYROBP")
)

plot_marker_dotplot_focal(
  sobj, lym_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters (focal only)",
  filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_LYM_focalOnly.pdf"),
  width = 10, height = 6
)

plot_marker_dotplot_vs_ref(
  sobj, lym_panel, populations = c("B CELLS","T CELLS"),
  title = "Lymphocyte subtype markers -- LYM subclusters vs reference",
  filename = file.path(out_plot_dir, "01_dotplot_lymphocyte_markers_LYM_vsRef.pdf"),
  width = 10, height = 6
)

# 4. IFN signature -- OPC + OLIGO + IMMUNE together, one DotPlot -------------
# cell_id_subcluster2 already carries the population prefix (e.g. "OPC|OPC-IFN|6"). Rather than just ordering clusters back-to-back, pull DotPlot's own long-format data (id/features.plot/pct.exp/avg.exp.scaled) and rebuild the plot with facet_grid(cell_type ~ ., scales = "free", space = "free") so OPC/OLIGO/IMMUNE render as visually separated row-blocks -- the candidate IFN clusters (OPC-IFN/Oligo-IFN/MIMS-IFN) compare directly against their own population's other clusters, with a clear visual break between populations.
ifn_panel <- list(
  IFN_signature = c(
    "ISG15", "MX1", "MX2", "IFIT1", "IFIT2", "IFIT3", "IFITM1", "IFITM3", "RSAD2",
    "OAS1", "OAS2", "OAS3", "OASL", "STAT1", "STAT2", "IRF1", "IRF7", "JAK1", "JAK2",
    "GBP1", "GBP5", "EPSTI1", "PARP14", "XAF1", "APOL1", "APOL2", "APOL3", "NLRC5", "B2M"
  )
)

ifn_genes_present <- intersect(ifn_panel$IFN_signature, rownames(sobj))
missing_ifn <- setdiff(ifn_panel$IFN_signature, ifn_genes_present)
if (length(missing_ifn) > 0) message("Missing from object, dropped: ", paste(missing_ifn, collapse = ", "))

sobj_ifn <- subset(sobj, cells = colnames(sobj)[sobj$cell_id %in% c("OPC", "OLIGO", "IMMUNE")])
Idents(sobj_ifn) <- "cell_id_subcluster2"

dp_ifn <- DotPlot(sobj_ifn, features = ifn_genes_present)

id_to_pop_ifn <- sobj_ifn@meta.data %>%
  distinct(cell_id_subcluster2, cell_id) %>%
  mutate(cluster_num = suppressWarnings(as.integer(str_extract(cell_id_subcluster2, "[0-9]+$")))) %>%
  rename(id = cell_id_subcluster2, cell_type = cell_id) %>%
  mutate(cell_type = factor(cell_type, levels = c("OPC", "OLIGO", "IMMUNE"))) %>%
  arrange(cell_type, cluster_num)

df_plot_ifn <- dp_ifn$data %>%
  left_join(id_to_pop_ifn, by = "id") %>%
  mutate(id = factor(id, levels = id_to_pop_ifn$id))

p_ifn <- ggplot(df_plot_ifn, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, col = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  facet_grid(cell_type ~ ., scales = "free", space = "free") +
  theme_cowplot() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(hjust = 1, angle = 90),
        strip.text.y = element_text(angle = 0)) +
  scale_color_gradient(low = "lightgrey", high = "blue") +
  labs(title = "IFN signature across OPC / OLIGO / IMMUNE subclusters",
       x = "IFN signature genes", y = "cell_id_subcluster2")

ggsave(plot = p_ifn, filename = file.path(out_plot_dir, "01_dotplot_IFNsignature_OPC_OLIGO_IMMUNE.pdf"),
       width = 16, height = 8)

# 5. OPALIN / oligodendrocyte maturation (Jakel 2019, Castelo-Branco lab) ----
# focal-vs-reference calibrated DotPlot: OLIGO cells shown at subcluster resolution, everything else pooled as "Ref_<cell_id>" -- run twice, once per fine-label column. Literature check: OPALIN marks the mature "Oligo6" state in Jakel et al. 2019 Nature (human MS white matter oligodendrocyte heterogeneity, Castelo-Branco lab); KLK6 marks "Oligo5"; RASGRF1 marks "Oligo1". Full Oligo1-6 subtype assignment would need the paper's supplementary marker table for high confidence -- this panel is a maturation/heterogeneity screen, not a claim of exact subtype identity.
opalin_panel <- list(
  oligo_maturation = c(
    "PDGFRA", "CSPG4", "GPR17", "BCAS1", "OLIG2", "SOX10", "MYRF", "PLP1", "MBP",
    "MOBP", "MAG", "MOG", "OPALIN", "TPPP", "CNP", "GFAP", "DPP10",
    "KLK6", "RASGRF1", "RBFOX1", "ANLN"
  )
)

plot_focal_vs_ref_dotplot <- function(obj,
                                      focal_pop,
                                      focal_label_col,
                                      panel,
                                      out_file,
                                      title) {
  obj$plot_id <- if_else(obj$cell_id == focal_pop,
                          as.character(obj[[focal_label_col, drop = TRUE]]),
                          paste0("Ref_", obj$cell_id))

  focal_ids <- sort(unique(obj$plot_id[obj$cell_id == focal_pop]))
  ref_ids   <- sort(unique(obj$plot_id[obj$cell_id != focal_pop]))
  id_levels <- c(focal_ids, ref_ids)
  obj$plot_id <- factor(obj$plot_id, levels = id_levels)
  Idents(obj) <- "plot_id"

  panel_filt <- lapply(panel, function(gs) intersect(gs, rownames(obj)))
  missing <- setdiff(unlist(panel), unlist(panel_filt))
  if (length(missing) > 0) message("Missing from object, dropped: ", paste(missing, collapse = ", "))
  panel_filt <- panel_filt[lengths(panel_filt) > 0]

  sep_y <- length(focal_ids) + 0.5
  p <- DotPlot(obj, features = panel_filt, dot.scale = 6, cluster.idents = FALSE) +
    RotatedAxis() +
    geom_hline(yintercept = sep_y, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    labs(title = title,
         caption = paste0(focal_pop, " subclusters below the line | Ref_* = pooled other cell_id above")) +
    theme(strip.text = element_text(angle = 90),
          plot.caption = element_text(hjust = 0, colour = "grey40", size = 7))

  ggsave(filename = out_file, plot = p, width = 14, height = length(id_levels) * 0.4 + 3.5)
  return(p)
}

plot_focal_vs_ref_dotplot(
  sobj, "OLIGO", "cell_id_subcluster", opalin_panel,
  out_file = file.path(out_plot_dir, "01_dotplot_OPALIN_oligoMaturation_calibrated_subclusterLabel.pdf"),
  title = "OLIGO subclusters vs reference cell types -- OPALIN/maturation panel"
)

plot_focal_vs_ref_dotplot(
  sobj, "OLIGO", "cell_id_subcluster2", opalin_panel,
  out_file = file.path(out_plot_dir, "01_dotplot_OPALIN_oligoMaturation_calibrated_subclusterLabel2.pdf"),
  title = "OLIGO subclusters (with res0.4 id) vs reference cell types -- OPALIN/maturation panel"
)
