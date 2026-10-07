# AIM ---------------------------------------------------------------------
# this is the equivalent of 11_quick_check_EBV_test_rawBarcodes_AbsintaBrainWMCX.R but for the Visium brain dataset processed with spaceranger (host + viral reference)
# no postQC object is available for this dataset, therefore the cross-check against postQC is removed
# the data are treated as a non-spatial (sn-like) dataset: no spatial coordinates are used, and the "called cells" are the barcodes in spaceranger's filtered_feature_bc_matrix

# libraries -----------------------------------------------------------------
library(Seurat)
library(tidyverse)

# specify the version of Seurat Assay -------------------------------------
options(Seurat.object.assay.version = "v5")

# parameters ------------------------------------------------------------------
dir_spaceranger <- "/beegfs/scratch/ric.absinta/ric.absinta/pipeline/snakemake_spaceranger_tests_HPC/visium_brain_out/viral/results/spaceranger/merged"
# list.dirs keeps only the sample folders (the <sample>_finished.log files are skipped)
sample_ids <- list.dirs(dir_spaceranger, recursive = FALSE, full.names = FALSE)
message("samples found: ", length(sample_ids))

out_table_dir <- "../../out/table/analysis_R45_pixi"
# dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
out_prefix <- "13_quick_check_EBV_test_rawBarcodes_VisiumBrain"

# EBV/viral marker panel -- identical (all 86 annotated genes) to scr/test_EBV/07_quick_check_EBV_test_Spinal.R
list_markers_viral <- list(
  EBV_LATENT = c("EBNA-1", "EBNA-2", "EBNA-3A", "EBNA-3B-EBNA-3C", "EBNA-LP", "LMP-1", "RPMS1"),
  EBV_LYTIC_IE = c("BZLF1", "BRLF1"),
  EBV_LYTIC_E = c("BMRF1", "BALF2", "BSLF2-BMLF1", "BALF5", "BHRF1"),
  EBV_LYTIC_L = c("BLLF1", "BcLF1", "BLRF2", "BALF4"),
  EBV_OTHER = c(
    "BNLF2b", "BNLF2a", "BARF1", "BALF1", "BARF0", "BALF3", "A73", "LF2", "LF1", "BILF1", "LF3",
    "BILF2", "BdRF1", "BVRF2", "BVLF1", "BVRF1", "BXRF1", "BXLF2", "BXLF1", "BTRF1", "BcRF1",
    "BDLF3", "BDLF2", "BDLF1", "BGLF2", "BGLF1", "BDLF4", "BDLF3.5", "BGRF1-BDRF1", "BBLF1",
    "BGLF5", "BGLF4", "BGLF3.5", "BGLF3", "BBRF3", "BBLF2-BBLF3", "BBRF2", "BBRF1", "BBLF4",
    "BKRF4", "BKRF3", "BKRF2", "BRRF2", "BRRF1", "BZLF2", "BLLF2", "BLRF1", "BLLF3", "BSRF1",
    "BSLF1", "BMRF2", "BaRF1", "BORF2", "BORF1", "BPLF1", "BOLF1", "BFRF3", "BFRF2", "BFRF1",
    "BFRF1A", "BFLF2", "BFLF1", "BHLF1", "BWRF1", "BCRF1",
    "rna-NC-007605.1:6956..7128", "rna-NC-007605.1:6629..6795", "BNRF1"
  )
)
viral_genes <- unique(unlist(list_markers_viral))
stopifnot(length(viral_genes) == 86)

# processing --------------------------------------------------------------
# per sample: load the raw (full-barcode) matrix, characterize every barcode, then trim to viral genes and keep only nonzero hits
# sample_id <- "MA7678_1"

# each sample returns a list: a one-row summary (always) + its viral hits (possibly empty)
list_res <- map(sample_ids, function(sample_id) {
  raw_dir <- file.path(dir_spaceranger, sample_id, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    message("sample: ", sample_id, "  -- no raw_feature_bc_matrix, skipping")
    return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = NA_integer_, n_genes = NA_integer_,
                                 n_viral_genes_present = NA_integer_, n_viral_hits = 0L, total_viral_UMI = 0, n_called_cell_hits = 0L),
                hits = tibble()))
  }
  message("sample: ", sample_id, "  loading raw matrix...")

  # genes x ALL barcodes, full transcriptome
  raw_mat <- Read10X(data.dir = raw_dir)
  # CreateSeuratObject sanitizes "_" -> "-" in feature names internally (Assay5 constructor)
  # Read10X's raw output still has underscores (BSLF2_BMLF1, EBNA-3B_EBNA-3C, BGRF1_BDRF1, BBLF2_BBLF3, both rna-NC_007605.1:... entries), which wouldn't match list_markers_viral's hyphenated names. Building the object also gets nCount_RNA nFeature_RNA computed for free (min.cells=0/min.features=0 defaults keep every raw barcode).
  # this is also in keeping with the processing in 07_quick_check_EBV_test.R
  obj <- CreateSeuratObject(counts = raw_mat, project = sample_id)
  rm(raw_mat); gc()
  message("  ", ncol(obj), " raw barcodes x ", nrow(obj), " genes")

  # percent.mt still needs its own call
  obj$pct_mt_full <- PercentageFeatureSet(obj, pattern = "^MT-")

  # officially called cells for this sample (spaceranger's own cell calling)
  called_barcodes_file <- file.path(dir_spaceranger, sample_id, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
  called_barcodes <- if (file.exists(called_barcodes_file)) {
    readLines(gzfile(called_barcodes_file))
  } else {
    character(0)
  }

  # now trim to the viral panel and pull out nonzero entries (sparse-native, no densifying)
  n_genes_total <- nrow(obj)
  genes_present <- intersect(viral_genes, rownames(obj))
  # check if there is any gene missing
  setdiff(viral_genes, rownames(obj))

  viral_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[genes_present, , drop = FALSE]
  meta <- obj@meta.data
  rm(obj); gc()

  n_barcodes <- nrow(meta)
  nz <- summary(viral_counts)  # i/j/x for every nonzero entry, sparse-native
  if (nrow(nz) == 0) {
    message("  no nonzero viral-panel entries in this sample's raw matrix")
    return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = n_barcodes, n_genes = n_genes_total,
                                 n_viral_genes_present = length(genes_present), n_viral_hits = 0L, total_viral_UMI = 0, n_called_cell_hits = 0L),
                hits = tibble()))
  }

  out <- tibble(
    sample_id = sample_id,
    gene = rownames(viral_counts)[nz$i],
    barcode = colnames(viral_counts)[nz$j],
    count = nz$x,
    nCount_RNA_full = meta$nCount_RNA[nz$j],
    nFeature_RNA_full = meta$nFeature_RNA[nz$j],
    pct_mt_full = meta$pct_mt_full[nz$j],
    is_called_tissue = colnames(viral_counts)[nz$j] %in% called_barcodes
  )
  message("  ", nrow(out), " nonzero viral-panel (gene, barcode) hit(s), total UMI: ", sum(out$count))
  return(list(summary = tibble(sample_id = sample_id, n_raw_barcodes = n_barcodes, n_genes = n_genes_total,
                               n_viral_genes_present = length(genes_present), n_viral_hits = nrow(out),
                               total_viral_UMI = sum(out$count), n_called_cell_hits = sum(out$is_called_tissue)),
              hits = out))
})

# per-sample summary (same numbers as printed in the log), written even when there are no hits
df_summary <- map_dfr(list_res, "summary")
df_hits <- map_dfr(list_res, "hits")
df_summary
write_tsv(df_summary, file.path(out_table_dir, paste0(out_prefix, "_perSampleSummary.tsv")))

if (nrow(df_hits) == 0) {
  message("no barcode in the entire raw dataset (all droplets, all samples) has a nonzero viral-panel count")
  df_hits_full <- df_hits
} else {
  df_hits_full <- df_hits %>% arrange(desc(count))
  message(nrow(df_hits_full), " viral-positive (gene, barcode) hit(s) found across all raw barcodes, ",
          "all samples -- ", sum(df_hits_full$is_called_tissue), " in spaceranger's called-cell list")
}
df_hits_full

write_tsv(df_hits_full, file.path(out_table_dir, paste0(out_prefix, "_viralPositiveBarcodes_raw.tsv")))
