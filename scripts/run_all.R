############################################################
# Run all analysis steps
# Execute from the repository root after placing input data in data/.
############################################################

scripts <- file.path("scripts", c(
  "01_prepare_10x_data.R",
  "02_qc_normalization_umap_clustering.R",
  "03_celltype_annotation.R",
  "04_exploratory_celltype_DE.R",
  "05_progenitor_volcano_ISG_heatmap.R",
  "06_ZIKV_IFNb_correlation_sensitivity.R",
  "07_make_supplementary_tables.R",
  "08_transcriptional_response_summary.R",
  "09_bulk_RNAseq_contextual_GSE97919.R"
))

for (script in scripts) {
  message("
==============================")
  message("Running: ", script)
  message("==============================")
  source(script, local = FALSE)
}
