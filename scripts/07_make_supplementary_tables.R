############################################################
# Zika scRNA-seq project – STEP 7
# Supplementary Tables 2 and 3
#
# Purpose:
#   Generate reproducibility tables documented in the analysis:
#     - Supplementary Table 2: sample mapping and dataset-use labels
#     - Supplementary Table 3: cell-type marker genes and curated ISG panel
#
# Analysis notes:
#   - GSE238140 is the primary fetal neural scRNA-seq dataset.
#   - GSE230571 is used only as a contextual Vero/moDC comparison.
#   - GSE97919 is included as contextual bulk RNA-seq information.
#   - The Vero/moDC and bulk datasets are not used as validation of fetal
#     neural lineage-specific responses.
#   - FSS13025 is annotated as a Cambodian 2010 Asian-lineage isolate.
#
# Inputs:
#   Optional:
#     results/dc_step2_norm_qc.rds
#
# Outputs:
#   supplementary_tables/Supplementary_Table_2_sample_mapping.csv
#   supplementary_tables/Supplementary_Table_3_gene_lists.csv
#   supplementary_tables/Supplementary_Tables_2_3.xlsx  if writexl is installed
#   session_info/STEP7_sessionInfo.txt
############################################################


## ==========================================================
## 0. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP7] Project root: ", project_root)

required_packages <- c("dplyr", "tibble", "readr")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP7] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 7."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})


## ==========================================================
## 1. Directories
## ==========================================================

results_dir <- file.path(project_root, "results")
supp_dir <- file.path(project_root, "supplementary_tables")
session_dir <- file.path(project_root, "session_info")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 2. Supplementary Table 2 – sample mapping
## ==========================================================

## -------------------------------
## 2A. GSE238140 fetal neural scRNA-seq
## -------------------------------

supp_table2_brain <- tibble(
  dataset = "Fetal neural scRNA-seq",
  geo_accession = "GSE238140",
  geo_sample_id = c("GSM7659279", "GSM7659280", "GSM7659281", "GSM7659282"),
  cell_system = "Human fetal brain-derived neural tissue",
  condition = c("Mock", "ZIKV-BR", "ZIKV-FSS/FSS13025", "IFNβ"),
  analysis_label = c("Mock", "ZIKV_BR", "ZIKV_FSS", "IFNb"),
  study_use = c(
    "Primary fetal neural dataset",
    "Primary fetal neural dataset",
    "Primary fetal neural dataset",
    "Primary fetal neural dataset"
  ),
  interpretation_scope = c(
    "Reference/control condition",
    "Exploratory viral exposure condition",
    "Exploratory viral exposure condition",
    "Cytokine stimulation reference condition"
  ),
  notes = c(
    "Mock control",
    "Brazilian ZIKV exposure condition, referred to as ZIKV-BR",
    "FSS13025 isolate; Cambodian 2010 Asian-lineage isolate, referred to as ZIKV-FSS",
    "Recombinant interferon-beta stimulation, referred to as IFNβ"
  )
)


## -------------------------------
## 2B. GSE230571 Vero/moDC scRNA-seq
## -------------------------------

dc_rds_path <- file.path(results_dir, "dc_step2_norm_qc.rds")

if (file.exists(dc_rds_path)) {

  message("[STEP7] Reading dc object for sample IDs: ", dc_rds_path)

  dc <- readRDS(dc_rds_path)

  if (!inherits(dc, "Seurat")) {
    warning("[STEP7] dc_step2_norm_qc.rds exists but is not a Seurat object. Creating summary-level GSE230571 rows.")

    supp_table2_dc <- tibble(
      dataset = "Vero/moDC scRNA-seq",
      geo_accession = "GSE230571",
      geo_sample_id = c("ZIKV-Vero-cells", "p22086 panel"),
      cell_system = c("Vero cells", "Primary human monocyte-derived dendritic cells"),
      condition = c("ZIKV-exposed", "moDC condition panel"),
      analysis_label = c("ZIKV_Vero_cells", "moDC_panel"),
      study_use = "Contextual comparison only",
      interpretation_scope = "Broad IFN/ISG-associated contextual comparison; not evidence for fetal neural lineage-specific responses",
      notes = "Used only for contextual cross-system comparison"
    )

  } else {

    dc_meta <- dc@meta.data %>%
      as.data.frame() %>%
      tibble::rownames_to_column("cell_barcode")

    if (!"sample_id" %in% colnames(dc_meta)) {
      warning("[STEP7] dc object does not contain sample_id metadata. Creating summary-level GSE230571 rows.")

      supp_table2_dc <- tibble(
        dataset = "Vero/moDC scRNA-seq",
        geo_accession = "GSE230571",
        geo_sample_id = c("ZIKV-Vero-cells", "p22086 panel"),
        cell_system = c("Vero cells", "Primary human monocyte-derived dendritic cells"),
        condition = c("ZIKV-exposed", "moDC condition panel"),
        analysis_label = c("ZIKV_Vero_cells", "moDC_panel"),
        study_use = "Contextual comparison only",
        interpretation_scope = "Broad IFN/ISG-associated contextual comparison; not evidence for fetal neural lineage-specific responses",
        notes = "Used only for contextual cross-system comparison"
      )

    } else {

      supp_table2_dc <- dc_meta %>%
        distinct(sample_id) %>%
        mutate(
          dataset = "Vero/moDC scRNA-seq",
          geo_accession = "GSE230571",
          geo_sample_id = sample_id,
          cell_system = case_when(
            grepl("ZIKV-Vero-cells", sample_id) ~ "Vero cells",
            grepl("p22086", sample_id) ~ "Primary human monocyte-derived dendritic cells",
            TRUE ~ "Unknown"
          ),
          condition = case_when(
            grepl("ZIKV-Vero-cells", sample_id) ~ "ZIKV-exposed",
            grepl("p22086", sample_id) ~ "moDC condition panel",
            TRUE ~ "Unknown"
          ),
          analysis_label = case_when(
            grepl("ZIKV-Vero-cells", sample_id) ~ "ZIKV_Vero_cells",
            grepl("p22086", sample_id) ~ "moDC_panel",
            TRUE ~ "Unknown"
          ),
          study_use = "Contextual comparison only",
          interpretation_scope = "Broad IFN/ISG-associated contextual comparison; not evidence for fetal neural lineage-specific responses",
          notes = "Used only for contextual cross-system comparison"
        ) %>%
        select(
          dataset,
          geo_accession,
          geo_sample_id,
          cell_system,
          condition,
          analysis_label,
          study_use,
          interpretation_scope,
          notes
        )
    }
  }

} else {

  message("[STEP7] dc_step2_norm_qc.rds not found. Creating summary-level GSE230571 rows.")

  supp_table2_dc <- tibble(
    dataset = "Vero/moDC scRNA-seq",
    geo_accession = "GSE230571",
    geo_sample_id = c("ZIKV-Vero-cells", "p22086 panel"),
    cell_system = c("Vero cells", "Primary human monocyte-derived dendritic cells"),
    condition = c("ZIKV-exposed", "moDC condition panel"),
    analysis_label = c("ZIKV_Vero_cells", "moDC_panel"),
    study_use = "Contextual comparison only",
    interpretation_scope = "Broad IFN/ISG-associated contextual comparison; not evidence for fetal neural lineage-specific responses",
    notes = "Used only for contextual cross-system comparison"
  )
}


## -------------------------------
## 2C. GSE97919 bulk RNA-seq
## -------------------------------

supp_table2_bulk <- tibble(
  dataset = "Bulk RNA-seq",
  geo_accession = "GSE97919",
  geo_sample_id = "GSM-level samples available in GEO metadata",
  cell_system = "Human cerebral organoid / neural bulk RNA-seq samples",
  condition = "ZIKV-exposed/control",
  analysis_label = "bulk_ZIKV_vs_control",
  study_use = "Contextual bulk RNA-seq analysis only",
  interpretation_scope = "Contextual support for broad interferon-associated transcriptional activation; not lineage-resolved evidence",
  notes = "Used only as contextual bulk transcriptomic comparison because experimental system and modality differ from fetal neural scRNA-seq"
)


## -------------------------------
## 2D. Combine and save
## -------------------------------

supp_table2 <- bind_rows(
  supp_table2_brain,
  supp_table2_dc,
  supp_table2_bulk
)

supp_table2_path <- file.path(supp_dir, "Supplementary_Table_2_sample_mapping.csv")
readr::write_csv(supp_table2, supp_table2_path)
message("[STEP7] Saved: ", supp_table2_path)


## ==========================================================
## 3. Supplementary Table 3 – marker genes and ISG panel
## ==========================================================

supp_table3 <- tibble(
  gene_list_category = c(
    rep("Cell-type marker", 7),
    "Contextual progenitor subtype marker",
    "Contextual progenitor subtype marker",
    "Contextual progenitor subtype marker",
    "ISG signature"
  ),
  cell_type_or_signature = c(
    "Progenitor-enriched radial glia-like",
    "Cycling NPC",
    "Excitatory neuron",
    "Inhibitory neuron",
    "Astrocyte",
    "Microglia",
    "OPC",
    "Apical radial glia context",
    "Outer/basal radial glia context",
    "Intermediate progenitor context",
    "Canonical ISG panel"
  ),
  genes = c(
    "SOX2, PAX6, NES",
    "MKI67, TOP2A, HMGB2",
    "SLC17A7, NEUROD6",
    "GAD1, GAD2, DLX1",
    "GFAP, AQP4",
    "C1QA, C1QB, TYROBP",
    "PDGFRA, OLIG1, OLIG2",
    "NOTCH1, HES1, HES5",
    "HOPX, PTPRZ1, FAM107A",
    "EOMES",
    "IFITM1, MX1, OAS1, OAS2, OAS3, IFIT1, IFI6, ISG15, STAT1, IRF7, RSAD2"
  ),
  use_in_this_study = c(
    "Cell-type annotation",
    "Cell-type annotation",
    "Cell-type annotation",
    "Cell-type annotation",
    "Cell-type annotation",
    "Cell-type annotation",
    "Cell-type annotation",
    "Contextual marker set used to support cautious annotation discussion",
    "Contextual marker set used to support cautious annotation discussion",
    "Contextual marker set used to support cautious annotation discussion",
    "ISG heatmaps, ISG summaries, AddModuleScore, and ZIKV–IFNβ correlation sensitivity analyses"
  ),
  notes = c(
    "Interpreted cautiously as progenitor-enriched radial glia-like because SOX2, PAX6, and NES are shared across multiple progenitor states",
    "Cycling progenitor-associated markers",
    "Excitatory neuron-associated markers",
    "Inhibitory neuron-associated markers",
    "Astrocyte-associated markers",
    "Microglia-associated markers",
    "Oligodendrocyte precursor-associated markers",
    "Not used to claim a resolved apical radial glia subtype; included to address marker-resolution limitations",
    "Not used to claim a resolved outer/basal radial glia subtype; included to address marker-resolution limitations",
    "Not used to claim a resolved intermediate progenitor subtype; included to address marker-resolution limitations",
    "Curated canonical interferon-stimulated gene panel"
  )
)

supp_table3_path <- file.path(supp_dir, "Supplementary_Table_3_gene_lists.csv")
readr::write_csv(supp_table3, supp_table3_path)
message("[STEP7] Saved: ", supp_table3_path)


## ==========================================================
## 4. Optional Excel output
## ==========================================================

if (requireNamespace("writexl", quietly = TRUE)) {

  xlsx_path <- file.path(supp_dir, "Supplementary_Tables_2_3.xlsx")

  writexl::write_xlsx(
    list(
      "Supp_Table_2_sample_mapping" = supp_table2,
      "Supp_Table_3_gene_lists" = supp_table3
    ),
    path = xlsx_path
  )

  message("[STEP7] Saved Excel file: ", xlsx_path)

} else {
  message("[STEP7] Package 'writexl' is not installed. CSV files were saved only.")
  message("[STEP7] Optional Excel output can be enabled with: install.packages('writexl')")
}


## ==========================================================
## 5. Save reproducibility information
## ==========================================================

interpretation_note <- c(
  "STEP 7 supplementary table interpretation note",
  "",
  "Supplementary Table 2 provides dataset/sample mapping and analysis labels.",
  "Supplementary Table 3 provides cell-type marker genes and the curated canonical ISG panel.",
  "The fetal neural dataset GSE238140 is the primary dataset.",
  "GSE230571 and GSE97919 are contextual comparison datasets only.",
  "They are not used as independent evidence for fetal neural lineage-specific responses.",
  "FSS13025 is annotated as a Cambodian 2010 Asian-lineage isolate and referred to as ZIKV-FSS."
)

writeLines(
  interpretation_note,
  file.path(results_dir, "step7_supplementary_tables_interpretation_note.txt")
)

sink(file.path(session_dir, "STEP7_sessionInfo.txt"))
sessionInfo()
sink()

message("\n[STEP7] Supplementary Table 2 preview:")
print(supp_table2)

message("\n[STEP7] Supplementary Table 3 preview:")
print(supp_table3)

message("\n[STEP7] Done.")

############################################################
# End of STEP 7
############################################################
