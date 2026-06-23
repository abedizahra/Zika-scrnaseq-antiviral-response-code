############################################################
# Zika scRNA-seq project – STEP 4
# Exploratory cell-type-level differential expression
#
# Purpose:
#   Run Seurat::FindMarkers within each annotated fetal neural
#   cell type for:
#     - ZIKV_BR  vs Mock
#     - ZIKV_FSS vs Mock
#     - IFNb     vs Mock
#
# Analysis notes:
#   - This dataset has one biological sample per fetal neural condition.
#   - Cells are NOT biological replicates.
#   - Cell-level Wilcoxon tests are used only as exploratory ranking and
#     descriptive summaries.
#   - P values/FDR values should not be interpreted as formal
#     condition-level inference.
#   - DEG counts are descriptive summaries of observed transcriptional
#     patterns, not statistically powered comparisons among conditions.
#
# Inputs:
#   results/brain_step3_annotated.rds
#
# Outputs:
#   results/DEG_sc/*.csv
#   results/step4_analysis_parameters.txt
#   results/step4_interpretation_note.txt
#   figures/step4/main/*.tiff
#   figures/step4/supplementary/*.tiff
#   session_info/STEP4_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

analysis_parameters <- list(
  seed = 1234,
  test_use = "wilcox",
  min_cells_per_group = 20,
  logfc_threshold = 0.25,
  min_pct = 0.10,
  adjusted_p_cutoff_for_summary = 0.05,
  logfc_cutoff_for_summary = 0.25
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP4] Project root: ", project_root)

required_packages <- c(
  "Seurat",
  "dplyr",
  "tibble",
  "tidyr",
  "purrr",
  "ggplot2",
  "readr",
  "scales",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP4] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 4."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(patchwork)
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)


## ==========================================================
## 2. Directories
## ==========================================================

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
session_dir <- file.path(project_root, "session_info")

deg_dir <- file.path(results_dir, "DEG_sc")
fig_step4_main <- file.path(figures_dir, "step4", "main")
fig_step4_supp <- file.path(figures_dir, "step4", "supplementary")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step4_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step4_supp, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Plot theme and save helper
## ==========================================================

theme_journal <- theme_classic(base_size = 13) +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x  = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.55, "cm"),
    axis.line = element_line(linewidth = 0.8, color = "black")
  )

save_tiff_step4 <- function(plot_obj, out_path, width = 6, height = 5, dpi = 600) {
  if (has_ragg) {
    ggplot2::ggsave(
      filename = out_path,
      plot = plot_obj,
      device = ragg::agg_tiff,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw"
    )
  } else {
    ggplot2::ggsave(
      filename = out_path,
      plot = plot_obj,
      device = "tiff",
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      compression = "lzw"
    )
  }
  message("[STEP4] Saved: ", out_path)
}


## ==========================================================
## 4. Load annotated STEP 3 object
## ==========================================================

brain_annotated_path <- file.path(results_dir, "brain_step3_annotated.rds")

if (!file.exists(brain_annotated_path)) {
  stop(
    "[STEP4] Cannot find: ", brain_annotated_path,
    "\nRun STEP 3 first and make sure it saves results/brain_step3_annotated.rds."
  )
}

brain <- readRDS(brain_annotated_path)

if (!inherits(brain, "Seurat")) {
  stop("[STEP4] The loaded object is not a Seurat object.")
}

required_metadata <- c("cell_type", "condition")
missing_metadata <- setdiff(required_metadata, colnames(brain@meta.data))

if (length(missing_metadata) > 0) {
  stop(
    "[STEP4] Missing required metadata column(s): ",
    paste(missing_metadata, collapse = ", "),
    "\nRun the corrected STEP 3 annotation script before STEP 4."
  )
}

message("[STEP4] Brain cells: ", ncol(brain))
message("[STEP4] Cell types:")
print(table(brain$cell_type))
message("[STEP4] Conditions:")
print(table(brain$condition))


## ==========================================================
## 5. Choose assay robustly
## ==========================================================

# Note:
# Do NOT use names(Assays(brain)).
# In many Seurat objects, names(Assays(brain)) returns NULL even when
# assays exist. Use Assays(brain) directly.

assays_available <- tryCatch(
  as.character(Seurat::Assays(brain)),
  error = function(e) character(0)
)

message("[STEP4] Assays available: ", paste(assays_available, collapse = ", "))

if ("SCT" %in% assays_available) {
  assay_to_use <- "SCT"
} else if ("RNA" %in% assays_available) {
  assay_to_use <- "RNA"
} else {
  assay_to_use <- DefaultAssay(brain)
  warning(
    "[STEP4] Neither SCT nor RNA was found by Assays(). ",
    "Using DefaultAssay(brain) = ", assay_to_use
  )
}

DefaultAssay(brain) <- assay_to_use
message("[STEP4] Assay selected for exploratory DE: ", assay_to_use)

# Prepare assay for FindMarkers.
# If SCT exists, use PrepSCTFindMarkers when possible.
# If only RNA exists, ensure log-normalized data are present.
if (assay_to_use == "SCT") {
  brain <- tryCatch(
    {
      message("[STEP4] Running PrepSCTFindMarkers() for SCT assay...")
      PrepSCTFindMarkers(brain)
    },
    error = function(e) {
      warning(
        "[STEP4] PrepSCTFindMarkers failed; continuing with SCT assay as-is. Message: ",
        conditionMessage(e)
      )
      brain
    }
  )
} else if (assay_to_use == "RNA") {
  message("[STEP4] Using RNA assay. Checking/creating normalized data layer if needed...")
  has_data_layer <- TRUE

  invisible(tryCatch(
    {
      tmp <- SeuratObject::LayerData(brain, assay = "RNA", layer = "data")
      if (ncol(tmp) == 0 || nrow(tmp) == 0) has_data_layer <<- FALSE
    },
    error = function(e) {
      has_data_layer <<- FALSE
    }
  ))

  if (!has_data_layer) {
    message("[STEP4] RNA data layer not found. Running NormalizeData()...")
    brain <- NormalizeData(
      brain,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
}


## ==========================================================
## 6. Helper: exploratory DE per cell type
## ==========================================================

run_scDE_celltype <- function(brain_obj,
                              cell_type_name,
                              condition_A = "Mock",
                              condition_B = "ZIKV_BR",
                              assay = assay_to_use,
                              min_cells_per_group = analysis_parameters$min_cells_per_group,
                              test_use = analysis_parameters$test_use,
                              logfc_threshold = analysis_parameters$logfc_threshold,
                              min_pct = analysis_parameters$min_pct) {

  message(
    "\n[STEP4] Exploratory DE: ",
    cell_type_name, " | ", condition_B, " vs ", condition_A
  )

  cells_ct <- WhichCells(brain_obj, expression = cell_type == cell_type_name)

  if (length(cells_ct) == 0) {
    warning("[STEP4] No cells found for cell type: ", cell_type_name)
    return(NULL)
  }

  cond_vec <- brain_obj$condition[cells_ct]
  tab_cond <- table(cond_vec)
  print(tab_cond)

  if (!all(c(condition_A, condition_B) %in% names(tab_cond))) {
    warning(
      "[STEP4] Missing one or both conditions for ",
      cell_type_name, ": ", condition_B, " vs ", condition_A
    )
    return(NULL)
  }

  if (any(tab_cond[c(condition_A, condition_B)] < min_cells_per_group)) {
    warning(
      "[STEP4] Too few cells for ", cell_type_name,
      " in ", condition_B, " vs ", condition_A,
      ". Need at least ", min_cells_per_group, " cells per condition."
    )
    return(NULL)
  }

  Idents(brain_obj) <- brain_obj$condition

  markers <- tryCatch(
    {
      FindMarkers(
        brain_obj,
        ident.1 = condition_B,
        ident.2 = condition_A,
        assay = assay,
        test.use = test_use,
        logfc.threshold = logfc_threshold,
        min.pct = min_pct,
        cells = cells_ct
      )
    },
    error = function(e) {
      warning(
        "[STEP4] FindMarkers failed for ",
        cell_type_name, " | ", condition_B, " vs ", condition_A,
        ". Message: ", conditionMessage(e)
      )
      return(NULL)
    }
  )

  if (is.null(markers) || nrow(markers) == 0) {
    warning(
      "[STEP4] No markers returned for ",
      cell_type_name, " | ", condition_B, " vs ", condition_A
    )
    return(NULL)
  }

  markers <- markers %>%
    rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC)) %>%
    mutate(
      cell_type = cell_type_name,
      contrast = paste(condition_B, "vs", condition_A),
      condition_A = condition_A,
      condition_B = condition_B,
      assay_used = assay,
      test_used = test_use,
      interpretation = "Exploratory cell-level DE only; cells are not biological replicates.",
      rank = row_number()
    )

  safe_cell_type <- gsub("[^A-Za-z0-9_]+", "_", cell_type_name)
  out_csv <- file.path(
    deg_dir,
    paste0("scDE_", safe_cell_type, "_", condition_B, "_vs_", condition_A, ".csv")
  )

  write.csv(markers, out_csv, row.names = FALSE)
  message("[STEP4] Saved DE table: ", out_csv)

  return(markers)
}


## ==========================================================
## 7. Run selected contrasts
## ==========================================================

all_cell_types <- as.character(sort(unique(brain$cell_type)))

preferred_cell_types <- c(
  "Radial_Glia",
  "Cycling_NPC",
  "Excitatory_Neuron",
  "Inhibitory_Neuron",
  "Astrocyte",
  "Microglia",
  "OPC"
)

cell_types_of_interest <- preferred_cell_types[
  preferred_cell_types %in% all_cell_types
]

if (length(cell_types_of_interest) == 0) {
  warning(
    "[STEP4] None of the preferred cell type names were found. ",
    "Using all non-Unknown cell types present in the object."
  )
  cell_types_of_interest <- setdiff(all_cell_types, "Unknown")
}

message("[STEP4] Cell types selected for exploratory DE:")
print(cell_types_of_interest)

contrasts <- tibble(
  condition_A = c("Mock", "Mock", "Mock"),
  condition_B = c("ZIKV_BR", "ZIKV_FSS", "IFNb")
)

all_scDE_results <- list()

for (ct in cell_types_of_interest) {
  for (i in seq_len(nrow(contrasts))) {
    A <- contrasts$condition_A[i]
    B <- contrasts$condition_B[i]
    res_name <- paste(ct, B, "vs", A, sep = "_")

    all_scDE_results[[res_name]] <- run_scDE_celltype(
      brain_obj = brain,
      cell_type_name = ct,
      condition_A = A,
      condition_B = B,
      assay = assay_to_use
    )
  }
}


## ==========================================================
## 8. Combine and summarize results
## ==========================================================

combined_deg <- all_scDE_results %>%
  purrr::discard(is.null) %>%
  dplyr::bind_rows()

if (nrow(combined_deg) == 0) {
  stop("[STEP4] No exploratory DE tables were produced.")
}

combined_path <- file.path(deg_dir, "scDE_all_contrasts_long.csv")
write.csv(combined_deg, combined_path, row.names = FALSE)
message("[STEP4] Saved combined DE table: ", combined_path)

p_adj_cutoff <- analysis_parameters$adjusted_p_cutoff_for_summary
logfc_cutoff <- analysis_parameters$logfc_cutoff_for_summary

deg_counts <- combined_deg %>%
  tibble::as_tibble() %>%
  dplyr::filter(
    p_val_adj < p_adj_cutoff,
    abs(avg_log2FC) >= logfc_cutoff
  ) %>%
  dplyr::count(
    cell_type,
    contrast,
    name = "n_genes_meeting_exploratory_criteria"
  ) %>%
  dplyr::arrange(cell_type, contrast)

deg_counts_long_path <- file.path(deg_dir, "scDE_DEG_counts_long.csv")
write.csv(deg_counts, deg_counts_long_path, row.names = FALSE)
message("[STEP4] Saved exploratory DEG count summary: ", deg_counts_long_path)

deg_counts_wide <- deg_counts %>%
  tidyr::pivot_wider(
    names_from = contrast,
    values_from = n_genes_meeting_exploratory_criteria,
    values_fill = 0
  )

deg_counts_wide_path <- file.path(deg_dir, "scDE_DEG_counts_wide.csv")
write.csv(deg_counts_wide, deg_counts_wide_path, row.names = FALSE)
message("[STEP4] Saved wide exploratory DEG count summary: ", deg_counts_wide_path)


## ==========================================================
## 9. Main barplot: exploratory gene counts
## ==========================================================

deg_counts$cell_type <- factor(deg_counts$cell_type, levels = cell_types_of_interest)

p_deg_counts <- ggplot(
  deg_counts,
  aes(
    x = cell_type,
    y = n_genes_meeting_exploratory_criteria,
    fill = contrast
  )
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  labs(
    title = "",
    x = "Cell type",
    y = "Genes meeting exploratory DE criteria",
    fill = "Contrast"
  ) +
  theme_journal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"))

save_tiff_step4(
  p_deg_counts,
  out_path = file.path(fig_step4_main, "scDE_exploratory_gene_counts_by_celltype_contrast.tiff"),
  width = 9,
  height = 5.5,
  dpi = 600
)


## ==========================================================
## 10. Supplementary volcano plots
## ==========================================================

plot_volcano <- function(df, cell_type_label, contrast_label) {
  df <- df %>%
    mutate(
      neg_log10_padj = -log10(p_val_adj + 1e-300),
      direction = case_when(
        p_val_adj < p_adj_cutoff & avg_log2FC >= logfc_cutoff ~ "Up",
        p_val_adj < p_adj_cutoff & avg_log2FC <= -logfc_cutoff ~ "Down",
        TRUE ~ "Not significant"
      )
    )

  ggplot(df, aes(x = avg_log2FC, y = neg_log10_padj, color = direction)) +
    geom_point(alpha = 0.70, size = 0.8) +
    geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(p_adj_cutoff),
      linetype = "dashed",
      linewidth = 0.4,
      color = "black"
    ) +
    scale_color_manual(
      values = c(
        "Down" = "#24b96c",
        "Not significant" = "grey75",
        "Up" = "#e76e5b"
      )
    ) +
    labs(
      title = paste0(cell_type_label, " (", contrast_label, ")"),
      x = "avg_log2FC",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_journal +
    theme(
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
      legend.position = "top",
      legend.direction = "horizontal"
    )
}

split_list <- combined_deg %>%
  filter(is.finite(p_val_adj), is.finite(avg_log2FC)) %>%
  split(list(.$cell_type, .$contrast), drop = TRUE)

volcano_plots <- list()

for (nm in names(split_list)) {
  df_sub <- split_list[[nm]]
  this_ct <- unique(df_sub$cell_type)
  this_contrast <- unique(df_sub$contrast)

  if (length(this_ct) != 1 || length(this_contrast) != 1) next

  p_vol <- plot_volcano(df_sub, this_ct, this_contrast)

  safe_ct <- gsub("[^A-Za-z0-9_]+", "_", this_ct)
  safe_contrast <- gsub("[^A-Za-z0-9_]+", "_", this_contrast)

  out_file <- file.path(
    fig_step4_supp,
    paste0("volcano_", safe_ct, "_", safe_contrast, ".tiff")
  )

  save_tiff_step4(
    p_vol,
    out_path = out_file,
    width = 5.5,
    height = 5.0,
    dpi = 600
  )

  volcano_plots[[paste(this_ct, this_contrast, sep = " | ")]] <- p_vol
}

if (length(volcano_plots) > 0) {
  n_cols <- 3
  n_rows <- ceiling(length(volcano_plots) / n_cols)

  volcano_grid_plot <- patchwork::wrap_plots(volcano_plots, ncol = n_cols)

  save_tiff_step4(
    volcano_grid_plot,
    out_path = file.path(fig_step4_main, "Supplementary_volcano_all_celltypes_all_contrasts.tiff"),
    width = 12,
    height = max(8, 4 * n_rows),
    dpi = 600
  )
}


## ==========================================================
## 11. Reproducibility files
## ==========================================================

parameters_path <- file.path(results_dir, "step4_analysis_parameters.txt")
capture.output(
  {
    cat("STEP 4 exploratory DE analysis parameters\n")
    cat("Assay used:", assay_to_use, "\n")
    cat("Available assays:", paste(assays_available, collapse = ", "), "\n\n")
    str(analysis_parameters)
  },
  file = parameters_path
)
message("[STEP4] Saved parameters: ", parameters_path)

interpretation_note <- c(
  "STEP 4 interpretation note",
  "",
  "This analysis uses Seurat::FindMarkers within annotated fetal neural cell types.",
  "The fetal neural dataset contains one biological sample per experimental condition.",
  "Cells are not biological replicates.",
  "Therefore, Wilcoxon tests, adjusted P values, and DEG counts are used only as exploratory cell-level descriptive summaries.",
  "They should not be interpreted as formal condition-level inference or statistically powered evidence of strain-dependent effects.",
  "Genes are described as meeting exploratory differential expression criteria when they satisfy the specified log2 fold-change, detection, and adjusted P-value thresholds.",
  "",
  paste0("Assay used for exploratory DE: ", assay_to_use)
)

interpretation_note_path <- file.path(results_dir, "step4_interpretation_note.txt")
writeLines(interpretation_note, interpretation_note_path)
message("[STEP4] Saved interpretation note: ", interpretation_note_path)

sink(file.path(session_dir, "STEP4_sessionInfo.txt"))
sessionInfo()
sink()
message("[STEP4] Saved session information.")

message("[STEP4] STEP 4 completed successfully.")

############################################################
# End of STEP 4
############################################################
