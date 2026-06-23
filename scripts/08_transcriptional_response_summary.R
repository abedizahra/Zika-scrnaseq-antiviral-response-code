############################################################
# Zika scRNA-seq project – STEP 8
# Transcriptional response magnitude and ISG summary
#
# Purpose:
#   Summarize transcriptional response magnitude in
#   progenitor-enriched radial glia-like cells and cycling NPCs:
#     - Number of genes meeting exploratory DE criteria, split into
#       upregulated and downregulated genes
#     - Mean log2FC of the curated core ISG panel
#     - Two-panel summary figure for downstream figure assembly
#
# Notes:
#   - STEP 4 DE tables were generated using cell-level Wilcoxon tests.
#   - The fetal neural dataset has one biological sample per condition.
#   - Cells are not biological replicates.
#   - DEG counts are descriptive summaries only and should not be interpreted
#     as formal condition-level inference.
#   - Mean ISG log2FC summarizes transcriptional IFN/ISG responsiveness and
#     does not measure IFNβ production, viral RNA burden, infection frequency,
#     susceptibility, antiviral efficacy, or cell-intrinsic antiviral initiation.
#
# Inputs:
#   results/DEG_sc/scDE_Radial_Glia_ZIKV_BR_vs_Mock.csv
#   results/DEG_sc/scDE_Radial_Glia_ZIKV_FSS_vs_Mock.csv
#   results/DEG_sc/scDE_Radial_Glia_IFNb_vs_Mock.csv
#   results/DEG_sc/scDE_Cycling_NPC_ZIKV_BR_vs_Mock.csv
#   results/DEG_sc/scDE_Cycling_NPC_ZIKV_FSS_vs_Mock.csv
#   results/DEG_sc/scDE_Cycling_NPC_IFNb_vs_Mock.csv
#
# Outputs:
#   results/step8_DEG_ISG_summary.csv
#   results/step8_DEG_ISG_summary_long_direction.csv
#   figures/step8/main/Figure6_AB_transcriptional_response_summary.tiff
#   figures/step8/main/step8_DEG_burden_up_down.tiff
#   figures/step8/main/step8_mean_core_ISG_log2FC.tiff
#   results/step8_interpretation_note.txt
#   session_info/STEP8_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

analysis_parameters <- list(
  seed = 1234,
  logfc_cutoff = 0.25,
  adjusted_p_cutoff = 0.05,
  core_isg = c(
    "IFITM1", "MX1", "OAS1", "OAS2", "OAS3",
    "IFIT1", "IFI6", "ISG15", "STAT1", "IRF7", "RSAD2"
  )
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP8] Project root: ", project_root)

required_packages <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "purrr",
  "tidyr",
  "tibble",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP8] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 8."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(scales)
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)


## ==========================================================
## 2. Directories
## ==========================================================

results_dir <- file.path(project_root, "results")
deg_dir <- file.path(results_dir, "DEG_sc")
figures_dir <- file.path(project_root, "figures")
session_dir <- file.path(project_root, "session_info")

fig_step8_main <- file.path(figures_dir, "step8", "main")
fig_step8_supp <- file.path(figures_dir, "step8", "supplementary")

if (!dir.exists(results_dir)) {
  stop("[STEP8] results/ directory not found. Run previous steps first.")
}

if (!dir.exists(deg_dir)) {
  stop("[STEP8] results/DEG_sc/ directory not found. Run STEP 4 first.")
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step8_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step8_supp, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Plot theme and save helper
## ==========================================================

step8_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0.5),
      axis.title.x = element_text(face = "bold", size = base_size, color = "black"),
      axis.title.y = element_text(face = "bold", size = base_size, color = "black"),
      axis.text.x = element_text(face = "bold", size = base_size - 2, color = "black"),
      axis.text.y = element_text(face = "bold", size = base_size - 2, color = "black"),
      legend.text = element_text(face = "bold", size = base_size - 2),
      legend.title = element_blank(),
      legend.key.size = unit(0.55, "cm"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

save_plot_tiff <- function(plot_obj, out_path, width, height, dpi = 600) {
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
  message("[STEP8] Saved: ", out_path)
}


## ==========================================================
## 4. Input file table
## ==========================================================

files_info <- tibble(
  cell_type = rep(c("Radial glia-like", "Cycling NPC"), each = 3),
  cell_type_file = rep(c("Radial_Glia", "Cycling_NPC"), each = 3),
  condition = rep(c("ZIKV_BR", "ZIKV_FSS", "IFNb"), times = 2),
  condition_label = rep(c("ZIKV-BR", "ZIKV-FSS", "IFNβ"), times = 2),
  filename = c(
    "scDE_Radial_Glia_ZIKV_BR_vs_Mock.csv",
    "scDE_Radial_Glia_ZIKV_FSS_vs_Mock.csv",
    "scDE_Radial_Glia_IFNb_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_BR_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_FSS_vs_Mock.csv",
    "scDE_Cycling_NPC_IFNb_vs_Mock.csv"
  )
)

missing_files <- files_info$filename[
  !file.exists(file.path(deg_dir, files_info$filename))
]

if (length(missing_files) > 0) {
  stop(
    "[STEP8] Missing required STEP 4 DE file(s):\n",
    paste(file.path(deg_dir, missing_files), collapse = "\n"),
    "\nRun STEP 4 successfully before STEP 8."
  )
}


## ==========================================================
## 5. Summarize DE tables
## ==========================================================

lfc_thresh <- analysis_parameters$logfc_cutoff
padj_thresh <- analysis_parameters$adjusted_p_cutoff
core_isg <- analysis_parameters$core_isg

summarize_de_table <- function(filename, cell_type, condition, condition_label) {
  full_path <- file.path(deg_dir, filename)
  de <- readr::read_csv(full_path, show_col_types = FALSE)

  required_cols <- c("gene", "avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(de))

  if (length(missing_cols) > 0) {
    stop(
      "[STEP8] Missing required columns in ", full_path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }

  de_sig <- de %>%
    filter(
      p_val_adj < padj_thresh,
      abs(avg_log2FC) >= lfc_thresh
    )

  isg_sub <- de %>%
    filter(gene %in% core_isg)

  missing_isg <- setdiff(core_isg, isg_sub$gene)

  tibble(
    cell_type = cell_type,
    condition = condition,
    condition_label = condition_label,
    n_up = sum(de_sig$avg_log2FC > 0, na.rm = TRUE),
    n_down = sum(de_sig$avg_log2FC < 0, na.rm = TRUE),
    n_genes_meeting_exploratory_criteria = nrow(de_sig),
    mean_ISG_log2FC = mean(isg_sub$avg_log2FC, na.rm = TRUE),
    n_ISGs_present = nrow(isg_sub),
    n_ISGs_missing = length(missing_isg),
    missing_ISGs = paste(missing_isg, collapse = ";"),
    interpretation = "Descriptive exploratory summary only; cells are not biological replicates."
  )
}

summary_df <- purrr::pmap_dfr(
  list(
    files_info$filename,
    files_info$cell_type,
    files_info$condition,
    files_info$condition_label
  ),
  summarize_de_table
)

summary_path <- file.path(results_dir, "step8_DEG_ISG_summary.csv")
readr::write_csv(summary_df, summary_path)
message("[STEP8] Saved summary table: ", summary_path)

summary_long <- summary_df %>%
  select(cell_type, condition, condition_label, n_up, n_down) %>%
  tidyr::pivot_longer(
    cols = c(n_up, n_down),
    names_to = "direction",
    values_to = "count"
  ) %>%
  mutate(
    direction = recode(
      direction,
      n_up = "Upregulated",
      n_down = "Downregulated"
    ),
    condition_label = factor(condition_label, levels = c("ZIKV-BR", "ZIKV-FSS", "IFNβ")),
    cell_type = factor(cell_type, levels = c("Radial glia-like", "Cycling NPC"))
  )

summary_long_path <- file.path(results_dir, "step8_DEG_ISG_summary_long_direction.csv")
readr::write_csv(summary_long, summary_long_path)
message("[STEP8] Saved long summary table: ", summary_long_path)


## ==========================================================
## 6. Panel A – exploratory DEG burden
## ==========================================================

p_deg <- ggplot(
  summary_long,
  aes(x = condition_label, y = count, fill = direction)
) +
  geom_col(width = 0.62) +
  facet_wrap(~ cell_type, nrow = 1) +
  scale_fill_manual(
    values = c(
      "Upregulated" = "#e76e5b",
      "Downregulated" = "#24b96c"
    )
  ) +
  labs(
    x = "Condition",
    y = "Number of genes",
    title = "Genes meeting exploratory DE criteria"
  ) +
  step8_theme(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")
  )

save_plot_tiff(
  p_deg,
  out_path = file.path(fig_step8_main, "step8_DEG_burden_up_down.tiff"),
  width = 8.5,
  height = 4.8
)


## ==========================================================
## 7. Panel B – mean ISG log2FC
## ==========================================================

isg_plot_df <- summary_df %>%
  mutate(
    condition_label = factor(condition_label, levels = c("ZIKV-BR", "ZIKV-FSS", "IFNβ")),
    cell_type = factor(cell_type, levels = c("Radial glia-like", "Cycling NPC"))
  )

p_isg <- ggplot(
  isg_plot_df,
  aes(x = condition_label, y = mean_ISG_log2FC, fill = cell_type)
) +
  geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.45) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  scale_fill_manual(
    values = c(
      "Radial glia-like" = "#7BC67B",
      "Cycling NPC" = "#B58AD9"
    )
  ) +
  labs(
    x = "Condition",
    y = "Mean log2FC of core ISGs",
    title = "Core ISG transcriptional responsiveness"
  ) +
  step8_theme(base_size = 13) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold")
  )

save_plot_tiff(
  p_isg,
  out_path = file.path(fig_step8_main, "step8_mean_core_ISG_log2FC.tiff"),
  width = 8.0,
  height = 4.8
)


## ==========================================================
## 8. Assemble Figure 6 panels A–B
## ==========================================================

figure6_ab <- p_deg / p_isg +
  patchwork::plot_layout(heights = c(1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

save_plot_tiff(
  figure6_ab,
  out_path = file.path(fig_step8_main, "Figure6_AB_transcriptional_response_summary.tiff"),
  width = 9,
  height = 8.5
)

## ==========================================================
## 9. Panels C and D – contextual Vero/moDC ISG analysis
## ==========================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(patchwork)

results_dir <- file.path(getwd(), "results")
figures_dir <- file.path(getwd(), "figures")
fig_step9_main <- file.path(figures_dir, "step9", "main")

dir.create(fig_step9_main, recursive = TRUE, showWarnings = FALSE)

dc_path <- file.path(results_dir, "dc_step2_norm_qc.rds")

if (!file.exists(dc_path)) {
  stop("Cannot find results/dc_step2_norm_qc.rds")
}

dc <- readRDS(dc_path)

if (!inherits(dc, "Seurat")) {
  stop("dc_step2_norm_qc.rds is not a Seurat object.")
}

## ----------------------------------------------------------
## 10. define labels for Vero vs moDC
## ----------------------------------------------------------

meta_cols <- colnames(dc@meta.data)

if (!"condition_group" %in% meta_cols) {
  if ("sample_id" %in% meta_cols) {
    dc$condition_group <- dplyr::case_when(
      grepl("ZIKV-Vero-cells", dc$sample_id) ~ "ZIKV_Vero_cells",
      grepl("p22086", dc$sample_id) ~ "moDC_panel",
      TRUE ~ "Unknown"
    )
  } else {
    stop("Neither 'condition_group' nor 'sample_id' exists in dc metadata.")
  }
}

dc$condition_group_plot <- dplyr::case_when(
  dc$condition_group %in% c("ZIKV_Vero_cells", "ZIKV-Vero-cells") ~ "ZIKV-exposed Vero cells",
  dc$condition_group %in% c("p22086_panel", "moDC_panel", "p22086") ~ "Primary moDCs",
  grepl("Vero", dc$condition_group) ~ "ZIKV-exposed Vero cells",
  grepl("p22086|moDC|DC", dc$condition_group) ~ "Primary moDCs",
  TRUE ~ "Unknown"
)

dc$condition_group_plot <- factor(
  dc$condition_group_plot,
  levels = c("ZIKV-exposed Vero cells", "Primary moDCs", "Unknown")
)

print(table(dc$condition_group_plot, useNA = "ifany"))

## ----------------------------------------------------------
## 11. choose assay
## ----------------------------------------------------------

available_assays <- tryCatch(
  as.character(Seurat::Assays(dc)),
  error = function(e) character(0)
)

if ("SCT" %in% available_assays) {
  assay_to_use <- "SCT"
} else if ("RNA" %in% available_assays) {
  assay_to_use <- "RNA"
} else {
  assay_to_use <- DefaultAssay(dc)
}

DefaultAssay(dc) <- assay_to_use

## if RNA is used, make sure normalized data exist
if (assay_to_use == "RNA") {
  has_data_layer <- TRUE
  
  invisible(tryCatch(
    {
      tmp <- SeuratObject::LayerData(dc, assay = "RNA", layer = "data")
      if (ncol(tmp) == 0 || nrow(tmp) == 0) has_data_layer <<- FALSE
    },
    error = function(e) {
      has_data_layer <<- FALSE
    }
  ))
  
  if (!has_data_layer) {
    dc <- NormalizeData(
      dc,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
}

## ----------------------------------------------------------
## 12. ISG module score
## ----------------------------------------------------------

core_isg <- c(
  "IFITM1", "MX1", "OAS1", "OAS2", "OAS3",
  "IFIT1", "IFI6", "ISG15", "STAT1", "IRF7", "RSAD2"
)

genes_available <- rownames(dc[[assay_to_use]])
core_isg_present <- intersect(core_isg, genes_available)
core_isg_missing <- setdiff(core_isg, core_isg_present)

print(core_isg_present)
print(core_isg_missing)

if (length(core_isg_present) < 3) {
  stop("Too few ISG genes found in the selected assay.")
}

set.seed(1234)

dc <- AddModuleScore(
  object = dc,
  features = list(core_isg_present),
  name = "ISGscore",
  assay = assay_to_use
)

## summary table
dc_isg_summary <- dc@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(condition_group_plot) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_ISGscore = mean(ISGscore1, na.rm = TRUE),
    median_ISGscore = median(ISGscore1, na.rm = TRUE),
    sd_ISGscore = sd(ISGscore1, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_csv(
  dc_isg_summary,
  file.path(results_dir, "dc_ISGscore_summary_by_condition_group.csv")
)

## ----------------------------------------------------------
## 13. plot theme helper
## ----------------------------------------------------------

journal_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", color = "black"),
      axis.text.x = element_text(size = base_size - 2, face = "bold", color = "black"),
      axis.text.y = element_text(size = base_size - 2, face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2, face = "bold"),
      legend.key.size = unit(0.55, "cm"),
      axis.line = element_line(linewidth = 0.8, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

save_plot_tiff <- function(plot_obj, out_path, width, height, dpi = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
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
}

## ----------------------------------------------------------
## 14. Panel C – violin plot
## ----------------------------------------------------------

meta_plot <- dc@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell_id") %>%
  dplyr::filter(!is.na(condition_group_plot))

p_C <- ggplot(
  meta_plot,
  aes(x = condition_group_plot, y = ISGscore1, fill = condition_group_plot)
) +
  geom_violin(width = 0.78, scale = "width", trim = TRUE, alpha = 0.85) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.35) +
  scale_fill_manual(
    values = c(
      "ZIKV-exposed Vero cells" = "#9ecae1",
      "Primary moDCs" = "#fdae6b",
      "Unknown" = "grey80"
    )
  ) +
  labs(
    x = NULL,
    y = "ISG module score",
    title = "Core ISG program in Vero/moDC dataset"
  ) +
  journal_theme(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1, face = "bold")
  )

save_plot_tiff(
  p_C,
  file.path(fig_step9_main, "Figure6_panel_C_dc_ISGscore_violin.tiff"),
  width = 5.2,
  height = 4.4
)

## ----------------------------------------------------------
## 15. Panel D – UMAP colored by ISG score
## ----------------------------------------------------------

reduction_names <- names(dc@reductions)

if (!"umap" %in% reduction_names) {
  stop("UMAP reduction not found in dc object.")
}

p_D <- FeaturePlot(
  dc,
  features = "ISGscore1",
  reduction = "umap",
  order = TRUE
) +
  scale_colour_viridis_c(option = "plasma") +
  labs(
    title = "ISG module score in Vero/moDC cells"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

save_plot_tiff(
  p_D,
  file.path(fig_step9_main, "Figure6_panel_D_dc_ISGscore_UMAP.tiff"),
  width = 5.8,
  height = 4.6
)


figure6 <- (p_deg | p_isg) / (p_C | p_D) +
  patchwork::plot_layout(heights = c(1, 1.05), widths = c(1, 1)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 18))

save_plot_tiff(
  figure6,
  file.path(fig_step9_main, "Figure6_contextual_ISG_comparison.tiff"),
  width = 13,
  height = 10
)
## ==========================================================
## 16. Reproducibility files
## ==========================================================

parameters_path <- file.path(results_dir, "step8_analysis_parameters.txt")
capture.output(
  {
    cat("STEP 8 transcriptional response magnitude and ISG summary\n\n")
    str(analysis_parameters)
  },
  file = parameters_path
)
message("[STEP8] Saved parameters: ", parameters_path)

interpretation_note <- c(
  "STEP 8 interpretation note",
  "",
  "This step summarizes STEP 4 exploratory cell-level DE tables.",
  "The fetal neural dataset contains one biological sample per experimental condition.",
  "Cells are not biological replicates.",
  "Therefore, DEG burden values are descriptive summaries only.",
  "They should not be interpreted as statistically powered condition-level differences or definitive strain-dependent effects.",
  "Mean ISG log2FC summarizes transcriptional IFN/ISG responsiveness.",
  "It does not measure IFNβ production, viral RNA burden, infection frequency, susceptibility, antiviral efficacy, or cell-intrinsic antiviral initiation."
)

interpretation_note_path <- file.path(results_dir, "step8_interpretation_note.txt")
writeLines(interpretation_note, interpretation_note_path)
message("[STEP8] Saved interpretation note: ", interpretation_note_path)

sink(file.path(session_dir, "STEP8_sessionInfo.txt"))
sessionInfo()
sink()
message("[STEP8] Saved session information.")

message("[STEP8] STEP 8 completed successfully.")

############################################################
# End of STEP 8
############################################################
