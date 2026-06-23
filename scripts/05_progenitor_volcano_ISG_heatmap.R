############################################################
# Zika scRNA-seq project – STEP 5
# Progenitor-focused exploratory volcano plots and ISG heatmap
#
# Purpose:
#   Generate descriptive figures from STEP 4
#   exploratory cell-type-level DE summaries:
#     - Volcano plots for progenitor-enriched radial glia-like cells
#       and cycling NPCs
#     - Core ISG log2FC heatmap across the same contrasts
#     - Combined multi-panel figure
#
# Analysis notes:
#   - STEP 4 DE tables were generated from cell-level Wilcoxon tests.
#   - The fetal neural dataset has one biological sample per condition.
#   - Cells are NOT biological replicates.
#   - These plots summarize exploratory transcriptional patterns only.
#   - P values/FDR values and DEG counts should not be interpreted as
#     formal condition-level statistical inference.
#   - ISG heatmap values represent transcriptional IFN/ISG
#     responsiveness and do not imply IFNβ production, viral burden,
#     infection frequency, susceptibility, or cell-intrinsic antiviral
#     initiation.
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
#   figures/step5/main/Figure4_progenitor_exploratory_patterns.tiff
#   figures/step5/main/volcano_6panel_progenitors.tiff
#   figures/step5/main/core_ISG_heatmap_progenitors.tiff
#   figures/step5/supplementary/volcano_*.tiff
#   results/step5_core_ISG_log2FC_long.csv
#   results/step5_core_ISG_log2FC_wide.csv
#   results/step5_interpretation_note.txt
#   session_info/STEP5_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

analysis_parameters <- list(
  seed = 1234,
  logfc_threshold = 0.25,
  adjusted_p_cutoff = 0.05,
  volcano_y_cap = 60,
  dpi = 600,
  core_isg = c(
    "IFITM1", "MX1", "OAS1", "OAS2", "OAS3",
    "IFIT1", "IFI6", "ISG15", "STAT1", "IRF7", "RSAD2"
  )
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP5] Project root: ", project_root)

required_packages <- c(
  "dplyr",
  "readr",
  "ggplot2",
  "patchwork",
  "tidyr",
  "purrr",
  "tibble",
  "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP5] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 5."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(purrr)
  library(tibble)
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

fig_step5_main <- file.path(figures_dir, "step5", "main")
fig_step5_supp <- file.path(figures_dir, "step5", "supplementary")

if (!dir.exists(results_dir)) {
  stop("[STEP5] results/ directory not found. Run previous steps first.")
}

if (!dir.exists(deg_dir)) {
  stop("[STEP5] results/DEG_sc/ directory not found. Run STEP 4 first.")
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step5_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step5_supp, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Plot theme and save helper
## ==========================================================

theme_journal <- function(base_size = 13) {
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

save_plot_tiff <- function(plot_obj, out_path, width, height, dpi = analysis_parameters$dpi) {
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
  message("[STEP5] Saved: ", out_path)
}


## ==========================================================
## 4. Input table
## ==========================================================

input_files <- tibble(
  panel_id = c(
    "RG_ZIKV_BR", "RG_ZIKV_FSS", "RG_IFNb",
    "NPC_ZIKV_BR", "NPC_ZIKV_FSS", "NPC_IFNb"
  ),
  cell_type = c(
    "Radial_Glia", "Radial_Glia", "Radial_Glia",
    "Cycling_NPC", "Cycling_NPC", "Cycling_NPC"
  ),
  cell_type_label = c(
    "Radial glia-like", "Radial glia-like", "Radial glia-like",
    "Cycling NPC", "Cycling NPC", "Cycling NPC"
  ),
  condition = c("ZIKV_BR", "ZIKV_FSS", "IFNb", "ZIKV_BR", "ZIKV_FSS", "IFNb"),
  contrast_label = c(
    "ZIKV-BR vs Mock", "ZIKV-FSS vs Mock", "IFNβ vs Mock",
    "ZIKV-BR vs Mock", "ZIKV-FSS vs Mock", "IFNβ vs Mock"
  ),
  filename = c(
    "scDE_Radial_Glia_ZIKV_BR_vs_Mock.csv",
    "scDE_Radial_Glia_ZIKV_FSS_vs_Mock.csv",
    "scDE_Radial_Glia_IFNb_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_BR_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_FSS_vs_Mock.csv",
    "scDE_Cycling_NPC_IFNb_vs_Mock.csv"
  )
)

missing_files <- input_files$filename[
  !file.exists(file.path(deg_dir, input_files$filename))
]

if (length(missing_files) > 0) {
  stop(
    "[STEP5] Missing required STEP 4 DE file(s):\n",
    paste(file.path(deg_dir, missing_files), collapse = "\n"),
    "\nRun STEP 4 successfully before STEP 5."
  )
}


## ==========================================================
## 5. Volcano plot helper
## ==========================================================

make_volcano <- function(csv_file,
                         title,
                         lfc_thresh = analysis_parameters$logfc_threshold,
                         padj_thresh = analysis_parameters$adjusted_p_cutoff,
                         y_cap = analysis_parameters$volcano_y_cap) {
  
  full_path <- file.path(deg_dir, csv_file)
  
  de <- readr::read_csv(full_path, show_col_types = FALSE)
  
  required_cols <- c("gene", "avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(de))
  
  if (length(missing_cols) > 0) {
    stop(
      "[STEP5] Missing required columns in ", full_path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  de <- de %>%
    mutate(
      logFC = avg_log2FC,
      padj = p_val_adj,
      neg_log10_padj = -log10(padj + 1e-300),
      neg_log10_padj = pmin(neg_log10_padj, y_cap),
      direction = case_when(
        padj < padj_thresh & logFC >= lfc_thresh ~ "Up",
        padj < padj_thresh & logFC <= -lfc_thresh ~ "Down",
        TRUE ~ "Not significant"
      )
    )
  
  x_lim <- quantile(de$logFC, probs = c(0.01, 0.99), na.rm = TRUE)
  max_abs <- max(abs(x_lim), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  
  ggplot(de, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(aes(color = direction), alpha = 0.80, size = 1.2) +
    geom_vline(
      xintercept = c(-lfc_thresh, lfc_thresh),
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    geom_hline(
      yintercept = -log10(padj_thresh),
      linetype = "dashed",
      color = "black",
      linewidth = 0.4
    ) +
    scale_color_manual(
      values = c(
        "Down" = "#24b96c",
        "Not significant" = "grey78",
        "Up" = "#e76e5b"
      ),
      breaks = c("Down", "Not significant", "Up")
    ) +
    coord_cartesian(xlim = c(-max_abs, max_abs)) +
    labs(
      x = "Log2 fold change",
      y = "-Log10 adjusted p-value",
      title = title
    ) +
    theme_journal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
}


## ==========================================================
## 6. Build and save volcano plots
## ==========================================================

volcano_plots <- purrr::pmap(
  list(input_files$filename, input_files$cell_type_label, input_files$contrast_label),
  function(filename, cell_type_label, contrast_label) {
    make_volcano(
      csv_file = filename,
      title = paste0(cell_type_label, " – ", contrast_label)
    )
  }
)

names(volcano_plots) <- input_files$panel_id

for (i in seq_along(volcano_plots)) {
  out_path <- file.path(
    fig_step5_supp,
    paste0("volcano_", input_files$panel_id[i], ".tiff")
  )
  
  save_plot_tiff(
    volcano_plots[[i]],
    out_path = out_path,
    width = 4.2,
    height = 4.0
  )
}

# Robust patchwork assembly.
# Using wrap_plots() avoids layout errors such as:
# "Need 3 panels, but together nrow and ncol only provide 2."
volcano_grid <- patchwork::wrap_plots(
  volcano_plots[
    c(
      "RG_ZIKV_BR", "RG_ZIKV_FSS", "RG_IFNb",
      "NPC_ZIKV_BR", "NPC_ZIKV_FSS", "NPC_IFNb"
    )
  ],
  ncol = 3,
  byrow = TRUE
)

save_plot_tiff(
  volcano_grid,
  out_path = file.path(fig_step5_main, "volcano_6panel_progenitors.tiff"),
  width = 13,
  height = 8
)


## ==========================================================
## 7. Core ISG heatmap
## ==========================================================

core_isg <- analysis_parameters$core_isg

read_isg_logfc <- function(csv_file, panel_id, cell_type_label, contrast_label) {
  full_path <- file.path(deg_dir, csv_file)
  
  de <- readr::read_csv(full_path, show_col_types = FALSE)
  
  tibble(gene = core_isg) %>%
    left_join(
      de %>%
        select(gene, avg_log2FC, p_val_adj),
      by = "gene"
    ) %>%
    mutate(
      panel_id = panel_id,
      cell_type = cell_type_label,
      contrast = contrast_label,
      missing_from_de_table = is.na(avg_log2FC)
    )
}

isg_long <- purrr::pmap_dfr(
  list(
    input_files$filename,
    input_files$panel_id,
    input_files$cell_type_label,
    input_files$contrast_label
  ),
  read_isg_logfc
)

isg_long_path <- file.path(results_dir, "step5_core_ISG_log2FC_long.csv")
readr::write_csv(isg_long, isg_long_path)
message("[STEP5] Saved ISG long table: ", isg_long_path)

isg_wide <- isg_long %>%
  select(gene, panel_id, avg_log2FC) %>%
  tidyr::pivot_wider(names_from = panel_id, values_from = avg_log2FC) %>%
  arrange(match(gene, core_isg))

isg_wide_path <- file.path(results_dir, "step5_core_ISG_log2FC_wide.csv")
readr::write_csv(isg_wide, isg_wide_path)
message("[STEP5] Saved ISG wide table: ", isg_wide_path)

# Row scaling for visualization only.
# Missing values are kept as missing and displayed separately via na.value.
isg_heatmap_df <- isg_wide %>%
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "panel_id",
    values_to = "avg_log2FC"
  ) %>%
  group_by(gene) %>%
  mutate(
    row_mean = mean(avg_log2FC, na.rm = TRUE),
    row_sd = sd(avg_log2FC, na.rm = TRUE),
    scaled_log2FC = case_when(
      is.na(avg_log2FC) ~ NA_real_,
      is.na(row_sd) | row_sd == 0 ~ 0,
      TRUE ~ (avg_log2FC - row_mean) / row_sd
    )
  ) %>%
  ungroup() %>%
  left_join(
    input_files %>% select(panel_id, cell_type_label, contrast_label),
    by = "panel_id"
  ) %>%
  mutate(
    gene = factor(gene, levels = rev(core_isg)),
    panel_id = factor(panel_id, levels = input_files$panel_id),
    panel_label = paste(cell_type_label, contrast_label, sep = "\n")
  )

panel_labels <- input_files %>%
  mutate(panel_label = paste(cell_type_label, contrast_label, sep = "\n")) %>%
  pull(panel_label)

names(panel_labels) <- input_files$panel_id

p_isg_heatmap <- ggplot(
  isg_heatmap_df,
  aes(x = panel_id, y = gene, fill = scaled_log2FC)
) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_x_discrete(labels = panel_labels) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey80",
    name = "Row-scaled\nlog2FC"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Canonical ISG transcriptional responsiveness"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(face = "bold", color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold"),
    panel.grid = element_blank()
  )

save_plot_tiff(
  p_isg_heatmap,
  out_path = file.path(fig_step5_main, "core_ISG_heatmap_progenitors.tiff"),
  width = 8.5,
  height = 5.0
)


## ==========================================================
## 8. Combined main figure
## ==========================================================

# Robust combined figure assembly.
# Do not use the `/` patchwork operator here because nested layouts can fail
# in some R/patchwork versions. wrap_plots(list(...), ncol = 1) is more stable.
combined_figure <- patchwork::wrap_plots(
  list(volcano_grid, p_isg_heatmap),
  ncol = 1,
  heights = c(2.1, 1.15)
)

save_plot_tiff(
  combined_figure,
  out_path = file.path(fig_step5_main, "Figure4_progenitor_exploratory_patterns.tiff"),
  width = 13,
  height = 13
)


## ==========================================================
## 9. Reproducibility files
## ==========================================================

parameters_path <- file.path(results_dir, "step5_analysis_parameters.txt")
capture.output(
  {
    cat("STEP 5 progenitor volcano and ISG heatmap parameters\n\n")
    str(analysis_parameters)
  },
  file = parameters_path
)
message("[STEP5] Saved parameters: ", parameters_path)

interpretation_note <- c(
  "STEP 5 interpretation note",
  "",
  "This step visualizes exploratory cell-level DE summaries generated in STEP 4.",
  "The fetal neural dataset contains one biological sample per experimental condition.",
  "Cells are not biological replicates.",
  "Therefore, the volcano plots and ISG heatmap should be interpreted as descriptive visual summaries only.",
  "The ISG heatmap summarizes transcriptional IFN/ISG responsiveness.",
  "It does not measure IFNβ production, viral RNA burden, infection frequency, susceptibility, or cell-intrinsic antiviral initiation.",
  "Missing ISG log2FC values are retained as missing values and displayed with the heatmap missing-value color rather than being imputed as zero."
)

interpretation_note_path <- file.path(results_dir, "step5_interpretation_note.txt")
writeLines(interpretation_note, interpretation_note_path)
message("[STEP5] Saved interpretation note: ", interpretation_note_path)

sink(file.path(session_dir, "STEP5_sessionInfo.txt"))
sessionInfo()
sink()
message("[STEP5] Saved session information.")

message("[STEP5] STEP 5 completed successfully.")

############################################################
# End of STEP 5
############################################################
