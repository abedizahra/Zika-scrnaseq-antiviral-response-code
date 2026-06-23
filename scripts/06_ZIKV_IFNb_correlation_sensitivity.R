############################################################
# Zika scRNA-seq project – STEP 6
# ZIKV–IFNβ transcriptional similarity sensitivity analysis
#
# Purpose:
#   Compare ZIKV-associated and IFNβ-associated log2FC patterns
#   in progenitor-enriched radial glia-like cells and cycling NPCs.
#
# Correlation sensitivity analyses:
#   1. Specify the primary gene universe used for correlation.
#   2. Add sensitivity analyses using:
#        - primary exploratory DE-table genes
#        - expanded all-expressed/tested genes from no-threshold FindMarkers
#        - Pearson correlation
#        - Spearman correlation
#        - canonical ISGs only
#        - non-ISG genes
#        - genes excluding the canonical ISG panel
#   3. Interpret correlations descriptively, not as functional equivalence.
#
# Analysis notes:
#   - The fetal neural dataset has one biological sample per condition.
#   - Cells are NOT biological replicates.
#   - Correlations are descriptive summaries of log2FC alignment.
#   - High ZIKV–IFNβ correlation can be inflated by gene-universe choice,
#     especially when using filtered exploratory DE tables enriched for ISGs.
#   - These analyses do not establish functional equivalence between
#     virus-associated and cytokine-stimulated programs.
#
# Inputs:
#   results/brain_step3_annotated.rds
#   results/DEG_sc/scDE_*.csv from STEP 4
#
# Outputs:
#   results/Supplementary_Table_1_ZIKV_IFNb_correlation_sensitivity.csv
#   results/correlation_sensitivity/*.csv
#   figures/step6/main/Figure5_ZIKV_IFNb_transcriptional_similarity.tiff
#   figures/step6/main/ZIKV_IFNb_primary_correlation_scatterplots.tiff
#   figures/step6/main/ZIKV_IFNb_correlation_sensitivity_summary.tiff
#   figures/step6/supplementary/*.tiff
#   results/step6_interpretation_note.txt
#   session_info/STEP6_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

analysis_parameters <- list(
  seed = 1234,
  primary_gene_universe = "primary_exploratory_DE_tables_from_STEP4",
  expanded_gene_universe = "FindMarkers logfc.threshold = 0, min.pct = 0.01",
  test_use = "wilcox",
  min_cells_per_group = 20,
  expanded_logfc_threshold = 0,
  expanded_min_pct = 0.01,
  core_isg = c(
    "IFITM1", "MX1", "OAS1", "OAS2", "OAS3",
    "IFIT1", "IFI6", "ISG15", "STAT1", "IRF7", "RSAD2"
  )
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP6] Project root: ", project_root)

required_packages <- c(
  "Seurat",
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
    "[STEP6] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 6."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
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

fig_step6_main <- file.path(figures_dir, "step6", "main")
fig_step6_supp <- file.path(figures_dir, "step6", "supplementary")
sensitivity_dir <- file.path(results_dir, "correlation_sensitivity")

if (!dir.exists(results_dir)) {
  stop("[STEP6] results/ directory not found. Run previous steps first.")
}

if (!dir.exists(deg_dir)) {
  stop("[STEP6] results/DEG_sc/ directory not found. Run STEP 4 first.")
}

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step6_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step6_supp, recursive = TRUE, showWarnings = FALSE)
dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Plot theme and save helper
## ==========================================================

step6_theme <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.y = element_text(size = base_size, face = "bold", color = "black"),
      axis.text = element_text(size = base_size - 2, face = "bold", color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 2, face = "bold"),
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
  message("[STEP6] Saved: ", out_path)
}


## ==========================================================
## 4. Core functions
## ==========================================================

core_isg <- analysis_parameters$core_isg

safe_cor <- function(x, y, method = c("pearson", "spearman")) {
  method <- match.arg(method)
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = method))
}

read_primary_de <- function(filename, prefix) {
  full_path <- file.path(deg_dir, filename)
  
  if (!file.exists(full_path)) {
    stop("[STEP6] Missing DE file: ", full_path)
  }
  
  de <- readr::read_csv(full_path, show_col_types = FALSE)
  
  required_cols <- c("gene", "avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(de))
  
  if (length(missing_cols) > 0) {
    stop(
      "[STEP6] Missing required columns in ", full_path, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  out <- de %>%
    dplyr::select(gene, avg_log2FC, p_val_adj)
  
  colnames(out) <- c(
    "gene",
    paste0("logFC_", prefix),
    paste0("padj_", prefix)
  )
  
  out
}

build_pair <- function(file_zikv,
                       file_ifnb,
                       zikv_label,
                       cell_type_label,
                       gene_universe_label) {
  de_zikv <- read_primary_de(file_zikv, zikv_label)
  de_ifnb <- read_primary_de(file_ifnb, "IFNb")

  full_join(de_zikv, de_ifnb, by = "gene") %>%
    mutate(
      cell_type = cell_type_label,
      zikv_condition = zikv_label,
      gene_universe = gene_universe_label,
      is_ISG = gene %in% core_isg
    )
}

summarize_correlations <- function(df, xcol, ycol, gene_universe_label) {
  df_complete <- df %>%
    filter(!is.na(.data[[xcol]]), !is.na(.data[[ycol]]))

  make_one <- function(dat, gene_class_label) {
    tibble(
      gene_universe = gene_universe_label,
      gene_class = gene_class_label,
      n_genes = nrow(dat),
      pearson_r = safe_cor(dat[[xcol]], dat[[ycol]], method = "pearson"),
      spearman_r = safe_cor(dat[[xcol]], dat[[ycol]], method = "spearman")
    )
  }

  bind_rows(
    make_one(df_complete, "all_genes_in_universe"),
    make_one(df_complete %>% filter(is_ISG), "canonical_ISGs_only"),
    make_one(df_complete %>% filter(!is_ISG), "non_ISG_genes"),
    make_one(df_complete %>% filter(!gene %in% core_isg), "excluding_canonical_ISGs")
  )
}


## ==========================================================
## 5. Primary correlation from STEP 4 exploratory DE tables
## ==========================================================

message("[STEP6] Building primary correlation datasets from STEP 4 DE tables...")

primary_pairs_info <- tibble(
  comparison = c("RG_BR", "RG_FSS", "NPC_BR", "NPC_FSS"),
  cell_type = c("Radial glia-like", "Radial glia-like", "Cycling NPC", "Cycling NPC"),
  zikv_condition = c("BR", "FSS", "BR", "FSS"),
  file_zikv = c(
    "scDE_Radial_Glia_ZIKV_BR_vs_Mock.csv",
    "scDE_Radial_Glia_ZIKV_FSS_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_BR_vs_Mock.csv",
    "scDE_Cycling_NPC_ZIKV_FSS_vs_Mock.csv"
  ),
  file_ifnb = c(
    "scDE_Radial_Glia_IFNb_vs_Mock.csv",
    "scDE_Radial_Glia_IFNb_vs_Mock.csv",
    "scDE_Cycling_NPC_IFNb_vs_Mock.csv",
    "scDE_Cycling_NPC_IFNb_vs_Mock.csv"
  )
)

primary_pairs <- purrr::pmap(
  list(
    primary_pairs_info$file_zikv,
    primary_pairs_info$file_ifnb,
    primary_pairs_info$zikv_condition,
    primary_pairs_info$cell_type
  ),
  function(file_zikv, file_ifnb, zikv_condition, cell_type) {
    build_pair(
      file_zikv = file_zikv,
      file_ifnb = file_ifnb,
      zikv_label = zikv_condition,
      cell_type_label = cell_type,
      gene_universe_label = "primary_exploratory_DE_tables"
    )
  }
)

names(primary_pairs) <- primary_pairs_info$comparison

primary_corr_summary <- purrr::imap_dfr(primary_pairs, function(df, nm) {
  zikv_label <- unique(df$zikv_condition)
  ycol <- paste0("logFC_", zikv_label)
  xcol <- "logFC_IFNb"

  summarize_correlations(
    df = df,
    xcol = xcol,
    ycol = ycol,
    gene_universe_label = "primary_exploratory_DE_tables"
  ) %>%
    mutate(
      comparison = nm,
      cell_type = unique(df$cell_type),
      zikv_condition = zikv_label
    )
})

primary_pairs_long <- purrr::imap_dfr(primary_pairs, function(df, nm) {
  df %>% mutate(comparison = nm)
})

primary_pairs_path <- file.path(sensitivity_dir, "primary_DE_table_correlation_input_long.csv")
readr::write_csv(primary_pairs_long, primary_pairs_path)
message("[STEP6] Saved primary correlation input: ", primary_pairs_path)


## ==========================================================
## 6. Expanded all-expressed/tested gene sensitivity analysis
## ==========================================================

brain_rds_path <- file.path(results_dir, "brain_step3_annotated.rds")

if (!file.exists(brain_rds_path)) {
  stop("[STEP6] Cannot find ", brain_rds_path, ". Run STEP 3 first.")
}

brain <- readRDS(brain_rds_path)

if (!inherits(brain, "Seurat")) {
  stop("[STEP6] brain_step3_annotated.rds is not a Seurat object.")
}

required_metadata <- c("cell_type", "condition")
missing_metadata <- setdiff(required_metadata, colnames(brain@meta.data))

if (length(missing_metadata) > 0) {
  stop(
    "[STEP6] Missing required metadata column(s): ",
    paste(missing_metadata, collapse = ", ")
  )
}

assays_available <- tryCatch(
  as.character(Seurat::Assays(brain)),
  error = function(e) character(0)
)

message("[STEP6] Assays available: ", paste(assays_available, collapse = ", "))

if ("SCT" %in% assays_available) {
  assay_to_use <- "SCT"
} else if ("RNA" %in% assays_available) {
  assay_to_use <- "RNA"
} else {
  assay_to_use <- DefaultAssay(brain)
  warning(
    "[STEP6] Neither SCT nor RNA was found by Assays(). ",
    "Using DefaultAssay(brain) = ", assay_to_use
  )
}

DefaultAssay(brain) <- assay_to_use
message("[STEP6] Assay selected for expanded sensitivity FindMarkers: ", assay_to_use)

if (assay_to_use == "SCT") {
  brain <- tryCatch(
    {
      message("[STEP6] Running PrepSCTFindMarkers()...")
      PrepSCTFindMarkers(brain)
    },
    error = function(e) {
      warning(
        "[STEP6] PrepSCTFindMarkers failed; continuing as-is. Message: ",
        conditionMessage(e)
      )
      brain
    }
  )
} else if (assay_to_use == "RNA") {
  message("[STEP6] Using RNA assay. Checking/creating normalized data layer if needed...")

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
    message("[STEP6] RNA data layer not found. Running NormalizeData()...")
    brain <- NormalizeData(
      brain,
      assay = "RNA",
      normalization.method = "LogNormalize",
      scale.factor = 10000,
      verbose = FALSE
    )
  }
}

run_all_expressed_fc <- function(brain_obj,
                                 cell_type_name,
                                 condition_B,
                                 condition_A = "Mock",
                                 assay = assay_to_use,
                                 min_cells_per_group = analysis_parameters$min_cells_per_group) {

  message(
    "[STEP6] Expanded sensitivity FindMarkers: ",
    cell_type_name, " | ", condition_B, " vs ", condition_A
  )

  cells_ct <- WhichCells(brain_obj, expression = cell_type == cell_type_name)

  if (length(cells_ct) == 0) {
    warning("[STEP6] No cells for ", cell_type_name)
    return(NULL)
  }

  cond_vec <- brain_obj$condition[cells_ct]
  tab_cond <- table(cond_vec)
  print(tab_cond)

  if (!all(c(condition_A, condition_B) %in% names(tab_cond))) {
    warning("[STEP6] Missing condition for ", cell_type_name, ": ", condition_B)
    return(NULL)
  }

  if (any(tab_cond[c(condition_A, condition_B)] < min_cells_per_group)) {
    warning("[STEP6] Not enough cells for ", cell_type_name, ": ", condition_B)
    return(NULL)
  }

  Idents(brain_obj) <- brain_obj$condition

  fm <- tryCatch(
    {
      FindMarkers(
        brain_obj,
        ident.1 = condition_B,
        ident.2 = condition_A,
        assay = assay,
        test.use = analysis_parameters$test_use,
        logfc.threshold = analysis_parameters$expanded_logfc_threshold,
        min.pct = analysis_parameters$expanded_min_pct,
        cells = cells_ct
      )
    },
    error = function(e) {
      warning(
        "[STEP6] Expanded FindMarkers failed for ",
        cell_type_name, " | ", condition_B, " vs ", condition_A,
        ". Message: ", conditionMessage(e)
      )
      NULL
    }
  )

  if (is.null(fm) || nrow(fm) == 0) {
    return(NULL)
  }

  fm %>%
    rownames_to_column("gene") %>%
    transmute(
      gene = gene,
      avg_log2FC = avg_log2FC,
      p_val_adj = p_val_adj,
      cell_type = cell_type_name,
      condition_B = condition_B,
      condition_A = condition_A,
      assay_used = assay,
      gene_universe = "expanded_all_expressed_tested_genes"
    )
}

expanded_files_info <- tibble(
  comparison = c(
    "RG_BR", "RG_FSS", "RG_IFNb",
    "NPC_BR", "NPC_FSS", "NPC_IFNb"
  ),
  cell_type = c(
    "Radial_Glia", "Radial_Glia", "Radial_Glia",
    "Cycling_NPC", "Cycling_NPC", "Cycling_NPC"
  ),
  condition_B = c("ZIKV_BR", "ZIKV_FSS", "IFNb", "ZIKV_BR", "ZIKV_FSS", "IFNb")
)

expanded_fc_long <- purrr::pmap_dfr(
  list(expanded_files_info$cell_type, expanded_files_info$condition_B),
  function(cell_type, condition_B) {
    run_all_expressed_fc(
      brain_obj = brain,
      cell_type_name = cell_type,
      condition_B = condition_B,
      condition_A = "Mock",
      assay = assay_to_use
    )
  }
)


if (nrow(expanded_fc_long) == 0) {
  stop("[STEP6] Expanded sensitivity FindMarkers returned no results.")
}

expanded_fc_path <- file.path(sensitivity_dir, "expanded_all_expressed_tested_FindMarkers_long.csv")
readr::write_csv(expanded_fc_long, expanded_fc_path)
message("[STEP6] Saved expanded all-expressed FindMarkers table: ", expanded_fc_path)

make_expanded_pair <- function(cell_type_name,
                               cell_type_label,
                               zikv_condition,
                               zikv_label) {

  zikv_df <- expanded_fc_long %>%
    filter(cell_type == cell_type_name, condition_B == zikv_condition) %>%
    select(gene, logFC_ZIKV = avg_log2FC)

  ifnb_df <- expanded_fc_long %>%
    filter(cell_type == cell_type_name, condition_B == "IFNb") %>%
    select(gene, logFC_IFNb = avg_log2FC)

  full_join(zikv_df, ifnb_df, by = "gene") %>%
    mutate(
      cell_type = cell_type_label,
      zikv_condition = zikv_label,
      gene_universe = "expanded_all_expressed_tested_genes",
      is_ISG = gene %in% core_isg
    )
}

expanded_pairs <- list(
  RG_BR = make_expanded_pair("Radial_Glia", "Radial glia-like", "ZIKV_BR", "BR"),
  RG_FSS = make_expanded_pair("Radial_Glia", "Radial glia-like", "ZIKV_FSS", "FSS"),
  NPC_BR = make_expanded_pair("Cycling_NPC", "Cycling NPC", "ZIKV_BR", "BR"),
  NPC_FSS = make_expanded_pair("Cycling_NPC", "Cycling NPC", "ZIKV_FSS", "FSS")
)

expanded_corr_summary <- purrr::imap_dfr(expanded_pairs, function(df, nm) {
  summarize_correlations(
    df = df,
    xcol = "logFC_IFNb",
    ycol = "logFC_ZIKV",
    gene_universe_label = "expanded_all_expressed_tested_genes"
  ) %>%
    mutate(
      comparison = nm,
      cell_type = unique(df$cell_type),
      zikv_condition = unique(df$zikv_condition)
    )
})

expanded_pairs_long <- purrr::imap_dfr(expanded_pairs, function(df, nm) {
  df %>% mutate(comparison = nm)
})

expanded_pairs_path <- file.path(sensitivity_dir, "expanded_correlation_input_long.csv")
readr::write_csv(expanded_pairs_long, expanded_pairs_path)
message("[STEP6] Saved expanded correlation input: ", expanded_pairs_path)


## ==========================================================
## 7. Save full Supplementary Table 1
## ==========================================================

correlation_summary <- bind_rows(
  primary_corr_summary,
  expanded_corr_summary
) %>%
  select(
    comparison,
    cell_type,
    zikv_condition,
    gene_universe,
    gene_class,
    n_genes,
    pearson_r,
    spearman_r
  ) %>%
  arrange(cell_type, zikv_condition, gene_universe, gene_class)

supp_table1_path <- file.path(
  results_dir,
  "Supplementary_Table_1_ZIKV_IFNb_correlation_sensitivity.csv"
)

readr::write_csv(correlation_summary, supp_table1_path)
message("[STEP6] Saved Supplementary Table 1: ", supp_table1_path)


## ==========================================================
## 8. Manuscript-style main Figure 5
## ==========================================================

# Panel A: Pearson correlation barplot using primary DE-table genes, all genes.
panelA_df <- primary_corr_summary %>%
  filter(
    gene_universe == "primary_exploratory_DE_tables",
    gene_class == "all_genes_in_universe"
  ) %>%
  mutate(
    strain_bar = case_when(
      comparison %in% c("RG_BR", "NPC_BR") ~ "BR",
      comparison %in% c("RG_FSS", "NPC_FSS") ~ "FSS",
      TRUE ~ zikv_condition
    ),
    cell_type_bar = case_when(
      comparison %in% c("RG_BR", "RG_FSS") ~ "Radial glia",
      comparison %in% c("NPC_BR", "NPC_FSS") ~ "Cycling NPC",
      TRUE ~ cell_type
    ),
    strain_bar = factor(strain_bar, levels = c("BR", "FSS")),
    cell_type_bar = factor(cell_type_bar, levels = c("Cycling NPC", "Radial glia"))
  )

p_A <- ggplot(panelA_df, aes(x = strain_bar, y = pearson_r, fill = cell_type_bar)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.62) +
  coord_cartesian(ylim = c(0, 1.05)) +
  scale_fill_manual(
    values = c("Cycling NPC" = "#B58AD9", "Radial glia" = "#7BC67B")
  ) +
  labs(
    x = "ZIKV strain",
    y = "Pearson correlation\n(ZIKV vs IFNβ)",
    title = "Pearson correlation between ZIKV- and IFNβ-associated responses"
  ) +
  step6_theme(base_size = 12) +
  theme(legend.position = "right")


make_scatter_manuscript <- function(df, panel_title) {
  zikv_label <- unique(df$zikv_condition)
  ycol <- paste0("logFC_", zikv_label)
  
  df_plot <- df %>%
    filter(!is.na(logFC_IFNb), !is.na(.data[[ycol]])) %>%
    mutate(
      is_ISG_label = ifelse(is_ISG, "Canonical ISG", "Other gene")
    )
  
  r_val <- safe_cor(df_plot$logFC_IFNb, df_plot[[ycol]], method = "pearson")
  
  ggplot(df_plot, aes(x = logFC_IFNb, y = .data[[ycol]])) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.35) +
    geom_point(aes(color = is_ISG_label), alpha = 0.75, size = 1.2) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.45) +
    scale_color_manual(
      values = c("Canonical ISG" = "#CC2B3C", "Other gene" = "grey70")
    ) +
    labs(
      x = "IFNβ log2FC",
      y = paste0("ZIKV-", zikv_label, " log2FC"),
      title = paste0(panel_title, " (r = ", sprintf("%.2f", r_val), ")")
    ) +
    step6_theme(base_size = 11) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 11, face = "bold", hjust = 0.5)
    )
}

p_B <- make_scatter_manuscript(
  primary_pairs$RG_BR,
  "Radial glia-like – ZIKV-BR vs IFNβ"
)

p_C <- make_scatter_manuscript(
  primary_pairs$RG_FSS,
  "Radial glia-like – ZIKV-FSS vs IFNβ"
)

p_D <- make_scatter_manuscript(
  primary_pairs$NPC_BR,
  "Cycling NPC – ZIKV-BR vs IFNβ"
)

p_E <- make_scatter_manuscript(
  primary_pairs$NPC_FSS,
  "Cycling NPC – ZIKV-FSS vs IFNβ"
)

scatter_grid <- (p_B | p_C) / (p_D | p_E)

figure5 <- p_A / scatter_grid +
  patchwork::plot_layout(heights = c(0.85, 2.15)) +
  patchwork::plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

save_plot_tiff(
  figure5,
  out_path = file.path(fig_step6_main, "Figure5_ZIKV_IFNb_transcriptional_similarity.tiff"),
  width = 11,
  height = 10
)

save_plot_tiff(
  scatter_grid,
  out_path = file.path(fig_step6_main, "Figure5_scatter_panels_B_to_E.tiff"),
  width = 11,
  height = 7.5
)

save_plot_tiff(
  p_A,
  out_path = file.path(fig_step6_main, "Figure5_panel_A_Pearson_barplot.tiff"),
  width = 6,
  height = 4
)


## ==========================================================
## 9. Correlation sensitivity summary plot
## Supplementary figure only, not main Figure 5
## ==========================================================

corr_plot_df <- correlation_summary %>%
  mutate(
    comparison_label = case_when(
      comparison == "RG_BR" ~ "Radial glia-like\nZIKV-BR",
      comparison == "RG_FSS" ~ "Radial glia-like\nZIKV-FSS",
      comparison == "NPC_BR" ~ "Cycling NPC\nZIKV-BR",
      comparison == "NPC_FSS" ~ "Cycling NPC\nZIKV-FSS",
      TRUE ~ comparison
    ),
    gene_universe = factor(
      gene_universe,
      levels = c(
        "primary_exploratory_DE_tables",
        "expanded_all_expressed_tested_genes"
      ),
      labels = c(
        "Primary DE-table genes",
        "Expanded all-expressed/tested genes"
      )
    ),
    gene_class = factor(
      gene_class,
      levels = c(
        "all_genes_in_universe",
        "canonical_ISGs_only",
        "non_ISG_genes",
        "excluding_canonical_ISGs"
      ),
      labels = c(
        "All genes in universe",
        "Canonical ISGs only",
        "Non-ISG genes",
        "Excluding canonical ISGs"
      )
    )
  )

p_pearson <- ggplot(
  corr_plot_df,
  aes(x = comparison_label, y = pearson_r, fill = gene_class)
) +
  geom_col(position = position_dodge(width = 0.80), width = 0.70) +
  facet_wrap(~ gene_universe, nrow = 1) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    x = NULL,
    y = "Pearson r",
    fill = NULL,
    title = "ZIKV–IFNβ log2FC correlation sensitivity: Pearson"
  ) +
  step6_theme(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "bottom"
  )

p_spearman <- ggplot(
  corr_plot_df,
  aes(x = comparison_label, y = spearman_r, fill = gene_class)
) +
  geom_col(position = position_dodge(width = 0.80), width = 0.70) +
  facet_wrap(~ gene_universe, nrow = 1) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(
    x = NULL,
    y = "Spearman ρ",
    fill = NULL,
    title = "ZIKV–IFNβ log2FC correlation sensitivity: Spearman"
  ) +
  step6_theme(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, face = "bold"),
    legend.position = "bottom"
  )

sensitivity_plot <- p_pearson / p_spearman

save_plot_tiff(
  sensitivity_plot,
  out_path = file.path(fig_step6_supp, "ZIKV_IFNb_correlation_sensitivity_summary.tiff"),
  width = 14,
  height = 10
)


## ==========================================================
## 11. Reproducibility files
## ==========================================================

parameters_path <- file.path(results_dir, "step6_analysis_parameters.txt")
capture.output(
  {
    cat("STEP 6 ZIKV–IFNβ transcriptional similarity sensitivity analysis\n\n")
    cat("Assay used for expanded sensitivity:", assay_to_use, "\n")
    cat("Available assays:", paste(assays_available, collapse = ", "), "\n\n")
    str(analysis_parameters)
  },
  file = parameters_path
)
message("[STEP6] Saved parameters: ", parameters_path)

interpretation_note <- c(
  "STEP 6 interpretation note",
  "",
  "This step compares ZIKV-associated and IFNβ-associated log2FC vectors descriptively.",
  "The fetal neural dataset contains one biological sample per experimental condition.",
  "Cells are not biological replicates.",
  "Primary correlations are calculated using genes present in the STEP 4 exploratory DE tables.",
  "Because filtered DE-table gene universes can enrich for shared ISG induction, correlation values may be inflated.",
  "Sensitivity analyses therefore include Pearson and Spearman correlations across canonical ISGs, non-ISG genes, genes excluding the canonical ISG panel, and an expanded all-expressed/tested gene universe.",
  "These correlations describe transcriptional alignment and IFN/ISG responsiveness.",
  "They do not demonstrate functional equivalence between virus-associated and IFNβ-stimulated programs."
)

interpretation_note_path <- file.path(results_dir, "step6_interpretation_note.txt")
writeLines(interpretation_note, interpretation_note_path)
message("[STEP6] Saved interpretation note: ", interpretation_note_path)

sink(file.path(session_dir, "STEP6_sessionInfo.txt"))
sessionInfo()
sink()
message("[STEP6] Saved session information.")

message("[STEP6] STEP 6 completed successfully.")

############################################################
# End of STEP 6
############################################################
