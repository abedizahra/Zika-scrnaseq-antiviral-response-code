############################################################
# Zika scRNA-seq project – STEP 3
# Cell-type annotation and pseudobulk preparation
#
# Purpose:
#   - Annotate fetal neural scRNA-seq clusters using canonical markers
#   - Explicitly document the radial glia annotation limitation specified
#     by analysiss: this cluster is treated as a progenitor-enriched
#     radial glia-like population, not as resolved aRG/oRG/IP subtypes
#   - Summarize cell-type composition descriptively
#   - Construct pseudobulk count profiles by sample_id and cell_type
#
# Analysis notes:
#   - No hard-coded setwd() path is used.
#   - The current working directory is treated as the project root.
#   - This script expects STEP 2 output:
#       results/brain_step2_SCT_qc.rds
#   - Pseudobulk profiles are constructed for descriptive aggregation only.
#     Pseudobulk differential expression is not performed here because the
#     primary fetal neural dataset has one biological sample per condition.
#   - Cell numbers are not treated as biological replicates.
#   - UMAP coordinates are dimensionless embedding units.
############################################################


## ==========================================================
## 0. Project root, parameters, and packages
## ==========================================================

set.seed(1234)

project_root <- getwd()
message("[STEP3] Project root: ", project_root)

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "tibble",
  "tidyr",
  "readr",
  "scales",
  "patchwork"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP3] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running STEP 3."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
})

has_ragg <- requireNamespace("ragg", quietly = TRUE)

message("[STEP3] R version: ", R.version.string)
message("[STEP3] Seurat version: ", as.character(utils::packageVersion("Seurat")))


## ==========================================================
## 1. Directories
## ==========================================================

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
tables_dir  <- file.path(project_root, "tables")
session_dir <- file.path(project_root, "session_info")

fig_step3_dir <- file.path(figures_dir, "step3")
fig_main_dir  <- file.path(figures_dir, "main")
fig_supp_dir  <- file.path(figures_dir, "supplementary")

dir.create(results_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step3_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_main_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir,  recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 2. Plot theme and save helper
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
    legend.key.size = unit(0.6, "cm"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

save_tiff <- function(plot_obj, out_path, width = 8, height = 6, dpi = 600) {
  if (has_ragg) {
    ggsave(
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
    ggsave(
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
  message("[STEP3] Saved: ", out_path)
}


## ==========================================================
## 3. Load STEP 2 fetal neural object
## ==========================================================

brain_rds_path <- file.path(results_dir, "brain_step2_SCT_qc.rds")

if (!file.exists(brain_rds_path)) {
  stop(
    "[STEP3] Cannot find STEP 2 output: ", brain_rds_path,
    "\nRun STEP 2 first."
  )
}

brain <- readRDS(brain_rds_path)

if (!inherits(brain, "Seurat")) {
  stop("[STEP3] Loaded object is not a Seurat object: ", brain_rds_path)
}

required_meta <- c("sample_id", "condition", "virus_strain", "treatment", "seurat_clusters")
missing_meta <- setdiff(required_meta, colnames(brain@meta.data))

if (length(missing_meta) > 0) {
  stop(
    "[STEP3] Missing required metadata column(s) in brain object: ",
    paste(missing_meta, collapse = ", "),
    "\nCheck STEP 2 output."
  )
}

message("[STEP3] Brain cells after STEP 2 QC: ", ncol(brain))
message("[STEP3] Conditions:")
print(table(brain$condition))
message("[STEP3] Clusters:")
print(table(brain$seurat_clusters))


## ==========================================================
## 4. Marker definitions and annotation notes
## ==========================================================

# Primary markers used to assign broad cell types.
# Annotation note:
#   SOX2/PAX6/NES are not sufficient to distinguish apical radial glia,
#   outer/basal radial glia, and intermediate progenitors. Therefore, the
#   radial glia assignment is interpreted as a progenitor-enriched
#   radial glia-like population, not a fully resolved radial glia subtype.

celltype_marker_table <- tibble::tribble(
  ~marker_group, ~cell_type_short, ~cell_type_label, ~gene, ~interpretation_note,
  "Primary annotation", "Radial_Glia", "Progenitor-enriched radial glia-like", "SOX2", "Shared progenitor marker; does not resolve aRG/oRG/IP subtypes",
  "Primary annotation", "Radial_Glia", "Progenitor-enriched radial glia-like", "PAX6", "Shared progenitor marker; does not resolve aRG/oRG/IP subtypes",
  "Primary annotation", "Radial_Glia", "Progenitor-enriched radial glia-like", "NES",  "Shared progenitor marker; does not resolve aRG/oRG/IP subtypes",
  "Primary annotation", "Cycling_NPC", "Cycling NPC", "MKI67", "Cycling/proliferative neural progenitor marker",
  "Primary annotation", "Cycling_NPC", "Cycling NPC", "TOP2A", "Cycling/proliferative neural progenitor marker",
  "Primary annotation", "Cycling_NPC", "Cycling NPC", "HMGB2", "Cycling/proliferative neural progenitor marker",
  "Primary annotation", "Excitatory_Neuron", "Excitatory neuron", "SLC17A7", "Excitatory neuronal marker",
  "Primary annotation", "Excitatory_Neuron", "Excitatory neuron", "NEUROD6", "Excitatory neuronal marker",
  "Primary annotation", "Inhibitory_Neuron", "Inhibitory neuron", "GAD1", "Inhibitory neuronal marker",
  "Primary annotation", "Inhibitory_Neuron", "Inhibitory neuron", "GAD2", "Inhibitory neuronal marker",
  "Primary annotation", "Inhibitory_Neuron", "Inhibitory neuron", "DLX1", "Inhibitory lineage marker",
  "Primary annotation", "Astrocyte", "Astrocyte", "GFAP", "Astrocyte-associated marker",
  "Primary annotation", "Astrocyte", "Astrocyte", "AQP4", "Astrocyte-associated marker",
  "Primary annotation", "Microglia", "Microglia", "C1QA", "Microglial marker",
  "Primary annotation", "Microglia", "Microglia", "C1QB", "Microglial marker",
  "Primary annotation", "Microglia", "Microglia", "TYROBP", "Microglial marker",
  "Primary annotation", "OPC", "OPC", "PDGFRA", "Oligodendrocyte progenitor marker",
  "Primary annotation", "OPC", "OPC", "OLIG1", "Oligodendrocyte lineage marker",
  "Primary annotation", "OPC", "OPC", "OLIG2", "Oligodendrocyte lineage marker",
  "Subtype context only", "aRG_context", "Apical radial glia context marker", "NOTCH1", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "aRG_context", "Apical radial glia context marker", "HES1", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "aRG_context", "Apical radial glia context marker", "HES5", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "oRG_context", "Outer/basal radial glia context marker", "HOPX", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "oRG_context", "Outer/basal radial glia context marker", "PTPRZ1", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "oRG_context", "Outer/basal radial glia context marker", "FAM107A", "Shown only to document unresolved progenitor subtype context",
  "Subtype context only", "IP_context", "Intermediate progenitor context marker", "EOMES", "Shown only to document unresolved progenitor subtype context"
)

marker_table_path <- file.path(tables_dir, "Supplementary_Table_3_celltype_markers.csv")
readr::write_csv(celltype_marker_table, marker_table_path)
message("[STEP3] Saved marker table: ", marker_table_path)

DefaultAssay(brain) <- if ("SCT" %in% names(Assays(brain))) "SCT" else "RNA"
Idents(brain) <- "seurat_clusters"

marker_vec_all <- unique(celltype_marker_table$gene)
marker_vec_present <- marker_vec_all[marker_vec_all %in% rownames(brain)]
marker_vec_missing <- setdiff(marker_vec_all, marker_vec_present)

missing_marker_path <- file.path(results_dir, "step3_missing_marker_genes.csv")
readr::write_csv(
  tibble(gene = marker_vec_missing),
  missing_marker_path
)

if (length(marker_vec_missing) > 0) {
  warning(
    "[STEP3] Some marker genes were not found in the object and will be omitted from DotPlot: ",
    paste(marker_vec_missing, collapse = ", ")
  )
}

if (length(marker_vec_present) == 0) {
  stop("[STEP3] None of the specified marker genes were found in the object.")
}


## ==========================================================
## 5. DotPlot for marker-based annotation
## ==========================================================

p_brain_markers <- DotPlot(
  brain,
  features = marker_vec_present,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  RotatedAxis() +
  ggtitle("Fetal neural marker expression by Seurat cluster") +
  scale_color_gradientn(colors = c("grey90", "#4B9CD3", "#08306B")) +
  theme_journal +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold")
  )

save_tiff(
  p_brain_markers,
  file.path(fig_step3_dir, "brain_DotPlot_markers_by_cluster.tiff"),
  width = 13,
  height = 6
)


## ==========================================================
## 6. Cluster-to-cell-type mapping
## ==========================================================

# This mapping is intentionally explicit and saved for reproducibility.
# If the marker DotPlot indicates a different mapping for a rerun or updated
# object, edit this table and rerun the script.

cluster_to_celltype <- tibble(
  seurat_clusters = as.character(0:15),
  cell_type = c(
    "Radial_Glia",        # 0
    "Radial_Glia",        # 1
    "Cycling_NPC",        # 2
    "Cycling_NPC",        # 3
    "Excitatory_Neuron",  # 4
    "Excitatory_Neuron",  # 5
    "Excitatory_Neuron",  # 6
    "Inhibitory_Neuron",  # 7
    "Inhibitory_Neuron",  # 8
    "Astrocyte",          # 9
    "Inhibitory_Neuron",  # 10
    "Microglia",          # 11
    "Astrocyte",          # 12
    "Microglia",          # 13
    "OPC",                # 14
    "Unknown"             # 15
  )
) %>%
  mutate(
    cell_type_label = case_when(
      cell_type == "Radial_Glia" ~ "Progenitor-enriched radial glia-like",
      cell_type == "Cycling_NPC" ~ "Cycling NPC",
      cell_type == "Excitatory_Neuron" ~ "Excitatory neuron",
      cell_type == "Inhibitory_Neuron" ~ "Inhibitory neuron",
      cell_type == "Astrocyte" ~ "Astrocyte",
      cell_type == "Microglia" ~ "Microglia",
      cell_type == "OPC" ~ "OPC",
      TRUE ~ "Unknown"
    ),
    annotation_confidence = case_when(
      cell_type == "Radial_Glia" ~ "Broad progenitor annotation; aRG/oRG/IP subtypes not resolved",
      cell_type == "Unknown" ~ "Unassigned",
      TRUE ~ "Broad canonical marker-based annotation"
    )
  )

cluster_map_path <- file.path(results_dir, "step3_cluster_to_celltype_mapping.csv")
readr::write_csv(cluster_to_celltype, cluster_map_path)
message("[STEP3] Saved cluster mapping: ", cluster_map_path)

brain_meta <- brain@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(seurat_clusters = as.character(seurat_clusters)) %>%
  left_join(cluster_to_celltype, by = "seurat_clusters")

brain_meta$cell_type[is.na(brain_meta$cell_type)] <- "Unknown"
brain_meta$cell_type_label[is.na(brain_meta$cell_type_label)] <- "Unknown"
brain_meta$annotation_confidence[is.na(brain_meta$annotation_confidence)] <- "Unassigned"

rownames(brain_meta) <- brain_meta$cell_id
brain_meta$cell_id <- NULL
brain@meta.data <- brain_meta

cell_type_order <- c(
  "Radial_Glia",
  "Cycling_NPC",
  "Excitatory_Neuron",
  "Inhibitory_Neuron",
  "Astrocyte",
  "Microglia",
  "OPC",
  "Unknown"
)

brain$cell_type <- factor(brain$cell_type, levels = cell_type_order)

message("[STEP3] Cell-type composition:")
print(table(brain$cell_type))


## ==========================================================
## 7. Cell-type composition table and figures
## ==========================================================

celltype_cols <- c(
  "Radial_Glia"        = "#E64B35",
  "Cycling_NPC"        = "#4DBBD5",
  "Excitatory_Neuron"  = "#00A087",
  "Inhibitory_Neuron"  = "#3C5488",
  "Astrocyte"          = "#F39B7F",
  "Microglia"          = "#8491B4",
  "OPC"                = "#91D1C2",
  "Unknown"            = "#7F7F7F"
)

celltype_label_map <- c(
  "Radial_Glia"        = "Radial glia-like",
  "Cycling_NPC"        = "Cycling NPC",
  "Excitatory_Neuron"  = "Excitatory neuron",
  "Inhibitory_Neuron"  = "Inhibitory neuron",
  "Astrocyte"          = "Astrocyte",
  "Microglia"          = "Microglia",
  "OPC"                = "OPC",
  "Unknown"            = "Unknown"
)


celltype_comp <- brain@meta.data %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  mutate(
    condition = as.character(condition),
    cell_type = as.character(cell_type)
  ) %>%
  dplyr::count(condition, cell_type, name = "n_cells") %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(fraction_of_cells = n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cell_type = factor(cell_type, levels = cell_type_order),
    cell_type_label = celltype_label_map[as.character(cell_type)]
  )

celltype_comp_path <- file.path(results_dir, "step3_celltype_composition_by_condition.csv")
readr::write_csv(celltype_comp, celltype_comp_path)
message("[STEP3] Saved cell-type composition table: ", celltype_comp_path)

p_brain_umap_celltype <- DimPlot(
  brain,
  reduction = "umap",
  group.by = "cell_type",
  label = FALSE,
  pt.size = 0.2,
  shuffle = TRUE
) +
  scale_color_manual(
    values = celltype_cols,
    labels = celltype_label_map,
    drop = FALSE
  ) +
  ggtitle("Fetal neural UMAP colored by annotated cell type") +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cell type"
  ) +
  theme_journal +
  theme(legend.position = "right")

save_tiff(
  p_brain_umap_celltype,
  file.path(fig_step3_dir, "brain_UMAP_by_cell_type.tiff"),
  width = 7,
  height = 6
)

p_compbar <- ggplot(celltype_comp, aes(x = condition, y = fraction_of_cells, fill = cell_type)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.25) +
  scale_fill_manual(
    values = celltype_cols,
    labels = celltype_label_map,
    drop = FALSE
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Descriptive cell-type composition by condition",
    x = "Condition",
    y = "Fraction of cells",
    fill = "Cell type"
  ) +
  theme_journal +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    legend.position = "right"
  )

save_tiff(
  p_compbar,
  file.path(fig_step3_dir, "brain_celltype_composition_by_condition.tiff"),
  width = 8,
  height = 6
)


## ==========================================================
## 8. Save annotated object
## ==========================================================

brain_annotated_path <- file.path(results_dir, "brain_step3_annotated.rds")
saveRDS(brain, brain_annotated_path)
message("[STEP3] Saved annotated brain object: ", brain_annotated_path)


## ==========================================================
## 9. Pseudobulk construction
## ==========================================================

message("[STEP3] Constructing pseudobulk count profiles by sample_id and cell_type...")

DefaultAssay(brain) <- "RNA"

brain$group_id <- paste(brain$sample_id, brain$cell_type, sep = "__")

# Use counts. Seurat versions differ in whether AggregateExpression uses slot
# or layer, so try the v5-style layer argument first and fall back if needed.
agg <- tryCatch(
  {
    AggregateExpression(
      brain,
      assays = "RNA",
      layer = "counts",
      group.by = "group_id",
      return.seurat = FALSE
    )
  },
  error = function(e) {
    message("[STEP3] AggregateExpression(layer='counts') failed; retrying with slot='counts'.")
    AggregateExpression(
      brain,
      assays = "RNA",
      slot = "counts",
      group.by = "group_id",
      return.seurat = FALSE
    )
  }
)

pseudobulk_counts <- agg$RNA

pb_meta <- tibble(group_id = colnames(pseudobulk_counts)) %>%
  tidyr::separate(group_id, into = c("sample_id", "cell_type"), sep = "__", remove = FALSE) %>%
  left_join(
    brain@meta.data %>%
      as.data.frame() %>%
      select(sample_id, condition, virus_strain, treatment) %>%
      distinct(),
    by = "sample_id"
  ) %>%
  mutate(
    cell_type_label = celltype_label_map[cell_type],
    analysis_note = "Pseudobulk profiles are descriptive aggregation units; pseudobulk DE is not performed because each condition has one biological sample."
  )

pb_counts_path <- file.path(results_dir, "brain_pseudobulk_counts.rds")
pb_meta_path   <- file.path(results_dir, "brain_pseudobulk_metadata.csv")
pb_note_path   <- file.path(results_dir, "step3_pseudobulk_analysis_note.txt")

saveRDS(pseudobulk_counts, pb_counts_path)
readr::write_csv(pb_meta, pb_meta_path)

writeLines(
  c(
    "STEP 3 pseudobulk note",
    "Pseudobulk count profiles were generated by aggregating raw RNA counts by sample_id and annotated cell_type.",
    "Because the fetal neural dataset contains one biological sample per experimental condition, pseudobulk differential expression testing is not performed.",
    "Cell numbers are not treated as biological replicates.",
    "These profiles are used only for descriptive summaries and visualization support."
  ),
  con = pb_note_path
)

message("[STEP3] Pseudobulk matrix dimensions: ", paste(dim(pseudobulk_counts), collapse = " x "))
message("[STEP3] Saved pseudobulk counts: ", pb_counts_path)
message("[STEP3] Saved pseudobulk metadata: ", pb_meta_path)


## ==========================================================
## 10. Assemble Figure 2
## ==========================================================

fig2_plot <- (
  p_brain_markers + theme(plot.margin = margin(t = 10, r = 10, b = 20, l = 10))
) / (
  p_brain_umap_celltype + theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10)) |
    p_compbar + theme(plot.margin = margin(t = 10, r = 10, b = 30, l = 10))
) +
  plot_annotation(tag_levels = "A")

fig2_out_path <- file.path(fig_main_dir, "Figure2_fetal_neural_celltype_architecture.tiff")

save_tiff(
  fig2_plot,
  fig2_out_path,
  width = 13,
  height = 11
)

# Keep a backwards-compatible filename for older manuscript drafts.
fig2_legacy_path <- file.path(fig_main_dir, "Fig2_fetal_brain_celltypes_patchwork.tiff")
save_tiff(
  fig2_plot,
  fig2_legacy_path,
  width = 13,
  height = 11
)


## ==========================================================
## 11. Session information
## ==========================================================

sink(file.path(session_dir, "STEP3_sessionInfo.txt"))
sessionInfo()
sink()

message("[STEP3] Saved session information.")
message("[STEP3] STEP 3 completed successfully.")

############################################################
# End of STEP 3
############################################################
