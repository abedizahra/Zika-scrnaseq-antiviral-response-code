############################################################
# Zika scRNA-seq project – STEP 2
# Quality control, normalization, dimensionality reduction,
# UMAP embedding, and clustering
#
# Purpose:
#   This script takes the STEP 1 Seurat objects, performs QC,
#   normalization, PCA, UMAP, clustering, and saves processed
#   objects plus diagnostic figures.
#
# Reproducibility notes:
#   - No hard-coded setwd() path is used.
#   - Run from the repository root, or from inside the repository.
#   - If STEP 1 RDS files exist, they are loaded from results/.
#   - If STEP 1 RDS files do not exist, scripts/01_prepare_10x_data.R
#     is sourced automatically to create the input objects.
#   - Random seed and analysis parameters are fixed below.
#   - Session information and QC summaries are saved for traceability.
#
# Expected input options:
#   Preferred:
#     results/brain_step1_raw.rds
#     results/dc_step1_raw.rds
#
#   Fallback:
#     scripts/01_prepare_10x_data.R creates objects named `brain` and `dc`
#
# Outputs:
#   results/brain_step2_SCT_qc.rds
#   results/dc_step2_norm_qc.rds
#   results/step2_qc_summary.csv
#   figures/step2/*.tiff
#   figures/main/Figure1_global_transcriptional_structure.tiff
#   figures/supplementary/Supplementary_Figure_S1_QC_and_PCA_diagnostics.tiff
#   session_info/STEP2_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

analysis_parameters <- list(
  seed = 1234,
  brain_min_features = 200,
  brain_max_features = 6000,
  brain_min_counts = 500,
  brain_max_percent_mt = 12,
  brain_sct_vars_to_regress = c("percent.mt", "nCount_RNA"),
  brain_pca_dims = 1:30,
  brain_cluster_resolution = 0.4,
  dc_min_features = 200,
  dc_max_features = 6000,
  dc_min_counts = 500,
  dc_max_percent_mt = 10,
  dc_variable_features = 3000,
  dc_pca_dims = 1:30,
  dc_cluster_resolution = 0.3
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

find_project_root <- function(start_dir = getwd()) {
  # For local checking: use the current R working directory as the project root.
  # This avoids requiring README.md or a scripts/ folder before the repo is fully organized.
  normalizePath(start_dir, winslash = "/", mustWork = TRUE)
}

project_root <- find_project_root()
message("[STEP2] Project root: ", project_root)

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "tibble",
  "patchwork",
  "ragg",
  "future"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP2] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(patchwork)
  library(ragg)
  library(future)
})

## Avoid future memory limit errors during SCTransform.
plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

message("[STEP2] R version: ", R.version.string)
message("[STEP2] Seurat version: ", as.character(utils::packageVersion("Seurat")))


## ==========================================================
## 2. Output directories
## ==========================================================

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
session_dir <- file.path(project_root, "session_info")

fig_step2_dir <- file.path(figures_dir, "step2")
fig_main_dir  <- file.path(figures_dir, "main")
fig_supp_dir  <- file.path(figures_dir, "supplementary")

dir.create(results_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(fig_step2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_main_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp_dir,  recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Load STEP 1 objects
## ==========================================================

brain_step1_rds <- Sys.getenv(
  "BRAIN_STEP1_RDS",
  unset = file.path(results_dir, "brain_step1_raw.rds")
)

dc_step1_rds <- Sys.getenv(
  "DC_STEP1_RDS",
  unset = file.path(results_dir, "dc_step1_raw.rds")
)

step1_script_candidates <- c(
  file.path(project_root, "scripts", "01_prepare_10x_data.R"),
  file.path(project_root, "01_prepare_10x_data.R")
)

step1_script <- step1_script_candidates[file.exists(step1_script_candidates)][1]
if (length(step1_script) == 0 || is.na(step1_script)) {
  step1_script <- file.path(project_root, "scripts", "01_prepare_10x_data.R")
}

load_step1_objects <- function() {
  if (file.exists(brain_step1_rds) && file.exists(dc_step1_rds)) {
    message("[STEP2] Loading STEP 1 RDS files:")
    message("[STEP2]   ", brain_step1_rds)
    message("[STEP2]   ", dc_step1_rds)
    
    brain <- readRDS(brain_step1_rds)
    dc <- readRDS(dc_step1_rds)
    
    return(list(brain = brain, dc = dc))
  }
  
  if (exists("brain", envir = .GlobalEnv) &&
      inherits(get("brain", envir = .GlobalEnv), "Seurat") &&
      exists("dc", envir = .GlobalEnv) &&
      inherits(get("dc", envir = .GlobalEnv), "Seurat")) {
    message("[STEP2] Using existing STEP 1 objects from the current R session.")
    
    return(list(
      brain = get("brain", envir = .GlobalEnv),
      dc = get("dc", envir = .GlobalEnv)
    ))
  }
  
  if (!file.exists(step1_script)) {
    stop(
      "[STEP2] STEP 1 objects were not found and STEP 1 script is missing: ",
      step1_script,
      "\nRun scripts/01_prepare_10x_data.R first, or save STEP 1 objects as:\n",
      "  results/brain_step1_raw.rds\n",
      "  results/dc_step1_raw.rds"
    )
  }
  
  message("[STEP2] STEP 1 RDS files were not found. Sourcing STEP 1 script...")
  source(step1_script, local = .GlobalEnv)
  
  if (!exists("brain", envir = .GlobalEnv) ||
      !inherits(get("brain", envir = .GlobalEnv), "Seurat")) {
    stop("[STEP2] STEP 1 did not create a valid Seurat object named `brain`.")
  }
  
  if (!exists("dc", envir = .GlobalEnv) ||
      !inherits(get("dc", envir = .GlobalEnv), "Seurat")) {
    stop("[STEP2] STEP 1 did not create a valid Seurat object named `dc`.")
  }
  
  brain <- get("brain", envir = .GlobalEnv)
  dc <- get("dc", envir = .GlobalEnv)
  
  saveRDS(brain, brain_step1_rds)
  saveRDS(dc, dc_step1_rds)
  
  message("[STEP2] Saved STEP 1 brain object for future runs: ", brain_step1_rds)
  message("[STEP2] Saved STEP 1 Vero/moDC object for future runs: ", dc_step1_rds)
  
  return(list(brain = brain, dc = dc))
}

step1_objects <- load_step1_objects()
brain <- step1_objects$brain
dc <- step1_objects$dc

if (!inherits(brain, "Seurat")) {
  stop("[STEP2] `brain` is not a Seurat object.")
}

if (!inherits(dc, "Seurat")) {
  stop("[STEP2] `dc` is not a Seurat object.")
}

message("[STEP2] Brain cells before QC: ", ncol(brain))
message("[STEP2] Vero/moDC cells before QC: ", ncol(dc))


## ==========================================================
## 4. Plot theme and save helper
## ==========================================================

theme_journal <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 13, face = "bold", color = "black"),
    axis.title.y = element_text(size = 13, face = "bold", color = "black"),
    axis.text.x = element_text(size = 11, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11, face = "bold"),
    legend.key.size = unit(0.6, "cm"),
    axis.line = element_line(linewidth = 0.8, color = "black")
  )

save_tiff <- function(plot_obj, out_path, width = 8, height = 6, dpi = 600) {
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
  message("[STEP2] Saved: ", out_path)
}


## ==========================================================
## 5. Add sample metadata
## ==========================================================

## ------------------------------
## 5A. Fetal neural dataset
## ------------------------------

message("[STEP2] Adding metadata for the fetal neural dataset...")
message("[STEP2] Brain sample IDs:")
print(table(brain$sample_id))

sample_info_brain <- tibble(
  sample_id    = c("GSM7659279", "GSM7659280", "GSM7659281", "GSM7659282"),
  condition    = c("Mock",       "ZIKV_BR",    "ZIKV_FSS",   "IFNb"),
  virus_strain = c("None",       "Brazil",     "FSS13025",   "None"),
  treatment    = c("None",       "None",       "None",       "IFNb")
)

brain_meta <- brain@meta.data %>%
  rownames_to_column("cell_id") %>%
  left_join(sample_info_brain, by = "sample_id")

if (any(is.na(brain_meta$condition))) {
  missing_samples <- unique(brain_meta$sample_id[is.na(brain_meta$condition)])
  stop(
    "[STEP2] Some brain cells have missing condition labels. ",
    "Check sample_id mapping for: ",
    paste(missing_samples, collapse = ", ")
  )
}

rownames(brain_meta) <- brain_meta$cell_id
brain_meta$cell_id <- NULL
brain@meta.data <- brain_meta

brain$condition <- factor(
  brain$condition,
  levels = c("Mock", "ZIKV_BR", "ZIKV_FSS", "IFNb")
)

brain$virus_strain <- factor(
  brain$virus_strain,
  levels = c("None", "Brazil", "FSS13025")
)

brain$treatment <- factor(
  brain$treatment,
  levels = c("None", "IFNb")
)

message("[STEP2] Brain condition table:")
print(table(brain$condition))


## ------------------------------
## 5B. Vero/moDC contextual dataset
## ------------------------------

message("[STEP2] Adding metadata for the Vero/moDC contextual dataset...")
message("[STEP2] Vero/moDC sample IDs:")
print(table(dc$sample_id))

dc_meta <- dc@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    cell_system = case_when(
      grepl("ZIKV-Vero-cells", sample_id) ~ "Vero",
      grepl("p22086", sample_id)          ~ "moDC_panel",
      TRUE                                ~ "Unknown"
    ),
    condition_group = case_when(
      grepl("ZIKV-Vero-cells", sample_id) ~ "ZIKV_Vero_cells",
      grepl("p22086", sample_id)          ~ "moDC_panel",
      TRUE                                ~ "Unknown"
    )
  )

rownames(dc_meta) <- dc_meta$cell_id
dc_meta$cell_id <- NULL
dc@meta.data <- dc_meta

dc$cell_system <- factor(
  dc$cell_system,
  levels = c("Vero", "moDC_panel", "Unknown")
)

dc$condition_group <- factor(
  dc$condition_group,
  levels = c("ZIKV_Vero_cells", "moDC_panel", "Unknown")
)

message("[STEP2] Vero/moDC condition group table:")
print(table(dc$condition_group))


## ==========================================================
## 6. Quality control
## ==========================================================

qc_summary <- tibble(
  dataset = character(),
  stage = character(),
  cells = integer()
)

qc_summary <- bind_rows(
  qc_summary,
  tibble(dataset = "brain", stage = "before_qc", cells = ncol(brain)),
  tibble(dataset = "dc", stage = "before_qc", cells = ncol(dc))
)

## ------------------------------
## 6A. Fetal neural dataset QC
## ------------------------------

message("[STEP2] Running QC for the fetal neural dataset...")

DefaultAssay(brain) <- "RNA"
brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")

p_brain_qc_raw <- VlnPlot(
  brain,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("Fetal neural dataset QC before filtering") +
  theme_journal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_tiff(
  p_brain_qc_raw,
  file.path(fig_step2_dir, "brain_QC_violin_raw.tiff"),
  width = 12,
  height = 4
)

brain <- subset(
  brain,
  subset = nFeature_RNA > analysis_parameters$brain_min_features &
    nFeature_RNA < analysis_parameters$brain_max_features &
    nCount_RNA > analysis_parameters$brain_min_counts &
    percent.mt < analysis_parameters$brain_max_percent_mt
)

message("[STEP2] Brain cells after QC: ", ncol(brain))
message("[STEP2] Brain condition table after QC:")
print(table(brain$condition))

qc_summary <- bind_rows(
  qc_summary,
  tibble(dataset = "brain", stage = "after_qc", cells = ncol(brain))
)

p_brain_qc_post <- VlnPlot(
  brain,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("Fetal neural dataset QC after filtering") +
  theme_journal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_tiff(
  p_brain_qc_post,
  file.path(fig_step2_dir, "brain_QC_violin_post.tiff"),
  width = 12,
  height = 4
)


## ------------------------------
## 6B. Vero/moDC dataset QC
## ------------------------------

message("[STEP2] Running QC for the Vero/moDC contextual dataset...")

DefaultAssay(dc) <- "RNA"
dc[["percent.mt"]] <- PercentageFeatureSet(dc, pattern = "^MT-")

p_dc_qc_raw <- VlnPlot(
  dc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition_group",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("Vero/moDC dataset QC before filtering") +
  theme_journal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_tiff(
  p_dc_qc_raw,
  file.path(fig_step2_dir, "dc_QC_violin_raw.tiff"),
  width = 12,
  height = 4
)

dc <- subset(
  dc,
  subset = nFeature_RNA > analysis_parameters$dc_min_features &
    nFeature_RNA < analysis_parameters$dc_max_features &
    nCount_RNA > analysis_parameters$dc_min_counts &
    percent.mt < analysis_parameters$dc_max_percent_mt
)

message("[STEP2] Vero/moDC cells after QC: ", ncol(dc))
message("[STEP2] Vero/moDC condition group table after QC:")
print(table(dc$condition_group))

qc_summary <- bind_rows(
  qc_summary,
  tibble(dataset = "dc", stage = "after_qc", cells = ncol(dc))
)

p_dc_qc_post <- VlnPlot(
  dc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "condition_group",
  ncol = 3,
  pt.size = 0
) +
  ggtitle("Vero/moDC dataset QC after filtering") +
  theme_journal +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_tiff(
  p_dc_qc_post,
  file.path(fig_step2_dir, "dc_QC_violin_post.tiff"),
  width = 12,
  height = 4
)


## ==========================================================
## 7. Normalization, PCA, UMAP, and clustering
## ==========================================================

## ------------------------------
## 7A. Fetal neural dataset
## ------------------------------

message("[STEP2] Running SCTransform, PCA, UMAP, and clustering for the fetal neural dataset...")

DefaultAssay(brain) <- "RNA"

brain <- SCTransform(
  brain,
  vars.to.regress = analysis_parameters$brain_sct_vars_to_regress,
  verbose = FALSE
)

brain <- RunPCA(
  brain,
  assay = "SCT",
  verbose = FALSE
)

p_brain_elbow <- ElbowPlot(brain, ndims = 50) +
  ggtitle("Fetal neural dataset PCA elbow plot") +
  theme_journal

save_tiff(
  p_brain_elbow,
  file.path(fig_step2_dir, "brain_PCA_elbow.tiff"),
  width = 6,
  height = 5
)

brain <- RunUMAP(
  brain,
  dims = analysis_parameters$brain_pca_dims,
  reduction = "pca",
  seed.use = analysis_parameters$seed,
  verbose = FALSE
)

brain <- FindNeighbors(
  brain,
  dims = analysis_parameters$brain_pca_dims,
  verbose = FALSE
)

brain <- FindClusters(
  brain,
  resolution = analysis_parameters$brain_cluster_resolution,
  verbose = FALSE
)

message("[STEP2] Brain clusters:")
print(table(brain$seurat_clusters))


## ------------------------------
## 7B. Vero/moDC contextual dataset
## ------------------------------

message("[STEP2] Running log-normalization, PCA, UMAP, and clustering for the Vero/moDC dataset...")

DefaultAssay(dc) <- "RNA"

dc <- NormalizeData(
  dc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

dc <- FindVariableFeatures(
  dc,
  selection.method = "vst",
  nfeatures = analysis_parameters$dc_variable_features,
  verbose = FALSE
)

dc <- ScaleData(
  dc,
  verbose = FALSE
)

dc <- RunPCA(
  dc,
  npcs = 50,
  verbose = FALSE
)

p_dc_elbow <- ElbowPlot(dc, ndims = 50) +
  ggtitle("Vero/moDC dataset PCA elbow plot") +
  theme_journal

save_tiff(
  p_dc_elbow,
  file.path(fig_step2_dir, "dc_PCA_elbow.tiff"),
  width = 6,
  height = 5
)

dc <- RunUMAP(
  dc,
  dims = analysis_parameters$dc_pca_dims,
  reduction = "pca",
  seed.use = analysis_parameters$seed,
  verbose = FALSE
)

dc <- FindNeighbors(
  dc,
  dims = analysis_parameters$dc_pca_dims,
  verbose = FALSE
)

dc <- FindClusters(
  dc,
  resolution = analysis_parameters$dc_cluster_resolution,
  verbose = FALSE
)

message("[STEP2] Vero/moDC clusters:")
print(table(dc$seurat_clusters))


## ==========================================================
## 8. UMAP plots
## ==========================================================

## ------------------------------
## 8A. Fetal neural UMAP plots
## ------------------------------

p_brain_umap_condition <- DimPlot(
  brain,
  reduction = "umap",
  group.by = "condition",
  pt.size = 0.2,
  shuffle = TRUE
) +
  ggtitle("Fetal neural dataset by exposure condition") +
  theme_journal

p_brain_umap_cluster <- DimPlot(
  brain,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2,
  shuffle = TRUE
) +
  ggtitle("Fetal neural dataset by Seurat cluster") +
  theme_journal +
  NoLegend()

save_tiff(
  p_brain_umap_condition,
  file.path(fig_step2_dir, "brain_UMAP_by_condition.tiff"),
  width = 7,
  height = 6
)

save_tiff(
  p_brain_umap_cluster,
  file.path(fig_step2_dir, "brain_UMAP_by_cluster.tiff"),
  width = 7,
  height = 6
)


## ------------------------------
## 8B. Vero/moDC UMAP plots
## ------------------------------

p_dc_umap_condition <- DimPlot(
  dc,
  reduction = "umap",
  group.by = "condition_group",
  pt.size = 0.2,
  shuffle = TRUE
) +
  ggtitle("Vero/moDC dataset by condition group") +
  theme_journal

p_dc_umap_cluster <- DimPlot(
  dc,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.2,
  shuffle = TRUE
) +
  ggtitle("Vero/moDC dataset by Seurat cluster") +
  theme_journal +
  NoLegend()

save_tiff(
  p_dc_umap_condition,
  file.path(fig_step2_dir, "dc_UMAP_by_condition_group.tiff"),
  width = 7,
  height = 6
)

save_tiff(
  p_dc_umap_cluster,
  file.path(fig_step2_dir, "dc_UMAP_by_cluster.tiff"),
  width = 7,
  height = 6
)


## ==========================================================
## 9. Assemble main and supplementary figures
## ==========================================================

fig1_global_structure <-
  (p_brain_umap_condition | p_brain_umap_cluster) /
  (p_dc_umap_condition    | p_dc_umap_cluster) +
  plot_annotation(tag_levels = "A")

save_tiff(
  fig1_global_structure,
  file.path(fig_main_dir, "Figure1_global_transcriptional_structure.tiff"),
  width = 13,
  height = 11
)

supp_fig_s1 <-
  (p_brain_qc_raw / p_brain_qc_post / p_dc_qc_raw / p_dc_qc_post /
     (p_brain_elbow | p_dc_elbow)) +
  plot_annotation(tag_levels = "A")

save_tiff(
  supp_fig_s1,
  file.path(fig_supp_dir, "Supplementary_Figure_S1_QC_and_PCA_diagnostics.tiff"),
  width = 14,
  height = 22
)


## ==========================================================
## 10. Save processed objects and reproducibility files
## ==========================================================

brain_rds_path <- file.path(results_dir, "brain_step2_SCT_qc.rds")
dc_rds_path    <- file.path(results_dir, "dc_step2_norm_qc.rds")
qc_summary_path <- file.path(results_dir, "step2_qc_summary.csv")
parameters_path <- file.path(results_dir, "step2_analysis_parameters.txt")

saveRDS(brain, brain_rds_path)
saveRDS(dc, dc_rds_path)

write.csv(qc_summary, qc_summary_path, row.names = FALSE)

capture.output(
  str(analysis_parameters),
  file = parameters_path
)

message("[STEP2] Saved brain object: ", brain_rds_path)
message("[STEP2] Saved Vero/moDC object: ", dc_rds_path)
message("[STEP2] Saved QC summary: ", qc_summary_path)
message("[STEP2] Saved analysis parameters: ", parameters_path)


## ==========================================================
## 11. Save session information
## ==========================================================

sink(file.path(session_dir, "STEP2_sessionInfo.txt"))
sessionInfo()
sink()

message("[STEP2] Saved session information.")
message("[STEP2] STEP 2 completed successfully.")

############################################################
# End of STEP 2
############################################################
