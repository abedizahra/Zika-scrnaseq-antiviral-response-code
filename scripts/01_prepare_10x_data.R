############################################################
# Zika scRNA-seq project – STEP 1
# Prepare public 10X single-cell RNA-seq datasets
#
# Datasets:
#   - GSE238140: human fetal brain-derived neural tissue
#                exposed ex vivo to Mock, ZIKV-BR, ZIKV-FSS,
#                or IFNβ
#   - GSE230571: ZIKV-exposed Vero cells and primary human
#                monocyte-derived dendritic cells (moDCs)
#
# Purpose:
#   This script prepares GEO supplementary files into standard
#   10X-style sample directories and loads them into Seurat objects.
#
# Reproducibility notes:
#   - No hard-coded setwd() path is used.
#   - Run this script from the repository root.
#   - All paths are relative to the project root.
#   - Large raw GEO archives are expected under:
#       data/GSE238140/
#       data/GSE230571/
#   - Multi-modal 10X outputs are handled by retaining only
#     the Gene Expression matrix.
#   - Seurat v5 layers are collapsed using JoinLayers().
############################################################


## ==========================================================
## 0. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP1] Project root: ", project_root)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
})

message("[STEP1] R version: ", R.version.string)
message("[STEP1] Seurat version: ", as.character(utils::packageVersion("Seurat")))


## ==========================================================
## 1. Helper function: locate a 10X matrix directory
## ==========================================================

find_10x_matrix_dir <- function(sample_root) {
  candidate_dirs <- c(
    sample_root,
    list.dirs(sample_root, recursive = TRUE, full.names = TRUE)
  )
  
  has_matrix <- file.exists(file.path(candidate_dirs, "matrix.mtx")) |
    file.exists(file.path(candidate_dirs, "matrix.mtx.gz"))
  
  candidate_dirs <- candidate_dirs[has_matrix]
  
  if (length(candidate_dirs) == 0) {
    return(NA_character_)
  }
  
  filtered_idx <- grepl("filtered_feature", candidate_dirs, ignore.case = TRUE)
  
  if (any(filtered_idx)) {
    candidate_dirs <- candidate_dirs[filtered_idx]
  }
  
  return(candidate_dirs[1])
}


## ==========================================================
## 2. Helper function: load 10X samples into one Seurat object
## ==========================================================

load_sc10x_dataset <- function(root_path,
                               project_name,
                               min.cells = 3,
                               min.features = 200) {
  message("[", project_name, "] Scanning root path: ", root_path)
  
  sample_dirs <- list.dirs(root_path, recursive = FALSE, full.names = TRUE)
  
  if (length(sample_dirs) == 0) {
    stop("[", project_name, "] No sample directories found under: ", root_path)
  }
  
  sample_list <- list()
  
  for (sample_dir in sample_dirs) {
    data_dir <- find_10x_matrix_dir(sample_dir)
    
    if (is.na(data_dir)) {
      warning(
        "[", project_name, "] No matrix.mtx(.gz) found under sample directory: ",
        sample_dir,
        " -- skipping."
      )
      next
    }
    
    message(
      "[", project_name, "] Reading sample ",
      basename(sample_dir),
      " from ",
      data_dir
    )
    
    counts_raw <- Read10X(data.dir = data_dir)
    
    # Some 10X outputs contain multiple feature modalities, such as
    # Gene Expression and Antibody Capture. For this transcriptomic
    # analysis, only the Gene Expression matrix is retained.
    if (is.list(counts_raw)) {
      message(
        "[", project_name, "] Multi-modal 10X data detected: ",
        paste(names(counts_raw), collapse = ", "),
        ". Retaining Gene Expression matrix."
      )
      
      if ("Gene Expression" %in% names(counts_raw)) {
        counts <- counts_raw[["Gene Expression"]]
      } else {
        warning(
          "[", project_name,
          "] Gene Expression matrix not explicitly found; using the first matrix."
        )
        counts <- counts_raw[[1]]
      }
    } else {
      counts <- counts_raw
    }
    
    sample_id <- basename(sample_dir)
    
    seurat_obj <- CreateSeuratObject(
      counts = counts,
      project = project_name,
      min.cells = min.cells,
      min.features = min.features
    )
    
    seurat_obj$sample_id <- sample_id
    sample_list[[sample_id]] <- seurat_obj
  }
  
  if (length(sample_list) == 0) {
    stop("[", project_name, "] No valid 10X samples were loaded.")
  }
  
  if (length(sample_list) == 1) {
    combined <- sample_list[[1]]
  } else {
    combined <- merge(
      x = sample_list[[1]],
      y = sample_list[2:length(sample_list)],
      add.cell.ids = names(sample_list),
      project = project_name
    )
  }
  
  combined <- JoinLayers(combined)
  
  message(
    "[", project_name, "] Loaded ",
    ncol(combined),
    " cells from ",
    length(sample_list),
    " samples."
  )
  
  return(combined)
}


## ==========================================================
## 3. Prepare GSE238140 fetal neural dataset
## ==========================================================

# GSE238140 is the primary fetal neural single-cell RNA-seq dataset.
# It contains human fetal brain-derived neural tissue exposed ex vivo
# to Mock, ZIKV-BR, ZIKV-FSS, or IFNβ.

brain_dir <- file.path(project_root, "data", "GSE238140")

if (!dir.exists(brain_dir)) {
  stop("[GSE238140] Data directory not found: ", brain_dir)
}

level1_brain_dir <- file.path(brain_dir, "level1")

if (!dir.exists(level1_brain_dir)) {
  dir.create(level1_brain_dir, showWarnings = FALSE, recursive = TRUE)
}

brain_raw_file <- list.files(
  brain_dir,
  pattern = "^GSE238140_RAW",
  full.names = TRUE
)

if (length(brain_raw_file) == 0) {
  stop("[GSE238140] RAW archive not found in: ", brain_dir)
}

if (length(list.files(level1_brain_dir)) == 0) {
  message("[GSE238140] Extracting main RAW archive to level1 directory...")
  untar(brain_raw_file, exdir = level1_brain_dir)
} else {
  message("[GSE238140] level1 directory already populated; skipping main archive extraction.")
}

message("[GSE238140] level1 contents:")
print(list.files(level1_brain_dir))


## ----------------------------------------------------------
## 3A. Extract each GSM-level archive into raw10x/GSM*
## ----------------------------------------------------------

brain_raw10x_dir <- file.path(brain_dir, "raw10x")

if (dir.exists(brain_raw10x_dir)) {
  unlink(brain_raw10x_dir, recursive = TRUE, force = TRUE)
}

dir.create(brain_raw10x_dir, showWarnings = FALSE, recursive = TRUE)

inner_tars_brain <- list.files(
  level1_brain_dir,
  pattern = "\\.tar.gz$",
  full.names = TRUE
)

if (length(inner_tars_brain) == 0) {
  stop("[GSE238140] No GSM-level *.tar.gz files found in: ", level1_brain_dir)
}

for (tarfile in inner_tars_brain) {
  gsm_id <- sub("^((GSM[0-9]+)).*", "\\1", basename(tarfile))
  
  out_dir <- file.path(brain_raw10x_dir, gsm_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(
    "[GSE238140] Extracting sample archive: ",
    basename(tarfile),
    " -> ",
    out_dir
  )
  
  untar(tarfile, exdir = out_dir)
}

sample_dirs_brain <- list.dirs(
  brain_raw10x_dir,
  recursive = FALSE,
  full.names = TRUE
)

message("[GSE238140] Sample directories under raw10x:")
print(sample_dirs_brain)


## ==========================================================
## 4. Prepare GSE230571 Vero/moDC contextual dataset
## ==========================================================

# GSE230571 is used only as a contextual comparison across distinct
# cellular systems. It contains ZIKV-exposed Vero cells and primary
# human monocyte-derived dendritic cells.

dc_dir <- file.path(project_root, "data", "GSE230571")

if (!dir.exists(dc_dir)) {
  stop("[GSE230571] Data directory not found: ", dc_dir)
}

level1_dc_dir <- file.path(dc_dir, "level1")

if (!dir.exists(level1_dc_dir)) {
  dir.create(level1_dc_dir, showWarnings = FALSE, recursive = TRUE)
}

dc_raw_file <- list.files(
  dc_dir,
  pattern = "^GSE230571_RAW",
  full.names = TRUE
)

if (length(dc_raw_file) == 0) {
  stop("[GSE230571] RAW archive not found in: ", dc_dir)
}

if (length(list.files(level1_dc_dir)) == 0) {
  message("[GSE230571] Extracting main RAW archive to level1 directory...")
  untar(dc_raw_file, exdir = level1_dc_dir)
} else {
  message("[GSE230571] level1 directory already populated; skipping main archive extraction.")
}

message("[GSE230571] level1 contents:")
print(list.files(level1_dc_dir))


## ----------------------------------------------------------
## 4A. Reconstruct 10X-style folders for each sample
## ----------------------------------------------------------

all_files_dc <- list.files(level1_dc_dir, full.names = TRUE)

matrix_barcode_files <- all_files_dc[
  grepl("barcodes.tsv.gz|matrix.mtx.gz", all_files_dc)
]

if (length(matrix_barcode_files) == 0) {
  stop("[GSE230571] No barcode or matrix files found in: ", level1_dc_dir)
}

file_names_dc <- basename(matrix_barcode_files)

sample_tags <- unique(
  gsub("-(barcodes.tsv.gz|matrix.mtx.gz)$", "", file_names_dc)
)

message("[GSE230571] Sample tags detected:")
print(sample_tags)

dc_raw10x_dir <- file.path(dc_dir, "raw10x")

if (dir.exists(dc_raw10x_dir)) {
  unlink(dc_raw10x_dir, recursive = TRUE, force = TRUE)
}

dir.create(dc_raw10x_dir, showWarnings = FALSE, recursive = TRUE)

for (tag in sample_tags) {
  out_dir <- file.path(dc_raw10x_dir, tag)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  message("[GSE230571] Creating 10X-style sample directory: ", out_dir)
  
  these_files <- matrix_barcode_files[grepl(tag, matrix_barcode_files)]
  
  barcode_file <- these_files[grepl("barcodes.tsv.gz", these_files)]
  matrix_file  <- these_files[grepl("matrix.mtx.gz", these_files)]
  
  if (length(barcode_file) == 0 || length(matrix_file) == 0) {
    warning("[GSE230571] Missing barcode or matrix file for sample tag: ", tag)
    next
  }
  
  file.copy(
    barcode_file,
    file.path(out_dir, "barcodes.tsv.gz"),
    overwrite = TRUE
  )
  
  file.copy(
    matrix_file,
    file.path(out_dir, "matrix.mtx.gz"),
    overwrite = TRUE
  )
}

sample_dirs_dc <- list.dirs(
  dc_raw10x_dir,
  recursive = FALSE,
  full.names = TRUE
)

message("[GSE230571] Sample directories under raw10x:")
print(sample_dirs_dc)


## ----------------------------------------------------------
## 4B. Attach the appropriate feature file to each sample
## ----------------------------------------------------------

feature_files <- list.files(
  dc_dir,
  pattern = "features.tsv(\\.gz)?$",
  full.names = TRUE
)

message("[GSE230571] Feature files detected:")
print(feature_files)

if (length(feature_files) == 0) {
  stop("[GSE230571] No feature files found in: ", dc_dir)
}

for (sample_dir in sample_dirs_dc) {
  tag <- basename(sample_dir)
  
  if (grepl("ZIKV-Vero-cells", tag)) {
    feature_file <- feature_files[grepl("ZIKV-Vero-cells", feature_files)]
  } else if (grepl("p22086", tag)) {
    feature_file <- feature_files[grepl("p22086", feature_files)]
  } else {
    warning("[GSE230571] Could not assign feature file for sample: ", tag)
    next
  }
  
  if (length(feature_file) == 0) {
    warning("[GSE230571] No matching feature file found for sample: ", tag)
    next
  }
  
  if (length(feature_file) > 1) {
    warning(
      "[GSE230571] Multiple matching feature files found for sample: ",
      tag,
      ". Using the first match."
    )
    feature_file <- feature_file[1]
  }
  
  out_feature_file <- if (grepl("\\.gz$", feature_file)) {
    file.path(sample_dir, "features.tsv.gz")
  } else {
    file.path(sample_dir, "features.tsv")
  }
  
  file.copy(feature_file, out_feature_file, overwrite = TRUE)
}

if (length(sample_dirs_dc) > 0) {
  message("[GSE230571] Example contents of the first reconstructed sample directory:")
  print(list.files(sample_dirs_dc[1]))
}


## ==========================================================
## 5. Create Seurat objects
## ==========================================================

brain_raw10x_root <- file.path(brain_dir, "raw10x")
dc_raw10x_root    <- file.path(dc_dir, "raw10x")

brain <- load_sc10x_dataset(
  root_path = brain_raw10x_root,
  project_name = "ZIKV_fetal_neural"
)

dc <- load_sc10x_dataset(
  root_path = dc_raw10x_root,
  project_name = "ZIKV_Vero_moDC_context"
)

message("[STEP1] Final fetal neural cells: ", ncol(brain))
message("[STEP1] Final Vero/moDC cells: ", ncol(dc))


## ==========================================================
## 6. Save session information
## ==========================================================

session_dir <- file.path(project_root, "session_info")
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)

sink(file.path(session_dir, "STEP1_sessionInfo.txt"))
sessionInfo()
sink()

message("[STEP1] Session information saved to session_info/STEP1_sessionInfo.txt")


############################################################
# End of STEP 1
############################################################
