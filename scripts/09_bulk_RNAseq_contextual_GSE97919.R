############################################################
# Zika scRNA-seq project – STEP 9
# Contextual bulk RNA-seq analysis of GSE97919
#
# Purpose:
#   Analyze the GSE97919 bulk RNA-seq dataset as a contextual
#   comparison for IFN/ISG-associated transcriptional activation.
#
# Inputs:
#   One of the following under the project root:
#     - data/GSE97919/
#     - GSE97919_RAW_extracted/
#     - GSE97919_RAW.tar
#     - GSE97919_RAW_extracted/
#     - GSE97919_RAW.tar
#
# Outputs:
#   results/bulk_GSE97919/
#   figures/supplementary/Supplementary_Figure_S4_bulk_GSE97919_contextual_analysis.png
#   figures/supplementary/Supplementary_Figure_S4_bulk_GSE97919_contextual_analysis.pdf
#   session_info/STEP9_bulk_sessionInfo.txt
############################################################


## ==========================================================
## 0. Reproducibility settings
## ==========================================================

set.seed(1234)

core_isg <- c(
  "IFITM1", "MX1", "OAS1", "OAS2", "OAS3",
  "IFIT1", "IFI6", "ISG15", "STAT1", "IRF7", "RSAD2"
)


## ==========================================================
## 1. Project root and packages
## ==========================================================

project_root <- getwd()
message("[STEP9] Project root: ", project_root)

required_packages <- c(
  "data.table",
  "dplyr",
  "stringr",
  "tidyr",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "magick"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "[STEP9] Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(magick)
})

has_cairo <- requireNamespace("Cairo", quietly = TRUE)


## ==========================================================
## 2. Directories
## ==========================================================

results_dir <- file.path(project_root, "results")
figures_dir <- file.path(project_root, "figures")
supp_fig_dir <- file.path(figures_dir, "supplementary")
session_dir <- file.path(project_root, "session_info")

out_dir <- file.path(results_dir, "bulk_GSE97919")
extract_dir <- file.path(out_dir, "GSE97919_RAW_extracted")
panel_dir <- file.path(out_dir, "figure_panels")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(session_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)


## ==========================================================
## 3. Locate GSE97919 input files
## ==========================================================

count_file_pattern <- "\\.txt\\.gz$|\\.txt$|\\.tsv\\.gz$|\\.tsv$|\\.csv\\.gz$|\\.csv$"

is_valid_gse97919_count_file <- function(x) {
  bn <- basename(x)
  
  grepl("^GSM", bn, ignore.case = TRUE) &
    grepl("\\.(txt|tsv|csv)(\\.gz)?$", bn, ignore.case = TRUE) &
    !grepl("barcodes|features|matrix|genes", bn, ignore.case = TRUE)
}

count_search_dirs <- unique(c(
  Sys.getenv("GSE97919_DIR", unset = NA_character_),
  file.path(project_root, "data", "GSE97919"),
  file.path(project_root, "GSE97919_RAW_extracted"),
  file.path(project_root, "GSE97919_RAW_extracted")
))

count_search_dirs <- count_search_dirs[!is.na(count_search_dirs) & nzchar(count_search_dirs)]
count_search_dirs <- count_search_dirs[dir.exists(count_search_dirs)]

tar_search_dirs <- unique(c(
  Sys.getenv("GSE97919_DIR", unset = NA_character_),
  file.path(project_root, "data", "GSE97919"),
  project_root
))

tar_search_dirs <- tar_search_dirs[!is.na(tar_search_dirs) & nzchar(tar_search_dirs)]
tar_search_dirs <- tar_search_dirs[dir.exists(tar_search_dirs)]

message("[STEP9] Count-file search directories:")
print(count_search_dirs)

message("[STEP9] Archive search directories:")
print(tar_search_dirs)

files <- character(0)

for (dir_i in count_search_dirs) {
  files_i <- list.files(
    dir_i,
    pattern = count_file_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  files_i <- files_i[is_valid_gse97919_count_file(files_i)]
  files <- c(files, files_i)
}

files <- unique(files)

if (length(files) == 0) {
  tar_candidates <- character(0)
  
  for (dir_i in tar_search_dirs) {
    tar_i <- list.files(
      dir_i,
      pattern = "^GSE97919_RAW(\\.tar|\\.tar\\.gz|\\.tgz)?$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
    
    tar_candidates <- c(tar_candidates, tar_i)
  }
  
  tar_candidates <- unique(tar_candidates)
  
  if (length(tar_candidates) == 0) {
    stop(
      "[STEP9] Could not find GSE97919 input files.\n",
      "Expected extracted GSM*.txt.gz files or GSE97919_RAW.tar under:\n",
      "  data/GSE97919/\n",
      "  GSE97919_RAW_extracted/\n",
      "  GSE97919_RAW_extracted/\n",
      "  \n",
      "  project root"
    )
  }
  
  tar_path <- tar_candidates[1]
  message("[STEP9] Extracting archive: ", tar_path)
  untar(tar_path, exdir = extract_dir)
  
  files <- list.files(
    extract_dir,
    pattern = count_file_pattern,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  files <- files[is_valid_gse97919_count_file(files)]
}

if (length(files) == 0) {
  stop("[STEP9] No valid GSE97919 GSM-level count files found.")
}

files <- unique(normalizePath(files, winslash = "/", mustWork = TRUE))

# The same GSE97919 files may exist both in the project root and in .
# Keep one copy per filename to avoid duplicate sample IDs and duplicate columns.
files <- files[!duplicated(basename(files))]

message("[STEP9] Number of GSE97919 count files after filename deduplication: ", length(files))
message("[STEP9] First files:")
print(utils::head(basename(files), 10))


## ==========================================================
## 4. Parse sample metadata
## ==========================================================

parse_sample_info <- function(path) {
  base <- basename(path)
  base <- sub("\\.gz$", "", base, ignore.case = TRUE)
  base <- sub("\\.txt$|\\.tsv$|\\.csv$", "", base, ignore.case = TRUE)
  
  parts <- strsplit(base, "_")[[1]]
  gsm <- parts[1]
  label <- if (length(parts) > 1) paste(parts[-1], collapse = "_") else base
  
  condition <- dplyr::case_when(
    str_detect(label, regex("control|mock", ignore_case = TRUE)) ~ "Control",
    str_detect(label, regex("ifn|interferon", ignore_case = TRUE)) ~ "IFNb",
    str_detect(label, regex("zika|zikv", ignore_case = TRUE)) ~ "ZIKV",
    TRUE ~ "Unknown"
  )
  
  timepoint <- str_extract(label, "d\\d+")
  line <- str_extract(label, "Z\\d+")
  replicate <- ifelse(
    length(parts) >= 2 && str_detect(parts[2], "^\\d+$"),
    parts[2],
    NA_character_
  )
  
  data.frame(
    sample_id = gsm,
    label = label,
    condition = condition,
    timepoint = timepoint,
    line = line,
    replicate = replicate,
    file = path,
    stringsAsFactors = FALSE
  )
}

meta <- bind_rows(lapply(files, parse_sample_info)) %>%
  arrange(sample_id)

# If the same GSM appears more than once because data were found in two directories,
# retain the first instance.
if (anyDuplicated(meta$sample_id) > 0) {
  duplicated_samples <- unique(meta$sample_id[duplicated(meta$sample_id)])
  warning(
    "[STEP9] Duplicate sample IDs detected and deduplicated: ",
    paste(duplicated_samples, collapse = ", ")
  )
  meta <- meta %>%
    distinct(sample_id, .keep_all = TRUE)
}

if (any(meta$condition == "Unknown")) {
  warning(
    "[STEP9] Some samples have Unknown condition and will be excluded from DESeq2: ",
    paste(meta$sample_id[meta$condition == "Unknown"], collapse = ", ")
  )
}

meta <- meta %>%
  filter(condition != "Unknown")

if (nrow(meta) < 4) {
  stop("[STEP9] Too few samples with known conditions after metadata parsing.")
}

if (!all(c("Control", "ZIKV") %in% meta$condition)) {
  stop("[STEP9] Control and ZIKV samples are required for ZIKV vs Control DE analysis.")
}

if (all(is.na(meta$timepoint))) {
  meta$timepoint <- "all"
}

if (all(is.na(meta$line))) {
  meta$line <- "all"
}

meta$condition <- factor(meta$condition, levels = c("Control", "ZIKV", "IFNb"))
meta$condition <- droplevels(meta$condition)
meta$timepoint <- factor(meta$timepoint)
meta$line <- factor(meta$line)

metadata_path <- file.path(out_dir, "GSE97919_sample_metadata.tsv")
data.table::fwrite(meta, metadata_path, sep = "\t")
message("[STEP9] Saved metadata: ", metadata_path)

message("[STEP9] Condition table:")
print(table(meta$condition))


## ==========================================================
## 5. Read counts and build count matrix
## ==========================================================

read_counts_one <- function(path) {
  dt <- data.table::fread(path)
  
  if (ncol(dt) < 2) {
    stop("[STEP9] Count file has fewer than two columns: ", basename(path))
  }
  
  col_names <- colnames(dt)
  col_names_lower <- tolower(col_names)
  
  if ("geneid" %in% col_names_lower) {
    gene_col <- col_names[which(col_names_lower == "geneid")[1]]
  } else if ("gene" %in% col_names_lower) {
    gene_col <- col_names[which(col_names_lower == "gene")[1]]
  } else {
    gene_col <- col_names[1]
  }
  
  possible_count_cols <- which(
    col_names_lower %in% c(
      "count", "counts", "reads", "readcount", "read_count",
      "expected_count", "raw_count"
    )
  )
  
  if (length(possible_count_cols) >= 1) {
    count_col <- col_names[possible_count_cols[1]]
  } else {
    numeric_cols <- which(vapply(dt, is.numeric, logical(1)))
    numeric_cols <- setdiff(numeric_cols, which(col_names == gene_col))
    
    if (length(numeric_cols) == 0) {
      stop("[STEP9] No numeric count column detected in: ", basename(path))
    }
    
    count_col <- col_names[numeric_cols[length(numeric_cols)]]
  }
  
  data.frame(
    gene = as.character(dt[[gene_col]]),
    count = as.numeric(dt[[count_col]]),
    stringsAsFactors = FALSE
  )
}

count_list <- list()

for (i in seq_len(nrow(meta))) {
  sample_id <- meta$sample_id[i]
  file_i <- meta$file[i]
  
  message("[STEP9] Reading: ", basename(file_i))
  
  counts_i <- read_counts_one(file_i) %>%
    filter(!is.na(gene), nzchar(gene)) %>%
    group_by(gene) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")
  
  count_list[[sample_id]] <- counts_i
}

all_genes <- sort(unique(unlist(lapply(count_list, `[[`, "gene"))))

count_mat <- matrix(
  0,
  nrow = length(all_genes),
  ncol = nrow(meta),
  dimnames = list(all_genes, meta$sample_id)
)

for (sample_id in meta$sample_id) {
  counts_i <- count_list[[sample_id]]
  count_mat[counts_i$gene, sample_id] <- counts_i$count
}

count_mat <- count_mat[rowSums(count_mat) > 0, , drop = FALSE]

count_matrix_path <- file.path(out_dir, "GSE97919_counts_matrix.tsv")
count_matrix_df <- data.frame(
  gene = rownames(count_mat),
  count_mat,
  check.names = FALSE
)

data.table::fwrite(
  count_matrix_df,
  count_matrix_path,
  sep = "\t"
)
message("[STEP9] Saved count matrix: ", count_matrix_path)


## ==========================================================
## 6. DESeq2 and VST
## ==========================================================

col_data <- meta %>%
  select(sample_id, condition, timepoint, line, replicate, label) %>%
  tibble::column_to_rownames("sample_id")

count_mat <- count_mat[, rownames(col_data), drop = FALSE]

design_formula <- if (nlevels(col_data$line) > 1 && nlevels(col_data$timepoint) > 1) {
  ~ line + timepoint + condition
} else if (nlevels(col_data$timepoint) > 1) {
  ~ timepoint + condition
} else if (nlevels(col_data$line) > 1) {
  ~ line + condition
} else {
  ~ condition
}

message("[STEP9] DESeq2 design: ", paste(deparse(design_formula), collapse = ""))

dds <- DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = col_data,
  design = design_formula
)

dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)

saveRDS(dds, file.path(out_dir, "GSE97919_dds.rds"))
saveRDS(vsd, file.path(out_dir, "GSE97919_vsd.rds"))


## ==========================================================
## 7. Differential expression: ZIKV vs Control
## ==========================================================

res_zikv <- results(dds, contrast = c("condition", "ZIKV", "Control"))

res_zikv_df <- as.data.frame(res_zikv) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

de_path <- file.path(out_dir, "GSE97919_DE_ZIKV_vs_Control.tsv")
data.table::fwrite(res_zikv_df, de_path, sep = "\t")
message("[STEP9] Saved DE table: ", de_path)


## ==========================================================
## 8. Figure panel helpers
## ==========================================================

theme_journal <- theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.8, "lines")
  )

save_gg_panel <- function(plot_obj, out_path, width = 6.8, height = 5.6, dpi = 600) {
  ggplot2::ggsave(
    filename = out_path,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  if (!file.exists(out_path)) {
    stop("[STEP9] Failed to save plot: ", out_path)
  }
  message("[STEP9] Saved: ", out_path)
}

open_png_device <- function(out_path, width = 3400, height = 2200, res = 300) {
  if (has_cairo) {
    Cairo::CairoPNG(
      filename = out_path,
      width = width,
      height = height,
      units = "px",
      dpi = res,
      bg = "white"
    )
  } else {
    png(
      filename = out_path,
      width = width,
      height = height,
      res = res,
      bg = "white"
    )
  }
}

save_pheatmap_panel <- function(out_path, width, height, res, draw_expr) {
  open_png_device(out_path, width = width, height = height, res = res)
  ok <- FALSE
  
  tryCatch(
    {
      force(draw_expr)
      ok <- TRUE
    },
    error = function(e) {
      dev.off()
      stop(e)
    }
  )
  
  dev.off()
  
  if (!ok || !file.exists(out_path)) {
    stop("[STEP9] Failed to save heatmap: ", out_path)
  }
  
  message("[STEP9] Saved: ", out_path)
}


## ==========================================================
## 9. Panel A: PCA
## ==========================================================

pca_df <- plotPCA(vsd, intgroup = c("condition", "timepoint", "line"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_A <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = timepoint)) +
  geom_point(size = 4, alpha = 0.95) +
  xlab(paste0("PC1 (", percent_var[1], "%)")) +
  ylab(paste0("PC2 (", percent_var[2], "%)")) +
  ggtitle("PCA (VST)") +
  theme_journal

panel_A_path <- file.path(panel_dir, "S4A_PCA_VST.png")
save_gg_panel(p_A, panel_A_path)


## ==========================================================
## 10. Panel B: ISG heatmap
## ==========================================================

ann <- meta %>%
  select(sample_id, condition, timepoint, line) %>%
  as.data.frame()

rownames(ann) <- ann$sample_id
ann$sample_id <- NULL

col_order <- rownames(ann)[order(ann$condition, ann$timepoint, ann$line)]
ann_ordered <- ann[col_order, , drop = FALSE]

isg_present <- core_isg[core_isg %in% rownames(vsd)]

if (length(isg_present) < 3) {
  stop("[STEP9] Fewer than three core ISGs were found in the VST matrix.")
}

mat_isg <- assay(vsd)[isg_present, col_order, drop = FALSE]
mat_isg_z <- t(scale(t(mat_isg)))
mat_isg_z[is.na(mat_isg_z)] <- 0

panel_B_path <- file.path(panel_dir, "S4B_ISG_heatmap.png")

save_pheatmap_panel(
  panel_B_path,
  width = 3800,
  height = 1900,
  res = 300,
  draw_expr = {
    pheatmap::pheatmap(
      mat_isg_z,
      annotation_col = ann_ordered,
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_colnames = FALSE,
      fontsize = 22,
      fontsize_row = 22,
      main = "ISG panel (VST z-score)",
      border_color = NA
    )
  }
)


## ==========================================================
## 11. Panel C: Volcano
## ==========================================================

res_plot <- res_zikv_df %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    padj = as.numeric(padj),
    neglog10FDR = -log10(padj + 1e-300),
    significance = case_when(
      !is.na(padj) & padj < 0.05 & abs(log2FoldChange) >= 1 ~ "Significant",
      TRUE ~ "Not significant"
    )
  )

p_C <- ggplot(res_plot, aes(log2FoldChange, neglog10FDR)) +
  geom_point(aes(alpha = significance, size = significance), color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = 2, linewidth = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.8) +
  scale_alpha_manual(values = c("Significant" = 0.95, "Not significant" = 0.30)) +
  scale_size_manual(values = c("Significant" = 2.1, "Not significant" = 1.2)) +
  xlab("log2 fold change (ZIKV vs Control)") +
  ylab("-log10(FDR)") +
  ggtitle("Volcano: ZIKV vs Control") +
  theme_journal +
  theme(legend.position = "none")

panel_C_path <- file.path(panel_dir, "S4C_Volcano_ZIKV_vs_Control.png")
save_gg_panel(p_C, panel_C_path)


## ==========================================================
## 12. Panel D: sample correlation heatmap
## ==========================================================

normalized_log_counts <- log2(counts(dds, normalized = TRUE)[, col_order, drop = FALSE] + 1)
cor_mat <- cor(normalized_log_counts, method = "pearson")
cor_mat <- cor_mat[col_order, col_order]

panel_D_path <- file.path(panel_dir, "S4D_sample_correlation.png")

# This panel is saved in a wide format to preserve readability
# after it is merged into the 2x2 supplementary figure.
save_pheatmap_panel(
  panel_D_path,
  width = 4200,
  height = 2100,
  res = 300,
  draw_expr = {
    pheatmap::pheatmap(
      cor_mat,
      annotation_col = ann_ordered,
      annotation_row = ann_ordered,
      cluster_cols = FALSE,
      cluster_rows = FALSE,
      show_colnames = FALSE,
      show_rownames = TRUE,
      fontsize = 18,
      fontsize_row = 14,
      main = "Sample correlation (log2 normalized counts)",
      border_color = NA
    )
  }
)


## ==========================================================
## 13. Merge panels
## ==========================================================

required_panel_paths <- c(panel_A_path, panel_B_path, panel_C_path, panel_D_path)

missing_panels <- required_panel_paths[!file.exists(required_panel_paths)]

if (length(missing_panels) > 0) {
  stop(
    "[STEP9] The following panel file(s) are missing and cannot be merged:\n",
    paste(missing_panels, collapse = "\n")
  )
}

trim_white <- function(img) {
  magick::image_trim(img, fuzz = 10)
}

add_outside_label <- function(img, label, band_h = 140, label_size = 120) {
  info <- image_info(img)
  
  band <- image_blank(width = info$width, height = band_h, color = "white")
  band <- image_annotate(
    band,
    label,
    size = label_size,
    color = "black",
    gravity = "west",
    location = "+30+0"
  )
  
  image_append(c(band, img), stack = TRUE)
}

pad_to <- function(img, width, height) {
  image_extent(
    img,
    geometry = paste0(width, "x", height),
    gravity = "northwest",
    color = "white"
  )
}

pad_border <- function(img, pixels = 25) {
  image_border(img, "white", paste0(pixels, "x", pixels))
}

safe_image_read <- function(path) {
  if (!file.exists(path)) {
    stop("[STEP9] Missing image file: ", path)
  }
  magick::image_read(path)
}

img_A <- safe_image_read(panel_A_path) %>% trim_white() %>% add_outside_label("A") %>% pad_border()
img_B <- safe_image_read(panel_B_path) %>% trim_white() %>% add_outside_label("B") %>% pad_border()
img_C <- safe_image_read(panel_C_path) %>% trim_white() %>% add_outside_label("C") %>% pad_border()
img_D <- safe_image_read(panel_D_path) %>% trim_white() %>% add_outside_label("D") %>% pad_border()

left_width <- max(image_info(img_A)$width, image_info(img_C)$width)
right_width <- max(image_info(img_B)$width, image_info(img_D)$width)

top_height <- max(image_info(img_A)$height, image_info(img_B)$height)
bottom_height <- max(image_info(img_C)$height, image_info(img_D)$height)

img_A_box <- pad_to(img_A, left_width, top_height)
img_B_box <- pad_to(img_B, right_width, top_height)
img_C_box <- pad_to(img_C, left_width, bottom_height)
img_D_box <- pad_to(img_D, right_width, bottom_height)

top_row <- image_append(c(img_A_box, img_B_box), stack = FALSE)
bottom_row <- image_append(c(img_C_box, img_D_box), stack = FALSE)
supp_fig <- image_append(c(top_row, bottom_row), stack = TRUE)

supp_png <- file.path(
  supp_fig_dir,
  "Supplementary_Figure_S4_bulk_GSE97919_contextual_analysis.png"
)

supp_pdf <- file.path(
  supp_fig_dir,
  "Supplementary_Figure_S4_bulk_GSE97919_contextual_analysis.pdf"
)

image_write(supp_fig, supp_png, format = "png")
image_write(supp_fig, supp_pdf, format = "pdf")

message("[STEP9] Saved final supplementary figure:")
message("[STEP9]   ", supp_png)
message("[STEP9]   ", supp_pdf)


## ==========================================================
## 14. Notes and session information
## ==========================================================

writeLines(
  c(
    "STEP 9 bulk RNA-seq contextual analysis note",
    "",
    "GSE97919 was analyzed as a contextual bulk RNA-seq dataset.",
    "This analysis summarizes broad IFN/ISG-associated transcriptional activation.",
    "This analysis is interpreted only as a cross-system contextual comparison and not as evidence for fetal neural lineage-specific patterns."
  ),
  file.path(out_dir, "GSE97919_interpretation_note.txt")
)

sink(file.path(session_dir, "STEP9_bulk_sessionInfo.txt"))
sessionInfo()
sink()

message("[STEP9] STEP 9 completed successfully.")

