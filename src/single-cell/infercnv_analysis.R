# inferCNV analysis pipeline

library(Seurat)
library(infercnv)
library(ggplot2)

# Load object saved from the KLa co-culture
srt.cancer <- readRDS("srt.cancer.rds")


# Colors and ordering
condition_order <- c(
  "KL1", "KL4", "KL12", "KL24", "KL36",
  "FK1", "FK4", "FK12", "FK24", "FK36",
  "MK1", "MK4", "MK12", "MK24", "MK36"
)

cols <- c(
  'KL1'     = '#9ecae1',
  'KL4'     = '#4292c6',
  'KL12'    = '#2271b5',
  'KL24'    = '#09519c',
  'KL36'    = '#08306b',
  
  'FK1'     = '#c7e9c0',
  'FK4'     = '#a1d99b',
  'FK12'    = '#74c476',
  'FK24'    = '#248b45',
  'FK36'    = '#00441b',
  
  'MK1'     = '#fb6a4a',
  'MK4'     = '#ef3b2c',
  'MK12'    = '#cb1c1d',
  'MK24'    = '#a51516',
  'MK36'    = '#67090d'
)

srt.cancer$condition <- factor(
  srt.cancer$condition,
  levels = condition_order
)

# Visualization

DimPlot(
  srt.cancer,
  reduction = "umap",
  group.by = "condition",
  cols = cols,
  pt.size = 2,
  shuffle = TRUE
)

# Function: run inferCNV
run_infercnv <- function(
    seurat_obj,
    ref_group = "KL1",
    out_dir
) {
  
  counts_matrix <- GetAssayData(
    seurat_obj,
    assay = "RNA",
    layer = "counts"
  )
  
  annot <- as.matrix(seurat_obj$condition)
  
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_matrix,
    annotations_file  = annot,
    gene_order_file   = "gene_order.tsv",
    ref_group_names   = ref_group
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    resume_mode = FALSE,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = TRUE,
    analysis_mode = "cells",
    tumor_subcluster_partition_method = "qnorm",
    denoise = FALSE,
    HMM = FALSE,
    HMM_type = "i6",
    num_threads = 10,
    output_format = "pdf",
    plot_steps = FALSE
  )
  
  return(infercnv_obj)
}

# Run inferCNV
# All cancer cells from co-cultures:
cells_infercnv_all <- run_infercnv(
  seurat_obj = srt.cancer,
  ref_group = "KL1",
  out_dir = "../ALL_with_ref"
)

# Only MK cells, to be able to then find CNV based clones:
Cells_MK <- subset(
  srt.cancer,
  subset = treatment == "KL+Macs" | condition == "KL1"
)

cells_infercnv_MK <- run_infercnv(
  seurat_obj = Cells_MK,
  ref_group = "KL1",
  out_dir = "../MK_with_ref"
)

# Only FK cells, to be able to then find CNV based clones:
Cells_FK <- subset(
  srt.cancer,
  subset = treatment == "KL+Fib" | condition == "KL1"
)

cells_infercnv_FK <- run_infercnv(
  seurat_obj = Cells_FK,
  ref_group = "KL1",
  out_dir = "../FK_with_ref"
)

# InferCNV subclones
# Assign the appropriate cells and infercnv object, for example: 
Cells = Cells_MK
cells_infercnv = cells_infercnv_MK

# Cutting the infercnv-based tree to 3 clones
k_obs_groups <- 3

obs_hcl <- hclust(
  dist(
    t(
      cells_infercnv@expr.data[
        ,
        unlist(cells_infercnv@observation_grouped_cell_indices)
      ]
    ),
    method = "euclidean"
  ),
  method = "ward.D2"
)

subclone <- cutree(
  tree = as.hclust(obs_hcl),
  k = k_obs_groups
)

# Add subclone metadata
Cells$subclone <- NA
Cells$subclone[names(subclone)] <- subclone
Cells$subclone <- factor(Cells$subclone)

# UMAP visualization
DimPlot(
  Cells,
  reduction = "umap",
  group.by = "subclone",
  pt.size = 2,
  label = TRUE,
  label.size = 6
)

DimPlot(
  Cells[, Cells$condition != "KL1"],
  reduction = "umap",
  group.by = "subclone",
  cols = c(
    "1" = "#F0FFFF",
    "2" = "#C1CDCD",
    "3" = "#838B8B"
  ),
  pt.size = 2,
  label = TRUE,
  label.size = 6
)

# Module plotting
module_features <- c(
  "Module.4",
  "Module.9",
  "Module.13",
  "Module.15"
)

for (mod in module_features) {
  
  print(
    VlnPlot(
      Cells[, Cells$condition != "KL1"],
      features = mod,
      split.by = "subclone",
      group.by = "condition"
    )
  )
  
}

# Subclone-condition table
subclone_condition_table <- table(
  Cells$condition,
  Cells$subclone
)

heatmap(scale(subclone_condition_table))

