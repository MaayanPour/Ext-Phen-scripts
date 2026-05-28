library(Seurat)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(data.table)
library(stringr)
library(nichenetr)
library(dplyr)
library(RITAN)
library(RITANdata)
library(tidyr)
library(ggplot2)
library(tibble)
library(broom)
library(GeneOverlap)


#########################
# Load data
#########################

setwd("/path/4T1_Coculture/")

data_dir <- "../raw_feature_bc_matrix"
umis <- Read10X(data.dir = data_dir)

# Create Seurat object
T4T1.srt <- CreateSeuratObject(
  counts = umis[["Gene Expression"]],
  min.cells = 3,
  min.features = 200
)

# Add CMO assay
T4T1.srt[["CMO"]] <- CreateAssayObject(
  counts = umis[["Multiplexing Capture"]]
)

#########################
# QC and filtering
#########################

T4T1.srt[["percent.mt"]] <- PercentageFeatureSet(
  T4T1.srt,
  pattern = "^mt-"
)

cells <- fread(
  "../assignment_confidence_table.csv",
  select = "Barcodes"
)

T4T1.srt <- subset(
  T4T1.srt,
  cells = cells$Barcodes
)

#########################
# Demultiplexing
#########################

DefaultAssay(T4T1.srt) <- "CMO"

cmo_features <- c(
  "CMO301","CMO302","CMO303",
  "CMO304","CMO305","CMO306",
  "CMO307","CMO308","CMO309"
)

T4T1.srt <- subset(T4T1.srt, features = cmo_features)
T4T1.srt <- subset(T4T1.srt, subset = nCount_CMO > 0)

T4T1.srt <- NormalizeData(
  T4T1.srt,
  assay = "CMO",
  normalization.method = "CLR"
)

T4T1.srt <- HTODemux(
  T4T1.srt,
  assay = "CMO",
  positive.quantile = 0.995
)

levels(T4T1.srt$hash.ID)[levels(T4T1.srt$hash.ID) == "Negative"] <- "Unassigned"
levels(T4T1.srt$hash.ID)[levels(T4T1.srt$hash.ID) == "Doublet"] <- "Multiplet"

T4T1.srt$tag <- T4T1.srt$hash.ID

#########################
# RNA QC
#########################

DefaultAssay(T4T1.srt) <- "RNA"

VlnPlot(
  T4T1.srt,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

T4T1.srt <- subset(
  T4T1.srt,
  subset =
    nFeature_RNA > 700 &
    nFeature_RNA < 7500 &
    percent.mt < 20
)

#########################
# Standard preprocessing
#########################

run_seurat <- function(obj, dims = 1:10, resolution = 0.5) {

  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000
  )

  obj <- ScaleData(
    obj,
    vars.to.regress = c("nCount_RNA", "percent.mt")
  )

  obj <- RunPCA(
    obj,
    features = VariableFeatures(object = obj)
  )

  obj <- FindNeighbors(obj, dims = dims)

  obj <- FindClusters(
    obj,
    resolution = resolution
  )

  obj <- RunUMAP(
    obj,
    dims = dims,
    reduction = "pca"
  )

  return(obj)
}

T4T1.srt <- run_seurat(T4T1.srt)

#########################
# Visualization
#########################

DimPlot(
  T4T1.srt,
  reduction = "umap",
  group.by = "tag",
  pt.size = 2,
  label = FALSE
)

FeaturePlot(
  T4T1.srt,
  features = c("Cd68", "Vim", "Fn1", "Epcam")
)

#########################
# Cell cycle scoring
#########################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <- convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

m.s.genes <- intersect(m.s.genes, rownames(T4T1.srt))
m.g2m.genes <- intersect(m.g2m.genes, rownames(T4T1.srt))

T4T1.srt <- CellCycleScoring(
  T4T1.srt,
  s.features = m.s.genes,
  g2m.features = m.g2m.genes,
  set.ident = FALSE
)

DimPlot(T4T1.srt, group.by = "Phase")

#########################
# Helper functions
#########################

subset_cells <- function(obj, tag, clusters) {

  obj[, obj$tag == tag &
        obj$seurat_clusters %in% clusters]
}

annotate <- function(obj, condition, treatment, time) {

  obj$condition <- condition
  obj$treatment <- treatment
  obj$time <- time

  return(obj)
}

#########################
# Merge and annotation
#########################

DimPlot(srt_4T1, group.by = "seurat_clusters")

FeaturePlot(
  srt_4T1,
  features = c("Cd68", "Vim", "Fn1", "Epcam")
)


# Cancer cells in isolation
BT1  <- subset_cells(srt_4T1, "CMO301", "0")
BT12 <- subset_cells(srt_4T1, "CMO302", "0")
BT24 <- subset_cells(srt_4T1, "CMO303", "0")

# Fibroblast coculture
FT1  <- subset_cells(srt_4T1, "CMO304", c("0","1"))
FT12 <- subset_cells(srt_4T1, "CMO305", c("0","1"))
FT24 <- subset_cells(srt_4T1, "CMO306", c("0","1"))

# Macrophage coculture
MT1  <- subset_cells(srt_4T1, "CMO307", c("0","2"))
MT12 <- subset_cells(srt_4T1, "CMO308", c("0","2"))
MT24 <- subset_cells(srt_4T1, "CMO309", c("0","2"))

# Naive controls
Macs <- subset_cells(srt_4T1, "CMO310", "2")
Fibs <- subset_cells(srt_4T1, "CMO311", "1")

#########################
# Annotation
#########################

# Cancer cells in isolation
BT1  <- annotate(BT1,  "BT1",  "BT",       "week 1")
BT12 <- annotate(BT12, "BT12", "BT",       "week 12")
BT24 <- annotate(BT24, "BT24", "BT",       "week 24")

# Fibroblast coculture
FT1  <- annotate(FT1,  "FT1",  "BT+Fib",   "week 1")
FT12 <- annotate(FT12, "FT12", "BT+Fib",   "week 12")
FT24 <- annotate(FT24, "FT24", "BT+Fib",   "week 24")

# Macrophage coculture
MT1  <- annotate(MT1,  "MT1",  "BT+Macs",  "week 1")
MT12 <- annotate(MT12, "MT12", "BT+Macs",  "week 12")
MT24 <- annotate(MT24, "MT24", "BT+Macs",  "week 24")

# Naive controls
Macs <- annotate(Macs, "Macs", "Macs", "Naive")
Fibs <- annotate(Fibs, "Fibs", "Fibs", "Naive")

#########################
# Merge all cells
#########################

srt_4T1_all <- merge(
  BT1,
  y = list(
    BT12, BT24,
    FT1, FT12, FT24,
    MT1, MT12, MT24,
    Macs, Fibs
  )
)

#########################
# Subset cancer cells
#########################

srt_4T1_cancer <- subset(
  srt_4T1_all,
  subset =
    seurat_clusters == "0"
)


cols <- c(
  'BT1'  = '#9ecae1',
  'BT12' = '#2271b5',
  'BT24' = '#09519c',
  
  'FT1'  = '#c7e9c0',
  'FT12' = '#74c476',
  'FT24' = '#248b45',
  
  'MT1'  = '#fb6a4a',
  'MT12' = '#cb1c1d',
  'MT24' = '#a51516'
)

DimPlot(
  srt_4T1_cancer,
  reduction = "umap",
  group.by = "condition",
  cols = cols,
  pt.size = 2
)

#########################
# Hotspot modules
#########################

load("modules_4T1.RData")

srt_4T1_cancer <- AddModuleScore(
  object = srt_4T1_cancer,
  features = modules_hotspot,
  ctrl = 5,
  name = names(hotspot)
)


saveRDS(
  srt_4T1_cancer,
  file = "4T1_cancer.rds"
)

#########################
# Average module scores
#########################

module_cols <- paste0(
  "Module.",
  1:length(modules_4T1)
)

srt.modules <- CreateSeuratObject(
  counts = t(
    as.matrix(
      srt_4T1_cancer@meta.data[, module_cols]
    )
  )
)

srt.modules@meta.data <- srt_4T1_cancer@meta.data

cond_levels <- c(
  "BT1", "BT12", "BT24",
  "FT1", "FT12", "FT24",
  "MT1", "MT12", "MT24"
)

srt.modules$condition <- factor(
  srt.modules$condition,
  levels = cond_levels
)

srt.modules_avg <- AverageExpression(
  srt.modules,
  group.by = "condition",
  return.seurat = TRUE
)

#########################
# Heatmap
#########################

palette <- colorRampPalette(
  rev(brewer.pal(9, "RdBu"))
)(100)

dittoHeatmap(
  srt.modules_avg,
  genes = module_cols,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  heatmap.colors = palette,
  annot.by = "orig.ident",
  annotation_colors = list(orig.ident = cols)
)



