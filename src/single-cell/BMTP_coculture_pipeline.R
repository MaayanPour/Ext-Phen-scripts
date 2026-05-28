#########################
# Load data
#########################

setwd("/path/BMTP_Coculture/")

data_dir <- "../raw_feature_bc_matrix"
umis <- Read10X(data.dir = data_dir)

# Create Seurat object
BMTP.srt <- CreateSeuratObject(
  counts = umis[["Gene Expression"]],
  min.cells = 3,
  min.features = 200
)

# Add CMO assay
BMTP.srt[["CMO"]] <- CreateAssayObject(
  counts = umis[["Multiplexing Capture"]]
)

#########################
# QC and filtering
#########################

BMTP.srt[["percent.mt"]] <- PercentageFeatureSet(
  BMTP.srt,
  pattern = "^mt-"
)

cells <- fread(
  "../assignment_confidence_table.csv",
  select = "Barcodes"
)

BMTP.srt <- subset(
  BMTP.srt,
  cells = cells$Barcodes
)

#########################
# Demultiplexing
#########################

DefaultAssay(BMTP.srt) <- "CMO"

cmo_features <- c(
  "CMO301","CMO302","CMO303",
  "CMO304","CMO305","CMO306",
  "CMO307","CMO308","CMO309",
  "CMO310","CMO312"
)

BMTP.srt <- subset(BMTP.srt, features = cmo_features)
BMTP.srt <- subset(BMTP.srt, subset = nCount_CMO > 0)

BMTP.srt <- NormalizeData(
  BMTP.srt,
  assay = "CMO",
  normalization.method = "CLR"
)

BMTP.srt <- HTODemux(
  BMTP.srt,
  assay = "CMO",
  positive.quantile = 0.995
)

levels(BMTP.srt$hash.ID)[levels(BMTP.srt$hash.ID) == "Negative"] <- "Unassigned"
levels(BMTP.srt$hash.ID)[levels(BMTP.srt$hash.ID) == "Doublet"] <- "Multiplet"

BMTP.srt$tag <- BMTP.srt$hash.ID

#########################
# RNA QC
#########################

DefaultAssay(BMTP.srt) <- "RNA"

VlnPlot(
  BMTP.srt,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

BMTP.srt <- subset(
  BMTP.srt,
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

BMTP.srt <- run_seurat(BMTP.srt)

#########################
# Visualization
#########################

DimPlot(
  BMTP.srt,
  reduction = "umap",
  group.by = "tag",
  pt.size = 2,
  label = FALSE
)

FeaturePlot(
  BMTP.srt,
  features = c("Cd68", "Vim", "Fn1", "Epcam")
)

#########################
# Cell cycle scoring
#########################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <- convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

m.s.genes <- intersect(m.s.genes, rownames(BMTP.srt))
m.g2m.genes <- intersect(m.g2m.genes, rownames(BMTP.srt))

BMTP.srt <- CellCycleScoring(
  BMTP.srt,
  s.features = m.s.genes,
  g2m.features = m.g2m.genes,
  set.ident = FALSE
)

DimPlot(BMTP.srt, group.by = "Phase")

#########################
# Helper functions
#########################

subset_cells <- function(obj, tag, clusters) {
  
  obj[, obj$tag == tag &
        obj$seurat_clusters %in% clusters]
}

annotate <- function(obj, condition, treatment, week) {
  
  obj$condition <- condition
  obj$treatment <- treatment
  obj$week <- week
  
  return(obj)
}

#########################
# Merge and annotation
#########################

DimPlot(BMTP.srt, group.by = "seurat_clusters")

FeaturePlot(
  BMTP.srt,
  features = c("Cd68", "Vim", "Fn1", "Epcam")
)

#########################
# Assign metadata
#########################

# Cancer cells in isolation
PB1  <- subset_cells(BMTP.srt, "CMO301", "0")
PB12 <- subset_cells(BMTP.srt, "CMO302", "0")
PB24 <- subset_cells(BMTP.srt, "CMO303", "0")

# Fibroblast coculture
FP1  <- subset_cells(BMTP.srt, "CMO304", c("0","1"))
FP12 <- subset_cells(BMTP.srt, "CMO305", c("0","1"))
FP24 <- subset_cells(BMTP.srt, "CMO306", c("0","1"))

# Macrophage coculture
MP1  <- subset_cells(BMTP.srt, "CMO307", c("0","2"))
MP12 <- subset_cells(BMTP.srt, "CMO308", c("0","2"))
MP24 <- subset_cells(BMTP.srt, "CMO309", c("0","2"))

# Naive controls
Macs <- subset_cells(BMTP.srt, "CMO310", "2")
Fibs <- subset_cells(BMTP.srt, "CMO312", "1")

#########################
# Annotation
#########################

# Cancer cells in isolation
PB1  <- annotate(PB1,  "PB1",  "PB",       "week 1")
PB12 <- annotate(PB12, "PB12", "PB",       "week 12")
PB24 <- annotate(PB24, "PB24", "PB",       "week 24")

# Fibroblast coculture
FP1  <- annotate(FP1,  "FP1",  "PB+Fib",   "week 1")
FP12 <- annotate(FP12, "FP12", "PB+Fib",   "week 12")
FP24 <- annotate(FP24, "FP24", "PB+Fib",   "week 24")

# Macrophage coculture
MP1  <- annotate(MP1,  "MP1",  "PB+Macs",  "week 1")
MP12 <- annotate(MP12, "MP12", "PB+Macs",  "week 12")
MP24 <- annotate(MP24, "MP24", "PB+Macs",  "week 24")

# Naive controls
Macs <- annotate(Macs, "Macs", "Macs", "Naive")
Fibs <- annotate(Fibs, "Fibs", "Fibs", "Naive")

#########################
# Merge all cells
#########################

srt_PBMC_all <- merge(
  PB1,
  y = list(
    PB12, PB24,
    FP1, FP12, FP24,
    MP1, MP12, MP24,
    Macs, Fibs
  )
)

#########################
# Subset cancer cells
#########################

srt_PBMC_cancer <- subset(
  srt_PBMC_all,
  subset =
    seurat_clusters == "0"
)

#########################
# Visualization
#########################

cols <- c(
  'PB1'  = '#9ecae1',
  'PB12' = '#2271b5',
  'PB24' = '#09519c',
  
  'FP1'  = '#c7e9c0',
  'FP12' = '#74c476',
  'FP24' = '#248b45',
  
  'MP1'  = '#fb6a4a',
  'MP12' = '#cb1c1d',
  'MP24' = '#a51516'
)

DimPlot(
  srt_PBMC_cancer,
  reduction = "umap",
  group.by = "condition",
  cols = cols,
  pt.size = 2
)

#########################
# Hotspot modules
#########################

load("modules_BMTP.RData")

srt_PBMC_cancer <- AddModuleScore(
  object = srt_PBMC_cancer,
  features = modules_hotspot,
  ctrl = 5,
  name = names(hotspot)
)

saveRDS(
  srt_PBMC_cancer,
  file = "BMTP_cancer.rds"
)

#########################
# Average module scores
#########################

module_cols <- paste0(
  "Module.",
  1:length(modules_BMTP)
)

srt.modules <- CreateSeuratObject(
  counts = t(
    as.matrix(
      srt_PBMC_cancer@meta.data[, module_cols]
    )
  )
)

srt.modules@meta.data <- srt_PBMC_cancer@meta.data

cond_levels <- c(
  "PB1", "PB12", "PB24",
  "FP1", "FP12", "FP24",
  "MP1", "MP12", "MP24"
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