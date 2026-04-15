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
library(org.Mm.eg.db)



setwd("/path/Kla_Coculture/")
# Load data

data_dir <- "../raw_feature_bc_matrix"
umis <- Read10X(data.dir = data_dir)

# Create Seurat object
Kla.srt <- CreateSeuratObject(
  counts = umis[["Gene Expression"]],
  min.cells = 3,
  min.features = 200
)

# Add CMO assay
Kla.srt[["CMO"]] <- CreateAssayObject(counts = umis[["Multiplexing Capture"]])

# Percent mitochondrial genes 
Kla.srt[["percent.mt"]] <- PercentageFeatureSet(Kla.srt, pattern = "^mt-")

# Filtering barcodes
cells <- fread(
  "../assignment_confidence_table.csv",
  select = "Barcodes"
)

Kla.srt <- subset(Kla.srt, cells = cells$Barcodes)

# Demultiplexing
DefaultAssay(Kla.srt) <- "CMO"

# Keep only relevant CMO features
#for reaction 1:
cmo_features <- c("CMO301","CMO302","CMO303","CMO304","CMO305","CMO306","CMO307","CMO308","CMO309")
#for reaction 2:
#cmo_features <- c("CMO302","CMO303","CMO304","CMO305","CMO306","CMO307","CMO308","CMO309","CMO310","CMO311","CMO312")

Kla.srt <- subset(Kla.srt, features = cmo_features)
Kla.srt <- subset(Kla.srt, subset = nCount_CMO > 0)

# Normalize and demultiplex
Kla.srt <- NormalizeData(Kla.srt, assay = "CMO", normalization.method = "CLR")
Kla.srt <- HTODemux(Kla.srt, assay = "CMO", positive.quantile = 0.995)

# Clean labels
levels(Kla.srt$hash.ID)[levels(Kla.srt$hash.ID) == "Negative"] <- "Unassigned"
levels(Kla.srt$hash.ID)[levels(Kla.srt$hash.ID) == "Doublet"] <- "Multiplet"

# Store tag
Kla.srt$tag <- Kla.srt$hash.ID

# RNA QC
DefaultAssay(Kla.srt) <- "RNA"

VlnPlot(Kla.srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3, pt.size = 0)

Kla.srt <- subset(
  Kla.srt,
  subset = nFeature_RNA > 700 &
    nFeature_RNA < 7500 &
    percent.mt < 20
)

# Normalize data
Kla.srt <- NormalizeData(Kla.srt, normalization.method = "LogNormalize", scale.factor = 10000)
# Find variable genes
Kla.srt <- FindVariableFeatures(Kla.srt, selection.method = "vst", nfeatures = 2000)
# scale data 
Kla.srt <- ScaleData(Kla.srt, vars.to.regress = c("nCount_RNA", "percent.mt"))
# PCA
Kla.srt <- RunPCA(Kla.srt, features = VariableFeatures(object = Kla_R1))

Kla.srt <- FindNeighbors(Kla.srt, dims = 1:10)
Kla.srt <- FindClusters(Kla.srt, resolution = 0.5)
Kla.srt <- RunUMAP(Kla.srt, dims = 1:10,  reduction = "pca")

DimPlot(kcl.srt, reduction = "umap", group.by = 'tag',cols = c('Multiplet' = 'lightgrey',
                                                               'Unassigned' = 'darkgrey',
                                                               'CMO302'= '#c7e9c0', 'CMO303'= '#74c476','CMO304'= '#248b45', 'CMO305'= '#00441b',
                                                               'CMO306'= '#fb6a4a', 'CMO307'= '#cb1c1d', 'CMO308'= '#a51516','CMO309'= '#67090d', 'CMO312' = '#9ecae1',
                                                               'CMO310' = "purple", 'CMO311' = 'pink'),pt.size = 2, label =F)


#########################
# Adding cell cycle score
#########################

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Convert human to mouse
m.s.genes  <- convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

m.s.genes  <- intersect(m.s.genes, rownames(Kla.srt))
m.g2m.genes <- intersect(m.g2m.genes, rownames(Kla.srt))

# Cell cycle scoring
Kla.srt <- CellCycleScoring(
  Kla.srt,
  s.features = m.s.genes,
  g2m.features = m.g2m.genes,
  set.ident = FALSE   
)

DimPlot(Kla.srt, group.by = "Phase")

#########################
# Merging two reactions
#########################


# for reaction 1
rec1 = Kla.srt
rec1$reaction = "1"

# for reaction 2
rec2 = Kla.srt
rec2$reaction = "2"


# running for Kla.srt as rec1:
subset_cells <- function(obj, tag, clusters) {
  obj[, obj$tag == tag & obj$seurat_clusters %in% clusters]
}

annotate <- function(obj, condition, treatment, week) {
  obj$condition <- condition
  obj$treatment <- treatment
  obj$week <- week
  obj
}

# Visualization
DimPlot(rec1, group.by = "seurat_clusters")
FeaturePlot(rec1, features = c("Cd68", "Vim", "Fn1", "Epcam"))


# rec1 subsets
KL1 <- subset_cells(rec1, "CMO301", "0")
FK1 <- subset_cells(rec1, "CMO302", c("0","2"))
MK1 <- subset_cells(rec1, "CMO303", c("0","1"))

KL4 <- subset_cells(rec1, "CMO304", "0")
FK4 <- subset_cells(rec1, "CMO305", c("0","2"))
MK4 <- subset_cells(rec1, "CMO306", c("0","1"))

KL12 <- subset_cells(rec1, "CMO307", "0")
FK12 <- subset_cells(rec1, "CMO308", c("0","2"))
MK12 <- subset_cells(rec1, "CMO309", c("0","1"))


# rec1 annotation

KL1 <- annotate(KL, "KL", "KL", "week 1")
FK1 <- annotate(FK1,    "FK1",    "KL+Fib",  "week 1")
MK1 <- annotate(MK1,    "MK1",    "KL+Macs", "week 1")

KL4 <- annotate(KL4, "KL4", "KL", "week 4")
FK4 <- annotate(FK4, "FK4", "KL+Fib", "week 4")
MK4 <- annotate(MK4, "MK4", "KL+Macs", "week 4")

KL12 <- annotate(KL12, "KL12", "KL", "week 12")
FK12 <- annotate(FK12, "FK12", "KL+Fib", "week 12")
MK12 <- annotate(MK12, "MK12", "KL+Macs", "week 12")

# running for Kla.srt as rec2:
rec2$reaction <- "2"

DimPlot(rec2, group.by = "seurat_clusters")
FeaturePlot(rec2, features = c("Cd68", "Vim", "Fn1", "Epcam"))

KL1 <- subset_cells(rec2, "CMO301", "0")

KL24 <- subset_cells(rec2, "CMO305", "0")
FK24 <- subset_cells(rec2, "CMO306", c("0","2"))
MK24 <- subset_cells(rec2, "CMO307", c("0","1"))

KL36 <- subset_cells(rec2, "CMO308", "0")
FK36 <- subset_cells(rec2, "CMO309", c("0","2"))
MK36 <- subset_cells(rec2, "CMO310", c("0","1"))

Macs <- subset_cells(rec2, "CMO311", "1")
Fibs <- subset_cells(rec2, "CMO312", "2")


# rec2 annotation
KL1 <- annotate(KL1, "KL1", "KL", "week 1")

KL24 <- annotate(KL24, "KL24", "KL", "week 24")
FK24 <- annotate(FK24, "FK24", "KL+Fib", "week 24")
MK24 <- annotate(MK24, "MK24", "KL+Macs", "week 24")

KL36 <- annotate(KL36, "KL36", "KL", "week 36")
FK36 <- annotate(FK36, "FK36", "KL+Fib", "week 36")
MK36 <- annotate(MK36, "MK36", "KL+Macs", "week 36")

Macs <- annotate(Macs, "Macs", "Macs", "NA")
Fibs <- annotate(Fibs, "Fibs", "Fibs", "NA")


# Final merge
srt <- merge(
  KL1,
  y = list(
    KL1, FK1, MK1, KL4, FK4, MK4, KL12, FK12, MK12,
    KL24, FK24, MK24, KL36, FK36, MK36, Macs, Fibs
  )
)


# Standard merging
srt <- NormalizeData(srt, normalization.method = "LogNormalize", scale.factor = 10000)
srt <- FindVariableFeatures(srt, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(srt), 10)

srt <- ScaleData(srt, features = rownames(srt))
srt <- RunPCA(srt, features = VariableFeatures(srt))

srt <- FindNeighbors(srt, dims = 1:15)
srt <- FindClusters(srt, resolution = 0.2)
srt <- RunUMAP(srt, dims = 1:10)


# Visualization of merged Seurat
cols = c('KL1'= '#9ecae1', 'KL4'= '#4292c6','KL12'= '#2271b5', 'KL24'= '#09519c', 'KL36'= '#08306b', 
         'FK1'= '#c7e9c0','FK4'= '#a1d99b', 'FK12'= '#74c476','FK24'= '#248b45', 'FK36'= '#00441b',
         'MK1'= '#fb6a4a', 'MK4'= '#ef3b2c','MK12'= '#cb1c1d', 'MK24'= '#a51516','MK36'= '#67090d')

DimPlot(srt, reduction = "umap", pt.size = 1, label = TRUE)

DimPlot(srt, reduction = "umap", group.by = "reaction", pt.size = 2)
DimPlot(srt, reduction = "umap", group.by = "condition",cols = cols, pt.size = 2)
DimPlot(srt, reduction = "umap", group.by = "treatment", pt.size = 2)
DimPlot(srt, reduction = "umap", group.by = "week", pt.size = 1.5)
DimPlot(srt, reduction = "umap", group.by = "Phase", pt.size = 1.5)

#####################
# Subset cancer cells
#####################

srt.cancer <- subset(
  srt,
  subset = seurat_clusters %in% c("0", "1") &
    !treatment %in% c("Fibs", "Macs")
)
srt.cancer$condition = factor(srt.cancer$condition, levels = c("KL1", "KL4", "KL12", "KL24", "KL36","FK1", "FK4", "FK12", "FK24", "FK36", "MK1", "MK4", "MK12", "MK24", "MK36"))
srt.cancer$cell_type = "Cancer"

# Removing zero genes
srt.cancer <- subset(
  srt.cancer,
  features = rownames(GetAssayData(srt.cancer))[rowSums(GetAssayData(srt.cancer)) > 0]
)

#####
# DEG
#####

Idents(object = srt.cancer) <- srt.cancer@meta.data$treatment
srt.cancer.markers <- FindAllMarkers(srt.cancer, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- srt.cancer.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
srt.cancer$treatment = factor(srt.cancer$treatment, levels = c("KL", "KL+Macs","KL+Fib"))

palette <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)


dittoHeatmap(srt.cancer, 
             genes = top10$gene,
             assay = "scaled", 
             cluster_cols = F, 
             cluster_rows = FALSE,
             scale = "row",
             heatmap.colors = palette,
             annot.by = c("condition", "treatment"),
             annotation_colors = list(condition = cols,
                                      treatment = c("KL" = "#2271b5", "KL+Macs" = "brown3","KL+Fib" = "#74c476")))


# Averaging per condition:
srt.cancer_avg = AverageExpression(srt.cancer,group.by = "condition", return.seurat = T)
srt.cancer_avg$orig.ident = factor(srt.cancer_avg$orig.ident, levels = c("KL1", "KL4", "KL12", "KL24", "KL36","FK1", "FK4", "FK12", "FK24", "FK36", "MK1", "MK4", "MK12", "MK24", "MK36"))
dittoHeatmap(srt.cancer_avg, 
             genes = top10$gene,
             cluster_cols = F, 
             cluster_rows = FALSE,
             scale = "row",
             heatmap.colors = palette,
             annot.by = c("orig.ident"),
             annotation_colors = list(orig.ident = cols
             )
)

#################################################
# Loading, naming and scoring of Hotspot modules
#################################################

# Loading table 1, Hotspot modules 
hotspot = read.csv("../data/modules_table_hotspot.csv", sep="\t", header=TRUE)
modules_hotspot <- hotspot %>%
  select(starts_with("Module.")) %>%
  lapply(na.omit)

# Choose modules + resources
modules <- modules_hotspot
resources <- "MSigDB_Hallmarks"

e <- term_enrichment_by_subset(
  modules,
  q_value_threshold = 0.05,
  resources = resources,
  all_symbols = cached_coding_genes
)

plot_df <- as.data.frame(e) %>%
  pivot_longer(
    cols = starts_with("Module"),
    names_to = "Module",
    values_to = "score"
  ) %>%
  mutate(
    term = gsub("^MSigDB_Hallmarks\\.", "", name),
    term = gsub("^HALLMARK_", "", term)
  ) %>%
  filter(score > 1.4) 

ggplot(plot_df, aes(x = term, y = Module, fill = score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(score, 1)), size = 3) +
  scale_fill_gradient(low = "white", high = "#4A708B",
                      name = "-log10(q)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Pathway Enrichment (Hallmarks)",
    x = "Pathway",
    y = "Module"
  )


# Scoring the cancer cells for the module levels

srt.cancer <- AddModuleScore(
  object = srt.cancer,
  features = modules_hotspot,
  ctrl = 5,
  name = names(hotspot)
)

# Seurat for module scores
module_cols <- paste0("Module.", 1:15)
cond_levels <- c(
  "KL1","KL4","KL12","KL24","KL36",
  "FK1","FK4","FK12","FK24","FK36",
  "MK1","MK4","MK12","MK24","MK36"
)

#saving the cancer cells object
saveRDS(srt.cancer, file = "cancer.srt.rds")

srt.cancer.matrix <- as.matrix(srt.cancer@meta.data[, module_cols])
srt.modules <- CreateSeuratObject(counts = t(srt.cancer.matrix))
srt.modules@meta.data <- srt.cancer@meta.data

srt.modules$condition <- factor(srt.modules$condition, levels = cond_levels)
srt.modules_avg <- AverageExpression(
  srt.modules,
  group.by = "condition",
  return.seurat = TRUE
)

srt.modules_avg$orig.ident <- factor(
  srt.modules_avg$orig.ident,
  levels = cond_levels
)

dittoHeatmap(
  srt.modules,
  genes = module_cols,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  annot.by = "condition",
  heatmap.colors = palette,
  annotation_colors = list(condition = cols)
)

dittoHeatmap(
  srt.modules_avg, 
  genes = module_cols,
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  scale = "row",
  heatmap.colors = palette,
  annot.by = "orig.ident",
  annotation_colors = list(condition = cols)
)



#######################
#Human-Mouse conversion
#######################

convertMouseGeneList <- function(x){
  
  require("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
convertHumanGeneList <- function(x){
  
  library("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  
  genesV2 <- getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                    values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  humanx <- unique(genesV2[, 2])
  
  return(humanx)
}


####################################
#Overlapping genes with Gavish et al
####################################

# Convert mouse modules to human
Gavish = as.list(readxl::read_xlsx("../data/Gavish_MP_list.xlsx"))

hotspot_modules_human <- lapply(modules_hotspot, function(genes) {
  
  genes <- unique(genes)
  convertMouseGeneList(genes)
})
# Cleaning empty overlaps
hotspot_modules_human <- hotspot_modules_human[
  sapply(hotspot_modules_human, function(x) length(x) > 5)
]

Hyp.obj <- newGOM(hotspot_modules_human, Gavish)
pvals <- -log10(getMatrix(Hyp.obj, name = "pval"))

Heatmap(
  pvals,
  name = "-log10(p)",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)

################
#Chi square test
################

# Extract matrix from Seurat object
module_scores <- GetAssayData(srt.modules_avg)

# Convert format
long_df <- module_scores %>%
  as.matrix() %>%
  t() %>%                      # rows = conditions
  as.data.frame() %>%
  rownames_to_column("Condition") %>%
  pivot_longer(-Condition, names_to = "Module", values_to = "Score") %>%
  mutate(
    Group = gsub("[0-9]+", "", Condition),
    Time  = as.numeric(gsub("[^0-9]", "", Condition))
  )

ggplot(long_df, aes(Time, Score, color = Group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.2) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_minimal()

results <- long_df %>%
  group_by(Module) %>%
  group_modify(~ {
    fit <- lm(Score ~ Group * Time, data = .x)
    tidy_fit <- tidy(fit)
    
    tibble(
      p_time = tidy_fit$p.value[tidy_fit$term == "Time"],
      p_interaction = min(tidy_fit$p.value[grepl("Group.*:Time", tidy_fit$term)], na.rm = TRUE)
    )
  }) %>%
  ungroup()

condition_time_effects <- long_df %>%
  group_by(Module, Group) %>%
  summarise(
    p_time = coef(summary(lm(Score ~ Time)))[ "Time", "Pr(>|t|)" ],
    .groups = "drop"
  )

#############
#Linear Modal
#############

mat <- as.matrix(GetAssayData(srt.modules_avg))

modules <- rownames(mat)

times <- c(1,4,12,24,36)
conds <- c("KL","FK","MK")

time_vec <- rep(times, times = length(conds))
cond_vec <- rep(conds, each = length(times))

lm_results <- t(apply(mat, 1, function(y) {
  
  df <- data.frame(
    exp = as.numeric(y),
    time = time_vec,
    condition = factor(cond_vec, levels = conds)
  )
  
  fit <- lm(exp ~ time * condition, data = df)
  coefs <- summary(fit)$coefficients
  
  time_p <- coefs["time", "Pr(>|t|)"]
  time_est <- coefs["time", "Estimate"]
  
  int_rows <- grep("time:condition", rownames(coefs))
  
  int_est <- coefs[int_rows, "Estimate"]
  int_p   <- coefs[int_rows, "Pr(>|t|)"]
  
  names(int_est) <- rownames(coefs)[int_rows]
  names(int_p) <- rownames(coefs)[int_rows]
  
  c(
    Time_KL = -log10(time_p) * sign(time_est),
    FK_vs_time = -log10(int_p["time:conditionFK"]) * sign(int_est["time:conditionFK"]),
    MK_vs_time = -log10(int_p["time:conditionMK"]) * sign(int_est["time:conditionMK"])
  )
}))


colnames(lm_results) <- c("Time_KL", "FK_vs_time", "MK_vs_time")
rownames(lm_results) <- modules

Heatmap(
  lm_results,
  name = "-log10(p) * direction",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100)
)



#####
#Macs
#####

Macs = srt[,srt$seurat_clusters == '2'|srt$seurat_clusters == '4']
Macs = Macs[,Macs$condition == 'Macs'|Macs$condition == 'MK1'|Macs$condition == 'MK4'
            |Macs$condition == 'MK12'|Macs$condition == 'MK24'
            |Macs$condition == 'MK36']
Macs = Macs[,Macs$condition != 'KL1']

Macs <- NormalizeData(Macs, normalization.method = "LogNormalize", scale.factor = 10000)
Macs <- FindVariableFeatures(Macs, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Macs), 10)
top10

all.genes <- rownames(Macs)
Macs <- ScaleData(Macs, features = all.genes)
Macs <- RunPCA(Macs, features = VariableFeatures(object = Macs))

Macs$cell_type = "Macs"

#saving the Macrophage object
saveRDS(Macs, file = "Macs.rds")

#####
#Fibs
#####

Fibs = srt.all[,srt.all$seurat_clusters == '3']
Fibs = Fibs[,Fibs$condition == 'Fibs'|Fibs$condition == 'FK1'|Fibs$condition == 'FK4'
            |Fibs$condition == 'FK12'|Fibs$condition == 'FK24'
            |Fibs$condition == 'FK36']

Fibs <- NormalizeData(Fibs, normalization.method = "LogNormalize", scale.factor = 10000)
Fibs <- FindVariableFeatures(Fibs, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(Fibs), 10)
top10

all.genes <- rownames(Fibs)
Fibs <- ScaleData(Fibs, features = all.genes)
Fibs <- RunPCA(Fibs, features = VariableFeatures(object = Fibs))

Fibs$cell_type = "Fibs"

#saving the Fibroblast object
saveRDS(Fibs, file = "Fibs.rds")
