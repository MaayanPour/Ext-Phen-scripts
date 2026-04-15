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


# Load the saved Macs and Fibs data
Macs = readRDS("Macs.rds") 
Fibs = readRDS("Fibs.rds")
# Load states from Gavish et al, converted to mouse genes
load("../data/Gavish_Mac_mouse.RData")
load("../data/Gavish_Fib_mouse.RData")


# Function to analyze PCs across conditions
analyze_PC_trend <- function(seurat_obj, PC, conditions, cols, y_pos = NULL) {
  
  # Extract PCA 
  df <- data.frame(
    PC_value = Embeddings(seurat_obj, "pca")[, PC],
    condition = seurat_obj$condition
  )
  
  df$condition <- factor(df$condition, levels = conditions, ordered = TRUE)
  df$condition_num <- as.numeric(df$condition)
  
  # Statistics
  cor_test <- cor.test(df$condition_num, df$PC_value, method = "spearman")
  pval <- cor_test$p.value
  fit <- lm(PC_value ~ condition_num, data = df)
  r2 <- summary(fit)$r.squared
  
  
  if (is.null(y_pos)) {
    y_pos <- max(df$PC_value, na.rm = TRUE)
  }
  
  # Plot
  library(ggplot2)
  
  p <- ggplot(df, aes(x = condition, y = PC_value, fill = condition)) +
    geom_violin(scale = "width") +
    geom_jitter(width = 0.1, size = 0.5) +
    scale_fill_manual(values = cols) +
    theme_classic() +
    annotate(
      "text",
      x = 2,
      y = y_pos,
      label = paste0(
        ", p = ", signif(pval, 3),
        "\nR² = ", round(r2, 3)
      ),
      size = 4
    )
  
  return(list(plot = p, stats = c(p = pval, r2 = r2)))
}

# Printing for Macs
cols_MK <- c( 'Macs' = "orange", 'MK1' = '#fb6a4a', 'MK12' = '#cb1c1d', 'MK24' = '#a51516', 'MK36' = '#67090d' ) 
conds_MK <- c("Macs", "MK1", "MK12", "MK24", "MK36")

res_PC2 <- analyze_PC_trend(Macs, "PC_2", conds_MK, cols_MK)
res_PC3 <- analyze_PC_trend(Macs, "PC_3", conds_MK, cols_MK)

res_PC2$plot + scale_y_reverse()
res_PC3$plot

# Printing for Fibs
cols_FK <- c('Fibs' = "#C9C973", 'FK1'= '#c7e9c0', 'FK12'= '#74c476', 'FK24'= '#248b45','FK36'= '#00441b')
conds_FK <- c("Fibs", "FK1", "FK12", "FK24", "FK36")

res_PC2_fibs <- analyze_PC_trend(Fibs, "PC_2", conds_FK, cols_FK)
res_PC3_fibs <- analyze_PC_trend(Fibs, "PC_3", conds_FK, cols_FK)

res_PC2_fibs$plot
res_PC5_fibs$plot

# Print top genes
get_top_genes_list <- function(seurat_obj, pcs = paste0("PC_", 1:5), nfeatures = 100) {
  
  pca_loadings <- seurat_obj[["pca"]]@feature.loadings
  top_genes_list <- list()
  
  for (pc in pcs) {
    
    pc_loadings <- pca_loadings[, pc]
    sorted_pc <- sort(pc_loadings, decreasing = TRUE)
    
    top_genes_list[[pc]] <- list(
      positive = names(sorted_pc[1:nfeatures]),
      negative = names(sorted_pc[(length(sorted_pc)):(length(sorted_pc)-nfeatures+1)])
    )
  }
  
  return(top_genes_list)
}

build_PC_list <- function(top_genes_list, pcs = paste0("PC_", 1:5), n_genes = 20) {
  
  PC_list <- list()
  
  for (pc in pcs) {
    PC_list[[paste0(pc, "_pos")]] <- top_genes_list[[pc]]$positive[1:n_genes]
    PC_list[[paste0(pc, "_neg")]] <- top_genes_list[[pc]]$negative[1:n_genes]
  }
  
  return(PC_list)
}

# Calculate overlap
run_overlap <- function(PC_list, reference, palette, title) {
  
  Hyp.obj <- newGOM(PC_list, reference)
  
  drawHeatmap(
    Hyp.obj,
    log.scale = TRUE,
    ncolused = 5,
    grid.col = palette,
    note.col = "black"
  )
  
  return(Hyp.obj)
}

top_genes_list_fibs <- get_top_genes_list(Fibs)
PC_list_fibs <- build_PC_list(top_genes_list_fibs)

Hyp.obj_fibs <- run_overlap(
  PC_list_fibs,
  Gavish_Fib_mouse,
  palette = "Greens",
  title = "Fibs"
)

top_genes_list_macs <- get_top_genes_list(Macs)
PC_list_macs <- build_PC_list(top_genes_list_macs)

Hyp.obj_macs <- run_overlap(
  PC_list_macs,
  Gavish_Mac_mouse,
  palette = "Reds",
  title = "Macs"
)


