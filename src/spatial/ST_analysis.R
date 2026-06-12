#!/usr/bin/env Rscript

# Summary pipeline for ST malignant-normal pairing and cell-state correlation heatmaps.
library(Seurat)
library(pheatmap)
library(infercnv)

BASE_DIR <- "../SpatialTranscriptomics"
GYN_DIR <- file.path(BASE_DIR, "GYN")
OUT_DIR <- file.path(GYN_DIR, "ST_results")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GAVISH_PATH <- file.path(BASE_DIR, "Gavish.RData")
GENE_ORDER_FILE <- file.path(BASE_DIR, "gene_ordering_file.txt")

ST_LIST_PATHS <- list(
  OVC = file.path(GYN_DIR, "OVC.list.rds"),
  BRC = file.path(GYN_DIR, "BRC.list.rds"),
  LUNG = file.path(BASE_DIR, "Lung.list.rds")
)

RUN_INFERCNV <- TRUE
INFERCNV_THREADS <- 10
# set the number of clones from cutting the infercnv based tree. possible to change per sample based on the cnv heatmap
INFERCNV_K <- list(
  OVC = 7,
  BRC = 7,
  LUNG = 7
)

# After inferCNV, set which clones are
# malignant/normal/mixed for each sample.this will be based on marker expression and correlation with cell type based annotation
# Example: CLONE_TO_TUMOR_NORMAL$OVC$ovca7 <- list(Malignant = c("3"), Normal = c("0"))

CLONE_TO_TUMOR_NORMAL <- list(
  OVC = list(default = list(Malignant = character(), Normal = c("0"), Mixed = character())),
  BRC = list(default = list(Malignant = character(), Normal = c("0"), Mixed = character())),
  LUNG = list(default = list(Malignant = character(), Normal = c("0"), Mixed = character()))
)

PAIR_DISTANCE_CUTOFF <- list(
  OVC = 600,
  BRC = 600,
  LUNG = 600
)

# Optional per-sample overrides for samples with lower distances.
PAIR_DISTANCE_CUTOFF_SAMPLE <- list(
  BRC = c(BRCA2 = 25, BRCA6 = 300, BRCA7 = 200)
)

N_RANDOM_PAIRINGS <- 500
SET_SEED <- 123

CELL_TYPES <- c(
  "Fibroblasts.Score",
  "Endothelial.cell.Score",
  "Macrophages.Score",
  "T.cell.Score",
  "B.cell.Score"
)

MODULES_GAVISH <- c(
  "MP1..Cell.Cycle...G2.M",
  "MP2..Cell.Cycle...G1.S",
  "MP3..Cell.Cycle.HMG.rich",
  "MP4..Chromatin",
  "MP5.Stress",
  "MP6.Hypoxia",
  "MP7.Stress..in.vitro.",
  "MP12.EMT.I",
  "MP13.EMT.II",
  "MP14.EMT.III",
  "MP18.Interferon.MHC.II..II.",
  "MP19.Epithelial.Senescence",
  "MP20.MYC",
  "MP22.Secreted.I",
  "MP23.Secreted.II",
  "MP31.Alveolar"
)

CELL_TYPE_ORDER <- c(
  "Fibroblasts.Score",
  "Endothelial.cell.Score",
  "B.cell.Score",
  "T.cell.Score",
  "Macrophages.Score"
)

##################
# helper functions
##################

load_st_lists <- function(paths) {
  lists <- lapply(paths, readRDS)
  names(lists) <- names(paths)
  list(
    OVC.st = lists$OVC,
    BRC.st = lists$BRC,
    LUNG.st = lists$LUNG
  )
}

add_gavish_scores_to_list <- function(st_list, Gavish, ctrl = 5) {
  for (sample_name in names(st_list)) {
    message("Adding Gavish scores to: ", sample_name)
    
    st <- st_list[[sample_name]]
    
    for (module_name in names(Gavish)) {
      if (module_name %in% colnames(st@meta.data)) {
        next
      }
      
      before_cols <- colnames(st@meta.data)
      
      st <- AddModuleScore(
        object = st,
        features = list(Gavish[[module_name]]),
        ctrl = ctrl,
        name = "Gavish_tmp"
      )
      
      new_col <- setdiff(colnames(st@meta.data), before_cols)
      names(st@meta.data)[names(st@meta.data) == new_col[1]] <- module_name
    }
    
    st_list[[sample_name]] <- st
  }
  
  return(st_list)
}


get_spatial_coords <- function(st) {
  if (all(c("x_coords", "y_coords") %in% colnames(st@meta.data))) {
    coord <- cbind(st$x_coords, st$y_coords)
    rownames(coord) <- colnames(st)
    return(coord)
  }

  if (length(st@images) == 0) {
    stop("No x_coords/y_coords columns and no spatial image coordinates found.")
  }

  coord <- st@images[[1]]@coordinates[, c("imagerow", "imagecol"), drop = FALSE]
  coord <- as.matrix(coord)
  coord <- coord[colnames(st), , drop = FALSE]
  colnames(coord) <- c("x", "y")
  coord
}

#########################
# infercnv and annotation
#########################

run_infercnv_subclones <- function(st, sample_id, cohort, k_obs_groups) {
  if (!requireNamespace("infercnv", quietly = TRUE)) {
    stop("The infercnv package is not installed/available.")
  }
  if (!file.exists(GENE_ORDER_FILE)) {
    stop("Missing gene ordering file: ", GENE_ORDER_FILE)
  }

  label_col <- if ("malignant" %in% colnames(st@meta.data)) "malignant" else "Tumor_Normal"
  if (!label_col %in% colnames(st@meta.data)) {
    stop("Need either a malignant or Tumor_Normal column before inferCNV.")
  }

  counts <- GetAssayData(st, assay = "Spatial", slot = "counts")
  annot <- as.matrix(st@meta.data[[label_col]])
  rownames(annot) <- colnames(st)

  out_dir <- file.path(OUT_DIR, "infercnv", cohort, sample_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cnv <- infercnv::CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annot,
    gene_order_file = GENE_ORDER_FILE,
    ref_group_names = NULL
  )

  cnv <- infercnv::run(
    cnv,
    resume_mode = FALSE,
    cutoff = 0.1,
    out_dir = out_dir,
    cluster_by_groups = FALSE,
    analysis_mode = "subclusters",
    tumor_subcluster_partition_method = "qnorm",
    denoise = FALSE,
    HMM = FALSE,
    HMM_type = "i6",
    num_threads = INFERCNV_THREADS,
    plot_steps = FALSE
  )

  obs_cells <- unlist(cnv@observation_grouped_cell_indices)
  obs_hcl <- hclust(dist(t(cnv@expr.data[, obs_cells, drop = FALSE]), "euclidean"), "ward.D2")
  subclone <- cutree(as.hclust(obs_hcl), k = k_obs_groups)

  st$subclone <- "0"
  st$subclone[names(subclone)] <- as.character(subclone)
  st$subclone <- factor(st$subclone)

  saveRDS(cnv, file.path(out_dir, "infercnv_object.rds"))
  st
}

get_clone_rule <- function(cohort, sample_id) {
  rules <- CLONE_TO_TUMOR_NORMAL[[cohort]]
  if (!is.null(rules[[sample_id]])) {
    return(rules[[sample_id]])
  }
  rules$default
}

annotate_tumor_normal_by_clone <- function(st, cohort, sample_id) {
  if (!"subclone" %in% colnames(st@meta.data)) {
    if (!"Tumor_Normal" %in% colnames(st@meta.data)) {
      stop("No subclone or Tumor_Normal annotation found for ", cohort, " / ", sample_id)
    }
    return(st)
  }

  rule <- get_clone_rule(cohort, sample_id)
  st$Tumor_Normal <- "Mixed"
  st$Tumor_Normal[as.character(st$subclone) %in% rule$Malignant] <- "Malignant"
  st$Tumor_Normal[as.character(st$subclone) %in% rule$Normal] <- "Normal"
  st$Tumor_Normal[as.character(st$subclone) %in% rule$Mixed] <- "Mixed"

  if (length(rule$Malignant) == 0 && "Tumor_Normal" %in% colnames(st@meta.data)) {
    message("No malignant clone rule set for ", cohort, " / ", sample_id,
            ". Check CLONE_TO_TUMOR_NORMAL before using pair results.")
  }

  st$Tumor_Normal <- factor(st$Tumor_Normal, levels = c("Malignant", "Normal", "Mixed"))
  st
}

######################
# tumor-normal pairing
######################

pair_malignant_normal <- function(st, distance_cutoff) {
  coord <- get_spatial_coords(st)
  malignant_cells <- colnames(st)[st$Tumor_Normal == "Malignant"]
  normal_cells <- colnames(st)[st$Tumor_Normal == "Normal"]

  if (length(malignant_cells) == 0 || length(normal_cells) == 0) {
    warning("Skipping pairing because malignant or normal cells are missing.")
    return(matrix(character(), ncol = 2, dimnames = list(NULL, c("cancer", "noncancer"))))
  }

  dist_mat <- as.matrix(dist(rbind(coord[malignant_cells, , drop = FALSE],
                                   coord[normal_cells, , drop = FALSE])))
  dist_mat <- dist_mat[malignant_cells, normal_cells, drop = FALSE]

  pair_distance <- matrix(NA_real_, nrow = nrow(dist_mat), ncol = 1,
                          dimnames = list(rownames(dist_mat), "distance"))
  pair_normal <- matrix(NA_character_, nrow = nrow(dist_mat), ncol = 1,
                        dimnames = list(rownames(dist_mat), "noncancer"))

  remaining <- dist_mat
  for (cell in rownames(dist_mat)) {
    if (ncol(remaining) == 0 || !cell %in% rownames(remaining)) {
      next
    }
    closest <- names(which.min(remaining[cell, ]))
    pair_distance[cell, "distance"] <- remaining[cell, closest]
    pair_normal[cell, "noncancer"] <- closest
    if (remaining[cell, closest] < distance_cutoff) {
      remaining <- remaining[, setdiff(colnames(remaining), closest), drop = FALSE]
    }
  }

  keep <- !is.na(pair_distance[, "distance"]) & pair_distance[, "distance"] < distance_cutoff
  pairs <- cbind(cancer = rownames(pair_distance)[keep], noncancer = pair_normal[keep, "noncancer"])
  pairs[!is.na(pairs[, "noncancer"]), , drop = FALSE]
}

annotate_pairs <- function(st, pairs) {
  st$group2 <- as.character(st$Tumor_Normal)
  if (nrow(pairs) > 0) {
    st$group2[pairs[, "cancer"]] <- "Malignant_pairs"
    st$group2[pairs[, "noncancer"]] <- "Normal_pairs"
  }
  st$group2 <- factor(st$group2, levels = c("Malignant", "Malignant_pairs", "Normal", "Normal_pairs", "Mixed"))
  st
}

##############
# correlations
##############

safe_scaled_vector <- function(st, feature, cells) {
  vals <- as.numeric(st[[feature]][cells, , drop = TRUE])
  if (all(is.na(vals)) || stats::sd(vals, na.rm = TRUE) == 0) {
    return(rep(NA_real_, length(vals)))
  }
  as.numeric(scale(vals))
}

stars_from_pvals <- function(pval_mat) {
  stars <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat))
  stars[pval_mat < 0.001] <- "***"
  stars[pval_mat < 0.01 & pval_mat >= 0.001] <- "**"
  stars[pval_mat < 0.05 & pval_mat >= 0.01] <- "*"
  stars
}

correlate_pairs_vs_random <- function(st, pairs, cancer_features, normal_features, n_iter) {
  cor_mat <- matrix(NA_real_, nrow = length(cancer_features), ncol = length(normal_features),
                    dimnames = list(cancer_features, normal_features))
  pval_mat <- cor_mat
  random_mean_mat <- cor_mat

  if (nrow(pairs) < 3) {
    warning("Fewer than 3 pairs; returning NA correlations.")
    return(list(cor = cor_mat, pval = pval_mat, random_mean = random_mean_mat, stars = stars_from_pvals(pval_mat)))
  }

  for (state_cancer in cancer_features) {
    for (state_normal in normal_features) {
      cancer_vec <- safe_scaled_vector(st, state_cancer, pairs[, "cancer"])
      normal_vec <- safe_scaled_vector(st, state_normal, pairs[, "noncancer"])

      if (all(is.na(cancer_vec)) || all(is.na(normal_vec))) {
        next
      }

      obs_cor <- cor(cancer_vec, normal_vec, use = "pairwise.complete.obs")
      rnd_cor <- replicate(n_iter, {
        random_normals <- sample(pairs[, "noncancer"], length(pairs[, "noncancer"]), replace = FALSE)
        random_vec <- safe_scaled_vector(st, state_normal, random_normals)
        cor(cancer_vec, random_vec, use = "pairwise.complete.obs")
      })

      cor_mat[state_cancer, state_normal] <- obs_cor
      pval_mat[state_cancer, state_normal] <- mean(abs(rnd_cor) >= abs(obs_cor), na.rm = TRUE)
      random_mean_mat[state_cancer, state_normal] <- mean(rnd_cor, na.rm = TRUE)
    }
  }

  list(cor = cor_mat, pval = pval_mat, random_mean = random_mean_mat, stars = stars_from_pvals(pval_mat))
}

##########
# heatmaps
##########

zavit_order <- function(mat, row_shift = 90, col_shift = 0) {
  mat <- as.matrix(mat)
  mat <- mat[rowSums(is.na(mat)) < ncol(mat), colSums(is.na(mat)) < nrow(mat), drop = FALSE]

  row_scaled <- t(scale(t(mat)))
  row_scaled[is.na(row_scaled)] <- 0
  row_pca <- prcomp(row_scaled)
  row_angle <- (atan2(row_pca$x[, 1], row_pca$x[, 2]) * 180 / pi + row_shift) %% 360
  row_order <- order(row_angle)

  col_scaled <- scale(mat)
  col_scaled[is.na(col_scaled)] <- 0
  col_pca <- prcomp(t(col_scaled))
  col_angle <- (atan2(col_pca$x[, 1], col_pca$x[, 2]) * 180 / pi + col_shift) %% 360
  col_order <- order(col_angle)

  list(
    matrix = col_scaled[row_order, col_order, drop = FALSE],
    row_order = rownames(mat)[row_order],
    col_order = colnames(mat)[col_order],
    row_angle = row_angle,
    col_angle = col_angle
  )
}

plot_correlation_heatmap <- function(cor_mat, stars_mat, cohort, file_prefix) {
  keep_rows <- rowSums(!is.na(cor_mat)) > 0
  keep_cols <- colSums(!is.na(cor_mat)) > 0
  cor_mat <- cor_mat[keep_rows, keep_cols, drop = FALSE]
  stars_mat <- stars_mat[keep_rows, keep_cols, drop = FALSE]

  if (length(intersect(CELL_TYPE_ORDER, colnames(cor_mat))) > 0) {
    ordered_cols <- c(intersect(CELL_TYPE_ORDER, colnames(cor_mat)), setdiff(colnames(cor_mat), CELL_TYPE_ORDER))
    cor_mat <- cor_mat[, ordered_cols, drop = FALSE]
    stars_mat <- stars_mat[, ordered_cols, drop = FALSE]
  }

  z <- zavit_order(t(cor_mat), row_shift = 90, col_shift = 0)
  stars_z <- t(stars_mat)[rownames(z$matrix), colnames(z$matrix), drop = FALSE]

  pdf(file.path(OUT_DIR, paste0(file_prefix, "_zavit_heatmap.pdf")), width = 10, height = 6)
  pheatmap(
    z$matrix,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = stars_z,
    color = colorRampPalette(c("blue", "white", "red"))(99),
    breaks = seq(-2, 2, length.out = 100),
    main = paste(cohort, "zavit-ordered paired correlations"),
    show_colnames = FALSE
  )
  dev.off()

  z
}

sample_distance_cutoff <- function(cohort, sample_id) {
  overrides <- PAIR_DISTANCE_CUTOFF_SAMPLE[[cohort]]
  if (!is.null(overrides) && sample_id %in% names(overrides)) {
    return(unname(overrides[[sample_id]]))
  }
  PAIR_DISTANCE_CUTOFF[[cohort]]
}

#######################
# run sample and cohort
#######################

run_one_sample <- function(st, cohort, sample_id, gavish_modules) {
  message("Running ", cohort, " / ", sample_id)

  st <- fix_module_names(st)

  if (RUN_INFERCNV) {
    st <- run_infercnv_subclones(st, sample_id, cohort, INFERCNV_K[[cohort]])
  }

  st <- annotate_tumor_normal_by_clone(st, cohort, sample_id)

  pairs <- pair_malignant_normal(st, sample_distance_cutoff(cohort, sample_id))
  st <- annotate_pairs(st, pairs)

  available_modules <- intersect(MODULES_GAVISH, colnames(st@meta.data))
  available_cell_types <- intersect(CELL_TYPES, colnames(st@meta.data))

  pair_stats <- correlate_pairs_vs_random(
    st = st,
    pairs = pairs,
    cancer_features = available_modules,
    normal_features = available_cell_types,
    n_iter = N_RANDOM_PAIRINGS
  )

  rownames(pair_stats$cor) <- paste0(rownames(pair_stats$cor), "_", sample_id)
  rownames(pair_stats$pval) <- paste0(rownames(pair_stats$pval), "_", sample_id)
  rownames(pair_stats$random_mean) <- paste0(rownames(pair_stats$random_mean), "_", sample_id)
  rownames(pair_stats$stars) <- paste0(rownames(pair_stats$stars), "_", sample_id)

  list(st = st, pairs = pairs, stats = pair_stats)
}

run_one_cohort <- function(st_list, cohort, gavish_modules) {
  results <- vector("list", length(st_list))
  names(results) <- names(st_list)

  for (sample_id in names(st_list)) {
    results[[sample_id]] <- run_one_sample(st_list[[sample_id]], cohort, sample_id, gavish_modules)
  }

  combined_cor <- do.call(rbind, lapply(results, function(x) x$stats$cor))
  combined_pval <- do.call(rbind, lapply(results, function(x) x$stats$pval))
  combined_random <- do.call(rbind, lapply(results, function(x) x$stats$random_mean))
  combined_stars <- stars_from_pvals(combined_pval)

  z <- plot_correlation_heatmap(
    cor_mat = combined_cor,
    stars_mat = combined_stars,
    cohort = cohort,
    file_prefix = cohort
  )

  cohort_result <- list(
    samples = results,
    combined_cor = combined_cor,
    combined_pval = combined_pval,
    combined_random_mean = combined_random,
    combined_stars = combined_stars,
    zavit = z
  )

  saveRDS(cohort_result, file.path(OUT_DIR, paste0(cohort, "_summary_result.rds")))
  cohort_result
}

#########################
# run analysis
#########################

set.seed(SET_SEED)

Gavish <- loadRData(GAVISH_PATH)
Gavish <- Gavish[intersect(names(Gavish), MODULES_GAVISH)]

st_objects <- load_st_lists(ST_LIST_PATHS)
OVC.st <- st_objects$OVC.st
BRC.st <- st_objects$BRC.st
LUNG.st <- st_objects$LUNG.st

OVC.st <- add_gavish_scores_to_list(OVC.st, Gavish)
BRC.st <- add_gavish_scores_to_list(BRC.st, Gavish)
LUNG.st <- add_gavish_scores_to_list(LUNG.st, Gavish)

all_results <- list(
  OVC = run_one_cohort(OVC.st, "OVC", Gavish),
  BRC = run_one_cohort(BRC.st, "BRC", Gavish),
  LUNG = run_one_cohort(LUNG.st, "LUNG", Gavish)
)

saveRDS(all_results, file.path(OUT_DIR, "all_cohorts_summary_result.rds"))

