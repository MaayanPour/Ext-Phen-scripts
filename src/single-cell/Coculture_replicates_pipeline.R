library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

#############
### Load data
#############

wk4  <- readRDS("srt.cancer.wk4.rds")
wk12 <- readRDS("srt.cancer.wk12.rds")
wk24 <- readRDS("srt.cancer.wk24.rds")

wk4$sample  <- paste0("wk4_",  wk4$condition)
wk12$sample <- paste0("wk12_", wk12$condition)
wk24$sample <- paste0("wk24_", wk24$condition)

############
### Merge
############

srt <- merge(
  x = wk4,
  y = list(wk12, wk24),
  project = "Coculture_Timecourse"
)

srt$week <- sub("_.*", "", srt$sample)
srt$group <- sub(".*_", "", srt$sample)
srt$replicate <- sub(".*([0-9]+)$", "\\1", srt$group)

###################
### Module score
###################

# Load modules
hotspot = read.csv("modules_table_hotspot.csv", sep="\t", header=TRUE)
modules_hotspot <- hotspot %>%
  select(starts_with("Module.")) %>%
  lapply(na.omit)
modules <- modules_hotspot

# Score seurat for modules
srt <- AddModuleScore(
  object = srt,
  features = modules,
  ctrl = 5,
  name = "Module."
)

##########################
### Average by replicate
##########################

meta <- srt@meta.data

df <- meta %>%
  dplyr::select(week, group, replicate, starts_with("Module.")) %>%
  group_by(week, group, replicate) %>%
  summarise(across(starts_with("Module."), mean), .groups = "drop")

###############
### Heatmap
###############

# for example- module 9
heat <- df %>%
  dplyr::select(week, group, Module.9) %>%
  tidyr::pivot_wider(names_from = week, values_from = Module.9)

heat <- heat %>%
  dplyr::select(group, wk4, wk12, wk24)

group_levels <- c(
  paste0("K", 1:4),
  paste0("FK", 1:4),
  paste0("MK", 1:4)
)

heat$group <- factor(heat$group, levels = group_levels)
heat <- heat[order(heat$group), ]

mat <- as.matrix(heat[, c("wk4", "wk12", "wk24")])
rownames(mat) <- heat$group

col_fun <- colorRamp2(
  c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)),
  c("blue", "white", "red")
)

Heatmap(
  mat,
  name = "Module.9",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  column_order = c("wk4", "wk12", "wk24"),
  col = col_fun
)

# selected modules
# for example- week 24:
modules_choose <- paste0("Module.", c(4,5,6,9,10,12,13,15))

module_order <- paste0("Module.", c(6,12,10,5,4,9,13,15))

df24 <- df %>%
  dplyr::filter(week == "wk24") %>%
  dplyr::select(group, all_of(module_order))

group_levels <- c(
  paste0("K", 1:4),
  paste0("FK", 1:4),
  paste0("MK", 1:4)
)

df24$group <- factor(df24$group, levels = group_levels)
df24 <- df24[order(df24$group), ]

mat24 <- as.matrix(df24[, module_order])
rownames(mat24) <- df24$group

mat24_t <- t(mat24)   # <-- flip orientation (modules = rows)

col_fun <- circlize::colorRamp2(
  c(min(mat24_t, na.rm = TRUE), 0, max(mat24_t, na.rm = TRUE)),
  c("blue", "white", "red")
)

ComplexHeatmap::Heatmap(
  mat24_t,
  name = "wk24_modules",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  col = col_fun,
  row_names_side = "left",
  column_names_rot = 45
)

##############
### Boxplots
#############
# for example- module 9
module <- "Module.9"

df$condition <- case_when(
  grepl("^K", df$group)  ~ "KL",
  grepl("^FK", df$group) ~ "FK",
  grepl("^MK", df$group) ~ "MK"
)

df$week <- factor(df$week, levels = c("wk4", "wk12", "wk24"))

df$cond_week <- paste0(
  df$condition,
  gsub("wk", "", df$week)
)

cols <- c(
  "KL4"  = "#4292c6",
  "KL12" = "#2271b5",
  "KL24" = "#09519c",
  "FK4"  = "#c7e9c0",
  "FK12" = "#74c476",
  "FK24" = "#248b45",
  "MK4"  = "#fb6a4a",
  "MK12" = "#cb1c1d",
  "MK24" = "#a51516"
)

ggplot(df, aes(x = week, y = .data[[module]], color = cond_week)) +
  
  geom_boxplot(aes(fill = cond_week),
               position = position_dodge(width = 0.8),
               width = 0.6,
               outlier.shape = NA,
               alpha = 0.7) +
  
  geom_point(position = position_jitterdodge(
    jitter.width = 0.15,
    dodge.width = 0.8
  ),
  size = 2) +
  
  theme_classic() +
  
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  
  labs(x = "Time", y = module)
