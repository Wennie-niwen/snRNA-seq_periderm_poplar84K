# ============================================================
# Script: 15_snRNAseq_scoring_analysis.R
# Description: Cell cycle phase visualization, G2M ridge plot, 
#              and Heterochromatin module scoring.
# Associated Figures: Figure 4B, 4C, 5E, 6D
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggridges)
  library(stringr)
  library(SCP) # Required for CaseMatch
  library(grid)
  library(tidydr)
  library(SeuratObject)
  library(tidyverse)
  library(data.table)
  library(purrr)
})



set.seed(123)

# ===== Setup directories =====
data_dir <- "data"
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ===== Load Seurat object =====
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# ============================================================
# 1. Cell Cycle Phase UMAP: Figure 4B
# ============================================================
cat("Plotting Cell Cycle Phase UMAP...\n")
print(table(seurat_object$Phase))

phase_colors <- c("G1" = '#625D9E', "S" = '#E5D2DD', "G2M" = "#57C3F3")

p_phase <- DimPlot(seurat_object, group.by = "Phase", reduction = "umap", label = FALSE) +
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"), ends = "last", type = "closed")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  theme(plot.margin = margin(5.5, 15, 5.5, 5.5),
        panel.grid = element_blank(),
        aspect.ratio = 1) +
  scale_color_manual(values = phase_colors, breaks = c("G1", "S", "G2M")) +
  scale_fill_manual(values = phase_colors)

ggsave(file.path(plot_dir, "UMAP_Phase.pdf"), p_phase, width = 5, height = 3)

# ============================================================
# 2. G2M Score Ridge Plot (Phellogen vs VC): Figure 4C
# ============================================================
cat("Plotting G2M Score Ridge Plot...\n")
# Assign groups
seurat_object$group3 <- "Other"
is_cluster3 <- as.character(seurat_object$seurat_clusters) == "3"
seurat_object$group3[seurat_object$orig.ident %in% c("EPI", "PER") & is_cluster3] <- "Phellogen"
seurat_object$group3[seurat_object$orig.ident %in% "PVC" & is_cluster3] <- "VascularCambium"
seurat_object$group3 <- factor(seurat_object$group3, levels = c("Phellogen", "VascularCambium", "Other"))

plot_data <- FetchData(seurat_object, vars = c("G2M.Score", "group3"))

p_ridge <- ggplot(plot_data, aes(x = G2M.Score, y = group3, fill = group3)) +
  geom_density_ridges(scale = 5, rel_min_height = 0.02, bandwidth = 0.03,
                      alpha = 0.85, color = "white", linewidth = 0.3) +
  coord_cartesian(xlim = c(-0.1, 0.5)) +
  scale_fill_manual(values = c(
    "Phellogen"       = "#F0AC54",
    "VascularCambium" = "#609895",
    "Other"           = "lightgrey"
  )) +
  scale_y_discrete(expand = expansion(add = c(0.1, 0.3))) +   
  labs(x = "G2M.Score", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y     = element_text(size = 11, color = "black"),
    axis.line.y     = element_blank(),
    axis.ticks.y    = element_blank()
  )

ggsave(file.path(plot_dir, "Ridge_G2M_Score.pdf"), p_ridge, width = 5, height = 3)

# ============================================================
# 3. Heterochromatin Score: Figure 5E
# ============================================================
cat("Calculating Heterochromatin Score...\n")
hetero_genes <- read.csv(file.path(data_dir, "hetero_marker_genes.csv"))$gene

# Match genes using SCP::CaseMatch
hetero_genes <- SCP::CaseMatch(search = hetero_genes, match = rownames(seurat_object))

# Add module score
seurat_object <- AddModuleScore(
  object = seurat_object,
  features = list(hetero_genes),
  name = "Heterochromatin",
  ctrl = 105,          
  seed = 123
)

# Rename clusters to cell types
if (all(levels(Idents(seurat_object)) %in% as.character(0:9))) {
  seurat_object <- RenameIdents(seurat_object,
                                "0" = "Phloem parenchyma",
                                "1" = "Phellem/cortex collenchyma",
                                "2" = "Epidermis/cortex collenchyma",
                                "3" = "Vascular cambium/phellogen",
                                "4" = "Cortex parenchyma",
                                "5" = "Companion cell",
                                "6" = "Phloem fiber",
                                "7" = "Phelloderm",
                                "8" = "Sieve element",
                                "9" = "Unknown"
  )
}
seurat_object$cell_type <- Idents(seurat_object)
seurat_object$cell_type <- factor(
  seurat_object$cell_type,
  levels = c("Phloem parenchyma", "Phellem/cortex collenchyma",
             "Epidermis/cortex collenchyma", "Vascular cambium/phellogen",
             "Cortex parenchyma", "Companion cell", "Phloem fiber",
             "Phelloderm", "Sieve element", "Unknown")
)

my_colors <- c(
  "Phloem parenchyma"            = "#A6CEE3",
  "Phellem/cortex collenchyma"   = "#1F78B4",
  "Epidermis/cortex collenchyma" = "#B2DF8A",
  "Vascular cambium/phellogen"   = "#33A02C",
  "Cortex parenchyma"            = "#FDBF6F",
  "Companion cell"               = "#FF7F00",
  "Phloem fiber"                 = "#FB9A99",
  "Phelloderm"                   = "#E31A1C",
  "Sieve element"                = "#CAB2D6",
  "Unknown"                      = "#6A3D9A"
)

p_hetero <- VlnPlot(seurat_object,
                    features = "Heterochromatin1",
                    group.by = "cell_type",
                    cols = my_colors,
                    pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +        
  labs(x = "", y = "Heterochromatin score", title = "Heterochromatin score")  

ggsave(file.path(plot_dir, "Vln_Heterochromatin_score.pdf"), p_hetero, width = 6, height = 4)

cat("All scoring analyses completed successfully!\n")



# ============================================================
# 4. Cortex parenchyma, collenchyma and subepidermis Score: Figure 6D
# ============================================================

set.seed(11)


# Create subdirectories for organized output
dirs <- c(
  file.path(out_dir, "signature_data"),
  file.path(plot_dir, "signature_plots")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Target cells for signature scoring
TARGET_CELLS <- c("Cortex parenchyma", "Cortex collenchyma (PER)", "Cortex collenchyma (EPI)", "Subepidermis")

# Filtering thresholds
PADJ_CUTOFF <- 0.05; MIN_LOG2FC <- 0.25; MIN_PCT1 <- 0.10; MIN_PCT_DIFF <- 0.10
TOP_N <- 30; MIN_SIGNATURE_GENES <- 3
ASSAY_USE <- "RNA"; N_BIN <- 24; CTRL_N <- 100; SCORE_SEED <- 717

# Colors
CELL_COLORS <- c(
  "Cortex parenchyma" = "#FDBF6F",
  "Cortex collenchyma (PER)" = "#ACC864",
  "Cortex collenchyma (EPI)" = "#8BA446",
  "Subepidermis" = "#3C6831"
)

other_cell_colors <- c(
  "Phloem parenchyma" = "#A6CEE3",
  "Phellem/cortex collenchyma" = "#1F78B4",
  "Epidermis/cortex collenchyma" = "#B2DF8A",
  "Vascular cambium/phellogen" = "#33A02C",
  "Cortex parenchyma" = "#FDBF6F",
  "Companion cell" = "#FF7F00",
  "Phloem fiber" = "#FB9A99",
  "Phelloderm" = "#E31A1C",
  "Sieve element" = "#CAB2D6",
  "Unknown" = "#6A3D9A"
)

# ============================================================
# 2. Extract Signature Genes from HVG Table
# ============================================================
cat("Loading HVG table...\n")
hvg <- fread(file.path(data_dir, "HVG.csv"), header = TRUE, data.table = FALSE, check.names = FALSE)

# Clean and standardize HVG data
hvg <- hvg %>%
  mutate(
    gene = trimws(as.character(gene)), 
    cluster = trimws(as.character(cluster)),
    # Assume pct.1 and pct.2 might be percentages
    pct.1 = ifelse(max(pct.1, na.rm = TRUE) > 1, pct.1 / 100, pct.1),
    pct.2 = ifelse(max(pct.2, na.rm = TRUE) > 1, pct.2 / 100, pct.2),
    pct_diff = pct.1 - pct.2
  ) %>%
  filter(!is.na(gene), gene != "", !is.na(cluster), cluster != "")

# Match clusters to target cells (handles prefixes)
match_cluster_target <- function(cluster_name) {
  matched <- TARGET_CELLS[vapply(TARGET_CELLS, function(target) {
    cluster_name == target || endsWith(cluster_name, paste0("_", target))
  }, logical(1))]
  if (length(matched) == 1) matched else NA_character_
}
hvg$signature_cell <- vapply(hvg$cluster, match_cluster_target, FUN.VALUE = character(1))

# Filter genes based on thresholds
hvg_filtered <- hvg %>% 
  filter(!is.na(signature_cell)) %>%
  filter(p_val_adj < PADJ_CUTOFF, avg_log2FC >= MIN_LOG2FC, 
         pct.1 >= MIN_PCT1, pct_diff >= MIN_PCT_DIFF) %>%
  arrange(signature_cell, p_val_adj, desc(avg_log2FC), desc(pct_diff)) %>%
  distinct(signature_cell, gene, .keep_all = TRUE) %>%
  group_by(signature_cell) %>%
  slice_head(n = TOP_N) %>%
  ungroup()

write_csv(hvg_filtered, file.path(out_dir, "signature_data", "filtered_signature_genes.csv"))

# ============================================================
# 3. Load Seurat Object and Prepare Metadata
# ============================================================
cat("Loading Seurat object...\n")
seurat_obj <- readRDS(file.path(data_dir, "Periderm.RDS"))

# Update meta.data with target cell tags
seurat_obj@meta.data <- seurat_obj@meta.data %>%
  mutate(
    target_cell = case_when(
      cell_type == "Cortex parenchyma" ~ "Cortex parenchyma",
      subcell_type == "Cortex parenchyma" ~ "Cortex parenchyma",
      subcell_type %in% c("Cortex collenchyma (PER)", "Cortex collenchyma (EPI)", "Subepidermis") ~ subcell_type,
      TRUE ~ NA_character_
    ),
    plot_group = case_when(
      !is.na(target_cell) ~ target_cell,
      !is.na(subcell_type) & subcell_type != "" ~ subcell_type,
      !is.na(cell_type) & cell_type != "" ~ cell_type,
      TRUE ~ "Unknown"
    )
  )

other_groups <- sort(setdiff(unique(seurat_obj$plot_group), TARGET_CELLS))
seurat_obj$plot_group <- factor(seurat_obj$plot_group, levels = c(TARGET_CELLS, other_groups))

# Handle Default Assay and Layers (Seurat v5 compatible)
DefaultAssay(seurat_obj) <- ASSAY_USE
if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  seurat_obj <- tryCatch(JoinLayers(seurat_obj, assay = ASSAY_USE), error = function(e) seurat_obj)
}

expr_mat <- GetAssayData(seurat_obj, assay = ASSAY_USE, layer = "data")
feature_names <- rownames(expr_mat)

# ============================================================
# 4. Resolve Signature Gene Names
# ============================================================
resolve_signature_features <- function(df, feature_names) {
  df %>% mutate(
    feature = case_when(
      gene %in% feature_names ~ gene,
      !is.na(Ptr_gene_id) & Ptr_gene_id %in% feature_names ~ Ptr_gene_id,
      !is.na(Ath_gene_id) & Ath_gene_id %in% feature_names ~ Ath_gene_id,
      !is.na(name) & name != "NA" & name %in% feature_names ~ name,
      TRUE ~ NA_character_
    )
  )
}

signature_mapping <- resolve_signature_features(hvg_filtered, feature_names)
signature_mapping_valid <- signature_mapping %>% filter(!is.na(feature)) %>% distinct(signature_cell, feature, .keep_all = TRUE)

# Create a named list of features for scoring
signature_list <- setNames(
  lapply(TARGET_CELLS, function(ct) {
    signature_mapping_valid %>% filter(signature_cell == ct) %>% pull(feature) %>% unique()
  }), 
  TARGET_CELLS
)

# ============================================================
# 5. Calculate AddModuleScore
# ============================================================
cat("Calculating module scores...\n")
# Use temp names to avoid underscores in cell type names causing issues
score_temp_cols <- paste0("HVGsig_", seq_along(TARGET_CELLS))
seurat_obj <- AddModuleScore(
  object = seurat_obj, features = signature_list, pool = feature_names,
  nbin = N_BIN, ctrl = CTRL_N, assay = ASSAY_USE,
  name = "HVGsig_", seed = SCORE_SEED
)

# Rename score columns for clarity
score_names <- paste0("Score_", make.names(TARGET_CELLS))
for (i in seq_along(TARGET_CELLS)) {
  seurat_obj@meta.data[[score_names[i]]] <- seurat_obj@meta.data[[score_temp_cols[i]]]
}
seurat_obj@meta.data[, score_temp_cols] <- NULL

# ============================================================
# 6. Visualize Scores (Violin Plots)
# ============================================================
cat("Generating violin plots...\n")

for (i in seq_along(TARGET_CELLS)) {
  ct <- TARGET_CELLS[i]
  score_col <- score_names[i]
  
  # Plot 1: All cells grouped by cell type
  p_all <- ggplot(seurat_obj@meta.data, aes(x = plot_group, y = .data[[score_col]], fill = plot_group)) +
    geom_violin(color = "grey35", linewidth = 0.35, trim = TRUE, scale = "width") +
    geom_boxplot(width = 0.12, fill = "white", color = "grey25", linewidth = 0.3, outlier.shape = NA) +
    scale_fill_manual(values = c(CELL_COLORS, other_cell_colors), drop = FALSE) +
    labs(title = paste0(ct, " Signature Score"), x = NULL, y = "AddModuleScore") +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      legend.position = "none"
    )
  
  safe_name <- gsub("[^A-Za-z0-9]+", "_", ct)
  ggsave(file.path(plot_dir, "signature_plots", paste0(safe_name, "_all_cells.pdf")), p_all, width = 6, height = 3)
}

cat("Signature scoring completed successfully!\n")