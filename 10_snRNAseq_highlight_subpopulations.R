# ============================================================
# Script: 10_highlight_subpopulations.R
# Description: Highlight specific cell subpopulations on UMAP
# Associated Figures: Figure 3A, 4A, 5A, 6A, 3H
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(grid)
  library(tidydr)
  library(dplyr)
})

# ===== Directories =====
data_dir <- "data"
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ===== Load object; remove unresolved mixed-label cells and fill NAs =====
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))
seurat_object$subcell_type <- coalesce(seurat_object$subcell_type,
                                       seurat_object$cell_type)

# ============================================================
# Universal function: highlight one group on a grey UMAP
# ============================================================
plot_highlight <- function(obj, cell_names, highlight_name, highlight_color) {
  hl_list <- list(cell_names)
  names(hl_list) <- highlight_name
  color_dict <- c(highlight_color)
  names(color_dict) <- highlight_name
  
  p <- DimPlot(obj,
               cells.highlight = hl_list,
               cols.highlight  = color_dict,
               cols            = "lightgrey",
               sizes.highlight = 1,
               pt.size         = 0.5,
               label           = FALSE) +
    theme_void() +
    theme_dr(xlength = 0.2, ylength = 0.2,
             arrow = grid::arrow(length = unit(0.1, "inches"),
                                 ends = "last", type = "closed")) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "") +
    theme(plot.margin  = margin(5.5, 15, 5.5, 5.5),
          panel.grid   = element_blank(),
          panel.border = element_blank(),
          aspect.ratio = 1) +
    scale_color_manual(
      values = c(color_dict, "Unselected" = "lightgrey"),
      labels = c(setNames(highlight_name, highlight_name), "Unselected" = "other")
    )
  return(p)
}

# ============================================================
# 1. Highlight Phellem: Figure 3A
# ============================================================
cat("Plotting Phellem...\n")
cells_phellem <- rownames(seurat_object@meta.data)[
  seurat_object$subcell_type %in% "Phellem"]
p1 <- plot_highlight(seurat_object, cells_phellem, "Phellem", "#D57B70")
ggsave(file.path(plot_dir, "UMAP_highlight_phellem.pdf"), p1, width = 5, height = 4)

# ============================================================
# 2. Highlight Phellogen (cluster 3 in EPI/PER): Figure 4A
# ============================================================
cat("Plotting Phellogen...\n")
cells_phellogen <- rownames(seurat_object@meta.data)[
  seurat_object$seurat_clusters %in% 3 &
    seurat_object$orig.ident %in% c("EPI", "PER")]
p2 <- plot_highlight(seurat_object, cells_phellogen, "Phellogen", "#F0AC54")
ggsave(file.path(plot_dir, "UMAP_highlight_phellogen.pdf"), p2, width = 5, height = 4)

# ============================================================
# 3. Highlight Vascular cambium (cluster 3 in PVC): Figure 4A
# ============================================================
cat("Plotting Vascular Cambium...\n")
cells_vascular <- rownames(seurat_object@meta.data)[
  seurat_object$seurat_clusters == 3 &
    seurat_object$orig.ident %in% "PVC"]
p3 <- plot_highlight(seurat_object, cells_vascular, "Vascular cambium", "#609895")
ggsave(file.path(plot_dir, "UMAP_highlight_vascularcambium.pdf"), p3, width = 5, height = 4)

# ============================================================
# 4. Highlight Phelloderm (cluster 7 in EPI/PER): Figure 5A
# ============================================================
cat("Plotting Phelloderm...\n")
cells_phelloderm <- rownames(seurat_object@meta.data)[
  seurat_object$seurat_clusters == 7 &
    seurat_object$orig.ident %in% c("EPI", "PER")]
p4 <- plot_highlight(seurat_object, cells_phelloderm, "Phelloderm", "#1f90b4")
ggsave(file.path(plot_dir, "UMAP_highlight_phelloderm.pdf"), p4, width = 5, height = 4)

# ============================================================
# Define common theme for multi-color plots to avoid redundancy
# ============================================================
my_theme <- theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = "last", type = "closed")) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  theme(plot.margin  = margin(5.5, 15, 5.5, 5.5),
        panel.grid   = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1)

# ============================================================
# 5. Highlight cortex subtypes (multi-color): Figure 6A
# ============================================================
cat("Plotting Cortex Subtypes...\n")
collen1_cells    <- rownames(seurat_object@meta.data)[
  seurat_object$subcell_type %in% "Cortex collenchyma (EPI)"]
collen2_cells    <- rownames(seurat_object@meta.data)[
  seurat_object$subcell_type %in% "Cortex collenchyma (PER)"]
parenchyma_cells <- rownames(seurat_object@meta.data)[
  seurat_object$cell_type %in% "Cortex parenchyma"]

highlight_list_multi <- list(
  "Cortex collenchyma (EPI)" = collen1_cells,
  "Cortex collenchyma (PER)" = collen2_cells,
  "Cortex parenchyma"       = parenchyma_cells
)
highlight_colors_multi <- c(
  "Cortex collenchyma (EPI)" = "#ACC763",
  "Cortex collenchyma (PER)" = "#89A346",
  "Cortex parenchyma"        = "#F6B86A"
)

p_multi <- DimPlot(seurat_object,
                   cells.highlight = highlight_list_multi,
                   cols.highlight  = highlight_colors_multi,
                   cols            = "lightgrey",
                   sizes.highlight = 1,
                   pt.size         = 0.5,
                   label           = FALSE) +
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 5), ncol = 1)) +
  scale_color_manual(
    values = c(highlight_colors_multi, "Unselected" = "lightgrey"),
    breaks = c("Cortex collenchyma (EPI)", "Cortex collenchyma (PER)",
               "Cortex parenchyma")
  )
ggsave(file.path(plot_dir, "UMAP_highlight_cortex_subtypes.pdf"),
       p_multi, width = 6, height = 5)

# ============================================================
# 6. Highlight Subepidermis: Figure 6A
# ============================================================
cat("Plotting Subepidermis...\n")
cells_subepidermis <- rownames(seurat_object@meta.data)[
  seurat_object$subcell_type %in% "Subepidermis"]
p6 <- plot_highlight(seurat_object, cells_subepidermis, "Subepidermis", "#3D6831")
ggsave(file.path(plot_dir, "UMAP_highlight_subepidermis.pdf"), p6, width = 5, height = 4)

# ============================================================
# 7. Dying (RNS3-positive) vs normal epidermis: Figure 3H
# ============================================================
cat("Plotting dying vs normal epidermis...\n")
epiderm <- subset(seurat_object, subcell_type == "Epidermis")
# Note: LayerData is compatible with Seurat v5
gene_expr <- LayerData(epiderm, assay = "RNA", layer = "data")["PAxG08Gg0069900", ]
epiderm$RNS3_status <- ifelse(gene_expr > 0, "RNS3_positive", "RNS3_negative")
print(table(epiderm$RNS3_status))

highlight_list <- list(
  "Dying epidermal cells"  = rownames(epiderm@meta.data)[epiderm$RNS3_status == "RNS3_positive"],
  "Normal epidermal cells" = rownames(epiderm@meta.data)[epiderm$RNS3_status == "RNS3_negative"]
)
cols_rns3 <- c(
  "Dying epidermal cells"  = "#6E511E",
  "Normal epidermal cells" = "#8FC320"
)

p7 <- DimPlot(seurat_object,
              cells.highlight = highlight_list,
              cols.highlight  = cols_rns3,
              cols            = "lightgrey",
              sizes.highlight = 1,
              pt.size         = 0.5,
              label           = FALSE) +
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_color_manual(
    values = c(cols_rns3, "Unselected" = "lightgrey"),
    labels = c("Dying epidermal cells" = "Dying epidermal cells",
               "Normal epidermal cells" = "Normal epidermal cells",
               "Unselected" = "other")
  )
ggsave(file.path(plot_dir, "UMAP_highlight_dyingEpidermis.pdf"),
       p7, width = 5, height = 4)

cat("All highlight UMAP plots completed.\n")