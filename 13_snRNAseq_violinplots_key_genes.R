# ============================================================
# Script: 13_snRNAseq_violinplots_key_genes.R
# Description: Generate violin plots for specific key genes 
#              in Epidermis vs Phellem, and VC vs Phellogen.
# Associated Figures: Figure 3F, 4H
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# Setup directories
data_dir <- "data"
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load Seurat object once
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# ============================================================
# 1. Epidermis vs Phellem (ABCG transporters): Figure 3F
# ============================================================
cat("\n--- Processing Epidermis vs Phellem ---\n")
seurat_sub1 <- subset(seurat_object, subcell_type %in% c("Epidermis", "Phellem"))
seurat_sub1$comparison <- factor(seurat_sub1$subcell_type, levels = c("Epidermis", "Phellem"))
print(table(seurat_sub1$comparison))

all_genes_1 <- c(
  'PAxG14Gg0063700' = 'ABCG2',
  'PAxG09Ag0041200' = 'ABCG15'
)
features_1 <- intersect(names(all_genes_1), rownames(seurat_sub1))

Idents(seurat_sub1) <- "comparison"
# Note: ident.1 is Epidermis, ident.2 is Phellem (capitalized to match factor levels)
markers_1 <- FindMarkers(seurat_sub1, ident.1 = "Epidermis", ident.2 = "Phellem",
                         features = features_1, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)
markers_1$gene <- rownames(markers_1)
markers_1$symbol <- all_genes_1[markers_1$gene]

# Keep specific symbols for plotting
keep_symbols_1 <- c("ABCG2", "ABCG15")
markers_deg_sub_1 <- markers_1[markers_1$symbol %in% keep_symbols_1, ]
markers_deg_sub_1 <- markers_deg_sub_1[match(keep_symbols_1, markers_deg_sub_1$symbol), ]
markers_deg_sub_1 <- na.omit(markers_deg_sub_1)

features_deg_1 <- markers_deg_sub_1$gene
cat("Genes to plot:", paste(markers_deg_sub_1$symbol, collapse = ", "), "\n")

my_colors_1 <- c("Epidermis" = "#8EC21F", "Phellem" = "#C5756C")
my_ymax <- 5 

p_list_1 <- lapply(features_deg_1, function(g) {
  fc  <- markers_deg_sub_1$avg_log2FC[markers_deg_sub_1$gene == g]
  sym <- all_genes_1[g]
  title_txt <- paste0(sym, "\nlogFC = ", sprintf("%.2f", fc))
  
  p <- VlnPlot(seurat_sub1, features = g, group.by = "comparison", cols = my_colors_1, pt.size = 0.15) +
    coord_cartesian(ylim = c(0, my_ymax)) +  
    NoLegend() +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 11, face = "bold", color = "black"),
      axis.text.x = element_text(size = 10, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 9, color = "black"),
      axis.line   = element_line(color = "black"),
      axis.ticks  = element_line(color = "black"),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      p$layers[[i]]$aes_params$alpha <- 0.15
      p$layers[[i]]$aes_params$size  <- 0.2
    }
  }
  p
})

combined_1 <- wrap_plots(p_list_1, ncol = 2)
ggsave(file.path(plot_dir, "VlnPlot_Epidermis_vs_Phellem_ABCG.pdf"), plot = combined_1, width = 4, height = 2)

# ============================================================
# 2. Vascular Cambium vs Phellogen: Figure 4H
# ============================================================
cat("\n--- Processing VC vs Phellogen ---\n")
# Subset cluster 3 from the original seurat_object
seurat_sub2 <- subset(seurat_object, seurat_clusters %in% c("3"))
seurat_sub2$group <- ifelse(seurat_sub2@meta.data$orig.ident %in% c("EPI", "PER"), "Phellogen", "VascularCambium")
seurat_sub2$group <- factor(seurat_sub2$group, levels = c("VascularCambium", "Phellogen"))
print(table(seurat_sub2$group))

all_genes_2 <- c(
  'PAxG10Ag0173900' = 'GPAT5',
  'PAxG07Gg0074800' = 'CYP86B1',
  'PAxG10Ag0194700' = 'PIP2;5',
  'PAxG10Gg0197600' = 'PIP2;2'
)
features_2 <- intersect(names(all_genes_2), rownames(seurat_sub2))

Idents(seurat_sub2) <- "group"
markers_2 <- FindMarkers(seurat_sub2, ident.1 = "VascularCambium", ident.2 = "Phellogen",
                         features = features_2, logfc.threshold = 0, min.pct = 0, only.pos = FALSE)
markers_2$gene <- rownames(markers_2)
markers_2$symbol <- all_genes_2[markers_2$gene]

keep_symbols_2 <- c("CYP86B1", "GPAT5", "PIP2;2", "PIP2;5")
markers_deg_sub_2 <- markers_2[markers_2$symbol %in% keep_symbols_2, ]
markers_deg_sub_2 <- markers_deg_sub_2[match(keep_symbols_2, markers_deg_sub_2$symbol), ]
markers_deg_sub_2 <- na.omit(markers_deg_sub_2)

features_deg_2 <- markers_deg_sub_2$gene
cat("Genes to plot:", paste(markers_deg_sub_2$symbol, collapse = ", "), "\n")

my_colors_2 <- c("VascularCambium" = "#609895", "Phellogen" = "#F0AC54")

p_list_2 <- lapply(features_deg_2, function(g) {
  fc  <- markers_deg_sub_2$avg_log2FC[markers_deg_sub_2$gene == g]
  sym <- all_genes_2[g]
  title_txt <- paste0(sym, "\nlogFC = ", sprintf("%.2f", fc))
  
  p <- VlnPlot(seurat_sub2, features = g, group.by = "group", cols = my_colors_2, pt.size = 0.15) +
    NoLegend() +
    labs(x = NULL, y = NULL, title = title_txt) +
    theme(
      plot.title  = element_text(hjust = 0.5, size = 10, color = "black"),
      axis.text.x = element_text(size = 9, angle = 30, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 8, color = "black"),
      axis.line   = element_line(color = "black"),
      axis.ticks  = element_line(color = "black"),
      plot.margin = margin(4, 4, 4, 4)
    )
  
  for (i in seq_along(p$layers)) {
    if (inherits(p$layers[[i]]$geom, "GeomPoint")) {
      p$layers[[i]]$aes_params$alpha <- 0.15
      p$layers[[i]]$aes_params$size  <- 0.2
    }
  }
  p
})

combined_2 <- wrap_plots(p_list_2, ncol = 2)
ggsave(file.path(plot_dir, "VlnPlot_VC_vs_Phellogen_4genes.pdf"), plot = combined_2, width = 3, height = 3.5)

cat("\nAll violin plots saved successfully!\n")