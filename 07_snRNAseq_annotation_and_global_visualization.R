# ============================================================
# Script: 07_snRNAseq_annotation_and_global_visualization.R
# Description: snRNAseq clusters annotation and global visualization of UMAP
# Associated Figures: Figure 2C, 2D and 2E
# ============================================================

library(Seurat)
library(Rcpp)
library(dplyr)
library(patchwork)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)
library(harmony)
library(tidyverse)
library(pheatmap)
library(ComplexHeatmap)

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")

cat("=== Object Structure ===\n")
cat("Seurat Version:", as.character(packageVersion("Seurat")), "\n")
cat("Assay Class:", class(seurat_object@assays$RNA), "\n")
cat("Number of cells:", ncol(seurat_object), "\n")
cat("Number of genes:", nrow(seurat_object), "\n\n")

cat("=== Data Layers ===\n")
cat("counts dimensions:", dim(seurat_object@assays$RNA@counts), "\n")
cat("data dimensions:",   dim(seurat_object@assays$RNA@data),   "\n")
cat("Number of all-zero columns in counts:", sum(colSums(seurat_object@assays$RNA@counts) == 0), "\n\n")

cat("=== Cell Name Check ===\n")
cn <- colnames(seurat_object)
cat("First 5 cell names:", head(cn, 5), "\n")
cat("Any duplicated cell names:", any(duplicated(cn)), "\n")
cat("Any NA in cell names:",  any(is.na(cn)), "\n")
cat("Any empty strings in cell names:", any(cn == ""), "\n")
cat("Total number of cells:", length(cn), "\n\n")

cat("=== Meta.data and Cell Names Alignment ===\n")
cat("meta.data rows:", nrow(seurat_object@meta.data), "\n")
cat("rownames(meta) matches colnames:",
    identical(rownames(seurat_object@meta.data), cn), "\n\n")

cat("=== Idents ===\n")
print(table(Idents(seurat_object)))
cat("Idents length:", length(Idents(seurat_object)), "\n")
cat("Any NA in Idents:", any(is.na(Idents(seurat_object))), "\n")


#findallmarkers------------
markers <- FindAllMarkers(seurat_object,
                          min.pct = 0.25,
                          logfc.threshold = 0.58) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10) %>%
  ungroup()
head(markers)

# 2. Calculate average expression for each cluster using AggregateExpression
avg_exp <- AggregateExpression(
  object = seurat_object,
  assays = 'RNA',
  group.by = 'ident',
  slot = 'data',           # Use normalized data
  return.seurat = FALSE
)$RNA
top10_matrix <- avg_exp[markers$gene, ]

# Create annotation information
name <-0:9
annotation_col <- data.frame(paste("Cluster",name," "))
colnames(annotation_col) <- "Cluster"
rownames(annotation_col) <- colnames(top10_matrix)

##Color definition--------
my_heatmap_colors <- c(
  colorRampPalette(colors = c("#2166AC","#4393C3","#92C5DE","#D1E5F0"))(50),
  colorRampPalette(colors = c("#FDDBC7","#F4A582","#D6604D","#B2182B"))(50)
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

cluster_levels <- c(
  "Phloem parenchyma",
  "Phellem/cortex collenchyma",
  "Epidermis/cortex collenchyma",
  "Vascular cambium/phellogen",
  "Cortex parenchyma",
  "Companion cell",
  "Phloem fiber",
  "Phelloderm",
  "Sieve element",
  "Unknown"
)

# 1. heatmap: TOP 10 markers of cluster0-9--------
# Annotation column: Cluster = 0..9
annotation_col1 <- data.frame(
  Cluster = factor(0:9, levels = 0:9),
  row.names = colnames(top10_matrix)
)

# Rename my_colors according to cluster numbers (0->Phloem parenchyma color, ...)
my_colors_by_cluster <- setNames(
  unname(my_colors[cluster_levels]),
  as.character(0:9)
)

pheatmap(
  top10_matrix,
  scale = "row",
  cluster_cols = FALSE,
  cluster_rows = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  fontsize_row = 8,
  annotation_col = annotation_col1,
  annotation_colors = list(Cluster = my_colors_by_cluster),   # Use my_colors
  color = my_heatmap_colors,
  border_color = NA,
  main = "Top 10 Marker Genes",
  filename = "plots/pheatmap_top10_markers_byCluster.pdf",
  width = 3.5, height = 3
)


# 2. Dotplot: reported genes----

# Rename Idents
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

cell_levels <- c(
  "Phloem parenchyma",
  "Phellem/cortex collenchyma",
  "Epidermis/cortex collenchyma",
  "Vascular cambium/phellogen",
  "Cortex parenchyma",
  "Companion cell",
  "Phloem fiber",
  "Phelloderm",
  "Sieve element",
  "Unknown"
)
seurat_object$cell_type <- factor(as.character(seurat_object$cell_type),
                                  levels = cell_levels)

# marker list
gene_table <- data.frame(
  cluster = c(0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8),
  gene = c(
    "PAxG05Ag0151100", "PAxG01Ag0041100",
    "PAxG01Gg0309400", "PAxG13Gg0072100",
    "PAxG01Ag0047500", "PAxG06Gg0202700",
    "PAxG09Ag0004200", "PAxG05Ag0126600",
    "PAxG11Ag0130700", "PAxG02Gg0095000",
    "PAxG15Gg0111000", "PAxG13Ag0100300",
    "PAxG02Ag0237400", "PAxG06Ag0170100",
    "PAxG09Gg0010000", "PAxG14Gg0018000",
    "PAxG17Gg0195600", "PAxG02Ag0054200"
  ),
  name = c(
    "UMAMIT20", "PP2-B10",
    "ASFT",     "PRX52",
    "LTPG2",    "KCS10",
    "CKS1",     "AIL2",
    "PSBR",     "AED3",
    "PP2-A10",  "SUC2",
    "CESA4",    "CESA7A",
    "BGAL8",    "WOX4a",
    "SEOR1",    "CALS7"
  ),
  stringsAsFactors = FALSE
)

features <- gene_table$gene
id2sym <- setNames(gene_table$name, gene_table$gene)

missing <- setdiff(features, rownames(seurat_object))
if (length(missing) > 0) {
  warning("The following genes are not in seurat_object:\n",
          paste(missing, collapse = "\n"))
}

# Dotplot
p <- DotPlot(seurat_object,
             features = features,
             group.by = "cell_type") +
  theme_classic() +
  theme(
    # All axis text in black
    axis.text   = element_text(size = 8, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    
    # Axis titles
    axis.title  = element_blank(),
    
    # Axis lines + ticks in black
    axis.line   = element_line(color = "black"),
    axis.ticks  = element_line(color = "black"),
    
    # Legend text in black
    legend.text  = element_text(size = 8, color = "black"),
    legend.title = element_text(size = 8, color = "black")
  ) +
  scale_color_distiller(palette = "RdBu") +
  scale_x_discrete(labels = id2sym)

# save
ggsave(file.path(plot_dir, "DotPlot_selected_markers.pdf"),
       plot = p, width = 6, height = 3.5)

cat("Completed DotPlot\n")


# 3. UMAP with cell type--------

rm(list=ls())
library(harmony)
library(ggsci)
library(dplyr) 
library(future)
library(Seurat)
library(cowplot)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(qs)
library(SCP)
library(grid)
library(tidydr) 


plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")
table(Idents(seurat_object))
# 0     1     2     3     4     5     6     7     8     9 
# 18162 10187  3923  2571  1904  1880  1497  1452   261    72 

#p1----------------

# Rename cell_type
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
print(table(seurat_object$cell_type))

# Color definition------------
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

# Plotting
p1 <- DimPlot(seurat_object,
              group.by = "cell_type",
              reduction = "umap",
              label = FALSE) +
  theme_void() +
  theme_dr(xlength = 0.2, ylength = 0.2,
           arrow = grid::arrow(length = unit(0.1, "inches"),
                               ends = "last", type = "closed")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  labs(x = "UMAP 1", y = "UMAP 2", title = "") +
  theme(plot.margin = margin(5.5, 15, 5.5, 5.5),
        panel.grid = element_blank(),
        aspect.ratio = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors)

# Save PDF
ggsave(file.path(plot_dir, "UMAP_celltype_renamed.pdf"),
       plot = p1, width = 6, height = 4)