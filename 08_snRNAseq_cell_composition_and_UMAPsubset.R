# ============================================================
# Script: 08_snRNAseq_cell_composition_and_UMAPsubset.R
# Description: snRNAseq cell composition of three samples, cell numbers of each cluster, and UMAP subset of EPI and PER samples
# Associated Figures: Figure 2F, S4D, and 2G
# ============================================================

library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")

# 1. CellRatio_bySample----------
print(table(seurat_object$orig.ident))
if (!"cell_type" %in% colnames(seurat_object@meta.data)) {
  seurat_object$cell_type <- Idents(seurat_object)
}
print(table(seurat_object$cell_type))

# Cell type color palette
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

# Create table
sample_table <- as.data.frame(
  table(seurat_object@meta.data$orig.ident,
        seurat_object@meta.data$cell_type)
)
names(sample_table) <- c("Samples", "Celltype", "CellNumber")

# Fix sample order to EPI, PER, PVC
sample_table$Samples <- factor(sample_table$Samples,
                               levels = c("EPI", "PER", "PVC"))

# Fix cell type order based on my_colors
sample_table$Celltype <- factor(sample_table$Celltype,
                                levels = names(my_colors))

# p1: X-axis as Samples, Y-axis as Cell proportion
p1 <- ggplot(sample_table,
             aes(x = Samples, weight = CellNumber, fill = Celltype)) +
  geom_bar(position = "fill", width = 0.7,
           linewidth = 0,          # Set line width to 0
           colour = NA) +          # Set border color to NA
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(labels = scales::percent) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        axis.line.x = element_line(colour = "black"),
        axis.line.y = element_line(colour = "black"),
        plot.title = element_text(lineheight = .8, face = "bold",
                                  hjust = 0.5, size = 15),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14)) +
  labs(y = "Percentage", x = NULL, title = "Cell type composition")

ggsave(file.path(plot_dir, "CellRatio_bySample_barplot.pdf"),
       plot = p1, width = 6, height = 5)

# p2: X-axis as Cell proportion, Y-axis as Samples (Flipped)
p2 <- p1 + coord_flip()

ggsave(file.path(plot_dir, "CellRatio_bySample_barplot_flipped.pdf"),
       plot = p2, width = 6, height = 3)


# 2. Cellnumber_byCluster----------
library(ggplot2)
library(dplyr)

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")

# Extract cell type information from Seurat object and count numbers
cell_counts <- seurat_object@meta.data %>%
  count(cell_type)

# Color palette
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
p <- ggplot(cell_counts, aes(x = n, y = reorder(cell_type, n), fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.7) + 
  scale_fill_manual(values = my_colors) +    
  theme_classic() +
  labs(x = "Number of cells", y = "Cell type") +
  theme(
    axis.text.y = element_text(size = 12, color = "black"),
    axis.text.x = element_text(size = 12, color = "black"),
    legend.position = "none"                 
  )
print(p)

# save
ggsave(file.path(plot_dir, "Cell_Count_Barplot.pdf"), p, width = 6, height = 4)


# 3. UMAP subset by samples----------
library(Seurat)
library(ggplot2)
library(grid)
library(tidydr)
library(patchwork)

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")

# Keep only the four specific cell types
keep_types <- c(
  "Phellem/cortex collenchyma",
  "Epidermis/cortex collenchyma",
  "Vascular cambium/phellogen",
  "Phelloderm"
)

seurat_sub <- subset(seurat_object, subset = cell_type %in% keep_types)
print(table(seurat_sub$cell_type))

# Color palette
my_colors <- c(
  "Phellem/cortex collenchyma"   = "#1F78B4",
  "Epidermis/cortex collenchyma" = "#B2DF8A",
  "Vascular cambium/phellogen"   = "#33A02C",
  "Phelloderm"                   = "#E31A1C"
)

# Plotting function (Keep borders, capitalize axis titles)
plot_umap_subset <- function(obj, title, color_map = my_colors) {
  p <- DimPlot(obj,
               group.by = "cell_type",
               reduction = "umap",
               label = FALSE) +
    theme_bw() +                                   
    guides(color = guide_legend(override.aes = list(size = 5))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = title) +
    theme(panel.grid      = element_blank(),
          panel.border    = element_rect(color = "black", fill = NA,
                                         linewidth = 0.8),
          axis.text       = element_text(color = "black", size = 9),
          axis.title      = element_text(color = "black", size = 11),
          axis.ticks      = element_line(color = "black"),
          plot.margin     = margin(5.5, 15, 5.5, 5.5),
          aspect.ratio    = 1,
          plot.title      = element_text(hjust = 0.5, size = 14)) +
    scale_color_manual(values = color_map, drop = TRUE) +
    scale_fill_manual(values = color_map, drop = TRUE)
  return(p)
}

# subset EPI / PER 
epi <- subset(seurat_sub, subset = orig.ident == "EPI")
cat("EPI cell count:", ncol(epi), "\n")
print(table(epi$cell_type))

per <- subset(seurat_sub, subset = orig.ident == "PER")
cat("PER cell count:", ncol(per), "\n")
print(table(per$cell_type))

p_epi <- plot_umap_subset(epi, title = "EPI")
p_per <- plot_umap_subset(per, title = "PER")

ggsave(file.path(plot_dir, "UMAP_EPI_celltype.pdf"),
       plot = p_epi, width = 5, height = 3)
ggsave(file.path(plot_dir, "UMAP_PER_celltype.pdf"),
       plot = p_per, width = 5, height = 3)

# Combine + add outer border
# Hide individual legends, let patchwork collect them to the right
combined <- (p_epi + theme(legend.position = "none")) +
  (p_per + theme(legend.position = "none")) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right",
        legend.justification = "center")

pdf(file.path(plot_dir, "UMAP_EPI_PER_combined.pdf"),
    width = 8, height = 3.2)
print(combined)
grid.rect(x = 0.5, y = 0.5, width = 0.98, height = 0.96,
          gp = gpar(lwd = 1.5, col = "black", fill = NA))
dev.off()