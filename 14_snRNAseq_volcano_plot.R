# ============================================================
# Script: 14_snRNAseq_volcano_plot.R
# Description: Generate a volcano plot for DEGs between 
#              Phellogen and Vascular Cambium.
# Associated Figures: Figure 4E
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(Seurat)
})

set.seed(123)

# Setup directories
data_dir <- "data"
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# 1. Load data
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# 2. Find Markers
cat("Running FindMarkers...\n")
Idents(seurat_object) <- "subcell_type"
markers <- FindMarkers(
  seurat_object,
  ident.1 = "Phellogen",   
  ident.2 = "Vascular cambium",        
  min.pct = 0.01,            
  logfc.threshold = 0.58,   
  only.pos = FALSE          
)
markers$gene <- rownames(markers)
DEG <- subset(markers, p_val_adj < 0.05)

# 3. Prepare plotting data
data <- DEG %>% 
  mutate(change = as.factor(ifelse(p_val_adj < 0.05 & abs(avg_log2FC) > 0.58,
                                   ifelse(avg_log2FC > 0.58, 'Up', 'Down'), 'No change'))) 

# Note: Modify this list if you want to highlight different genes
core_gene <- c(
  'PAxG03Gg0096900',
  "PAxG01Ag0115000",
  'PAxG07Ag0008100',
  'PAxG05Gg0059700',
  'PAxG13Ag0072300'
)

# Gene symbol mapping
gene_name_map <- c(
  'PAxG03Gg0096900' = 'PXYa',
  "PAxG01Ag0115000" = 'PXYb',
  'PAxG07Ag0008100' = 'MYB84',
  'PAxG05Gg0059700' = 'MYB92',
  'PAxG13Ag0072300' = 'PRX52'
)

data <- data %>%
  mutate(gene_display = case_when(
    gene %in% names(gene_name_map) ~ gene_name_map[gene],
    TRUE ~ gene  
  ))

# 4. Filter extreme values for better visualization
# Removing genes with -log10(q) > 50 to prevent axis squishing
data_filtered <- data %>%
  mutate(neg_log10_q = -log10(p_val_adj)) %>%
  filter(neg_log10_q < 50)

# 5. Generate Volcano Plot
cat("Generating volcano plot...\n")
p <- ggplot(data_filtered, aes(x = avg_log2FC, y = neg_log10_q)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "black") + # Corrected from -1, 1 to match 0.58 threshold
  geom_point(aes(size = neg_log10_q, color = neg_log10_q)) +
  scale_color_gradientn(
    values = seq(0, 1, 0.2),
    colors = c("#39489f", "#39bbec", "#f9ed36", "#f38466", "#b81f25"),
    limits = c(0, 50)  
  ) +
  scale_size_continuous(range = c(0, 1), limits = c(0, 50)) +  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = 'right',
    legend.justification = c(0, 1)
  ) +
  guides(
    col = guide_colorbar(
      title = "-Log10(q-value)",
      ticks.colour = NA,
      reverse = TRUE,
      title.vjust = 0.8,
      barheight = 8,
      barwidth = 1
    ),
    size = "none"
  ) +
  geom_text_repel(
    data = filter(data_filtered, gene %in% core_gene),
    max.overlaps = 20,
    aes(label = gene_display),
    size = 3,
    color = 'black'
  ) +
  xlab("Log2FC") +
  ylab("-Log10(q-value)")

# 6. Save the plot
ggsave(file.path(plot_dir, "Volcano_Phellogen_vs_VC.pdf"), plot = p, width = 6, height = 3)

cat("Volcano plot saved successfully!\n")

