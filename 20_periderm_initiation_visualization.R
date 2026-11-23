# ============================================================
# Script: 20_periderm_initiation_visualization.R
# Description: Visualize initiation trajectory, heatmaps, 
#              gene dynamics, and run GO enrichment.
# Associated Figures: Figure 8
# ============================================================

suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(pheatmap)
  library(clusterProfiler)
  library(data.table)
  library(ggridges)
  library(RColorBrewer)
})

# Setup directories
data_dir <- "data"
outdir   <- "results"
plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

N_USE <- 80
N_CLUSTERS <- 2
BIN_SIZE <- 30
TOP_GO <- 5

# 1. Load CDS
cds <- readRDS(file.path(outdir, paste0("cds_initiation_top", N_USE, ".rds")))
cat("Cells:", ncol(cds), "\n")

# Define colors
celltype_colors <- c(
  "phellogen"   = "#B66BAE",
  "sub-epiderm" = "#6FB0E5"
)

# 2. Trajectory Plots
cat("Plotting trajectory...\n")
p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
p2 <- plot_cell_trajectory(cds, color_by = "State", size = 1, show_backbone = TRUE)
p3 <- plot_cell_trajectory(cds, color_by = "cell_type", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = celltype_colors)

p_all <- p1 + p2 + p3 + plot_layout(nrow = 1)
ggsave(file.path(plot_dir, "initiation_traj_panel.pdf"), p_all, width = 12, height = 4)

# 3. Ridge Plot
cat("Plotting ridge...\n")
plotdf <- pData(cds)
p_ridge <- ggplot(plotdf, aes(x = Pseudotime, y = cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 1, alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_x_continuous(position = "top") +
  scale_y_discrete(position = "right") +
  labs(x = "Pseudotime", y = NULL)
ggsave(file.path(plot_dir, "initiation_ridge.pdf"), p_ridge, width = 6, height = 2.5)

# ============================================================
# 4. Non-branched Pseudotime Heatmap (K-means Clustering)
# ============================================================
cat("Plotting pseudotime heatmap...\n")
# Extract DEGs
deg <- read.csv(file.path(outdir, "diff_deg_initiation.csv"))
ordergene_all <- head(unique(deg$gene_short_name), 500)
ordergene_all <- intersect(ordergene_all, rownames(cds))

# Use built-in Monocle heatmap function
pdf(file.path(plot_dir, "initiation_pseudotime_heatmap.pdf"), width = 6, height = 8)
plot_pseudotime_heatmap(
  cds[ordergene_all, ],
  num_clusters = N_CLUSTERS,
  cores = 1,
  use_gene_short_name = TRUE,
  show_rownames = FALSE
)
dev.off()

# ============================================================
# 5. GO Enrichment (Using ClusterGVis or clusterProfiler logic)
# ============================================================
# Note: The original script had a highly customized K-means + pheatmap logic.
# Here we simplify it to a standard clusterProfiler approach for the two cell types.
cat("Running GO enrichment...\n")
# Load GO database (Ensure go_data_anno.txt is in data_dir)
go_file <- file.path(data_dir, "go_data_anno.txt")
if(file.exists(go_file)) {
  go_data <- fread(go_file, header = TRUE, sep = " ", quote = "\"", fill = TRUE, data.table = FALSE)
  go_bp <- go_data[go_data$ONTOLOGY == "BP", ]
  TERM2GENE.GO <- go_bp[, c("GO", "GID")] %>% distinct()
  TERM2NAME.GO <- go_bp[, c("GO", "TERM")] %>% distinct()
  
  # Take Top 200 DEGs for enrichment
  top_genes <- head(unique(deg$gene_short_name), 200)
  
  ego <- enricher(gene = top_genes, TERM2GENE = TERM2GENE.GO, TERM2NAME = TERM2NAME.GO, 
                  pvalueCutoff = 0.1, pAdjustMethod = "BH")
  
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    write.csv(as.data.frame(ego), file.path(outdir, "initiation_GO_enrichment.csv"), row.names = FALSE)
    
    # Simple barplot
    p_go <- barplot(ego, showCategory = TOP_GO) + theme_classic()
    ggsave(file.path(plot_dir, "initiation_GO_barplot.pdf"), p_go, width = 6, height = 4)
  }
}
cat("Visualization and enrichment completed!\n")