# ============================================================
# Script: 18_periderm_differentiation_visualization.R
# Description: Visualize trajectory, pseudotime, and gene dynamics.
# Associated Figures: Figure 7
# ============================================================

suppressPackageStartupMessages({
  library(monocle)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggridges)
  library(RColorBrewer)
  library(viridis)
})

# Set up directories
outdir   <- "results"
plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load pre-computed CDS and BEAM results
N_USE <- 90
cds <- readRDS(file.path(outdir, paste0("cds_top", N_USE, ".rds")))
BEAM_res <- readRDS(file.path(outdir, "BEAM_res_all_candidates.rds"))

# Define cell type colors
celltype_colors <- c(
  "phellem"    = "#D57B70",
  "phelloderm" = "#1f90b4",
  "phellogen"  = "#F0AC54"
)

# 2. Plot basic trajectory (Pseudotime, State, Cell type)
cat("Plotting trajectory base plots...\n")
p1 <- plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
p3 <- plot_cell_trajectory(cds, color_by = "group", size = 1, show_backbone = TRUE) + 
  scale_color_manual(values = celltype_colors)

p_all <- p1 + p3 + plot_layout(nrow = 1)
ggsave(file.path(plot_dir, "traj_basic_pseudotime_celltype.pdf"), p_all, width = 8, height = 4)

# 3. Ridge plots (Cell density along Pseudotime)
cat("Plotting ridge plots...\n")
plotdf <- pData(cds)
p_ridge <- ggplot(plotdf, aes(x = Pseudotime, y = cell_type, fill = cell_type)) +
  geom_density_ridges(scale = 1, alpha = 0.8) +
  scale_fill_manual(values = celltype_colors) +
  theme_minimal() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  scale_x_continuous(position = "top") +
  scale_y_discrete(position = "right") +
  labs(x = "Pseudotime", y = NULL)

ggsave(file.path(plot_dir, "traj_ridge_celltype.pdf"), p_ridge, width = 6, height = 2.5)

# 4. Gene dynamics (Specific genes on trajectory)
# Genes from the branched heatmap (Cell fate 1, Pre-branch, Cell fate 2)
marker_genes <- c(
  # ---- Cell fate 1 (Top branch) ----
  "PAxG09Gg0096900" = "LACS2",
  "PAxG08Gg0138500" = "KCS11",
  "PAxG10Gg0211200" = "GELP86",
  "PAxG10Gg0188800" = "ABCG16",
  "PAxG10Gg0176300" = "GPAT5",
  "PAxG12Ag0045600" = "NAC058",
  
  # ---- Pre-branch (Middle, small cluster) ----
  "PAxG10Gg0068200" = "KCS11",
  "PAxG09Ag0004200" = "CKS1",
  "PAxG17Ag0151800" = "TUA2",
  "PAxG02Gg0222000" = "XTH16",
  
  # ---- Cell fate 2 (Bottom branch) ----
  "PAxG06Gg0076100" = "HTA12",
  "PAxG05Gg0028800" = "GAMMA-H2AX",
  "PAxG03Ag0108600" = "TUB8",
  "PAxG01Gg0008000" = "CMT3",
  "PAxG12Ag0066900" = "DEL1"
)

cds_genes <- rownames(cds)
marker_in <- names(marker_genes)[names(marker_genes) %in% cds_genes]

for (g in marker_in) {
  sym <- marker_genes[g]
  p_gene <- plot_cell_trajectory(cds, markers = g, use_color_gradient = TRUE, show_backbone = TRUE) +
    viridis::scale_color_viridis(option = "magma") +
    ggtitle(sym) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))
  
  ggsave(file.path(plot_dir, paste0("traj_gene_", sym, ".pdf")), p_gene, width = 4, height = 3.5)
}

# 5. Branched Heatmap (requires BEAM result)
cat("Plotting branched heatmap...\n")
BEAM_sig <- subset(BEAM_res, qval < 0.05)
if (nrow(BEAM_sig) >= 20) {
  pdf(file.path(plot_dir, "BEAM_branched_heatmap.pdf"), width = 6, height = 8)
  plot_genes_branched_heatmap(
    cds[rownames(BEAM_sig), ], 
    branch_point = 1,
    num_clusters = 6, 
    cores = 1,
    use_gene_short_name = TRUE, 
    show_rownames = FALSE
  )
  dev.off()
}

cat("Trajectory visualization completed!\n")