# ============================================================
# Script: 21_cortex_initiation_trajectory.R
# Description: Subset phellogen and cortex collenchyma cells, 
#              construct Monocle trajectory, and visualize.
# Associated Figures: Figure 9A
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(123)

# Setup directories
data_dir <- "data"
outdir   <- "results"
plot_dir <- "plots"

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Extract target cells and merge
# ============================================================
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# Subset target cells
phellogen <- subset(seurat_object, subcell_type %in% c("Phellogen") & orig.ident %in% c("EPI"))

# Note: Fixed the duplicate assignment in the original script. 
# Assuming you meant to extract immature and mature collenchyma separately.
collen_I <- subset(seurat_object, subcell_type %in% c("Cortex collenchyma (EPI)") & orig.ident %in% c("EPI", "PER"))
collen_M <- subset(seurat_object, subcell_type %in% c("Cortex collenchyma (PER)") & orig.ident %in% c("EPI", "PER"))

# Add source prefix to orig_cluster
phellogen$orig_cluster <- "phellogen"
collen_I$orig_cluster <- "collen_I"
collen_M$orig_cluster <- "collen_M"

# Merge
CI <- merge(phellogen, y = c(collen_I, collen_M),
            add.cell.ids = c("phellogen", "collen_I", "collen_M"))
CI <- JoinLayers(CI)
print(table(CI$orig_cluster, CI$orig.ident))

# Save intermediate Seurat object
saveRDS(CI, file.path(data_dir, "Cortex_initiation.RDS"))

# ============================================================
# 2. Construct Monocle CDS
# ============================================================
cat("Constructing Monocle CDS...\n")
exp <- CI@assays$RNA$counts

fdata <- data.frame(
  gene_short_name = row.names(CI), 
  row.names = row.names(CI)
)
pdata <- CI@meta.data

fd <- new("AnnotatedDataFrame", data = fdata)
pd <- new("AnnotatedDataFrame", data = pdata)

CDS <- newCellDataSet(cellData = exp, phenoData = pd, featureData = fd)

# Clean up memory
rm(CI, fd, pd, fdata, pdata, exp)
gc()

# Calculate size factors and dispersions
CDS <- estimateSizeFactors(CDS)
CDS <- estimateDispersions(CDS)
gc()

# ============================================================
# 3. Differential Analysis & Trajectory Construction
# ============================================================
cat("Detecting genes and running differentialGeneTest...\n")
CDS <- detectGenes(CDS, min_expr = 0.1)

# Filter low-expression genes
expressed_genes <- fData(CDS) %>%
  subset(num_cells_expressed >= 15) %>%
  rownames()

# Remove mitochondrial and chloroplast genes
expressed_genes_clean <- expressed_genes[!grepl("GMtg|GPtg", expressed_genes)]

# Differential analysis by cluster
clustering_DEG_genes <- differentialGeneTest(
  CDS[expressed_genes_clean, ], 
  fullModelFormulaStr = '~orig_cluster'
)

# Filter significant DEGs
clustering_DEG_genes_q <- subset(clustering_DEG_genes, qval < 0.000000001)
write.csv(clustering_DEG_genes, file.path(outdir, "12_Cortex_initiation_DEG.csv"), row.names = FALSE)

# Ordering genes based on significant DEGs
ordering_genes <- rownames(clustering_DEG_genes_q)

cat("Constructing trajectory (DDRTree)...\n")
CDS <- setOrderingFilter(CDS, ordering_genes = ordering_genes)
CDS <- reduceDimension(CDS, method = 'DDRTree', max_components = 2)
CDS <- orderCells(CDS)

# Save Monocle object
saveRDS(CDS, file.path(outdir, "12_Cortex_initiation_CDS.rds"))

# ============================================================
# 4. Visualization
# ============================================================
cat("Visualizing trajectory...\n")
# Rename clusters for better plotting
CDS$orig_cluster_new <- as.character(CDS$orig_cluster)
CDS$orig_cluster_new[CDS$orig_cluster_new == "collen_I"] <- "collenchyma"
CDS$orig_cluster_new[CDS$orig_cluster_new == "collen_M"] <- "collenchyma"
CDS$orig_cluster_new <- factor(CDS$orig_cluster_new)

# 1. Color by Pseudotime
p1 <- plot_cell_trajectory(CDS, color_by = "Pseudotime", size = 1, show_backbone = FALSE)
ggsave(file.path(plot_dir, "12_Cortex_initiation_Pseudotime.pdf"), p1, width = 4, height = 4)

# 2. Color by State
p2 <- plot_cell_trajectory(CDS, color_by = "State", size = 1, show_backbone = FALSE)
ggsave(file.path(plot_dir, "12_Cortex_initiation_State.pdf"), p2, width = 4, height = 4)

# 3. Color by orig_cluster_new
p3 <- plot_cell_trajectory(CDS, color_by = "orig_cluster_new", size = 1, show_backbone = FALSE)
ggsave(file.path(plot_dir, "12_Cortex_initiation_CellType.pdf"), p3, width = 4, height = 4)

# 4. Color by original orig_cluster
p4 <- plot_cell_trajectory(CDS, color_by = "orig_cluster", size = 1, show_backbone = FALSE)
ggsave(file.path(plot_dir, "12_Cortex_initiation_OrigCluster.pdf"), p4, width = 4, height = 4)

cat("Trajectory analysis completed successfully!\n")