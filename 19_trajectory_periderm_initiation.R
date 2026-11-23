# ============================================================
# Script: 19_periderm_initiation_trajectory_core.R
# Description: Subset phellogen and sub-epidermis cells, convert to 
#              Monocle CDS, and construct a non-branched trajectory.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
})

set.seed(123)

# Setup directories
data_dir <- "data"
outdir   <- "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

N_USE <- 80
MIN_CELLS_EXPRESSED <- 10

# ============================================================
# 1. Extract cells and merge
# ============================================================
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# Subset target cells
phellogen <- subset(seurat_object, subcell_type %in% c("Phellogen") & orig.ident %in% c("EPI"))
subepidermis <- subset(seurat_object, subcell_type %in% c("Subepidermis"))

phellogen$orig_cluster <- "phellogen"
subepidermis$orig_cluster <- "sub-epiderm"

# Merge
periderm <- merge(subepidermis, y = phellogen, add.cell.ids = c("sub-epiderm", "phellogen"))
periderm <- JoinLayers(periderm)
periderm$cell_type <- factor(periderm$orig_cluster)
periderm$group <- periderm$cell_type

saveRDS(periderm, file = file.path(data_dir, "Periderm_initiation.RDS"))

# ============================================================
# 2. Convert to Monocle CDS
# ============================================================
cat("Constructing Monocle CDS...\n")
DefaultAssay(periderm) <- 'RNA'
data <- periderm@assays$RNA$counts

pd <- new('AnnotatedDataFrame', data = periderm@meta.data)
fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(data, phenoData = pd, featureData = fd,
                      lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# ============================================================
# 3. Differential Analysis
# ============================================================
cat("Running differentialGeneTest...\n")
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= MIN_CELLS_EXPRESSED))
expressed_genes_clean <- expressed_genes[!grepl("GMtg|GPtg", expressed_genes)]

diff <- differentialGeneTest(cds[expressed_genes_clean, ], fullModelFormulaStr = "~cell_type", cores = 4)
diff <- na.omit(diff)
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval), ]
write.csv(deg, file.path(outdir, "diff_deg_initiation.csv"), row.names = FALSE)

# ============================================================
# 4. Construct Trajectory (DDRTree & orderCells)
# ============================================================
cat("Constructing trajectory...\n")
ordergene <- head(unique(as.character(deg$gene_short_name)), N_USE)
ordergene <- intersect(ordergene, rownames(cds))

cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', verbose = F)
cds <- orderCells(cds)

saveRDS(cds, file.path(outdir, paste0("cds_initiation_top", N_USE, ".rds")))
cat("Core trajectory analysis completed!\n")