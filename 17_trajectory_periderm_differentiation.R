# ============================================================
# Script: 17_trajectory_periderm_differentiation.R
# Description: Subset periderm cells, convert to Monocle CDS, 
#              and run differential & BEAM analysis.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle)
  library(dplyr)
  library(data.table)
})

set.seed(123)

# Set up directories
data_dir <- "data"
outdir   <- "results"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Parameters
N_USE <- 90
BRANCH_POINT <- 1
MIN_CELLS_EXPRESSED <- 10
BEAM_QVAL_CUTOFF <- 0.05

# ============================================================
# 1. Extract target cells and merge
# ============================================================
cat("Loading Seurat object...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))

# Subset the three cell types
phellem <- subset(seurat_object, subcell_type %in% c("Phellem") & orig.ident %in% c("PER"))
phellogen <- subset(seurat_object, subcell_type %in% c("Phellogen") & orig.ident %in% c("PER"))
phelloderm <- subset(seurat_object, cell_type %in% c("Phelloderm") & orig.ident %in% c("PER"))

# Merge
periderm <- merge(phellem, y = c(phellogen, phelloderm), add.cell.ids = c("phellem", "phellogen", "phelloderm"))
periderm <- JoinLayers(periderm)

# Define cell type from barcode prefix
periderm$cell_type <- sub("_.*", "", colnames(periderm))
periderm$cell_type <- factor(periderm$cell_type)
periderm$group <- periderm$cell_type

DefaultAssay(periderm) <- 'RNA'
saveRDS(periderm, file = file.path(data_dir, "periderm_differentiation.RDS"))

# ============================================================
# 2. Convert to Monocle CDS
# ============================================================
cat("Constructing Monocle CDS...\n")
data <- periderm@assays$RNA$counts
if (!inherits(data, "dgCMatrix")) data <- as(data, "CsparseMatrix")

pd <- new('AnnotatedDataFrame', data = periderm@meta.data)
fData <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
fd <- new('AnnotatedDataFrame', data = fData)

cds <- newCellDataSet(
  data,
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

# ============================================================
# 3. QC & DE analysis
# ============================================================
cat("QC and filtering...\n")
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
expressed_genes_clean <- expressed_genes[!grepl("GMtg|GPtg", expressed_genes)]

saveRDS(cds, file.path(outdir, "cds_base.rds"))

cat("Running differentialGeneTest...\n")
diff <- differentialGeneTest(cds[expressed_genes_clean, ], fullModelFormulaStr = "~cell_type", cores = 4)
diff <- na.omit(diff)
deg <- subset(diff, qval < 0.01)
deg <- deg[order(deg$qval), ]
write.csv(deg, file.path(outdir, "diff_deg.csv"), row.names = FALSE)

# ============================================================
# 4. Construct Trajectory (reduceDimension & orderCells)
# ============================================================
cat("Constructing trajectory...\n")
ordergene <- head(unique(as.character(deg$gene_short_name)), N_USE)
ordergene <- intersect(ordergene, rownames(cds))

cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, reduction_method = 'DDRTree', verbose = F)
cds <- orderCells(cds)

saveRDS(cds, file.path(outdir, paste0("cds_top", N_USE, ".rds")))

# ============================================================
# 5. BEAM Analysis (Branch Expression Analysis)
# ============================================================
cat("Running BEAM analysis (this may take a while)...\n")
# Get all eligible genes
fd <- fData(cds)
candidate_ids <- rownames(fd)
if ("num_cells_expressed" %in% colnames(fd)) {
  candidate_ids <- candidate_ids[fd$num_cells_expressed >= MIN_CELLS_EXPRESSED]
}
candidate_ids <- candidate_ids[!grepl("GMtg|GPtg", fd[candidate_ids, "gene_short_name"])]

BEAM_res <- tryCatch({
  BEAM(cds[candidate_ids, ], branch_point = BRANCH_POINT, cores = 1)
}, error = function(e) {
  stop("BEAM analysis failed: ", conditionMessage(e))
})

saveRDS(BEAM_res, file.path(outdir, "BEAM_res_all_candidates.rds"))
BEAM_sig <- subset(BEAM_res, qval < BEAM_QVAL_CUTOFF)
write.csv(BEAM_sig, file.path(outdir, "BEAM_significant_genes.csv"), row.names = FALSE)

cat("Core trajectory analysis completed!\n")