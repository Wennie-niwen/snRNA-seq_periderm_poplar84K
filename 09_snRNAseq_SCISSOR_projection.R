# ============================================================
# Script: 09_snRNAseq_SCISSOR_projection.R
# Description: Projection of bulk RNA-seq phenotypic data onto the snRNA-seq UMAP using Scissor
# Associated Figures: Figure S4E
# ============================================================

set.seed(11)
packageVersion("Scissor")

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(edgeR)
  library(Scissor); library(ggplot2); library(patchwork)
})

OUTDIR <- 'results'
DATA_DIR <- file.path(OUTDIR, "DATA"); PLOTS_DIR <- file.path(OUTDIR, "PLOTS")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

# 1. Read data needed------------
BULK_FILE <- "data/SI_raw_counts.txt"
seurat_object <- readRDS("data/Periderm.RDS")

UMAP_REDUCTION <- "umap"
SAMPLE_IDS <- c(paste0("SI_S1_R", 1:3), paste0("SI_S3b_R", 1:3))
SCISSOR_ALPHA <- NULL; SCISSOR_CUTOFF <- 0.2

DefaultAssay(seurat_object) <- "RNA"
print(seurat_object); print(table(seurat_object$cell_type))
seurat_object[["RNA"]] <- as(seurat_object[["RNA"]], Class = "Assay")
NEIGHBOR_REDUCTION <- Reductions(seurat_object)[1]
emb <- Embeddings(seurat_object, reduction = NEIGHBOR_REDUCTION)
dims_use <- seq_len(min(30, ncol(emb)))
seurat_object <- FindNeighbors(seurat_object, reduction = NEIGHBOR_REDUCTION, dims = dims_use, graph.name = c("RNA_nn", "RNA_snn"), verbose = FALSE)

bulk_table <- read.delim(BULK_FILE, header = TRUE, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
print(dim(bulk_table)); print(colnames(bulk_table))

sc_features <- rownames(seurat_object[["RNA"]])
id_candidates <- intersect(c("gene", "Ptr_gene_id", "Ath_gene_id", "name"), colnames(bulk_table))
overlap_number <- sapply(id_candidates, function(id_col) {
  ids <- trimws(as.character(bulk_table[[id_col]]))
  sum(!is.na(ids) & ids != "" & ids %in% sc_features)
})
id_check <- data.frame(ID_column = names(overlap_number), Overlap_genes = as.numeric(overlap_number))
print(id_check)
write.csv(id_check, file.path(DATA_DIR, "gene_ID_overlap.csv"), row.names = FALSE)
best_id_col <- names(which.max(overlap_number))
cat("\nSelected Gene ID column:", best_id_col, "\n")

gene_ids <- trimws(as.character(bulk_table[[best_id_col]]))
count_df <- lapply(bulk_table[, SAMPLE_IDS, drop = FALSE], function(x) as.numeric(as.character(x)))
count_mat <- do.call(cbind, count_df); colnames(count_mat) <- SAMPLE_IDS
valid_gene <- (!is.na(gene_ids) & gene_ids != "")
count_mat <- count_mat[valid_gene, , drop = FALSE]; gene_ids <- gene_ids[valid_gene]
rownames(count_mat) <- gene_ids; storage.mode(count_mat) <- "numeric"
count_mat <- rowsum(count_mat, group = gene_ids, reorder = FALSE)

group <- factor(ifelse(grepl("^SI_S1_", colnames(count_mat)), "SI_S1", "SI_S3b"), levels = c("SI_S1", "SI_S3b"))
print(table(group))
dge <- DGEList(counts = count_mat, group = group)
keep_genes <- filterByExpr(dge, group = group)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
bulk_dataset <- cpm(dge, log = TRUE, prior.count = 1)
print(dim(bulk_dataset))
write.csv(bulk_dataset, file.path(DATA_DIR, "Bulk_logCPM.csv"))

common_genes <- intersect(rownames(seurat_object[["RNA"]]), rownames(bulk_dataset))
cat("\nNumber of common genes:", length(common_genes), "\n")
bulk_dataset <- bulk_dataset[common_genes, , drop = FALSE]
sc_dataset <- subset(seurat_object, features = common_genes)
common_genes_final <- rownames(sc_dataset[["RNA"]])
bulk_dataset <- bulk_dataset[common_genes_final, , drop = FALSE]

phenotype <- as.integer(group == "SI_S3b"); names(phenotype) <- colnames(bulk_dataset)
tag <- c("SI_S1", "SI_S3b")
print(data.frame(Sample = names(phenotype), Group = group, Phenotype = phenotype))

infos <- Scissor::Scissor(bulk_dataset = bulk_dataset, sc_dataset = sc_dataset, phenotype = phenotype, tag = tag, alpha = SCISSOR_ALPHA, cutoff = SCISSOR_CUTOFF, family = "binomial", Save_file = file.path(DATA_DIR, "Scissor_inputs.RData"))

all_cells <- colnames(sc_dataset)
get_scissor_cells <- function(x, all_cells) {
  x <- unlist(x, use.names = FALSE); x <- as.character(x)
  idx <- suppressWarnings(as.integer(x)); valid_idx <- !is.na(idx) & idx >= 1 & idx <= length(all_cells)
  x[valid_idx] <- all_cells[idx[valid_idx]]; intersect(x, all_cells)
}
Scissor_pos <- get_scissor_cells(infos$Scissor_pos, all_cells)
Scissor_neg <- get_scissor_cells(infos$Scissor_neg, all_cells)
cat("\nNumber of Scissor+ cells:\n"); print(length(Scissor_pos))
cat("\nNumber of Scissor- cells:\n"); print(length(Scissor_neg))
cat("\nTotal number of cells:\n"); print(length(all_cells))

scissor_group <- setNames(rep("Other cells", length(all_cells)), all_cells)
scissor_group[Scissor_pos] <- "SI_S3b-associated"
scissor_group[Scissor_neg] <- "SI_S1-associated"
Scissor_final <- seurat_object
Scissor_final$Scissor_group <- factor(scissor_group[colnames(Scissor_final)], levels = c("Other cells", "SI_S1-associated", "SI_S3b-associated"))
cat("\nFinal Scissor classification:\n"); print(table(Scissor_final$Scissor_group))
cat("\nScissor classification across cell types:\n")
celltype_summary <- as.data.frame(table(Cell_type = Scissor_final$cell_type, Scissor_group = Scissor_final$Scissor_group))
print(celltype_summary)
write.csv(celltype_summary, file.path(DATA_DIR, "Scissor_celltype_summary.csv"), row.names = FALSE)

p_scissor <- DimPlot(Scissor_final, reduction = UMAP_REDUCTION, group.by = "Scissor_group",
                     cols = c("Other cells" = "grey80", "SI_S1-associated" = "#ABD58A", "SI_S3b-associated" = "#1F71AA"),
                     raster = FALSE, pt.size = 0.6, order = c("SI_S1-associated", "SI_S3b-associated")) +
  labs(x = "UMAP 1", y = "UMAP 2") + ggtitle("Scissor") + theme_classic() +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(PLOTS_DIR, "Scissor_UMAP.pdf"), p_scissor, width = 4, height = 3)

umap_embedding <- Embeddings(Scissor_final, reduction = UMAP_REDUCTION)
umap_embedding <- umap_embedding[colnames(Scissor_final), , drop = FALSE]
umap_df <- data.frame(Cell = rownames(umap_embedding), UMAP_1 = umap_embedding[, 1], UMAP_2 = umap_embedding[, 2],
                      Cell_type = Scissor_final$cell_type, Scissor_group = Scissor_final$Scissor_group)
write.csv(umap_df, file.path(DATA_DIR, "Scissor_UMAP_coordinates.csv"), row.names = FALSE)

saveRDS(infos, file.path(DATA_DIR, "Scissor_result.rds"))
saveRDS(Scissor_final, file.path(DATA_DIR, "Scissor_final_seurat.rds"))