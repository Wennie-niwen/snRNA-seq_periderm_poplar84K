# ============================================================
# Script: 05_bulkRNAseq_heatmap_key_genes_suberin_cutin_metabolism.R
# Description: Heatmap showing expression dynamics of cutin- and suberin-associated genes in epidermal and phellem cells.
# Associated Figures: Figure S8
# ============================================================

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(matrixStats)

data_dir  <- "data"
outdir    <- 'results'
outdir2   <- 'plots'

if(!dir.exists(outdir2)) dir.create(outdir2, recursive = TRUE)

# ============ 1. Read sample metadata and TPM matrix ============
# Read metadata defined manually to keep script independent
meta_sample <- data.frame(
  sample = c(
    "SI_S1_R1", "SI_S1_R2", "SI_S1_R3",
    "SI_S2_R1", "SI_S2_R2", "SI_S2_R3",
    "SI_S3a_R1", "SI_S3a_R2", "SI_S3a_R3",
    "SI_S3b_R1", "SI_S3b_R2", "SI_S3b_R3",
    "EI_S1_R1", "EI_S1_R2", "EI_S1_R3",
    "EI_S2_R1", "EI_S2_R2", "EI_S2_R3",
    "EI_S3_R1", "EI_S3_R2", "EI_S3_R3",
    "EI_S1_V_R1", "EI_S1_V_R2", "EI_S1_V_R3",
    "EI_S2_V_R1", "EI_S2_V_R2", "EI_S2_V_R3",
    "EI_S3_V_R1", "EI_S3_V_R2", "EI_S3_V_R3"
  ),
  group_2 = c(
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 6),
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 3),
    rep("Vascular", 9)
  ),
  stringsAsFactors = FALSE 
)

# Read TPM files directly from GEO
TPM1 <- read.table(file.path(data_dir, "SI_TPM_matrix.txt"), header = T)
TPM2 <- read.table(file.path(data_dir, "EI_TPM_matrix.txt"), header = T)
TPM <- merge(TPM1, TPM2, by="gene") %>% column_to_rownames(var = "gene")

# ============ 2. Define target genes in specific biological order ============
# Ordered from top to bottom as it appears in the bubble plot
# Note: Heatmap draws from bottom to top, so we reverse (rev) the list.
genes_used <- rev(c(
  # 1. ABCG transporter (Top)
  'PAxG14Gg0063700', 'PAxG18Ag0068300', 'PAxG18Gg0071100', 'PAxG09Ag0041200',
  
  # 2. Apoplastic lipid transfer proteins
  'PAxG02Gg0045400', 'PAxG02Ag0046700', 'PAxG01Ag0047500', 'PAxG01Gg0046800',
  
  # 3. Extracellular polymerisation
  'PAxG07Ag0035000', 'PAxG07Gg0035100', 'PAxG12Ag0061100', 'PAxG12Gg0067800',
  'PAxG01Ag0181400', 'PAxG01Gg0172600', 'PAxG04Gg0042500', 'PAxG13Gg0090700',
  'PAxG15Gg0135100', 'PAxG15Ag0134000',
  
  # 4. Fatty acid ω-oxidation
  'PAxG01Ag0193400', 'PAxG01Gg0190900', 'PAxG07Gg0074800', 'PAxG05Ag0075200',
  
  # 5. Fatty acyl elongase (FAE) complex
  'PAxG08Gg0138500', 'PAxG10Gg0113400', 'PAxG10Ag0109900', 'PAxG06Gg0202700', 'PAxG06Ag0207400',
  
  # 6. Glycerylation (Bottom)
  'PAxG06Gg0092700', 'PAxG14Ag0066000'
))

# Check if these genes are in the matrix and filter
genes_used <- intersect(genes_used, rownames(TPM))
TPM_sample <- TPM[genes_used, ]
cat("Genes used for plotting:", dim(TPM_sample), "\n")

# ============ 3. Process SI and EI matrices with Z-score ============
# SI matrix
TPM_sample1 <- TPM_sample[, 1:12]
meta_sample1 <- meta_sample[1:12, ]
meta_sample1$group_2 <- factor(meta_sample1$group_2, levels = c("Stage I", "Stage II", "Stage III"))
meta_sample_sorted1 <- meta_sample1[order(meta_sample1$group_2), ]
mat1 <- TPM_sample1[, meta_sample_sorted1$sample]

mat_z1 <- (mat1 - rowMeans(mat1)) / rowSds(as.matrix(mat1))
mat_z1 <- as.matrix(mat_z1)
mat_z1[is.na(mat_z1)] <- 0

# EI matrix
TPM_sample2 <- TPM_sample[, 13:30]
meta_sample2 <- meta_sample[13:30, ]
meta_sample2$group_2 <- factor(meta_sample2$group_2, levels = c("Vascular", "Stage I", "Stage II", "Stage III"))
meta_sample_sorted2 <- meta_sample2[order(meta_sample2$group_2), ]
mat2 <- TPM_sample2[, meta_sample_sorted2$sample]

mat_z2 <- (mat2 - rowMeans(mat2)) / rowSds(as.matrix(mat2))
mat_z2 <- as.matrix(mat_z2)
mat_z2[is.na(mat_z2)] <- 0

# ============ 4. Color configuration ============
ht_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
col_fun   <- colorRamp2(seq(-2, 2, length.out = 100), ht_colors)

group_colors <- c(
  "Vascular"  = "lightgrey",
  "Stage I"   = "#56B4E9",
  "Stage II"  = "#F79431",
  "Stage III" = "#D44637"
)

# ============ 5. Calculate group means ============
# SI means
mat_z1_mean <- do.call(cbind, lapply(c("Stage I", "Stage II", "Stage III"), function(stage) {
  cols <- which(meta_sample_sorted1$group_2 == stage)
  rowMeans(mat_z1[, cols, drop = FALSE])
}))
colnames(mat_z1_mean) <- c("SI_Stage I", "SI_Stage II", "SI_Stage III")

# EI means
mat_z2_mean <- do.call(cbind, lapply(c("Vascular", "Stage I", "Stage II", "Stage III"), function(stage) {
  cols <- which(meta_sample_sorted2$group_2 == stage)
  rowMeans(mat_z2[, cols, drop = FALSE])
}))
colnames(mat_z2_mean) <- c("EI_Vascular", "EI_Stage I", "EI_Stage II", "EI_Stage III")

# Mean matrices column annotations
col_group1_mean <- factor(c("Stage I", "Stage II", "Stage III"), levels = c("Stage I", "Stage II", "Stage III"))
col_ha1_mean <- HeatmapAnnotation(
  Group = col_group1_mean,
  col   = list(Group = group_colors[levels(col_group1_mean)]),
  show_annotation_name = FALSE,
  show_legend = FALSE 
)

col_group2_mean <- factor(c("Vascular", "Stage I", "Stage II", "Stage III"), levels = c("Vascular", "Stage I", "Stage II", "Stage III"))
col_ha2_mean <- HeatmapAnnotation(
  Group = col_group2_mean,
  col   = list(Group = group_colors[levels(col_group2_mean)]),
  show_annotation_name = FALSE,
  show_legend = TRUE
)

# Row annotation based on gene categories
row_category <- case_when(
  grepl("KCS6|KCS10|GPAT8|HTH|ABCG15|LTPG1|GELP67|PEL3", genes_used) ~ "Epidermis-specific",
  grepl("KCS11|GPAT6|CYP86B1|ABCG2|ABCG32|LTPG5|GELP51|GELP38", genes_used) ~ "Phellem-specific",
  TRUE ~ "Other"
)
row_category <- factor(row_category, levels = c("Epidermis-specific", "Phellem-specific"))

row_ha <- rowAnnotation(
  Tissue = row_category,
  col = list(Tissue = c("Epidermis-specific" = "#2166AC", "Phellem-specific" = "#B2182B")),
  show_annotation_name = FALSE
)

# ============ 6. Draw Heatmaps ============
ht1_mean <- Heatmap(
  mat_z1_mean,
  name             = "Z-score (SI)",
  col              = col_fun,
  top_annotation   = col_ha1_mean,
  cluster_rows     = FALSE, 
  show_row_dend    = FALSE,
  cluster_columns  = FALSE, 
  column_names_rot = 45,    
  column_names_gp  = gpar(fontsize = 11, fontface = "bold"),
  show_row_names   = FALSE, 
  row_gap          = unit(2, "mm"),
  border           = TRUE,
  rect_gp          = gpar(col = "white", lwd = 1) 
)

ht2_mean <- Heatmap(
  mat_z2_mean,
  name             = "Z-score (EI)",
  col              = col_fun,
  left_annotation  = row_ha, 
  top_annotation   = col_ha2_mean,
  cluster_rows     = FALSE,
  show_row_dend    = FALSE,
  cluster_columns  = FALSE,
  column_names_rot = 45,
  column_names_gp  = gpar(fontsize = 11, fontface = "bold"),
  show_row_names   = TRUE,
  row_names_gp     = gpar(fontsize = 10, fontface = "italic"),
  row_gap          = unit(2, "mm"),
  border           = TRUE,
  rect_gp          = gpar(col = "white", lwd = 1)
)

# ============ 7. Save ============
pdf(file.path(outdir2, 'Heatmap_Specific_Isoforms_SI_EI_Mean.pdf'), width = 7, height = 6)
draw(ht1_mean + ht2_mean, merge_legend = TRUE)
dev.off()

cat("Heatmap for specific isoforms saved successfully.\n")
