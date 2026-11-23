# ============================================================
# Script: 01_bulkRNAseq_TCSeqCluster.R
# Description: TCseq clustering of genes shared between EI and SI that were upregulated at each of the three developmental stages
# Associated Figures: Figure 1D
# ============================================================

library(tidyverse)
library(data.table)
library(matrixStats)

data_dir  <- "data"
outdir  <- 'results'
outdir2 <- 'plots'

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if(!dir.exists(outdir2)) dir.create(outdir2, recursive = TRUE)

## 1. Input data ----
# Create the metadata dataframe for bulk RNA-seq samples
metadata <- data.frame(
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
  group = c(
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 6),
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 3),
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 3)
  ),
  annotation = c(
    rep("epiderm", 3), rep("peridermInitiation", 3), 
    rep("peridermMaintanance", 3), rep("Stage3", 3),
    rep("epiderm", 3), rep("peridermInitiation", 3), 
    rep("peridermMaintanance", 3),
    rep("PVC0", 3), rep("PVC1", 3), rep("PVC2", 3)
  ),
  batch = c(
    rep("batch2", 12), rep("batch1", 18)
  ),
  group_2 = c(
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 6),
    rep("Stage I", 3), rep("Stage II", 3), rep("Stage III", 3),
    rep("Vascular", 9)
  ),
  stringsAsFactors = FALSE 
)

# read TPM (Paths kept exactly as you requested)
TPM1 <- read.table(file.path(data_dir, "SI_TPM_matrix.txt"), header = T)
TPM2 <- read.table(file.path(data_dir, "EI_TPM_matrix.txt"), header = T)
TPM <- merge(TPM1, TPM2, by="gene")
TPM <- TPM %>% column_to_rownames(var = "gene")
cat("TPM:", dim(TPM), "\n")

# Calculate mean by sample name groups
s_I   <- metadata$sample[metadata$group_2 == "Stage I"]
s_II <- metadata$sample[metadata$group_2 == "Stage II"]
s_III   <- metadata$sample[metadata$group_2 == "Stage III"]
s_V   <- metadata$sample[metadata$group_2 == "Vascular"]
cat("Stage I sample:",   paste(s_I,   collapse = ", "), "\n")
cat("Stage II sample:", paste(s_II, collapse = ", "), "\n")
cat("Stage III sample:",   paste(s_III,   collapse = ", "), "\n")
cat("Vascular sample:",   paste(s_V,   collapse = ", "), "\n")

mean_I   <- rowMeans(TPM[, s_I])
mean_II <- rowMeans(TPM[, s_II])
mean_III   <- rowMeans(TPM[, s_III])
mean_V   <- rowMeans(TPM[, s_V])
data <- data.frame(StageI   = mean_I,
                   StageII = mean_II,
                   StageIII   = mean_III,
                   Vascular = mean_V)

# Filter: At least one sample > 1
data2 <- data[apply(data, 1, function(row) any(row > 1)), ]
cat("Before filter: ", nrow(data), " After filter: ", nrow(data2), "\n")

write.table(data2,
            file.path(outdir, "TPM_group_mean1.txt"),
            quote = FALSE, row.names = TRUE, sep = "\t")

## 2. Normalization: Divide each gene by its maximum value across four time points -----
max_per_gene <- apply(data2, 1, max)
data_normalized <- sweep(data2, 1, max_per_gene, "/")

# Verify if the maximum value of each row is 1
print(apply(data_normalized, 1, max, na.rm = TRUE)[1:5])

# Save
data_normalized$gene <- rownames(data_normalized)
data_normalized <- data_normalized[, c(1:4)]

write.table(data_normalized,
            file.path(outdir, "TPM_group_mean1_normalized.txt"),
            quote = FALSE, row.names = T, sep = "\t")

grp_levels <- c("StageI", "StageII", "StageIII","Vascular")

# 2. Clustering (TCseq) ----------
library(TCseq)
set.seed(123)
data_normalized <- as.matrix(data_normalized)

# for (K in c(4, 6, 8, 10)) {
#   set.seed(123)
#   tcseq_cluster <- timeclust(data_normalized, algo = 'cm', k = K, standardize = F)
#   p <- timeclustplot(tcseq_cluster, value = 'raw', cols = 4)
#   ggsave(file.path(outdir2, paste0('TCseq_K', K, '.pdf')), p, width = 8, height = 6)
# }

cluster_num <- 8
tcseq_cluster <- timeclust(data_normalized, algo = 'cm', k = cluster_num, standardize = F)
clu <- data.frame(
  gene  = names(tcseq_cluster@cluster),
  GROUP = as.integer(tcseq_cluster@cluster)
)

# Write .kgg
fwrite(clu,
       file.path(outdir, paste0('TPM_group_mean1_normalized_K_G', cluster_num, '.kgg')),
       sep = '\t', quote = FALSE, col.names = TRUE)

# Write .CDT 
cdt <- data.frame(
  gene    = rownames(data_normalized),
  cluster = clu$GROUP,
  as.data.frame(data_normalized),
  check.names = FALSE
)
fwrite(cdt,
       file.path(outdir, paste0('TPM_group_mean1_normalized_K_G', cluster_num, '.CDT')),
       sep = '\t', quote = FALSE, col.names = TRUE)

#### 1. Read CDT and kgg for G8 ----
file <- fread(file.path(outdir, 'TPM_group_mean1_normalized_K_G8.CDT'))
file <- na.omit(file)
cluster_info <- fread(file.path(outdir, 'TPM_group_mean1_normalized_K_G8.kgg'))

cat("CDT dimensions:", dim(file), "\n")
cat("kgg dimensions:", dim(cluster_info), "\n")
print(table(cluster_info$GROUP))

#### 2. Merge ----
file <- merge(file, cluster_info, by = "gene", all.x = TRUE)
cat("Merged column names:\n"); print(colnames(file))

if ("GROUP" %in% colnames(file)) {
  colnames(file)[colnames(file) == "GROUP"] <- "cluster"
} else {
  colnames(file)[ncol(file)] <- "cluster"
}

write.table(file, file.path(outdir, 'file_G8_TPM.csv'),
            row.names = FALSE, sep = ',')

#### 3. Merge with TPM matrix ----
TPM <- read.table(file.path(outdir, 'TPM_group_mean1_normalized.txt'),
                  sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

data <- TPM[file$gene, ]
data$cluster <- file$cluster
data$gene    <- rownames(data)

write.table(data, file.path(outdir, 'TPM_mean1_G8.csv'),
            row.names = FALSE, sep = ',')

#### 4. Plot cluster mean trend ----
object <- aggregate(data[, grp_levels], list(cluster = data$cluster), FUN = mean)

pdata <- pivot_longer(object,
                      cols = all_of(grp_levels),
                      names_to = "time", values_to = "value")

pdata$cluster <- as.factor(pdata$cluster)
pdata$time    <- factor(pdata$time, levels = grp_levels)

# Calculate gene count per cluster
cluster_count <- table(data$cluster)
cat("Gene count per cluster:\n"); print(cluster_count)

# Construct legend labels with counts
cluster_labels <- paste0("C", names(cluster_count),
                         " (", as.integer(cluster_count), ")")
names(cluster_labels) <- names(cluster_count)

n_clu <- length(unique(pdata$cluster))
my_colors <- c("#EF1828","#56B4E9","#262A35","#9FA3A8",
               "#FFD401","#FE7501","#B4A7D6","#8BC34A")[1:n_clu]

p <- ggplot(pdata, aes(x = time, y = value,
                       color = cluster, group = cluster)) +
  geom_point(size = 2) +
  geom_smooth(se = FALSE, method = "loess", span = 1) +
  theme_classic() +
  scale_color_manual(values = my_colors,
                     labels = cluster_labels) +   
  theme(axis.title = element_text(color = "black"),
        axis.text  = element_text(color = "black")) +
  xlab("Stage") +
  ylab("Mean TPM") +
  labs(title = "Cluster scheme: G8",
       color = "Cluster")                        

ggsave(file.path(outdir2, 'Cluster_mean_G8.pdf'),
       plot = p, width = 3.5, height = 2)
cat("Scheme G8 completed\n")


# ============================================================
chosen     <- "G8"
grp_levels <- c("StageI", "StageII", "StageIII", "Vascular")

# ---- 1. Read the final selected cluster table ----
data2 <- read.table(file.path(outdir, 'TPM_group_mean1_normalized.txt'),
                    sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
data2$gene <- rownames(data2)

cluster_info <- fread(file.path(outdir, paste0('TPM_group_mean1_normalized_K_', chosen, '.kgg')))
data3 <- merge(data2, cluster_info, by = "gene", all.x = TRUE)
colnames(data3)[ncol(data3)] <- 'cluster'

# ---- 2. Prepare Sample-level TPM and Z-score matrices ----
# Reload raw TPM for absolute expression threshold checking
TPM_full <- merge(TPM1, TPM2, by="gene") %>% column_to_rownames(var = "gene")

# Map cluster to initial group name
# The stage-specific clusters should be defined manually based on the visualization in plots/Cluster_mean_G8.pdf.
# Note: 
# - Stage I: Cluster 6 (Peaks at Stage I)
# - Stage II: Cluster 2 (Contains stage II specific genes, but peaks in vascular. Filtered later against Vascular)
# - Stage III: Clusters 4 & 7 (Peaks and descending curve for Stage III)
# - Stable: All others
cluster2group <- c(
  "6" = "Stage I",
  "2" = "Stage II",
  "4" = "Stage III",
  "7" = "Stage III",
  "1" = "Stable",
  "5" = "Stable",
  "3" = "Stable",
  "8" = "Stable"
)
data3$group <- factor(cluster2group[as.character(data3$cluster)], 
                      levels = c("Stable", "Stage I", "Stage II", "Stage III"))

# Get EI sample subsets (assuming EI samples are in columns 13 to 30 of metadata)
meta_sample <- metadata
meta_sample2 <- meta_sample[c(13:30), ]
meta_sample2$group_2 <- factor(meta_sample2$group_2, 
                               levels = c("Vascular", "Stage I", "Stage II", "Stage III"))
meta_sample_sorted2 <- meta_sample2[order(meta_sample2$group_2), ]

# Extract EI TPM matrix for these samples
mat2 <- TPM_full[data3$gene, meta_sample_sorted2$sample]
mat2 <- as.matrix(mat2)

# Z-score calculation
mat_z2 <- (mat2 - rowMeans(mat2)) / rowSds(mat2)
mat_z2[is.na(mat_z2)] <- 0

# ---- 3. Filtering based on Vascular threshold ----
col_group2 <- factor(meta_sample_sorted2$group_2, 
                     levels = c("Vascular", "Stage I", "Stage II", "Stage III"))

idx_vasc <- which(col_group2 == "Vascular")
idx_s1   <- which(col_group2 == "Stage I")
idx_s2   <- which(col_group2 == "Stage II")
idx_s3   <- which(col_group2 == "Stage III")

# Group means (Z-score)
mean_vasc <- rowMeans(mat_z2[, idx_vasc, drop = FALSE])
mean_s1   <- rowMeans(mat_z2[, idx_s1,   drop = FALSE])
mean_s2   <- rowMeans(mat_z2[, idx_s2,   drop = FALSE])
mean_s3   <- rowMeans(mat_z2[, idx_s3,   drop = FALSE])

# Group means (TPM)
mean_tpm_s2 <- rowMeans(mat2[, idx_s2, drop = FALSE])
mean_tpm_s3 <- rowMeans(mat2[, idx_s3, drop = FALSE])

# Thresholds
delta   <- 0.25
tpm_cut <- 1

# Conditions for Stage II and Stage III
keep_s2 <- (mean_s2 - mean_vasc > delta) & (mean_s2 - mean_s1 > delta) & (mean_tpm_s2 > tpm_cut)
keep_s3 <- (mean_s3 - mean_vasc > delta) & (mean_s3 - mean_s1 > delta) & (mean_tpm_s3 > tpm_cut)

red_genes_s2 <- rownames(mat_z2)[keep_s2]
red_genes_s3 <- rownames(mat_z2)[keep_s3]

cat("Genes passing Stage II filter:",  length(red_genes_s2), "\n")
cat("Genes passing Stage III filter:", length(red_genes_s3), "\n")

# ---- 4. Replace original Stage II and Stage III ----
old_s2_genes <- rownames(mat_z2)[data3$group == "Stage II"]
old_s3_genes <- rownames(mat_z2)[data3$group == "Stage III"]

drop_s2 <- setdiff(old_s2_genes, red_genes_s2)
drop_s3 <- setdiff(old_s3_genes, red_genes_s3)
drop_genes <- union(drop_s2, drop_s3)

# Keep final genes
keep_final <- !(data3$gene %in% drop_genes)
data3_final <- data3[keep_final, ]

# Factorize
group_names <- c("Stable", "Stage I", "Stage II", "Stage III")
data3_final$group <- factor(data3_final$group, levels = group_names)

cat("\nFinal gene distribution after filtering:\n")
print(table(data3_final$group))

# ---- 5. Save final output tables ----
# Save the full list of genes and their group assignment
final_gene_table <- data3_final[, c("gene", "cluster", "group")]
write.table(final_gene_table, 
            file.path(outdir, 'Final_Gene_Group_Assignment.txt'), 
            quote = FALSE, row.names = FALSE, sep = "\t")

# Save the three specific lists as requested (for downstream usage if needed)
df_s1 <- data.frame(gene = data3_final$gene[data3_final$group == "Stage I"], stage = "Stage I")
df_s2 <- data.frame(gene = data3_final$gene[data3_final$group == "Stage II"], stage = "Stage II")
df_s3 <- data.frame(gene = data3_final$gene[data3_final$group == "Stage III"], stage = "Stage III")

write.csv(df_s1, file.path(outdir, "StageI_specific_genes.csv"), row.names = FALSE, quote = FALSE)
write.csv(df_s2, file.path(outdir, "StageII_specific_genes.csv"), row.names = FALSE, quote = FALSE)
write.csv(df_s3, file.path(outdir, "StageIII_specific_genes.csv"), row.names = FALSE, quote = FALSE)

# ---- 6. Plot group mean trend ----
# Prepare data for plotting
pdata_mean <- aggregate(data3_final[, grp_levels], 
                        list(group = data3_final$group), FUN = mean)

pdata <- pivot_longer(pdata_mean,
                      cols = all_of(grp_levels),
                      names_to = "time", values_to = "value")

pdata$group <- factor(pdata$group, levels = group_names)
pdata$time  <- factor(pdata$time,  levels = grp_levels)

# Check: should be 4 groups × 4 stages = 16 rows
cat("pdata rows:", nrow(pdata), "\n")

# Filter out Vascular for plotting
pdata <- pdata %>% filter(time != "Vascular")

# Update group names and colors to exclude Vascular
group_names_plot <- c("Stable", "Stage I", "Stage II", "Stage III")

# Colors and plotting
group_colors <- c(
  "Stable"   = "#9FA3A8",   
  "Stage I"  = "#56B4E9",   
  "Stage II" = "#F79431",   
  "Stage III"= "#D44637"    
)

p_group <- ggplot(pdata, aes(x = time, y = value,
                             color = group, group = group)) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 0.8) +
  theme_classic() +
  scale_color_manual(values = group_colors,
                     breaks = group_names) +
  theme(axis.title = element_text(color = "black"),
        axis.text  = element_text(color = "black")) +
  xlab("Stage") +
  ylab("Mean normalized expression") +
  labs(color = NULL)

ggsave(file.path(outdir2, paste0('Group_mean_', chosen, '_specific.pdf')),
       plot = p_group, width = 5.5, height = 3)
