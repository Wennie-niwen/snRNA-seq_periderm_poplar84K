# ============================================================
# Script: 06_snRNAseq_sample_correlation.R
# Description: The correlation of gene expression profiles between the two biological replicates for the EPI,
# PER, and PVC samples
# Associated Figures: Figure S4C
# ============================================================


library(ggplot2)
library(ggpubr)

plot_dir <- "plots"
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

seurat_object <- readRDS("data/Periderm.RDS")

# Extract metadata
# 1. Regex extraction: capture labels such as EPI-1, PER-2
# Logic: match EPI, PER, or PVC followed by "-" and a digit
seurat_object$Replicate <- sub(".*_(EPI-[12]|PER-[12]|PVC-[12])_.*", "\\1", colnames(seurat_object))

# 2. Check whether the six samples were successfully separated
table(seurat_object$Replicate)

per1_obj <- subset(seurat_object, subset = Replicate == "PVC-1")
per2_obj <- subset(seurat_object, subset = Replicate == "PVC-2")

per1_mean <- rowMeans(GetAssayData(per1_obj, layer = "data"))
per2_mean <- rowMeans(GetAssayData(per2_obj, layer = "data"))

df <- data.frame(Rep_1 = per1_mean, Rep_2 = per2_mean)

# Plot drawing
p <- ggplot(df, aes(x = Rep_1, y = Rep_2)) +
  geom_point(color = "#9B59B6", alpha = 0.3, size = 1) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 1), 
    axis.text = element_text(color = "black", size = 11),       
    axis.ticks = element_line(color = "black"),                
    axis.ticks.length = unit(0.2, "cm"),                      
    plot.title = element_blank()                              
  ) +
  scale_x_continuous(expand = c(0.02, 0)) + 
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(x = "PVC-1", y = "PVC-2")

print(p)

# 3. Save PDF
ggsave(file.path(plot_dir, "PVC_Correlation_Box.pdf"), p, width = 3, height = 3)
