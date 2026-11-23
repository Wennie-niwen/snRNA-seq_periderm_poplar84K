# ============================================================
# Script: 11_GO_enrichment_subpopulations.R
# Description: Identify DEGs for specific cell subpopulations, 
#              perform GO enrichment, and generate barplots.
# Associated Figures: Figure 3C, 3K, 4F, 4G, 5B
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(ggplot2)
  library(dplyr)
  library(forcats)
  library(Seurat)
})

# Setup directories
data_dir <- "data"
outdir   <- "result"
plot_dir <- "plot"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Load data
cat("Loading Seurat object and GO annotation...\n")
seurat_object <- readRDS(file.path(data_dir, "Periderm.RDS"))
go_data <- read.table(file.path(data_dir, "go_data_anno.txt"), header = T)
GOannotation <- split(go_data, with(go_data, ONTOLOGY, TERM))
term2gene <- GOannotation$BP[, c("GO", "GID")]
term2name <- GOannotation$BP[, c("GO", "TERM")]

# ============================================================
# Universal function for DEG & GO Enrichment
# ============================================================
run_go_analysis <- function(obj, ident_1, ident_2, group_by, fill_color, output_name, 
                            p_cutoff = 0.05, use_pvalue = FALSE, 
                            min_pct = 0.25, logfc = 0.58, only_pos = TRUE, 
                            top_n = 10, custom_go_ids = NULL, plot_width = 4, plot_height = 4) {
  
  cat("\n=== Processing:", ident_1, "vs", ident_2, "===\n")
  
  # 1. Find Markers
  Idents(obj) <- group_by
  markers <- FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2, 
                         min.pct = min_pct, logfc.threshold = logfc, only.pos = only_pos)
  markers$gene <- rownames(markers)
  
  # Filter by p-value or p-adjust
  if(use_pvalue) {
    markers <- subset(markers, p_val < p_cutoff)
  } else {
    markers <- subset(markers, p_val_adj < p_cutoff)
  }
  
  # 2. GO Enrichment
  ego <- enricher(gene = markers$gene, TERM2GENE = term2gene, TERM2NAME = term2name)
  ego_res <- ego@result
  
  # Save raw GO results
  write.csv(ego_res, file.path(outdir, paste0("GO_raw_", output_name, ".csv")), row.names = FALSE)
  
  # 3. Prepare Plot Data
  if(use_pvalue) {
    ego_res$`-log10(p)` <- -log10(ego_res$pvalue)
    sort_col <- "pvalue"
  } else {
    ego_res$`-log10(p)` <- -log10(ego_res$p.adjust)
    sort_col <- "p.adjust"
  }
  
  if(!is.null(custom_go_ids)) {
    # Use custom GO IDs
    plot_data <- ego_res %>% 
      filter(ID %in% custom_go_ids) %>% 
      arrange(!!sym(sort_col)) %>% 
      mutate(Description = fct_reorder(Description, `-log10(p)`))
  } else {
    # Use Top N
    plot_data <- ego_res %>% 
      arrange(!!sym(sort_col)) %>% 
      slice_head(n = top_n) %>% 
      mutate(Description = fct_reorder(Description, `-log10(p)`))
  }
  
  if(nrow(plot_data) == 0) {
    warning("No GO terms to plot for ", output_name)
    return(NULL)
  }
  
  # 4. Plot
  x_off <- max(plot_data$`-log10(p)`) * 0.05
  
  p <- ggplot(plot_data, aes(x = `-log10(p)`, y = Description)) +
    geom_col(fill = fill_color, width = 0.6) +
    geom_text(aes(x = x_off, y = Description, label = Description),
              hjust = 0, size = 3.5, color = "black") +
    scale_y_discrete(position = "right") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(x = ifelse(use_pvalue, "-log10(pvalue)", "-log10(p.adjust)"), y = NULL) +
    theme_bw() +
    theme(
      legend.position    = "none",
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x        = element_text(size = 10, color = "black"),
      axis.title.x       = element_text(size = 14, color = "black", face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_rect(color = "black", linewidth = 0.8),
      axis.ticks.x       = element_line(color = "black"),
      plot.margin        = margin(8, 12, 8, 8)
    )
  
  ggsave(file.path(plot_dir, paste0(output_name, ".pdf")), p, width = plot_width, height = plot_height)
  cat("Saved:", output_name, "\n")
}

# ============================================================
# 1. Phellem vs Other: Figure 3C
# ============================================================
seurat_object$comparison_1 <- ifelse(seurat_object$subcell_type %in% "Phellem", "Phellem", "Other")
run_go_analysis(seurat_object, "Phellem", "Other", "comparison_1", 
                fill_color = "#D57B70", output_name = "GO_Phellem_vs_Other", top_n = 10)

# ============================================================
# 2. Phellogen (vs Vascular Cambium): Figure 4F
# ============================================================
run_go_analysis(seurat_object, "Phellogen", "Vascular cambium", "subcell_type", 
                fill_color = "#F0AC54", output_name = "GO_Phellogen_vs_VC", 
                min_pct = 0.01, top_n = 10, plot_width = 3)

# ============================================================
# 3. Vascular Cambium (vs Phellogen): Figure 4G
# ============================================================
# Note: For reverse plot, you can add `+ scale_x_reverse()` inside the function, 
# but here we just create a standard plot.
run_go_analysis(seurat_object, "Vascular cambium", "Phellogen", "subcell_type", 
                fill_color = "#609895", output_name = "GO_VC_vs_Phellogen", 
                min_pct = 0.01, top_n = 10, plot_width = 3)

# ============================================================
# 4. Phelloderm vs Other: Figure 5B
# ============================================================
seurat_object$comparison_2 <- ifelse(seurat_object$cell_type %in% "Phellem" & seurat_object$orig.ident %in% c("EPI","PER"), "Phellogen", "Other")
run_go_analysis(seurat_object, "Phellogen", "Other", "comparison_2", 
                fill_color = "#1f90b4", output_name = "GO_Phelloderm_vs_Other", 
                top_n = 15, plot_width = 3)

# ============================================================
# 5. Dead Epidermis (RNS3 positive vs negative): Figure 3K
# ============================================================
cat("\n=== Processing Dead Epidermis ===\n")
epiderm <- subset(seurat_object, subcell_type %in% c("Epidermis"))
gene_expr <- LayerData(epiderm, assay = "RNA", layer = "data")["PAxG08Gg0069900", ]
epiderm$RNS3_status <- ifelse(gene_expr > 0, "RNS3_positive", "RNS3_negative")

my_go_ids <- c('GO:0005996', 'GO:0009698', 'GO:0019318', 'GO:0046173', 'GO:0009631',
               'GO:0044042', 'GO:0007033', 'GO:0071555', 'GO:0043067', 'GO:0006914')

# Note: For the dead epidermis, only down-regulated genes were selected (avg_log2FC < 0)
# We will handle this manually before calling the function
Idents(epiderm) <- epiderm$RNS3_status
markers <- FindMarkers(epiderm, ident.1 = "RNS3_positive", ident.2 = "RNS3_negative", 
                       min.pct = 0.01, logfc.threshold = 0.58)
markers$gene <- rownames(markers)
markers <- subset(markers, p_val < 0.05 & avg_log2FC < 0)

ego <- enricher(gene = markers$gene, TERM2GENE = term2gene, TERM2NAME = term2name)
ego_res <- ego@result
ego_res$`-log10(p)` <- -log10(ego_res$pvalue)

plot_data <- ego_res %>% 
  filter(ID %in% my_go_ids) %>% 
  arrange(pvalue) %>% 
  mutate(Description = fct_reorder(Description, `-log10(p)`))

x_off <- max(plot_data$`-log10(p)`) * 0.05

p <- ggplot(plot_data, aes(x = `-log10(p)`, y = Description)) +
  geom_col(fill = "#CCB674", width = 0.6) +
  geom_text(aes(x = x_off, y = Description, label = Description),
            hjust = 0, size = 3.5, color = "black") +
  scale_y_discrete(position = "right") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "-log10(pvalue)", y = NULL) +
  theme_bw() +
  theme(
    legend.position    = "none",
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    axis.title.y       = element_blank(),
    axis.text.x        = element_text(size = 10, color = "black"),
    axis.title.x       = element_text(size = 14, color = "black", face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", linewidth = 0.8),
    axis.ticks.x       = element_line(color = "black"),
    plot.margin        = margin(8, 12, 8, 8)
  )

ggsave(file.path(plot_dir, "GO_DeadEpidermis.pdf"), p, width = 3, height = 3)
cat("Saved: GO_DeadEpidermis.pdf\n")