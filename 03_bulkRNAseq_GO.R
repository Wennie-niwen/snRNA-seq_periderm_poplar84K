# ============================================================
# Script: 03_bulkRNAseq_GO.R
# Description: TCseq clustering of genes shared between EI and SI that were upregulated at each of the three developmental stages
# Associated Figures: Figure 1D
# ============================================================

library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(data.table)

data_dir  <- "data"
outdir    <- 'results'
outdir2   <- 'plots'

if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if(!dir.exists(outdir2)) dir.create(outdir2, recursive = TRUE)

# Read 84K poplar GO annotation
# Note: Please place "go_data_anno.txt" in the data directory
go_data <- read.table(file.path(data_dir, "go_data_anno.txt"), header = T)

# Remove any residual \r in the columns
go_data[] <- lapply(go_data, function(x) {
  if (is.character(x)) gsub("\r", "", x) else x
})

cat("GO column names:\n"); print(colnames(go_data))

go_bp <- go_data[go_data$ONTOLOGY == "BP", ]
term2gene <- go_bp[, c("GO", "GID")]
term2name <- go_bp[, c("GO", "TERM")]

file <- fread(file.path(outdir, 'Final_Gene_Group_Assignment.txt'))
print(table(file$group))

group_names <- c("Stage I", "Stage II", "Stage III")

for (case in group_names) {
  
  cat("\n===== GO Enrichment:", case, "=====\n")
  genelist <- file$gene[file$group == case]
  cat("Number of genes:", length(genelist), "\n")
  
  ego <- tryCatch({
    enricher(gene          = genelist,
             TERM2GENE     = term2gene,
             TERM2NAME     = term2name,
             pvalueCutoff  = 0.05,
             pAdjustMethod = "BH",
             minGSSize     = 5,
             maxGSSize     = 500)
  }, error = function(e) {
    cat("  [Error]", conditionMessage(e), "\n"); NULL
  })
  
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("  No significant enrichment for", case, "\n"); next
  }
  
  term <- as.data.frame(ego)
  saveRDS(term, file.path(outdir, paste0(case, '_terms_GO.rds')))
  write.csv(term, file.path(outdir, paste0(case, '_terms_GO.csv')),
            row.names = FALSE)
}

# Colors (consistent with the trend plot)
group_colors <- c(
  "Stage I"   = "#56B4E9",
  "Stage II"  = "#F79431",
  "Stage III" = "#F08A8A"
)

# ============================================================
# Manually specify the GO IDs to display for each group
# ============================================================
go_pick <- list(
  "Stage I" = c(
    'GO:0009992',
    'GO:0044550',
    'GO:0010143',
    'GO:0006833',
    'GO:0042335'
  ),
  "Stage II" = c(
    'GO:0043043',
    'GO:0042255',
    'GO:0070925',
    'GO:0031571',
    'GO:0044819'
  ),
  "Stage III" = c(
    'GO:0009698',
    'GO:0000278',
    'GO:0006631',
    'GO:0010345',
    'GO:0009809'
  )
)

# ============================================================
# Loop for plotting
# ============================================================
for (case in group_names) {
  
  rds_file <- file.path(outdir, paste0(case, '_terms_GO.rds'))
  if (!file.exists(rds_file)) {
    cat("Skipped (no results):", case, "\n"); next
  }
  
  term <- readRDS(rds_file)
  if (is.null(term) || nrow(term) == 0) {
    cat("Skipped (empty results):", case, "\n"); next
  }
  
  # Filter by the specified GO list
  pick <- go_pick[[case]]
  
  if (length(pick) == 0) {
    cat("Skipped (empty list):", case, "\n"); next
  }
  
  data_plot <- term %>%
    filter(ID %in% pick) %>%
    mutate(Description = fct_reorder(Description, -log10(p.adjust)))
  
  cat("\n", case, ": List has", length(pick), "GOs, found", nrow(data_plot), "in results\n")
  
  if (nrow(data_plot) == 0) {
    cat("  Skipped (None of the specified GOs were enriched)\n"); next
  }
  
  # Print the actual terms displayed for verification
  print(data_plot[, c("ID", "Description", "p.adjust")])
  
  # Check which IDs from the list were not found
  missing_ids <- setdiff(pick, data_plot$ID)
  if (length(missing_ids) > 0) {
    cat("  [WARNING] Missing GO IDs:", paste(missing_ids, collapse = ", "), "\n")
  }
  
  # x-axis range
  x_max <- max(-log10(data_plot$p.adjust))
  
  p_bar <- ggplot(data_plot,
                  aes(x = -log10(p.adjust),
                      y = Description,
                      fill = case)) +
    geom_col(width = 0.7, color = NA, linewidth = 0.3) +
    geom_text(aes(x = -log10(p.adjust) * 0.03,
                  y = Description,
                  label = Description),
              size = 5, hjust = 0, color = "black") +
    scale_fill_manual(values = setNames(group_colors[case], case),
                      name = NULL) +
    scale_x_continuous(limits = c(0, x_max),
                       expand = c(0, 0)) +
    labs(x = "-log10(p.adjust)", y = NULL, title = case) +
    theme_classic(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12, face = "bold"),
      axis.text.y     = element_blank(),
      axis.ticks.y    = element_blank(),
      axis.title.y    = element_blank(),
      axis.text.x     = element_text(colour = 'black', size = 9),
      axis.line.x     = element_line(colour = 'black', linewidth = 1),
      axis.ticks.x    = element_line(colour = 'black', linewidth = 1),
      axis.title.x    = element_text(colour = 'black', size = 15),
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 9),
      panel.grid      = element_blank(),
      plot.margin     = margin(5, 10, 5, 5)
    )
  
  # Fixed: Removed the 'GO_bar' folder path issue, saving directly to outdir2
  ggsave(file.path(outdir2, paste0(case, '_GO_bar.pdf')),
         plot = p_bar, width = 5, height = 3)
  
  cat("Completed:", case, "\n")
}
