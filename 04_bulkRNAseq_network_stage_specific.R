# ============================================================
# Script: 04_bulkRNAseq_network_stage_specific.R
# Description: . Gene Ontology (GO) enrichment networks of Bulk RNA-seq across the three periderm developmental stages
# Associated Figures: Figure S2
# ============================================================

library(ggplot2)
library(ggraph)
library(igraph)
library(dplyr)
library(patchwork) # 用于拼接多张图

data_dir  <- "data"
outdir2   <- "plots"

if(!dir.exists(outdir2)) dir.create(outdir2, recursive = TRUE)

# 1. Read network data
df <- read.csv(file.path(data_dir, "network_data_group.csv"), stringsAsFactors = FALSE)

# 2. Define a function to draw the circular network for a specific group
draw_network <- function(data, group_name, main_color) {
  # Create igraph object for this specific group
  g <- graph_from_data_frame(data[, c("Source", "Target")], directed = FALSE)
  
  # Assign node types
  V(g)$type <- ifelse(V(g)$name %in% data$Target, "Gene", "GO Term")
  
  # Plot
  p <- ggraph(g, layout = 'linear', circular = TRUE) + 
    geom_edge_arc(aes(color = from), curvature = 0.3, width = 0.8, alpha = 0.6) + 
    geom_node_point(aes(color = type, size = type)) +
    geom_node_text(aes(label = name), 
                   size = 3.5, 
                   nudge_x = 0.15,  
                   nudge_y = 0.15,  
                   repel = FALSE) + 
    
    # Use main_color for Genes and a neutral color for GO terms
    scale_color_manual(
      name = "Node Type",
      values = c("Gene" = main_color, "GO Term" = "#9FA3A8"), 
      breaks = c("Gene", "GO Term")
    ) +
    scale_size_manual(values = c("Gene" = 3, "GO Term" = 6), guide = "none") +
    coord_fixed(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5), clip = "off") +
    theme_void() +
    labs(title = group_name) +
    theme(
      legend.position = "none", # Hide legend for individual plots, we will make a global one
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = main_color)
    )
  
  return(p)
}

# 3. Generate the three plots with specific colors matching your previous style
p1 <- draw_network(df %>% filter(Group == "Stage I"), "Stage I", "#56B4E9")
p2 <- draw_network(df %>% filter(Group == "Stage II"), "Stage II", "#F79431")
p3 <- draw_network(df %>% filter(Group == "Stage III"), "Stage III", "#D44637")

# 4. Combine them side by side
# Using patchwork to combine. Adjust widths if one stage has too many nodes.
final_plot <- p1 + p2 + p3 + 
  plot_layout(ncol = 3) +
  plot_annotation(
    title = 'Stage-Specific Gene-Concept Networks',
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, face = "bold"))
  )

# 5. Save the plot
ggsave(file.path(outdir2, 'Network_Stage_Specific_Combined.pdf'), 
       plot = final_plot, width = 18, height = 6)

cat("Combined network plots saved successfully.\n")
