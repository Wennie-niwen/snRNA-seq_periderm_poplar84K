# ============================================================
# Script: 16_snRNAseq_network_plot.R
# # Description: Network plot showed function and related genes
# Associated Figures: Figure 6B
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(igraph)
  library(scatterpie)
  library(ggforce)
  library(RColorBrewer)
})

# Setup directories
plot_dir <- "plots"
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# ============================================================
# 1. Embed raw data directly
# ============================================================
cat("Loading embedded gene-pathway data...\n")
df_long <- tribble(
  ~Pathway, ~Gene, ~number,
  # Cell wall synthesis
  "Cell wall synthesis and modification", "XTH5", 20,
  "Cell wall synthesis and modification", "PAxG03Gg0161800", 20,
  "Cell wall synthesis and modification", "PAxG05Gg0179500", 20,
  "Cell wall synthesis and modification", "XTR6", 20,
  "Cell wall synthesis and modification", "GAE6", 20,
  
  # Protective substance
  "Protective substance biosynthesis", "KCS2", 20,
  "Protective substance biosynthesis", "KCS11", 20,
  "Protective substance biosynthesis", "KCS6", 20,
  "Protective substance biosynthesis", "PRX52", 20,
  "Protective substance biosynthesis", "GPAT8", 20,
  "Protective substance biosynthesis", "LP1", 20,
  
  # Photosynthesis
  "Photosynthesis", "PSBO1", 20,
  "Photosynthesis", "PSBO2", 20,
  "Photosynthesis", "LHCA3", 20
)

# numbers represent dot size, 

# ============================================================
# 2. Data preparation
# ============================================================
gene_path <- df_long %>%
  distinct(Pathway, Gene) %>%
  group_by(Gene) %>%
  mutate(k = n()) %>%
  ungroup() %>%
  mutate(prop = 1 / k)

pathways <- sort(unique(df_long$Pathway))

wide <- gene_path %>%
  select(Gene, Pathway, prop) %>%
  pivot_wider(names_from = Pathway, values_from = prop, values_fill = 0)

gene_size <- df_long %>%
  distinct(Gene, number)

nodes <- wide %>%
  left_join(gene_size, by = "Gene")

# Construct edges
edge_list <- df_long %>%
  distinct(Pathway, Gene) %>%
  inner_join(distinct(df_long, Pathway, Gene), by = "Pathway", suffix = c(".x", ".y")) %>%
  filter(Gene.x < Gene.y) %>%
  select(from = Gene.x, to = Gene.y) %>%
  distinct()

# Create igraph object
g <- graph_from_data_frame(edge_list, vertices = nodes$Gene, directed = FALSE)

# Layout
set.seed(116)
layout_mat <- layout_with_fr(g)
layout_df  <- as.data.frame(layout_mat)
colnames(layout_df) <- c("x", "y")
layout_df$Gene <- V(g)$name

nodes_coords <- nodes %>%
  left_join(layout_df, by = "Gene")

edges_plot <- edge_list %>%
  left_join(layout_df %>% rename(from_x = x, from_y = y), by = c("from" = "Gene")) %>%
  left_join(layout_df %>% rename(to_x = x,   to_y = y),   by = c("to"   = "Gene")) %>%
  select(from, to, from_x, from_y, to_x, to_y)

# ============================================================
# 3. Colors & Pie sizing
# ============================================================
n_path <- length(pathways)
if(!requireNamespace("RColorBrewer", quietly = TRUE)) {
  pal <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(min(12, n_path), "Set3"))(n_path)
} else {
  pal <- RColorBrewer::brewer.pal(min(max(3, n_path), 12), "Set3")
  if(n_path > length(pal)) pal <- colorRampPalette(pal)(n_path)
}
names(pal) <- pathways

pie_df <- nodes_coords %>%
  select(Gene, x, y, number, all_of(pathways))

max_radius <- 0.3
pie_df <- pie_df %>%
  mutate(r = sqrt(number)) %>%
  mutate(r = r / max(r) * max_radius)

# ============================================================
# 4. Plotting
# ============================================================
cat("Plotting pathway pie network...\n")
p <- ggplot() +
  geom_segment(
    data = edges_plot,
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
    color = "grey75", alpha = 0.6, linewidth = 0.4 # Changed 'size' to 'linewidth' to avoid deprecation warning
  ) +
  geom_scatterpie(
    data = pie_df,
    aes(x = x, y = y, r = r),
    cols = pathways,
    color = "lightgrey",  
    lwd = 0.3,            
    alpha = 1
  ) +
  scale_fill_manual(
    name = "Pathway",
    values = pal
  ) +
  theme_void() +
  coord_equal(clip = "off") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20) 
  ) +
  geom_text(data = pie_df, size = 3,
            aes(x = x, y = y - r - 0.02, label = Gene))

# Save plot
ggsave(file.path(plot_dir, "hvg_PG_circle.pdf"), p, height = 5, width = 8)

cat("Network plot saved successfully!\n")
