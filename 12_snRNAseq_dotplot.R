# ============================================================
# Script: 12_snRNAseq_dotplot.R
# # Description: DotPlots for the periderm manuscript:
#   1) Simple marker DotPlots for Phellem / Phellogen / Phelloderm
#   2) Epidermis RNS3+ vs RNS3- grouped DotPlot (Death / Cell Wall / Stress)
#   3) Epidermis vs Phellem comparison, vertical faceted DotPlot
#   4) Cortex marker DotPlot (Collenchyma/Parenchyma/Subepidermis)
# Associated Figures: Figure 3B, 3L, 4D, 5D, 6C, S8
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

base_dir <- "data"
plot_dir <- "plots"
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ===== Load object; remove unresolved mixed-label cells and fill NAs =====
seurat_object <- readRDS(file.path(base_dir, "Periderm.RDS"))
seurat_object <- subset(seurat_object,
                        subset = subcell_type != "Phellem/cortex collenchyma")
seurat_object$subcell_type <- coalesce(seurat_object$subcell_type,
                                       seurat_object$cell_type)

# ============================================================
# 1. Simple marker DotPlots for the three periderm layers: Figure 3B, 4D, 5D
# ============================================================
all_markers <- data.frame(
  set = c(
    rep("Phellem", 5),
    rep("Phellogen", 4),
    rep("Phelloderm", 5)
  ),
  gene = c(
    "PAxG07Ag0002300", "PAxG02Gg0045300", "PAxG05Ag0141800",
    "PAxG12Ag0045600", "PAxG05Gg0059700",
    "PAxG01Ag0165100", "PAxG14Ag0101600", "PAxG17Ag0151800", "PAxG09Ag0004200",
    "PAxG04Ag0120700", "PAxG19Ag0045900", "PAxG19Ag0119200",
    "PAxG01Ag0009000", "PAxG10Ag0065800"
  ),
  name = c(
    "CYP716A1", "LTPG15", "MYB93", "NAC058", "MYB92",
    "CYCA2;4", "CSLD5", "TUA2", "CKS1",
    "MET1", "LHP1", "DDM1", "CMT3", "H1.1"
  ),
  stringsAsFactors = FALSE
)

my_levels <- c(
  "Phloem parenchyma", "Phloem precursor", "Starch sheath", "Phellem",
  "Epidermis", "Subepidermis",
  "Cortex collenchyma (EPI)", "Cortex collenchyma (PER)",
  "Vascular Cambium", "Phellogen", "Cortex parenchyma",
  "Companion cell", "Phloem fiber", "Phelloderm",
  "Sieve element", "Unknown"
)
seurat_object$subcell_type <- factor(seurat_object$subcell_type,
                                     levels = rev(my_levels))

for (layer in unique(all_markers$set)) {
  mk <- all_markers[all_markers$set == layer, ]
  features <- mk$gene
  id2sym   <- setNames(mk$name, mk$gene)
  
  missing <- setdiff(features, rownames(seurat_object))
  if (length(missing) > 0)
    warning("Missing genes in ", layer, ": ", paste(missing, collapse = ", "))
  
  p <- DotPlot(seurat_object, features = features, group.by = "subcell_type") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1,
                                 face = "italic", color = "black", size = 10),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title  = element_blank(),
      legend.text  = element_text(size = 9, color = "black"),
      legend.title = element_text(size = 9, color = "black"),
      axis.line   = element_line(color = "black"),
      axis.ticks  = element_line(color = "black")
    ) +
    scale_color_distiller(palette = "RdBu") +
    scale_x_discrete(labels = id2sym)
  
  ggsave(file.path(plot_dir, paste0("DotPlot_", layer, ".pdf")),
         plot = p, width = 5, height = 3)
}

# ============================================================
# 2. Epidermis: RNS3+ vs RNS3- (Death / Cell Wall / Stress): Figure 3L
# ============================================================
epiderm <- subset(seurat_object, subcell_type == "Epidermis")

# RNS3 status from normalized data layer
gene_expr <- LayerData(epiderm, assay = "RNA", layer = "data")["PAxG08Gg0069900", ]
epiderm$RNS3_status <- ifelse(gene_expr > 0, "RNS3_positive", "RNS3_negative")
epiderm$RNS3_status <- factor(epiderm$RNS3_status,
                              levels = c("RNS3_negative", "RNS3_positive"))
print(table(epiderm$RNS3_status))

genes_death <- c("PAxG12Gg0072900", "PAxG06Ag0114800", "PAxG02Gg0157500",
                 "PAxG10Gg0136100", "PAxG17Gg0260800", "PAxG02Ag0036400",
                 "PAxG15Ag0085400")
genes_wall <- c("PAxG06Ag0207400", "PAxG10Ag0109900", "PAxG16Gg0001100",
                "PAxG02Gg0060200", "PAxG13Ag0105000", "PAxG05Ag0075200",
                "PAxG09Ag0035200", "PAxG14Ag0125900")
genes_stress <- c("PAxG11Gg0129000", "PAxG01Ag0099400", "PAxG06Gg0219000",
                  "PAxG06Gg0250500", "PAxG09Gg0131600", "PAxG08Gg0153100")
all_features <- c(genes_death, genes_wall, genes_stress)

p_dot <- DotPlot(epiderm, features = all_features, group.by = "RNS3_status")
df_dot <- p_dot$data
df_dot$Category <- rep(
  c("Death & Autophagy", "Cell Wall & Wax", "Stress Response"),
  times = c(length(genes_death), length(genes_wall), length(genes_stress))
)

gene_labels <- c("VPS34", "LOS4", "ATG6", "ATG8h", "ATG5", "ATG9b", "CEP1",
                 "KCS10", "KCS6", "GAUT4", "CESA9", "KCS19", "CYP86B1", "CYP86A1",
                 "CER1", "GSTU19", "CBF4", "PFK3", "HK5", "CBF2", "GSTU7")
df_dot$features.plot <- factor(df_dot$features.plot,
                               levels = all_features, labels = gene_labels)

p_advanced <- ggplot(df_dot, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "grey", high = "#6E511E") +
  scale_size_continuous(range = c(0, 6)) +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", color = "black", size = 11),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title  = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.spacing = unit(0.1, "lines")
  )
ggsave(file.path(plot_dir, "DotPlot_epidermis_RNS3.pdf"),
       plot = p_advanced, width = 9, height = 4)

# ============================================================
# 3. Epidermis vs Phellem, vertical faceted DotPlot: Figure S8
# ============================================================
seurat_sub <- subset(seurat_object, subcell_type %in% c("Epidermis", "Phellem"))
seurat_sub$comparison <- seurat_sub$subcell_type
seurat_sub$comparison <- factor(seurat_sub$comparison,
                                levels = c("Epidermis", "Phellem"))
print(table(seurat_sub$comparison))

FAE <- c("PAxG08Gg0138500", "PAxG10Gg0113400", "PAxG10Ag0109900",
         "PAxG06Gg0202700", "PAxG06Ag0207400")
Glycerylation <- c("PAxG06Gg0092700", "PAxG14Ag0066000")
oxidation <- c("PAxG07Gg0074800", "PAxG05Ag0075200",
               "PAxG01Ag0193400", "PAxG01Gg0190900")
ABCG <- c("PAxG14Gg0063700", "PAxG18Ag0068300",
          "PAxG18Gg0071100", "PAxG09Ag0041200")
transfer <- c("PAxG02Gg0045400", "PAxG02Ag0046700",
              "PAxG01Ag0047500", "PAxG01Gg0046800")
Extracellular <- c("PAxG07Ag0035000", "PAxG07Gg0035100", "PAxG12Ag0061100",
                   "PAxG12Gg0067800", "PAxG01Ag0181400", "PAxG01Gg0172600",
                   "PAxG04Gg0042500", "PAxG13Gg0090700", "PAxG15Gg0135100",
                   "PAxG15Ag0134000")
all_features <- c(FAE, Glycerylation, oxidation, ABCG, transfer, Extracellular)

p_dot <- DotPlot(seurat_sub, features = all_features, group.by = "comparison")
df_dot <- p_dot$data
df_dot$Category <- rep(
  c("Fatty acyl elongase (FAE) complex", "Glycerylation",
    "Fatty acid omega-oxidation", "ABCG transporter",
    "Apoplastic lipid transfer proteins", "Extracellular polymerisation"),
  times = c(length(FAE), length(Glycerylation), length(oxidation),
            length(ABCG), length(transfer), length(Extracellular))
)

gene_labels <- c(
  "PAxG08Gg0138500 KCS11", "PAxG10Gg0113400 KCS6",
  "PAxG10Ag0109900 KCS6", "PAxG06Gg0202700 KCS10",
  "PAxG06Ag0207400 KCS10", "PAxG06Gg0092700 GPAT6",
  "PAxG14Ag0066000 GPAT8", "PAxG07Gg0074800 CYP86B1",
  "PAxG05Ag0075200 CYP86B1", "PAxG01Ag0193400 HTH",
  "PAxG01Gg0190900 HTH", "PAxG14Gg0063700 ABCG2",
  "PAxG18Ag0068300 ABCG32", "PAxG18Gg0071100 ABCG32",
  "PAxG09Ag0041200 ABCG15", "PAxG02Gg0045400 LTPG5",
  "PAxG02Ag0046700 LTPG5", "PAxG01Ag0047500 LTPG1",
  "PAxG01Gg0046800 LTPG1", "PAxG07Ag0035000 GELP51",
  "PAxG07Gg0035100 GELP51", "PAxG12Ag0061100 GELP38",
  "PAxG12Gg0067800 GELP38", "PAxG01Ag0181400 GELP67",
  "PAxG01Gg0172600 GELP67", "PAxG04Gg0042500 GELP102",
  "PAxG13Gg0090700 GELP36", "PAxG15Gg0135100 PEL3",
  "PAxG15Ag0134000 PEL3"
)
df_dot$features.plot <- factor(df_dot$features.plot,
                               levels = all_features, labels = gene_labels)

p_advanced_vertical <- ggplot(df_dot, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient(low = "#2166AC", high = "#B2182B") +
  scale_size_continuous(range = c(0, 6)) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.y = element_text(hjust = 1, face = "italic",
                               color = "black", size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1,
                               color = "black", size = 12, face = "bold"),
    axis.title  = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text.y = element_text(size = 10, face = "bold", angle = 0),
    panel.grid.major = element_line(color = "grey90"),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "right"
  )
ggsave(file.path(plot_dir, "DotPlot_epidermis_vs_phellem_vertical.pdf"),
       plot = p_advanced_vertical, width = 4, height = 10)

# ============================================================
# 4. Cortex subtypes: Subepidermis / Collenchyma (EPI) / (PER) / Parenchyma: Figure 6C
# ============================================================
seurat_sub <- subset(seurat_object, subcell_type %in% c(
  "Subepidermis", "Cortex collenchyma (EPI)",
  "Cortex collenchyma (PER)", "Cortex parenchyma"))
print(table(seurat_sub$subcell_type))

# Marker groups (one set per facet panel)
CollenchymaEPI  <- c("PAxG14Ag0153600", "PAxG01Ag0047500", "PAxG06Gg0202700")
Subepidermis    <- c("PAxG16Gg0105100", "PAxG01Gg0290500", "PAxG10Ag0194700")
CollenchymaPER <- c("PAxG01Ag0312400", "PAxG05Ag0075200", "PAxG18Ag0068300")
parenchyma      <- c("PAxG11Gg0128700", "PAxG03Gg0056800", "PAxG06Ag0162100")
all_features <- c(CollenchymaEPI, Subepidermis, CollenchymaPER, parenchyma)

# Extract DotPlot data using raw gene IDs
p_dot <- DotPlot(seurat_sub, features = all_features, group.by = "subcell_type")
df_dot <- p_dot$data

# Category column for faceting
df_dot$Category <- rep(
  c("Collenchyma(EPI)", "Subepidermis", "Collenchyma(PER)", "Parenchyma"),
  times = c(length(CollenchymaEPI), length(Subepidermis),
            length(CollenchymaPER), length(parenchyma))
)

# Replace raw IDs with gene symbols (order must match all_features)
gene_labels <- c(
  "CER1", "LTPG2", "KCS10",
  "GPAT6", "CAD9", "PIP2;5",
  "RWP1", "CYP86B1", "ABCG32",
  "GSTU19", "TT7", "UGT78D2"
)
df_dot$features.plot <- factor(df_dot$features.plot,
                              levels = all_features, labels = gene_labels)

# Faceted DotPlot
p_advanced <- ggplot(df_dot, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_distiller(palette = "RdBu") +
  scale_size_continuous(range = c(0, 6)) +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,
                               face = "italic", color = "black", size = 11),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title  = element_blank(),
    strip.background = element_rect(fill = "grey90", color = "black"),
    strip.text = element_text(size = 10, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.spacing = unit(0.1, "lines")
  )
ggsave(file.path(plot_dir, "DotPlot_cortex_subtypes.pdf"),
       plot = p_advanced, width = 9, height = 4)
