# Periderm single-cell and bulk RNA-seq analysis

Code for the Populus periderm development study. Analysis of Populus periderm development, including
snRNA-seq clustering, subpopulation characterization, trajectory inference,
scRNA-bulk projection (Scissor), and bulk RNA-seq time-course analysis.

## Data availability

Raw sequencing reads and large processed files are deposited at GEO
(GSE353574; SRA: PRJNA1551124) and are not stored in this repository.
The small tables and intermediate objects required to reproduce the
downstream analyses are included under `data/` and are fetched automatically
when you clone the repository; keep them in the `data/` directory.

| File | Description |
|---|---|
| `Periderm.RDS` | Processed Seurat object of the periderm snRNA-seq dataset (also available on GEO) |
| `EI_raw_counts.txt` | Bulk raw count matrix (EI samples) |
| `SI_raw_counts.txt` | Bulk raw count matrix (SI samples) |
| `EI_TPM_matrix.txt` | Bulk TPM matrix (EI samples) |
| `SI_TPM_matrix.txt` | Bulk TPM matrix (SI samples) |
| `hvg.csv` | Highly variable genes used for cortex scoring analysis (`15_snRNAseq_scoring_analysis.R`) |
| `hetero_marker_genes.csv` | Heterochromatin marker gene set used for module scoring (`15_snRNAseq_scoring_analysis.R`) |
| `go_data_anno.txt` | GO annotation table used for enrichment |
| `network_data_group.csv` | Gene-group assignments for the network visualization (`04_bulkRNAseq_network_stage_specific.R`) |

## Requirements

- R 4.5

### Command-line tools
- FastQC 0.11.9
- TrimGalore 0.6.10
- HISAT2 2.2.1
- featureCounts 2.0.3
- 10x Genomics Cell Ranger 7.1.0

### R packages
- Seurat 5.5.1
- DoubletFinder 2.0.6
- harmony 2.0.5
- clustree 0.5.1
- clusterProfiler 4.18.3
- ClusterGVis 0.0.9
- Scissor 2.0.0
- Monocle2 2.28.0
- TCseq 1.32.0
- ComplexHeatmap 2.26.0
- edgeR 4.6.3
- dplyr, ggplot2, patchwork, tidyr, forcats, grid, tidydr, vroom

## Repository structure
- `data/` — input tables and objects listed above
- `result/` — exported result tables produced by the scripts
- `plots/` — output figures
- Numbered analysis scripts are stored at the repository root; bulk and snRNA-seq
  scripts are distinguished by their filename prefix and are listed in execution order.

## How to run
Scripts are numbered in execution order. Run them from the repository root;
they read inputs from `data/`, write result tables to `result/`, and figures
to `plots/`.

## Notes
- `Periderm.RDS` contains the annotated Seurat object; raw fastq do not need
  to be reprocessed to reproduce the figures.
- Figure panels correspond to scripts as listed in the Methods section of
  the manuscript.
