# Single-Cell Clustering Analysis Plan

> **Feature Branch**: `feature/single-cell-clustering`  
> **Goal**: Reproduce paper's 7 microglial subpopulations and analyze TE expression at single-cell level

## Objective

Perform single-cell clustering analysis to:
1. Reproduce the paper's finding of **7 microglial subpopulations**
2. Analyze **TE expression patterns** at single-cell resolution
3. Identify which subtypes are **affected in Alzheimer's Disease**
4. Generate publication-quality results for novel TE findings

## Background

**Paper (Alsema et al. 2020)** identified 7 microglial subpopulations:
- Homeostatic microglia
- Inflammatory/activated states
- Disease-associated microglia (DAM)
- Regional-specific subtypes

**Our advantage**: We have TE quantification data that the paper didn't analyze!

## Data Available

- **31 samples** with scTE quantification (SRR11271993 - SRR11272023)
- **~12,942 total cells** (~80-84 cells per sample)
- **Per-cell expression matrices** in `scTE_output/SRR*/SRR*.csv`
- **Cell barcodes**: 84 valid 10bp barcodes
- **Features**: ~15,661 genes + TEs combined

## Analysis Pipeline

### Step 1: Load Single-Cell Data ✓
- [x] scTE output format: rows = features, columns = cells
- [x] Each sample has separate CSV file
- [ ] Combine all 31 samples into one matrix
- [ ] Add sample metadata (AD, CTR, CTR+, MCI)

### Step 2: Quality Control 📋
- [ ] Filter low-quality cells (min UMI count, min genes detected)
- [ ] Filter low-expression features (min cells expressing)
- [ ] Calculate QC metrics:
  - Total UMI counts per cell
  - Number of genes detected per cell
  - Percentage of TE expression per cell
- [ ] Generate QC plots

### Step 3: Normalization & Dimensionality Reduction 📋
- [ ] Normalize counts (CPM or log-normalization)
- [ ] Select highly variable genes
- [ ] PCA (Principal Component Analysis)
- [ ] UMAP or t-SNE for visualization
- [ ] Determine optimal number of PCs

### Step 4: Clustering 📋
- [ ] Leiden or Louvain clustering
- [ ] Optimize resolution parameter
- [ ] Aim for ~7 clusters (matching paper)
- [ ] Validate clusters with known markers:
  - CX3CR1, P2RY12, TMEM119 (homeostatic)
  - CD68, APOE (activated)
  - TREM2, CD33 (DAM markers)

### Step 5: Cluster Annotation 📋
- [ ] Differential expression per cluster
- [ ] Identify marker genes for each cluster
- [ ] Compare with paper's annotations
- [ ] Name clusters based on biology

### Step 6: TE Expression Analysis 📋
- [ ] Extract TE expression per cluster
- [ ] Identify cluster-specific TE patterns
- [ ] Compare TE expression: AD vs Control
- [ ] Identify which clusters show TE dysregulation

### Step 7: Visualization 📋
- [ ] UMAP colored by cluster
- [ ] UMAP colored by condition (AD vs CTR)
- [ ] Heatmap of marker genes per cluster
- [ ] Heatmap of TE expression per cluster
- [ ] Violin plots for key genes/TEs
- [ ] Proportion plots (cluster abundance by condition)

### Step 8: Statistical Analysis 📋
- [ ] Test for differential abundance (clusters in AD vs CTR)
- [ ] Cell-type-specific differential expression
- [ ] TE family enrichment analysis per cluster
- [ ] Integration with paper's findings

## Tools & Methods

### Primary Tool: Scanpy (Python)
```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Recommended workflow
adata = sc.read_csv("combined_matrix.csv")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
```

### Alternative: Seurat (R)
```r
library(Seurat)
library(dplyr)
library(ggplot2)

# Standard Seurat workflow
seurat_obj <- CreateSeuratObject(counts = count_matrix)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
```

## Expected Outputs

### 1. Clustering Results
- `clustering/cell_metadata.csv` - Cell assignments to clusters
- `clustering/cluster_markers.csv` - Marker genes per cluster
- `clustering/te_markers.csv` - Marker TEs per cluster

### 2. Visualizations
- `clustering/umap_clusters.png` - UMAP colored by cluster
- `clustering/umap_condition.png` - UMAP colored by disease state
- `clustering/marker_heatmap.png` - Heatmap of top markers
- `clustering/te_heatmap.png` - TE expression heatmap
- `clustering/cluster_proportions.png` - Abundance by condition

### 3. Statistical Results
- `clustering/differential_abundance.csv` - Cluster composition changes
- `clustering/cluster_specific_de.csv` - DE per cluster
- `clustering/te_enrichment.csv` - TE family enrichment

### 4. Documentation
- `clustering/RESULTS.md` - Summary of findings
- `clustering/METHODS.md` - Detailed methods

## Success Criteria

1. ✅ **Reproduce paper's clusters**: Identify 7 (±2) microglial subpopulations
2. ✅ **Validate with markers**: Clusters express expected marker genes
3. ✅ **Novel TE findings**: Identify cluster-specific TE expression patterns
4. ✅ **AD-specific changes**: Show which clusters/TEs change in disease
5. ✅ **Publication-ready**: High-quality plots and comprehensive analysis

## Timeline

- **Day 1**: Data loading, QC, normalization, dimensionality reduction
- **Day 2**: Clustering, annotation, TE analysis, visualization

## Getting Started

```bash
# Create analysis directory
mkdir -p clustering

# Create Python script for clustering
touch scripts/single_cell_clustering.py

# Install required packages (if needed)
pip install scanpy leidenalg python-igraph
```

## Notes

- **Remember**: scTE output has features as rows, cells as columns (opposite of Seurat default)
- **Cell barcodes**: Each cell is identified by its 10bp barcode + sample ID
- **Sample metadata**: Use `pseudobulk/pseudobulk_sample_info.csv` for condition labels
- **TE format**: TEs are named as `NAME#CLASS#FAMILY` (e.g., `MER61B#LTR#ERV1`)

## References

- Alsema et al. (2020) Front Mol Neurosci - Original paper
- Scanpy documentation: https://scanpy.readthedocs.io/
- Seurat vignettes: https://satijalab.org/seurat/
- scTE paper: Jin et al. (2015) Nat Commun

---

**Status**: 🚧 In Progress  
**Branch**: `feature/single-cell-clustering`  
**Started**: December 1, 2024
