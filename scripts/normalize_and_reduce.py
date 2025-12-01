#!/usr/bin/env python3
"""
Normalization and Dimensionality Reduction

This script performs:
1. Normalization (CPM + log transformation)
2. Highly variable gene selection
3. PCA dimensionality reduction
4. UMAP visualization
5. Saves processed data for clustering

Input: clustering/filtered_adata.h5ad (12,515 cells × 20,048 features)
Output: clustering/normalized_adata.h5ad (ready for clustering)

Author: GitHub Copilot
Date: December 1, 2024
"""

import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set plotting parameters
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', figsize=(8, 6))
sns.set_style('whitegrid')

# File paths
FILTERED_ADATA = 'clustering/filtered_adata.h5ad'
NORMALIZED_ADATA = 'clustering/normalized_adata.h5ad'
PLOT_DIR = 'clustering/normalization_plots'

# Create output directory
Path(PLOT_DIR).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("NORMALIZATION AND DIMENSIONALITY REDUCTION")
print("=" * 80)

# Load filtered data
print("\n[1/8] Loading filtered data...")
adata = sc.read_h5ad(FILTERED_ADATA)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} features")

# Save raw counts for later use (e.g., differential expression)
print("\n[2/8] Storing raw counts...")
adata.layers['counts'] = adata.X.copy()
print("  Saved raw counts in adata.layers['counts']")

# Normalize to counts per million (CPM)
print("\n[3/8] Normalizing to 10,000 counts per cell...")
sc.pp.normalize_total(adata, target_sum=1e4)
print("  Normalization complete (CPM)")

# Log-transform (log1p: log(x + 1))
print("\n[4/8] Log-transforming data...")
sc.pp.log1p(adata)
print("  Log transformation complete")

# Store normalized data
adata.layers['log_normalized'] = adata.X.copy()

# Select highly variable genes (for PCA)
print("\n[5/8] Selecting highly variable features...")

# We'll focus on genes for dimensionality reduction (TEs are often sparse)
gene_mask = adata.var['feature_type'] == 'Gene'
adata_genes = adata[:, gene_mask].copy()

print(f"  Working with {adata_genes.n_vars} genes")

# Find highly variable genes
# Using 'seurat' flavor (doesn't require scikit-misc dependency)
sc.pp.highly_variable_genes(
    adata_genes,
    n_top_genes=2000,
    flavor='seurat',
    subset=False
)

n_hvg = adata_genes.var['highly_variable'].sum()
print(f"  Selected {n_hvg} highly variable genes")

# Plot highly variable genes
sc.pl.highly_variable_genes(adata_genes, show=False)
plt.savefig(f'{PLOT_DIR}/highly_variable_genes.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/highly_variable_genes.png")
plt.close()

# Transfer HVG info back to main adata
adata.var['highly_variable'] = False
adata.var.loc[gene_mask, 'highly_variable'] = adata_genes.var['highly_variable'].values

# PCA on highly variable genes
print("\n[6/8] Running PCA...")
sc.tl.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
print("  PCA complete (50 components)")

# Plot variance ratio
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Elbow plot
axes[0].plot(range(1, 51), adata.uns['pca']['variance_ratio'], 'o-')
axes[0].set_xlabel('Principal Component')
axes[0].set_ylabel('Variance Ratio')
axes[0].set_title('PCA Variance Ratio (Elbow Plot)')
axes[0].axvline(x=30, color='red', linestyle='--', label='30 PCs (suggested)')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

# Cumulative variance
cumvar = np.cumsum(adata.uns['pca']['variance_ratio'])
axes[1].plot(range(1, 51), cumvar, 'o-')
axes[1].set_xlabel('Principal Component')
axes[1].set_ylabel('Cumulative Variance Ratio')
axes[1].set_title('Cumulative Variance Explained')
axes[1].axhline(y=0.8, color='red', linestyle='--', label='80% variance')
axes[1].axvline(x=30, color='red', linestyle='--', alpha=0.5)
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/pca_variance.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/pca_variance.png")
plt.close()

# Determine number of PCs to use (80% variance or 30 PCs, whichever is available)
cumvar = np.cumsum(adata.uns['pca']['variance_ratio'])
n_pcs_80_idx = np.where(cumvar >= 0.8)[0]
if len(n_pcs_80_idx) > 0:
    n_pcs_80 = n_pcs_80_idx[0] + 1
    n_pcs = min(n_pcs_80, 30)
    print(f"  PCs explaining 80% variance: {n_pcs_80}")
else:
    # If 50 PCs don't reach 80%, use all 30 PCs
    n_pcs = 30
    print(f"  50 PCs explain {cumvar[-1]*100:.1f}% variance (< 80%)")
print(f"  Using {n_pcs} PCs for downstream analysis")

# Compute neighborhood graph
print("\n[7/8] Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=n_pcs, method='umap', metric='euclidean')
print(f"  Neighborhood graph computed (k={15}, n_pcs={n_pcs})")

# UMAP
print("\n[8/8] Computing UMAP embedding...")
sc.tl.umap(adata, min_dist=0.3, spread=1.0)
print("  UMAP complete")

# Plot UMAP
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# UMAP colored by condition
conditions = adata.obs['condition'].unique()
colors = plt.cm.tab10(np.linspace(0, 1, len(conditions)))
color_map = dict(zip(conditions, colors))

for condition in conditions:
    mask = adata.obs['condition'] == condition
    axes[0].scatter(
        adata.obsm['X_umap'][mask, 0],
        adata.obsm['X_umap'][mask, 1],
        c=[color_map[condition]],
        label=condition,
        s=5,
        alpha=0.5
    )

axes[0].set_xlabel('UMAP 1')
axes[0].set_ylabel('UMAP 2')
axes[0].set_title('UMAP colored by Condition')
axes[0].legend(markerscale=3, loc='best')
axes[0].set_aspect('equal')

# UMAP colored by total counts
scatter = axes[1].scatter(
    adata.obsm['X_umap'][:, 0],
    adata.obsm['X_umap'][:, 1],
    c=adata.obs['n_counts'],
    s=5,
    cmap='viridis',
    alpha=0.5
)
axes[1].set_xlabel('UMAP 1')
axes[1].set_ylabel('UMAP 2')
axes[1].set_title('UMAP colored by Total Counts')
axes[1].set_aspect('equal')
plt.colorbar(scatter, ax=axes[1], label='Total counts')

plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/umap_overview.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/umap_overview.png")
plt.close()

# Additional UMAP plots
fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# UMAP colored by genes detected
scatter = axes[0].scatter(
    adata.obsm['X_umap'][:, 0],
    adata.obsm['X_umap'][:, 1],
    c=adata.obs['n_genes'],
    s=5,
    cmap='magma',
    alpha=0.5
)
axes[0].set_xlabel('UMAP 1')
axes[0].set_ylabel('UMAP 2')
axes[0].set_title('UMAP colored by Genes Detected')
axes[0].set_aspect('equal')
plt.colorbar(scatter, ax=axes[0], label='Genes detected')

# UMAP colored by TE %
scatter = axes[1].scatter(
    adata.obsm['X_umap'][:, 0],
    adata.obsm['X_umap'][:, 1],
    c=adata.obs['pct_te'],
    s=5,
    cmap='coolwarm',
    alpha=0.5
)
axes[1].set_xlabel('UMAP 1')
axes[1].set_ylabel('UMAP 2')
axes[1].set_title('UMAP colored by % TE Expression')
axes[1].set_aspect('equal')
plt.colorbar(scatter, ax=axes[1], label='% TE')

plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/umap_qc_metrics.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/umap_qc_metrics.png")
plt.close()

# Save processed data
print("\n[9/9] Saving normalized data...")
adata.write(NORMALIZED_ADATA)
print(f"  Saved: {NORMALIZED_ADATA}")

# Save processing report
report_path = f'{PLOT_DIR}/normalization_report.txt'
with open(report_path, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("NORMALIZATION AND DIMENSIONALITY REDUCTION REPORT\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {pd.Timestamp.now()}\n\n")
    
    f.write("DATA:\n")
    f.write(f"  Cells: {adata.n_obs}\n")
    f.write(f"  Features: {adata.n_vars}\n")
    f.write(f"    Genes: {(adata.var['feature_type'] == 'Gene').sum()}\n")
    f.write(f"    TEs: {(adata.var['feature_type'] == 'TE').sum()}\n\n")
    
    f.write("NORMALIZATION:\n")
    f.write("  Method: CPM normalization (10,000 counts per cell)\n")
    f.write("  Transformation: log1p (natural log)\n")
    f.write("  Layers stored: 'counts' (raw), 'log_normalized' (processed)\n\n")
    
    f.write("HIGHLY VARIABLE GENES:\n")
    f.write(f"  Method: Seurat v3\n")
    f.write(f"  Selected: {n_hvg} genes\n\n")
    
    f.write("PCA:\n")
    f.write(f"  Components computed: 50\n")
    f.write(f"  PCs for downstream: {n_pcs} (explaining {cumvar[n_pcs-1]*100:.1f}% variance)\n")
    f.write(f"  Variance explained by PC1: {adata.uns['pca']['variance_ratio'][0]*100:.1f}%\n")
    f.write(f"  Variance explained by PC2: {adata.uns['pca']['variance_ratio'][1]*100:.1f}%\n\n")
    
    f.write("NEIGHBORHOOD GRAPH:\n")
    f.write(f"  Neighbors: 15\n")
    f.write(f"  PCs used: {n_pcs}\n")
    f.write(f"  Method: UMAP\n")
    f.write(f"  Metric: Euclidean\n\n")
    
    f.write("UMAP:\n")
    f.write(f"  Min distance: 0.3\n")
    f.write(f"  Spread: 1.0\n")
    f.write(f"  Computed from {n_pcs} PCs\n")

print(f"  Saved: {report_path}")

print("\n" + "=" * 80)
print("NORMALIZATION AND DIMENSIONALITY REDUCTION COMPLETE!")
print("=" * 80)
print(f"\nProcessed data saved to: {NORMALIZED_ADATA}")
print(f"Plots saved to: {PLOT_DIR}/")
print(f"\nData ready for clustering with {n_pcs} PCs")
print("\nNext step: Clustering (Step 4)")
