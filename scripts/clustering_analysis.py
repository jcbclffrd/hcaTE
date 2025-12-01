#!/usr/bin/env python3
"""
Clustering Analysis

This script performs Leiden clustering to identify microglial subpopulations:
1. Run Leiden clustering at multiple resolutions
2. Find optimal resolution (~7 clusters to match paper)
3. Annotate clusters with marker genes
4. Generate cluster visualization plots
5. Save clustered data

Input: clustering/normalized_adata.h5ad
Output: clustering/clustered_adata.h5ad

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
NORMALIZED_ADATA = 'clustering/normalized_adata.h5ad'
CLUSTERED_ADATA = 'clustering/clustered_adata.h5ad'
PLOT_DIR = 'clustering/cluster_plots'

# Create output directory
Path(PLOT_DIR).mkdir(parents=True, exist_ok=True)

# Microglia marker genes from paper (Alsema et al. 2020)
MICROGLIA_MARKERS = {
    'Core identity': ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'CSF1R'],
    'Homeostatic': ['P2RY12', 'TMEM119', 'GPR34', 'OLFML3', 'SALL1'],
    'Activated': ['CD68', 'CD86', 'ITGAX', 'MHC-II genes'],
    'Proliferation': ['MKI67', 'TOP2A', 'PCNA'],
    'Inflammation': ['IL1B', 'TNF', 'IL6', 'CXCL10'],
    'Phagocytosis': ['CD68', 'TREM2', 'APOE', 'LPL']
}

print("=" * 80)
print("CLUSTERING ANALYSIS")
print("=" * 80)

# Load normalized data
print("\n[1/7] Loading normalized data...")
adata = sc.read_h5ad(NORMALIZED_ADATA)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} features")

# Test multiple resolutions to find optimal
print("\n[2/7] Testing multiple resolutions...")
resolutions = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
n_clusters = []

for res in resolutions:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    n_clust = len(adata.obs[f'leiden_{res}'].unique())
    n_clusters.append(n_clust)
    print(f"  Resolution {res:.1f}: {n_clust} clusters")

# Find resolution closest to 7 clusters (matching paper)
target_clusters = 7
best_idx = np.argmin([abs(n - target_clusters) for n in n_clusters])
best_res = resolutions[best_idx]
best_n = n_clusters[best_idx]

print(f"\nOptimal resolution: {best_res} ({best_n} clusters, target={target_clusters})")
adata.obs['leiden'] = adata.obs[f'leiden_{best_res}'].copy()

# Plot resolution comparison
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(resolutions, n_clusters, 'o-', linewidth=2, markersize=8)
ax.axhline(y=target_clusters, color='red', linestyle='--', label=f'Target ({target_clusters} clusters)')
ax.axvline(x=best_res, color='green', linestyle='--', label=f'Selected (res={best_res})')
ax.set_xlabel('Resolution Parameter', fontsize=12)
ax.set_ylabel('Number of Clusters', fontsize=12)
ax.set_title('Leiden Clustering: Resolution vs Number of Clusters', fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()
plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/resolution_comparison.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/resolution_comparison.png")
plt.close()

# Visualize clusters on UMAP
print("\n[3/7] Visualizing clusters on UMAP...")

fig, axes = plt.subplots(1, 2, figsize=(16, 6))

# UMAP colored by cluster
sc.pl.umap(adata, color='leiden', legend_loc='on data', 
           title=f'Leiden Clustering (resolution={best_res})',
           ax=axes[0], show=False, frameon=False)

# UMAP colored by condition
sc.pl.umap(adata, color='condition', title='Condition',
           ax=axes[1], show=False, frameon=False)

plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/umap_clusters.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/umap_clusters.png")
plt.close()

# Cluster size and composition
print("\n[4/7] Analyzing cluster composition...")

cluster_stats = adata.obs.groupby('leiden').agg({
    'n_counts': ['mean', 'median'],
    'n_genes': ['mean', 'median'],
    'pct_te': ['mean', 'median'],
    'condition': lambda x: x.value_counts().to_dict()
}).round(1)

print("\nCluster statistics:")
print(f"  Total clusters: {best_n}")
for cluster in sorted(adata.obs['leiden'].unique()):
    n_cells = (adata.obs['leiden'] == cluster).sum()
    pct = 100 * n_cells / adata.n_obs
    print(f"  Cluster {cluster}: {n_cells} cells ({pct:.1f}%)")

# Plot cluster sizes
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Bar plot of cluster sizes
cluster_counts = adata.obs['leiden'].value_counts().sort_index()
axes[0].bar(range(len(cluster_counts)), cluster_counts.values)
axes[0].set_xlabel('Cluster')
axes[0].set_ylabel('Number of Cells')
axes[0].set_title('Cluster Sizes')
axes[0].set_xticks(range(len(cluster_counts)))
axes[0].set_xticklabels(cluster_counts.index)
axes[0].grid(True, alpha=0.3, axis='y')

# Stacked bar plot of conditions per cluster
condition_counts = pd.crosstab(adata.obs['leiden'], adata.obs['condition'])
condition_counts.plot(kind='bar', stacked=True, ax=axes[1])
axes[1].set_xlabel('Cluster')
axes[1].set_ylabel('Number of Cells')
axes[1].set_title('Cluster Composition by Condition')
axes[1].legend(title='Condition', bbox_to_anchor=(1.05, 1), loc='upper left')
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{PLOT_DIR}/cluster_composition.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/cluster_composition.png")
plt.close()

# Find marker genes for each cluster
print("\n[5/7] Finding cluster marker genes...")
print("  This may take a few minutes...")

# Use Wilcoxon rank-sum test (fast and robust)
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False, layer='log_normalized')

# Plot top marker genes
print("  Generating marker gene plots...")

# Heatmap of top markers
sc.pl.rank_genes_groups_heatmap(
    adata, n_genes=10, groupby='leiden',
    show_gene_labels=True, show=False, cmap='viridis',
    figsize=(12, 10)
)
plt.savefig(f'{PLOT_DIR}/marker_genes_heatmap.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/marker_genes_heatmap.png")
plt.close()

# Dotplot of top markers
sc.pl.rank_genes_groups_dotplot(
    adata, n_genes=5, groupby='leiden',
    show=False, figsize=(14, 8)
)
plt.savefig(f'{PLOT_DIR}/marker_genes_dotplot.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {PLOT_DIR}/marker_genes_dotplot.png")
plt.close()

# Check microglia marker expression
print("\n[6/7] Checking microglia marker expression...")

# Find which markers are present in the dataset
available_markers = {}
for category, markers in MICROGLIA_MARKERS.items():
    available = [m for m in markers if m in adata.var_names]
    if available:
        available_markers[category] = available

print(f"\nAvailable microglia markers:")
for category, markers in available_markers.items():
    print(f"  {category}: {', '.join(markers)}")

# Plot microglia markers on UMAP
if available_markers:
    all_markers = [m for markers in available_markers.values() for m in markers]
    n_markers = len(all_markers)
    n_cols = 4
    n_rows = (n_markers + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_markers > 1 else [axes]
    
    for idx, marker in enumerate(all_markers):
        sc.pl.umap(adata, color=marker, ax=axes[idx], show=False, 
                   title=marker, frameon=False, cmap='Reds')
    
    # Hide empty subplots
    for idx in range(n_markers, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    plt.savefig(f'{PLOT_DIR}/microglia_markers_umap.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {PLOT_DIR}/microglia_markers_umap.png")
    plt.close()
    
    # Dotplot of microglia markers by cluster
    sc.pl.dotplot(adata, var_names=all_markers, groupby='leiden',
                  show=False, figsize=(12, 6))
    plt.savefig(f'{PLOT_DIR}/microglia_markers_dotplot.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {PLOT_DIR}/microglia_markers_dotplot.png")
    plt.close()

# Save clustered data
print("\n[7/7] Saving clustered data...")
adata.write(CLUSTERED_ADATA)
print(f"  Saved: {CLUSTERED_ADATA}")

# Save clustering report
report_path = f'{PLOT_DIR}/clustering_report.txt'
with open(report_path, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("CLUSTERING ANALYSIS REPORT\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {pd.Timestamp.now()}\n\n")
    
    f.write("CLUSTERING PARAMETERS:\n")
    f.write(f"  Method: Leiden\n")
    f.write(f"  Resolution: {best_res}\n")
    f.write(f"  Number of clusters: {best_n}\n")
    f.write(f"  Target clusters: {target_clusters} (from paper)\n\n")
    
    f.write("CLUSTER SIZES:\n")
    for cluster in sorted(adata.obs['leiden'].unique()):
        n_cells = (adata.obs['leiden'] == cluster).sum()
        pct = 100 * n_cells / adata.n_obs
        f.write(f"  Cluster {cluster}: {n_cells} cells ({pct:.1f}%)\n")
    
    f.write("\nCLUSTER COMPOSITION BY CONDITION:\n")
    composition = pd.crosstab(adata.obs['leiden'], adata.obs['condition'])
    f.write(composition.to_string() + "\n\n")
    
    f.write("MARKER GENES (top 5 per cluster):\n")
    for cluster in sorted(adata.obs['leiden'].unique()):
        cluster_idx = int(cluster)
        markers = sc.get.rank_genes_groups_df(adata, group=cluster).head(5)
        f.write(f"\nCluster {cluster}:\n")
        for _, row in markers.iterrows():
            f.write(f"  {row['names']}: logFC={row['logfoldchanges']:.2f}, pval={row['pvals_adj']:.2e}\n")
    
    f.write("\nMICROGLIA MARKER EXPRESSION:\n")
    for category, markers in available_markers.items():
        f.write(f"\n{category}:\n")
        for marker in markers:
            mean_expr = adata[:, marker].X.mean()
            pct_expr = (adata[:, marker].X > 0).sum() / adata.n_obs * 100
            f.write(f"  {marker}: mean={mean_expr:.2f}, {pct_expr:.1f}% cells expressing\n")

print(f"  Saved: {report_path}")

# Save marker genes to CSV
print("\n[8/8] Saving marker gene tables...")
for cluster in sorted(adata.obs['leiden'].unique()):
    markers = sc.get.rank_genes_groups_df(adata, group=cluster)
    markers.to_csv(f'{PLOT_DIR}/markers_cluster_{cluster}.csv', index=False)
    print(f"  Saved: {PLOT_DIR}/markers_cluster_{cluster}.csv")

print("\n" + "=" * 80)
print("CLUSTERING ANALYSIS COMPLETE!")
print("=" * 80)
print(f"\nClustered data saved to: {CLUSTERED_ADATA}")
print(f"Plots and reports saved to: {PLOT_DIR}/")
print(f"\nIdentified {best_n} clusters at resolution {best_res}")
print("\nNext step: TE expression analysis per cluster (Step 5)")
