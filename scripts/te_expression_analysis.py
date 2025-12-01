#!/usr/bin/env python3
"""
TE Expression Analysis per Cluster

This script analyzes transposable element (TE) expression across microglial clusters:
1. Calculate TE expression statistics per cluster
2. Find differentially expressed TEs between clusters
3. Analyze TE families and classes per cluster
4. Visualize TE expression patterns
5. Compare with gene expression

Input: clustering/clustered_adata.h5ad
Output: clustering/te_analysis/

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
CLUSTERED_ADATA = 'clustering/clustered_adata.h5ad'
OUTPUT_DIR = 'clustering/te_analysis'

# Create output directory
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("TE EXPRESSION ANALYSIS PER CLUSTER")
print("=" * 80)

# Load clustered data
print("\n[1/10] Loading clustered data...")
adata = sc.read_h5ad(CLUSTERED_ADATA)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} features")

# Separate genes and TEs
print("\n[2/10] Separating genes and TEs...")
gene_mask = adata.var['feature_type'] == 'Gene'
te_mask = adata.var['feature_type'] == 'TE'

n_genes = gene_mask.sum()
n_tes = te_mask.sum()
print(f"  Genes: {n_genes}")
print(f"  TEs: {n_tes}")

# Create separate AnnData objects for genes and TEs
adata_genes = adata[:, gene_mask].copy()
adata_tes = adata[:, te_mask].copy()

# Calculate TE expression statistics per cluster
print("\n[3/10] Calculating TE expression statistics per cluster...")

te_stats_list = []
for cluster in sorted(adata.obs['leiden'].unique()):
    cluster_cells = adata.obs['leiden'] == cluster
    n_cells = cluster_cells.sum()
    
    # Overall TE metrics
    te_pct_mean = adata.obs.loc[cluster_cells, 'pct_te'].mean()
    te_pct_median = adata.obs.loc[cluster_cells, 'pct_te'].median()
    
    # TE counts
    te_counts = adata_tes[cluster_cells, :].X.sum(axis=1)
    te_counts_mean = np.mean(te_counts)
    te_counts_median = np.median(te_counts)
    
    # Number of TEs detected per cell
    n_tes_per_cell = (adata_tes[cluster_cells, :].X > 0).sum(axis=1)
    n_tes_mean = np.mean(n_tes_per_cell)
    n_tes_median = np.median(n_tes_per_cell)
    
    te_stats_list.append({
        'cluster': cluster,
        'n_cells': n_cells,
        'te_pct_mean': te_pct_mean,
        'te_pct_median': te_pct_median,
        'te_counts_mean': te_counts_mean,
        'te_counts_median': te_counts_median,
        'n_tes_mean': n_tes_mean,
        'n_tes_median': n_tes_median
    })

te_stats = pd.DataFrame(te_stats_list)
print("\nTE statistics per cluster:")
print(te_stats.to_string(index=False))

# Save statistics
te_stats.to_csv(f'{OUTPUT_DIR}/te_stats_per_cluster.csv', index=False)
print(f"\n  Saved: {OUTPUT_DIR}/te_stats_per_cluster.csv")

# Plot TE statistics per cluster
print("\n[4/10] Plotting TE statistics...")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# TE percentage by cluster
adata.obs.boxplot(column='pct_te', by='leiden', ax=axes[0, 0])
axes[0, 0].set_xlabel('Cluster')
axes[0, 0].set_ylabel('% TE expression')
axes[0, 0].set_title('TE Expression % by Cluster')
axes[0, 0].get_figure().suptitle('')

# Bar plot of mean TE %
axes[0, 1].bar(te_stats['cluster'].astype(str), te_stats['te_pct_mean'])
axes[0, 1].set_xlabel('Cluster')
axes[0, 1].set_ylabel('Mean % TE')
axes[0, 1].set_title('Mean TE Expression % by Cluster')
axes[0, 1].grid(True, alpha=0.3, axis='y')

# Number of TEs detected per cell
te_stats.plot(x='cluster', y='n_tes_mean', kind='bar', ax=axes[1, 0], legend=False)
axes[1, 0].set_xlabel('Cluster')
axes[1, 0].set_ylabel('Mean # TEs detected')
axes[1, 0].set_title('TEs Detected per Cell by Cluster')
axes[1, 0].grid(True, alpha=0.3, axis='y')

# TE counts by cluster
te_stats.plot(x='cluster', y='te_counts_mean', kind='bar', ax=axes[1, 1], legend=False, color='coral')
axes[1, 1].set_xlabel('Cluster')
axes[1, 1].set_ylabel('Mean TE counts')
axes[1, 1].set_title('TE Counts per Cell by Cluster')
axes[1, 1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/te_stats_by_cluster.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/te_stats_by_cluster.png")
plt.close()

# Analyze TE families and classes
print("\n[5/10] Analyzing TE families and classes...")

# Count TEs by class
te_class_counts = adata_tes.var['te_class'].value_counts()
print("\nTE classes in dataset:")
print(te_class_counts.to_string())

# Count TEs by family
te_family_counts = adata_tes.var['te_family'].value_counts().head(20)
print("\nTop 20 TE families:")
print(te_family_counts.to_string())

# Plot TE class distribution
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Pie chart of TE classes
axes[0].pie(te_class_counts.values, labels=te_class_counts.index, autopct='%1.1f%%')
axes[0].set_title('TE Class Distribution')

# Bar chart of top TE families
te_family_counts.plot(kind='barh', ax=axes[1])
axes[1].set_xlabel('Number of TEs')
axes[1].set_ylabel('TE Family')
axes[1].set_title('Top 20 TE Families')
axes[1].invert_yaxis()

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/te_class_family_distribution.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/te_class_family_distribution.png")
plt.close()

# Find differentially expressed TEs between clusters
print("\n[6/10] Finding differentially expressed TEs...")
print("  This may take several minutes...")

# Run differential expression on TEs only
sc.tl.rank_genes_groups(adata_tes, 'leiden', method='wilcoxon', use_raw=False)

# Save top DE TEs per cluster
print("\n[7/10] Saving differential TE expression results...")
for cluster in sorted(adata_tes.obs['leiden'].unique()):
    de_tes = sc.get.rank_genes_groups_df(adata_tes, group=cluster)
    de_tes_sig = de_tes[de_tes['pvals_adj'] < 0.05].copy()
    
    # Parse TE names to add class/family info
    te_info = []
    for te_name in de_tes_sig['names']:
        if te_name in adata_tes.var_names:
            te_class = adata_tes.var.loc[te_name, 'te_class']
            te_family = adata_tes.var.loc[te_name, 'te_family']
            te_info.append({'te_name': te_name, 'te_class': te_class, 'te_family': te_family})
        else:
            te_info.append({'te_name': te_name, 'te_class': 'Unknown', 'te_family': 'Unknown'})
    
    te_info_df = pd.DataFrame(te_info)
    de_tes_sig = pd.concat([de_tes_sig.reset_index(drop=True), te_info_df], axis=1)
    
    de_tes_sig.to_csv(f'{OUTPUT_DIR}/de_tes_cluster_{cluster}.csv', index=False)
    print(f"  Cluster {cluster}: {len(de_tes_sig)} significant TEs (padj < 0.05)")

# Plot top DE TEs
print("\n[8/10] Plotting differentially expressed TEs...")

# Heatmap of top DE TEs
sc.pl.rank_genes_groups_heatmap(
    adata_tes, n_genes=10, groupby='leiden',
    show_gene_labels=True, show=False, cmap='viridis',
    figsize=(12, 10)
)
plt.savefig(f'{OUTPUT_DIR}/de_tes_heatmap.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/de_tes_heatmap.png")
plt.close()

# Dotplot of top DE TEs
sc.pl.rank_genes_groups_dotplot(
    adata_tes, n_genes=5, groupby='leiden',
    show=False, figsize=(14, 8)
)
plt.savefig(f'{OUTPUT_DIR}/de_tes_dotplot.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/de_tes_dotplot.png")
plt.close()

# Analyze specific TE classes per cluster
print("\n[9/10] Analyzing TE class expression per cluster...")

# Calculate mean expression of each TE class per cluster
te_class_expr = []
for cluster in sorted(adata.obs['leiden'].unique()):
    cluster_cells = adata.obs['leiden'] == cluster
    
    for te_class in adata_tes.var['te_class'].unique():
        class_tes = adata_tes.var['te_class'] == te_class
        if class_tes.sum() > 0:
            class_expr = float(adata_tes[cluster_cells, :][:, class_tes].X.mean())
            te_class_expr.append({
                'cluster': str(cluster),  # Convert to string for pivot
                'te_class': te_class,
                'mean_expression': class_expr
            })

te_class_expr_df = pd.DataFrame(te_class_expr)

# Pivot for heatmap
te_class_pivot = te_class_expr_df.pivot(index='te_class', columns='cluster', values='mean_expression')

# Ensure all values are numeric
te_class_pivot = te_class_pivot.astype(float)

# Plot heatmap of TE class expression
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(te_class_pivot, cmap='YlOrRd', annot=True, fmt='.2f', ax=ax, cbar_kws={'label': 'Mean Expression'})
ax.set_xlabel('Cluster')
ax.set_ylabel('TE Class')
ax.set_title('Mean TE Class Expression by Cluster')
plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/te_class_expression_heatmap.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/te_class_expression_heatmap.png")
plt.close()

# Save TE class expression table
te_class_pivot.to_csv(f'{OUTPUT_DIR}/te_class_expression_by_cluster.csv')
print(f"  Saved: {OUTPUT_DIR}/te_class_expression_by_cluster.csv")

# Compare TE vs gene expression
print("\n[10/10] Comparing TE vs gene expression per cluster...")

comparison_stats = []
for cluster in sorted(adata.obs['leiden'].unique()):
    cluster_cells = adata.obs['leiden'] == cluster
    
    # Gene expression
    gene_expr = adata_genes[cluster_cells, :].X.sum(axis=1)
    gene_mean = np.mean(gene_expr)
    
    # TE expression
    te_expr = adata_tes[cluster_cells, :].X.sum(axis=1)
    te_mean = np.mean(te_expr)
    
    # Ratio
    te_gene_ratio = te_mean / gene_mean if gene_mean > 0 else 0
    
    comparison_stats.append({
        'cluster': cluster,
        'mean_gene_expr': gene_mean,
        'mean_te_expr': te_mean,
        'te_gene_ratio': te_gene_ratio
    })

comparison_df = pd.DataFrame(comparison_stats)
print("\nTE vs Gene expression comparison:")
print(comparison_df.to_string(index=False))

# Save comparison
comparison_df.to_csv(f'{OUTPUT_DIR}/te_vs_gene_expression.csv', index=False)
print(f"\n  Saved: {OUTPUT_DIR}/te_vs_gene_expression.csv")

# Plot TE vs gene expression
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Stacked bar chart
x = np.arange(len(comparison_df))
width = 0.6
axes[0].bar(x, comparison_df['mean_gene_expr'], width, label='Genes', color='steelblue')
axes[0].bar(x, comparison_df['mean_te_expr'], width, bottom=comparison_df['mean_gene_expr'], 
            label='TEs', color='coral')
axes[0].set_xlabel('Cluster')
axes[0].set_ylabel('Mean Expression')
axes[0].set_title('Gene vs TE Expression by Cluster')
axes[0].set_xticks(x)
axes[0].set_xticklabels(comparison_df['cluster'])
axes[0].legend()
axes[0].grid(True, alpha=0.3, axis='y')

# TE/Gene ratio
axes[1].bar(comparison_df['cluster'].astype(str), comparison_df['te_gene_ratio'], color='green')
axes[1].set_xlabel('Cluster')
axes[1].set_ylabel('TE / Gene Expression Ratio')
axes[1].set_title('TE to Gene Expression Ratio by Cluster')
axes[1].grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig(f'{OUTPUT_DIR}/te_vs_gene_comparison.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {OUTPUT_DIR}/te_vs_gene_comparison.png")
plt.close()

# Generate summary report
print("\n[11/11] Generating summary report...")
report_path = f'{OUTPUT_DIR}/te_analysis_report.txt'
with open(report_path, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("TE EXPRESSION ANALYSIS REPORT\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {pd.Timestamp.now()}\n\n")
    
    f.write("DATASET SUMMARY:\n")
    f.write(f"  Total cells: {adata.n_obs}\n")
    f.write(f"  Total features: {adata.n_vars}\n")
    f.write(f"  Genes: {n_genes}\n")
    f.write(f"  TEs: {n_tes}\n")
    f.write(f"  Clusters: {len(adata.obs['leiden'].unique())}\n\n")
    
    f.write("TE STATISTICS PER CLUSTER:\n")
    f.write(te_stats.to_string(index=False) + "\n\n")
    
    f.write("TE CLASS DISTRIBUTION:\n")
    f.write(te_class_counts.to_string() + "\n\n")
    
    f.write("TOP 10 TE FAMILIES:\n")
    f.write(te_family_counts.head(10).to_string() + "\n\n")
    
    f.write("TE VS GENE EXPRESSION:\n")
    f.write(comparison_df.to_string(index=False) + "\n\n")
    
    f.write("KEY FINDINGS:\n")
    f.write(f"  - Mean TE expression across all cells: {adata.obs['pct_te'].mean():.2f}%\n")
    f.write(f"  - Cluster with highest TE %: {te_stats.loc[te_stats['te_pct_mean'].idxmax(), 'cluster']}\n")
    f.write(f"  - Cluster with lowest TE %: {te_stats.loc[te_stats['te_pct_mean'].idxmin(), 'cluster']}\n")
    f.write(f"  - Most abundant TE class: {te_class_counts.index[0]}\n")
    f.write(f"  - Most abundant TE family: {te_family_counts.index[0]}\n")
    
    # Count significant TEs per cluster
    f.write("\nDIFFERENTIALLY EXPRESSED TEs (padj < 0.05):\n")
    for cluster in sorted(adata_tes.obs['leiden'].unique()):
        de_file = f'{OUTPUT_DIR}/de_tes_cluster_{cluster}.csv'
        if os.path.exists(de_file):
            n_sig = len(pd.read_csv(de_file))
            f.write(f"  Cluster {cluster}: {n_sig} significant TEs\n")

print(f"  Saved: {report_path}")

print("\n" + "=" * 80)
print("TE EXPRESSION ANALYSIS COMPLETE!")
print("=" * 80)
print(f"\nResults saved to: {OUTPUT_DIR}/")
print("\nKey outputs:")
print(f"  - TE statistics: {OUTPUT_DIR}/te_stats_per_cluster.csv")
print(f"  - DE TEs per cluster: {OUTPUT_DIR}/de_tes_cluster_*.csv")
print(f"  - TE class analysis: {OUTPUT_DIR}/te_class_expression_by_cluster.csv")
print(f"  - Summary report: {OUTPUT_DIR}/te_analysis_report.txt")
print("\nNext step: Create publication-quality visualizations (Step 6)")
