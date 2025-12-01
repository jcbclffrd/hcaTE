#!/usr/bin/env python3
"""
Quality Control Filtering for Single-Cell Data

This script performs QC filtering on the loaded single-cell data:
1. Filters out low-quality cells (low counts, low genes)
2. Filters out low-expression features (expressed in too few cells)
3. Generates QC plots to document filtering decisions
4. Saves filtered AnnData object for downstream analysis

Input: clustering/raw_adata.h5ad (12,942 cells × 76,393 features)
Output: clustering/filtered_adata.h5ad (high-quality cells and features)

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
RAW_ADATA = 'clustering/raw_adata.h5ad'
FILTERED_ADATA = 'clustering/filtered_adata.h5ad'
QC_PLOT_DIR = 'clustering/qc_plots'

# Create output directory
Path(QC_PLOT_DIR).mkdir(parents=True, exist_ok=True)

print("=" * 80)
print("QUALITY CONTROL FILTERING")
print("=" * 80)

# Load raw data
print("\n[1/6] Loading raw data...")
adata = sc.read_h5ad(RAW_ADATA)
print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} features")
print(f"  Genes: {(adata.var['feature_type'] == 'Gene').sum()}")
print(f"  TEs: {(adata.var['feature_type'] == 'TE').sum()}")

# Store initial counts
n_cells_initial = adata.n_obs
n_features_initial = adata.n_vars

# Visualize QC metrics before filtering
print("\n[2/6] Generating QC plots (before filtering)...")

# Plot 1: Distribution of counts per cell
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Counts per cell
axes[0, 0].hist(adata.obs['n_counts'], bins=50, edgecolor='black')
axes[0, 0].set_xlabel('Total counts per cell')
axes[0, 0].set_ylabel('Number of cells')
axes[0, 0].set_title('Distribution of Total Counts')
axes[0, 0].axvline(adata.obs['n_counts'].median(), color='red', linestyle='--', label='Median')
axes[0, 0].legend()

# Genes per cell
axes[0, 1].hist(adata.obs['n_genes'], bins=50, edgecolor='black')
axes[0, 1].set_xlabel('Genes detected per cell')
axes[0, 1].set_ylabel('Number of cells')
axes[0, 1].set_title('Distribution of Genes Detected')
axes[0, 1].axvline(adata.obs['n_genes'].median(), color='red', linestyle='--', label='Median')
axes[0, 1].legend()

# TE percentage
axes[1, 0].hist(adata.obs['pct_te'], bins=50, edgecolor='black')
axes[1, 0].set_xlabel('% TE expression')
axes[1, 0].set_ylabel('Number of cells')
axes[1, 0].set_title('Distribution of TE Expression %')
axes[1, 0].axvline(adata.obs['pct_te'].median(), color='red', linestyle='--', label='Median')
axes[1, 0].legend()

# Scatter: counts vs genes
axes[1, 1].scatter(adata.obs['n_counts'], adata.obs['n_genes'], alpha=0.3, s=1)
axes[1, 1].set_xlabel('Total counts')
axes[1, 1].set_ylabel('Genes detected')
axes[1, 1].set_title('Counts vs Genes Detected')

plt.tight_layout()
plt.savefig(f'{QC_PLOT_DIR}/qc_metrics_before_filtering.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {QC_PLOT_DIR}/qc_metrics_before_filtering.png")
plt.close()

# Calculate filtering thresholds
print("\n[3/6] Calculating filtering thresholds...")

# For cells: use MAD (Median Absolute Deviation) approach
# Remove cells with counts or genes > 3 MADs below median
median_counts = adata.obs['n_counts'].median()
mad_counts = np.median(np.abs(adata.obs['n_counts'] - median_counts))
min_counts = max(500, median_counts - 3 * mad_counts)  # At least 500 counts

median_genes = adata.obs['n_genes'].median()
mad_genes = np.median(np.abs(adata.obs['n_genes'] - median_genes))
min_genes = max(200, median_genes - 3 * mad_genes)  # At least 200 genes

# For features: keep features expressed in at least 10 cells (0.1% of cells)
min_cells_expressing = 10

print(f"\nFiltering thresholds:")
print(f"  Min counts per cell: {min_counts:.0f}")
print(f"  Min genes per cell: {min_genes:.0f}")
print(f"  Min cells expressing feature: {min_cells_expressing}")

# Apply filters
print("\n[4/6] Applying filters...")

# Filter cells
print(f"\nBefore cell filtering: {adata.n_obs} cells")
sc.pp.filter_cells(adata, min_counts=min_counts)
print(f"After min_counts filter: {adata.n_obs} cells ({n_cells_initial - adata.n_obs} removed)")

n_after_counts = adata.n_obs
sc.pp.filter_cells(adata, min_genes=min_genes)
print(f"After min_genes filter: {adata.n_obs} cells ({n_after_counts - adata.n_obs} removed)")

# Filter features
print(f"\nBefore feature filtering: {adata.n_vars} features")
sc.pp.filter_genes(adata, min_cells=min_cells_expressing)
print(f"After min_cells filter: {adata.n_vars} features ({n_features_initial - adata.n_vars} removed)")

# Summary statistics
n_cells_final = adata.n_obs
n_features_final = adata.n_vars
n_genes_final = (adata.var['feature_type'] == 'Gene').sum()
n_tes_final = (adata.var['feature_type'] == 'TE').sum()

print("\n" + "=" * 80)
print("FILTERING SUMMARY")
print("=" * 80)
print(f"Cells: {n_cells_initial} → {n_cells_final} ({100 * n_cells_final / n_cells_initial:.1f}% retained)")
print(f"Features: {n_features_initial} → {n_features_final} ({100 * n_features_final / n_features_initial:.1f}% retained)")
print(f"  Genes: {n_genes_final} ({100 * n_genes_final / n_features_final:.1f}% of features)")
print(f"  TEs: {n_tes_final} ({100 * n_tes_final / n_features_final:.1f}% of features)")

# Cells by condition (after filtering)
print("\nCells by condition (after filtering):")
print(adata.obs['condition'].value_counts().to_string())

# QC metrics after filtering
print("\nQC metrics after filtering:")
print(f"  Mean counts per cell: {adata.obs['n_counts'].mean():.1f}")
print(f"  Median counts per cell: {adata.obs['n_counts'].median():.1f}")
print(f"  Mean genes per cell: {adata.obs['n_genes'].mean():.1f}")
print(f"  Median genes per cell: {adata.obs['n_genes'].median():.1f}")
print(f"  Mean % TE: {adata.obs['pct_te'].mean():.2f}%")
print(f"  Median % TE: {adata.obs['pct_te'].median():.2f}%")

# Visualize QC metrics after filtering
print("\n[5/6] Generating QC plots (after filtering)...")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Counts per cell
axes[0, 0].hist(adata.obs['n_counts'], bins=50, edgecolor='black', color='green', alpha=0.7)
axes[0, 0].set_xlabel('Total counts per cell')
axes[0, 0].set_ylabel('Number of cells')
axes[0, 0].set_title('Distribution of Total Counts (After Filtering)')
axes[0, 0].axvline(adata.obs['n_counts'].median(), color='red', linestyle='--', label='Median')
axes[0, 0].legend()

# Genes per cell
axes[0, 1].hist(adata.obs['n_genes'], bins=50, edgecolor='black', color='green', alpha=0.7)
axes[0, 1].set_xlabel('Genes detected per cell')
axes[0, 1].set_ylabel('Number of cells')
axes[0, 1].set_title('Distribution of Genes Detected (After Filtering)')
axes[0, 1].axvline(adata.obs['n_genes'].median(), color='red', linestyle='--', label='Median')
axes[0, 1].legend()

# TE percentage
axes[1, 0].hist(adata.obs['pct_te'], bins=50, edgecolor='black', color='green', alpha=0.7)
axes[1, 0].set_xlabel('% TE expression')
axes[1, 0].set_ylabel('Number of cells')
axes[1, 0].set_title('Distribution of TE Expression % (After Filtering)')
axes[1, 0].axvline(adata.obs['pct_te'].median(), color='red', linestyle='--', label='Median')
axes[1, 0].legend()

# Scatter: counts vs genes
axes[1, 1].scatter(adata.obs['n_counts'], adata.obs['n_genes'], alpha=0.3, s=1, color='green')
axes[1, 1].set_xlabel('Total counts')
axes[1, 1].set_ylabel('Genes detected')
axes[1, 1].set_title('Counts vs Genes Detected (After Filtering)')

plt.tight_layout()
plt.savefig(f'{QC_PLOT_DIR}/qc_metrics_after_filtering.png', dpi=150, bbox_inches='tight')
print(f"  Saved: {QC_PLOT_DIR}/qc_metrics_after_filtering.png")
plt.close()

# Create violin plots by condition
print("\n[6/6] Creating QC plots by condition...")

# Only plot conditions with enough cells
condition_counts = adata.obs['condition'].value_counts()
plot_conditions = condition_counts[condition_counts >= 50].index.tolist()

if len(plot_conditions) > 1:
    # Subset to plottable conditions
    adata_subset = adata[adata.obs['condition'].isin(plot_conditions)].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    
    # Counts by condition
    adata_subset.obs.boxplot(column='n_counts', by='condition', ax=axes[0])
    axes[0].set_xlabel('Condition')
    axes[0].set_ylabel('Total counts')
    axes[0].set_title('Counts per Cell by Condition')
    axes[0].get_figure().suptitle('')  # Remove default title
    
    # Genes by condition
    adata_subset.obs.boxplot(column='n_genes', by='condition', ax=axes[1])
    axes[1].set_xlabel('Condition')
    axes[1].set_ylabel('Genes detected')
    axes[1].set_title('Genes Detected by Condition')
    axes[1].get_figure().suptitle('')
    
    # TE % by condition
    adata_subset.obs.boxplot(column='pct_te', by='condition', ax=axes[2])
    axes[2].set_xlabel('Condition')
    axes[2].set_ylabel('% TE expression')
    axes[2].set_title('TE Expression % by Condition')
    axes[2].get_figure().suptitle('')
    
    plt.tight_layout()
    plt.savefig(f'{QC_PLOT_DIR}/qc_metrics_by_condition.png', dpi=150, bbox_inches='tight')
    print(f"  Saved: {QC_PLOT_DIR}/qc_metrics_by_condition.png")
    plt.close()
else:
    print("  Skipping condition plots (not enough labeled cells)")

# Save filtered data
print("\n[7/7] Saving filtered data...")
adata.write(FILTERED_ADATA)
print(f"  Saved: {FILTERED_ADATA}")

# Save filtering report
report_path = f'{QC_PLOT_DIR}/filtering_report.txt'
with open(report_path, 'w') as f:
    f.write("=" * 80 + "\n")
    f.write("QUALITY CONTROL FILTERING REPORT\n")
    f.write("=" * 80 + "\n\n")
    f.write(f"Date: {pd.Timestamp.now()}\n\n")
    
    f.write("FILTERING THRESHOLDS:\n")
    f.write(f"  Min counts per cell: {min_counts:.0f}\n")
    f.write(f"  Min genes per cell: {min_genes:.0f}\n")
    f.write(f"  Min cells expressing feature: {min_cells_expressing}\n\n")
    
    f.write("FILTERING RESULTS:\n")
    f.write(f"  Cells: {n_cells_initial} → {n_cells_final} ({100 * n_cells_final / n_cells_initial:.1f}% retained)\n")
    f.write(f"  Features: {n_features_initial} → {n_features_final} ({100 * n_features_final / n_features_initial:.1f}% retained)\n")
    f.write(f"    Genes: {n_genes_final} ({100 * n_genes_final / n_features_final:.1f}% of features)\n")
    f.write(f"    TEs: {n_tes_final} ({100 * n_tes_final / n_features_final:.1f}% of features)\n\n")
    
    f.write("CELLS BY CONDITION (after filtering):\n")
    f.write(adata.obs['condition'].value_counts().to_string() + "\n\n")
    
    f.write("QC METRICS (after filtering):\n")
    f.write(f"  Mean counts per cell: {adata.obs['n_counts'].mean():.1f}\n")
    f.write(f"  Median counts per cell: {adata.obs['n_counts'].median():.1f}\n")
    f.write(f"  Mean genes per cell: {adata.obs['n_genes'].mean():.1f}\n")
    f.write(f"  Median genes per cell: {adata.obs['n_genes'].median():.1f}\n")
    f.write(f"  Mean % TE: {adata.obs['pct_te'].mean():.2f}%\n")
    f.write(f"  Median % TE: {adata.obs['pct_te'].median():.2f}%\n")

print(f"  Saved: {report_path}")

print("\n" + "=" * 80)
print("QC FILTERING COMPLETE!")
print("=" * 80)
print(f"\nFiltered data saved to: {FILTERED_ADATA}")
print(f"QC plots saved to: {QC_PLOT_DIR}/")
print("\nNext step: Normalization and dimensionality reduction (Step 3)")
