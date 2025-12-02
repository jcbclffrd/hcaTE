#!/usr/bin/env python3
"""
Differential TE Expression Analysis - 10x Genomics Data
Compare AD vs MCI donors for transposable element expression

Analysis:
1. Load scTE count matrices
2. Filter and normalize
3. Differential expression: Donor2018-135 (AD) vs Donor2019-010 (MCI)
4. Identify TEs enriched in AD
5. Visualize results
"""

import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sc.set_figure_params(dpi=150, frameon=False, figsize=(6, 6))
sns.set_style("whitegrid")

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = "/home/jacobc/hcaTE"
SCTE_DIR = f"{BASE_DIR}/10x_scTE_output"
OUTPUT_DIR = f"{BASE_DIR}/10x_analysis"
FIGURE_DIR = f"{OUTPUT_DIR}/figures"

# Donors
DONOR_MCI = "Donor2019-010"  # MCI (control)
DONOR_AD = "Donor2018-135"   # AD (case)

# QC thresholds
MIN_GENES = 200
MAX_GENES = 5000
MAX_MT_PERCENT = 20

# TE filtering
MIN_CELLS_EXPRESSING = 10  # TE must be detected in at least 10 cells
MIN_COUNTS = 3  # Minimum total counts for a TE

# ============================================================================
# FUNCTIONS
# ============================================================================

def setup_directories():
    """Create output directories."""
    Path(OUTPUT_DIR).mkdir(exist_ok=True)
    Path(FIGURE_DIR).mkdir(exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Figure directory: {FIGURE_DIR}")


def load_scte_data(donor_name):
    """Load scTE output CSV for a donor."""
    csv_file = f"{SCTE_DIR}/{donor_name}/{donor_name}.csv"
    print(f"  Loading {donor_name}...")
    
    # Load CSV (rows=cells, cols=features)
    df = pd.read_csv(csv_file, index_col=0)
    
    print(f"    {df.shape[0]} cells × {df.shape[1]} features")
    
    return df


def separate_genes_and_tes(df):
    """
    Separate gene and TE counts from scTE output.
    TEs have '#' in their names, genes don't.
    """
    # Identify TEs (contain '#') and genes (don't contain '#')
    te_mask = df.columns.str.contains('#', regex=False)
    
    genes = df.loc[:, ~te_mask]
    tes = df.loc[:, te_mask]
    
    print(f"    Genes: {genes.shape[1]}")
    print(f"    TEs: {tes.shape[1]}")
    
    return genes, tes


def create_anndata(genes, tes, donor_name, condition):
    """
    Create AnnData object with gene counts as main matrix.
    Store TE data separately for later analysis.
    """
    # Create AnnData with gene counts
    adata = sc.AnnData(genes)
    
    # Store TE DataFrame separately (will be processed later)
    adata.uns['TE_df'] = tes
    
    # Add metadata
    adata.obs['donor'] = donor_name
    adata.obs['condition'] = condition
    adata.obs['n_counts'] = np.sum(genes.values, axis=1)
    adata.obs['n_genes'] = np.sum(genes.values > 0, axis=1)
    adata.obs['n_te_counts'] = np.sum(tes.values, axis=1)
    adata.obs['n_tes'] = np.sum(tes.values > 0, axis=1)
    
    # Calculate MT percentage
    mt_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mt'] = np.sum(
        adata[:, mt_genes].X, axis=1
    ) / np.sum(adata.X, axis=1) * 100
    
    print(f"    Created AnnData: {adata.shape}")
    print(f"    Stored {tes.shape[1]} TEs for later analysis")
    
    return adata


def qc_filter(adata):
    """Apply QC filtering."""
    print(f"  Before filtering: {adata.n_obs} cells")
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    adata = adata[adata.obs['n_genes'] < MAX_GENES, :]
    adata = adata[adata.obs['percent_mt'] < MAX_MT_PERCENT, :]
    
    # Filter genes (detected in at least 3 cells)
    sc.pp.filter_genes(adata, min_cells=3)
    
    print(f"  After filtering: {adata.n_obs} cells, {adata.n_vars} genes")
    
    return adata


def filter_tes(te_df):
    """Filter TEs based on detection criteria."""
    print(f"  Before filtering: {te_df.shape[1]} TEs")
    
    # Calculate TE statistics
    n_cells = (te_df > 0).sum(axis=0)
    total_counts = te_df.sum(axis=0)
    
    # Filter TEs
    keep = (n_cells >= MIN_CELLS_EXPRESSING) & (total_counts >= MIN_COUNTS)
    te_filtered = te_df.loc[:, keep]
    
    print(f"  After filtering: {te_filtered.shape[1]} TEs")
    print(f"    Removed {te_df.shape[1] - te_filtered.shape[1]} low-count TEs")
    
    return te_filtered


def normalize_and_log(adata):
    """Normalize and log-transform."""
    # Normalize to counts per 10k
    sc.pp.normalize_total(adata, target_sum=1e4)
    
    # Log transform
    sc.pp.log1p(adata)
    
    print(f"  Normalized and log-transformed")


def create_te_anndata(adata_gene, te_df_filtered):
    """
    Create a separate AnnData for TE analysis,
    keeping only cells that passed gene QC.
    """
    # Keep only cells that are in the filtered gene AnnData
    common_cells = adata_gene.obs_names.intersection(te_df_filtered.index)
    te_filtered_cells = te_df_filtered.loc[common_cells, :]
    
    # Create AnnData for TEs
    adata_te = sc.AnnData(te_filtered_cells)
    
    # Copy metadata from gene AnnData for matching cells
    adata_te.obs = adata_gene[common_cells, :].obs.copy()
    
    print(f"  Created TE AnnData: {adata_te.shape}")
    print(f"    Cells: {adata_te.n_obs} (matched from gene AnnData)")
    
    return adata_te


def differential_expression_tes(adata_te):
    """
    Perform differential expression analysis on TEs.
    Compare AD vs MCI.
    """
    print("\nDifferential TE Expression Analysis")
    print("="*70)
    
    # Normalize TE counts
    sc.pp.normalize_total(adata_te, target_sum=1e4)
    sc.pp.log1p(adata_te)
    
    # Run Wilcoxon rank-sum test
    print("  Running Wilcoxon rank-sum test...")
    sc.tl.rank_genes_groups(
        adata_te,
        groupby='condition',
        groups=['AD'],
        reference='MCI',
        method='wilcoxon',
        use_raw=False
    )
    
    # Extract results
    result = sc.get.rank_genes_groups_df(adata_te, group='AD')
    
    # Add statistics
    result['abs_logfc'] = result['logfoldchanges'].abs()
    result['-log10_pval'] = -np.log10(result['pvals_adj'] + 1e-300)
    
    # Parse TE names
    result['TE_name'] = result['names']
    result['TE_class'] = result['names'].str.split('#').str[1]
    result['TE_family'] = result['names'].str.split('#').str[2]
    
    # Calculate mean expression in each group
    ad_cells = adata_te.obs['condition'] == 'AD'
    mci_cells = adata_te.obs['condition'] == 'MCI'
    
    result['mean_AD'] = [adata_te[ad_cells, te].X.mean() for te in result['names']]
    result['mean_MCI'] = [adata_te[mci_cells, te].X.mean() for te in result['names']]
    result['pct_AD'] = [(adata_te[ad_cells, te].X > 0).mean() * 100 for te in result['names']]
    result['pct_MCI'] = [(adata_te[mci_cells, te].X > 0).mean() * 100 for te in result['names']]
    
    # Sort by significance
    result = result.sort_values('pvals_adj')
    
    return result, adata_te


def plot_te_volcano(result, output_file):
    """Volcano plot of TE differential expression."""
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Significance thresholds
    fc_threshold = 0.5  # log fold change
    pval_threshold = 0.05
    
    # Classify TEs
    sig_up = (result['logfoldchanges'] > fc_threshold) & (result['pvals_adj'] < pval_threshold)
    sig_down = (result['logfoldchanges'] < -fc_threshold) & (result['pvals_adj'] < pval_threshold)
    not_sig = ~(sig_up | sig_down)
    
    # Plot
    ax.scatter(result.loc[not_sig, 'logfoldchanges'], 
               result.loc[not_sig, '-log10_pval'],
               c='gray', alpha=0.3, s=10, label='Not Sig')
    
    ax.scatter(result.loc[sig_down, 'logfoldchanges'],
               result.loc[sig_down, '-log10_pval'],
               c='blue', alpha=0.6, s=20, label='Down in AD')
    
    ax.scatter(result.loc[sig_up, 'logfoldchanges'],
               result.loc[sig_up, '-log10_pval'],
               c='red', alpha=0.6, s=20, label='Up in AD')
    
    # Add threshold lines
    ax.axhline(-np.log10(pval_threshold), color='black', linestyle='--', linewidth=0.5)
    ax.axvline(fc_threshold, color='black', linestyle='--', linewidth=0.5)
    ax.axvline(-fc_threshold, color='black', linestyle='--', linewidth=0.5)
    
    # Labels
    ax.set_xlabel('Log2 Fold Change (AD vs MCI)', fontsize=12)
    ax.set_ylabel('-Log10 Adjusted P-value', fontsize=12)
    ax.set_title('Differential TE Expression: AD vs MCI\n10x Genomics Data', fontsize=14)
    ax.legend()
    
    # Annotate top TEs
    top_tes = result.nsmallest(10, 'pvals_adj')
    for _, row in top_tes.iterrows():
        if abs(row['logfoldchanges']) > fc_threshold and row['pvals_adj'] < pval_threshold:
            ax.annotate(row['TE_family'], 
                       xy=(row['logfoldchanges'], row['-log10_pval']),
                       xytext=(5, 5), textcoords='offset points',
                       fontsize=8, alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved volcano plot: {output_file}")
    plt.close()


def plot_top_tes(result, adata_te, output_file, n_top=20):
    """Heatmap of top differentially expressed TEs."""
    # Get top TEs (by adjusted p-value)
    top_tes = result.nsmallest(n_top, 'pvals_adj')['names'].tolist()
    
    # Create plot data
    adata_subset = adata_te[:, top_tes]
    
    # Sort cells by condition
    adata_subset = adata_subset[adata_subset.obs['condition'].argsort(), :]
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Get expression matrix
    expr = pd.DataFrame(
        adata_subset.X,
        index=adata_subset.obs_names,
        columns=[c.split('#')[2] for c in adata_subset.var_names]  # Use TE family name
    )
    
    # Plot
    sns.heatmap(expr.T, cmap='RdBu_r', center=0, 
                cbar_kws={'label': 'Log(CPM+1)'},
                xticklabels=False, yticklabels=True,
                ax=ax)
    
    # Add condition labels
    condition_colors = {'AD': 'red', 'MCI': 'blue'}
    colors = [condition_colors[c] for c in adata_subset.obs['condition']]
    
    ax.set_xlabel('Cells (sorted by condition)', fontsize=12)
    ax.set_ylabel('Transposable Elements', fontsize=12)
    ax.set_title(f'Top {n_top} Differentially Expressed TEs\nAD vs MCI', fontsize=14)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved heatmap: {output_file}")
    plt.close()


def plot_te_class_summary(result, output_file):
    """Bar plot of significant TEs by class."""
    # Filter significant TEs
    sig_tes = result[result['pvals_adj'] < 0.05].copy()
    
    if len(sig_tes) == 0:
        print("  No significant TEs found!")
        return
    
    # Count by class and direction
    sig_tes['direction'] = sig_tes['logfoldchanges'].apply(
        lambda x: 'Up in AD' if x > 0 else 'Down in AD'
    )
    
    class_counts = sig_tes.groupby(['TE_class', 'direction']).size().unstack(fill_value=0)
    
    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    class_counts.plot(kind='barh', stacked=False, ax=ax, 
                      color=['blue', 'red'], alpha=0.7)
    
    ax.set_xlabel('Number of Significant TEs', fontsize=12)
    ax.set_ylabel('TE Class', fontsize=12)
    ax.set_title('Significant TEs by Class (FDR < 0.05)', fontsize=14)
    ax.legend(title='Direction')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved TE class summary: {output_file}")
    plt.close()


def main():
    """Main analysis pipeline."""
    print("\n" + "="*70)
    print("Differential TE Expression Analysis - 10x Genomics")
    print("AD vs MCI Comparison")
    print("="*70)
    
    # Setup
    setup_directories()
    
    # Load data
    print("\n1. Loading scTE data...")
    df_mci = load_scte_data(DONOR_MCI)
    df_ad = load_scte_data(DONOR_AD)
    
    # Separate genes and TEs
    print("\n2. Separating genes and TEs...")
    genes_mci, tes_mci = separate_genes_and_tes(df_mci)
    genes_ad, tes_ad = separate_genes_and_tes(df_ad)
    
    # Create AnnData objects
    print("\n3. Creating AnnData objects...")
    adata_mci = create_anndata(genes_mci, tes_mci, DONOR_MCI, 'MCI')
    adata_ad = create_anndata(genes_ad, tes_ad, DONOR_AD, 'AD')
    
    # Combine donors
    print("\n4. Combining donors...")
    # Make cell barcodes unique by adding donor prefix
    adata_mci.obs_names = [f"MCI_{bc}" for bc in adata_mci.obs_names]
    adata_ad.obs_names = [f"AD_{bc}" for bc in adata_ad.obs_names]
    
    # Also update TE DataFrames
    tes_mci.index = adata_mci.obs_names
    tes_ad.index = adata_ad.obs_names
    
    adata = sc.concat([adata_mci, adata_ad])
    print(f"  Combined: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # QC filtering
    print("\n5. QC filtering...")
    adata = qc_filter(adata)
    
    # Normalize genes
    print("\n6. Normalizing gene counts...")
    normalize_and_log(adata)
    
    # Filter TEs
    print("\n7. Filtering TEs...")
    tes_combined = pd.concat([tes_mci, tes_ad])
    tes_filtered = filter_tes(tes_combined)
    
    # Create TE AnnData
    print("\n8. Creating TE AnnData...")
    adata_te = create_te_anndata(adata, tes_filtered)
    
    # Differential expression
    print("\n9. Running differential expression...")
    result, adata_te_norm = differential_expression_tes(adata_te)
    
    # Save results
    output_csv = f"{OUTPUT_DIR}/TE_differential_AD_vs_MCI.csv"
    result.to_csv(output_csv, index=False)
    print(f"\n  Saved results: {output_csv}")
    
    # Summary statistics
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total TEs tested: {len(result)}")
    print(f"Significant TEs (FDR < 0.05): {(result['pvals_adj'] < 0.05).sum()}")
    print(f"  Up in AD: {((result['pvals_adj'] < 0.05) & (result['logfoldchanges'] > 0)).sum()}")
    print(f"  Down in AD: {((result['pvals_adj'] < 0.05) & (result['logfoldchanges'] < 0)).sum()}")
    
    # Top 10 most significant
    print("\nTop 10 Most Significant TEs:")
    print("-"*70)
    top10 = result.nsmallest(10, 'pvals_adj')
    for _, row in top10.iterrows():
        direction = "↑" if row['logfoldchanges'] > 0 else "↓"
        print(f"  {direction} {row['TE_family']:20s} (Class: {row['TE_class']:15s}) "
              f"FC={row['logfoldchanges']:6.2f}, FDR={row['pvals_adj']:.2e}")
    
    # TE class breakdown
    print("\nSignificant TEs by Class:")
    print("-"*70)
    sig_tes = result[result['pvals_adj'] < 0.05]
    class_counts = sig_tes['TE_class'].value_counts()
    for te_class, count in class_counts.items():
        print(f"  {te_class:20s}: {count} TEs")
    
    # Create plots
    print("\n10. Creating visualizations...")
    plot_te_volcano(result, f"{FIGURE_DIR}/TE_volcano_AD_vs_MCI.png")
    plot_top_tes(result, adata_te_norm, f"{FIGURE_DIR}/TE_heatmap_top20.png", n_top=20)
    plot_te_class_summary(result, f"{FIGURE_DIR}/TE_class_summary.png")
    
    print("\n" + "="*70)
    print("Analysis complete!")
    print(f"Results: {output_csv}")
    print(f"Figures: {FIGURE_DIR}/")
    print("="*70)


if __name__ == "__main__":
    main()
