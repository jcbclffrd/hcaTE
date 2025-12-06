#!/usr/bin/env python3
"""
Quick check: Is HSPA6 differentially expressed in 10x data?
Compare AD vs MCI for HSPA6 and other heat shock proteins
"""

import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================

BASE_DIR = "/home/jacobc/hcaTE"
SCTE_DIR = f"{BASE_DIR}/10x_scTE_output"

# Donors
DONOR_MCI = "Donor2019-010"  # MCI
DONOR_AD = "Donor2018-135"   # AD

# Genes to check
GENES_OF_INTEREST = [
    'HSPA6', 'HSPA1A', 'HSPA1B', 'HSPA8', 'HSP90AA1', 
    'DNAJB1', 'JUN', 'FOS', 'COA1'
]

# ============================================================================
# LOAD DATA
# ============================================================================

print("=" * 80)
print("HSPA6 Expression Check - 10x Genomics Data")
print("=" * 80)

print("\n1. Loading scTE output...")
df_mci = pd.read_csv(f"{SCTE_DIR}/{DONOR_MCI}/{DONOR_MCI}.csv", index_col=0)
df_ad = pd.read_csv(f"{SCTE_DIR}/{DONOR_AD}/{DONOR_AD}.csv", index_col=0)

# scTE format: rows=cells, columns=features. Need to transpose.
print(f"   Original shape - MCI: {df_mci.shape}, AD: {df_ad.shape}")
df_mci = df_mci.T
df_ad = df_ad.T

print(f"   MCI: {df_mci.shape[1]} cells, {df_mci.shape[0]} features")
print(f"   AD:  {df_ad.shape[1]} cells, {df_ad.shape[0]} features")

# ============================================================================
# SEPARATE GENES FROM TEs
# ============================================================================

print("\n2. Separating genes from TEs...")

# Genes don't have '#' in their names
genes_mci = df_mci[~df_mci.index.str.contains('#', regex=False)]
genes_ad = df_ad[~df_ad.index.str.contains('#', regex=False)]

print(f"   MCI genes: {genes_mci.shape[0]}")
print(f"   AD genes: {genes_ad.shape[0]}")

# ============================================================================
# CHECK FOR GENES OF INTEREST
# ============================================================================

print("\n3. Checking for genes of interest...")

available_genes = []
for gene in GENES_OF_INTEREST:
    in_mci = gene in genes_mci.index
    in_ad = gene in genes_ad.index
    both = in_mci and in_ad
    
    status = "✓ FOUND" if both else "✗ MISSING"
    print(f"   {gene:12s} - {status}")
    
    if both:
        available_genes.append(gene)

if not available_genes:
    print("\n⚠ WARNING: None of the genes of interest were found!")
    print("   Showing first 20 genes in dataset:")
    print(genes_mci.index[:20].tolist())
    exit(0)

# ============================================================================
# QUICK QC FILTER
# ============================================================================

print("\n4. Quick QC filtering...")

# Calculate QC metrics for MCI
mci_genes_per_cell = (genes_mci > 0).sum(axis=0)
mci_counts_per_cell = genes_mci.sum(axis=0)

# Calculate QC metrics for AD
ad_genes_per_cell = (genes_ad > 0).sum(axis=0)
ad_counts_per_cell = genes_ad.sum(axis=0)

# Filter cells: 200-5000 genes, >500 counts
mci_qc_pass = (mci_genes_per_cell >= 200) & (mci_genes_per_cell <= 5000) & (mci_counts_per_cell > 500)
ad_qc_pass = (ad_genes_per_cell >= 200) & (ad_genes_per_cell <= 5000) & (ad_counts_per_cell > 500)

genes_mci_filt = genes_mci.loc[:, mci_qc_pass]
genes_ad_filt = genes_ad.loc[:, ad_qc_pass]

print(f"   MCI: {genes_mci_filt.shape[1]} cells pass QC (from {genes_mci.shape[1]})")
print(f"   AD:  {genes_ad_filt.shape[1]} cells pass QC (from {genes_ad.shape[1]})")

# ============================================================================
# ANALYZE EXPRESSION FOR EACH GENE
# ============================================================================

print("\n5. Expression analysis for genes of interest:")
print("=" * 80)

results = []

for gene in available_genes:
    # Get expression values
    mci_expr = genes_mci_filt.loc[gene].values
    ad_expr = genes_ad_filt.loc[gene].values
    
    # Calculate statistics
    mci_mean = np.mean(mci_expr)
    ad_mean = np.mean(ad_expr)
    mci_pct = (mci_expr > 0).sum() / len(mci_expr) * 100
    ad_pct = (ad_expr > 0).sum() / len(ad_expr) * 100
    
    # Wilcoxon rank-sum test
    statistic, pval = stats.ranksums(ad_expr, mci_expr)
    
    # Log fold change (pseudo-count of 1)
    log2fc = np.log2((ad_mean + 1) / (mci_mean + 1))
    
    # Store results
    results.append({
        'gene': gene,
        'mean_MCI': mci_mean,
        'mean_AD': ad_mean,
        'pct_MCI': mci_pct,
        'pct_AD': ad_pct,
        'log2FC': log2fc,
        'fold_change': (ad_mean + 1) / (mci_mean + 1),
        'pvalue': pval,
        'direction': 'UP in AD' if log2fc > 0 else 'DOWN in AD'
    })
    
    # Print summary
    print(f"\n{gene}:")
    print(f"   Expression:")
    print(f"      MCI: mean={mci_mean:.2f}, detected in {mci_pct:.1f}% of cells")
    print(f"      AD:  mean={ad_mean:.2f}, detected in {ad_pct:.1f}% of cells")
    print(f"   Differential expression:")
    print(f"      Fold change: {(ad_mean + 1) / (mci_mean + 1):.2f}x")
    print(f"      Log2 FC: {log2fc:.3f}")
    print(f"      p-value: {pval:.2e}")
    print(f"      Direction: {results[-1]['direction']}")
    
    if pval < 0.05:
        sig = "*** SIGNIFICANT ***"
    else:
        sig = "(not significant)"
    print(f"      Status: {sig}")

# ============================================================================
# SUMMARY TABLE
# ============================================================================

print("\n" + "=" * 80)
print("SUMMARY TABLE")
print("=" * 80)

df_results = pd.DataFrame(results)
df_results = df_results.sort_values('pvalue')

print("\nGenes ranked by p-value:")
print(df_results.to_string(index=False))

# Check for significance
sig_genes = df_results[df_results['pvalue'] < 0.05]

print(f"\n{len(sig_genes)} out of {len(available_genes)} genes are significant (p < 0.05)")

if len(sig_genes) > 0:
    print("\nSignificant genes:")
    for _, row in sig_genes.iterrows():
        print(f"   {row['gene']:12s}: {row['fold_change']:.2f}x, p={row['pvalue']:.2e} ({row['direction']})")

# ============================================================================
# SAVE RESULTS
# ============================================================================

output_file = f"{BASE_DIR}/10x_analysis/gene_expression_check.csv"
df_results.to_csv(output_file, index=False)
print(f"\n✓ Results saved to: {output_file}")

print("\n" + "=" * 80)
print("Analysis complete!")
print("=" * 80)
