#!/usr/bin/env python3
"""
Compare our gene expression results with the paper's Table S10 and S11.
Performs left join to see which genes are in the paper's data and which are missing.
"""

import pandas as pd
import sys

def load_and_clean_paper_table(filepath, donor_name):
    """Load paper table and clean up column names."""
    # Skip the first row (table title) and use second row as header
    df = pd.read_csv(filepath, skiprows=[0])
    
    # Convert numeric columns
    numeric_cols = ['p-value', 'avg_logFC', 'pct.1', 'pct.2', 'adjusted p-value']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    
    # Add source column
    df['paper_source'] = f'Table S10/11 ({donor_name})'
    
    return df

def main():
    print("="*80)
    print("COMPARING OUR RESULTS WITH PAPER'S TABLE S10 & S11")
    print("="*80)
    
    # Load our results
    print("\n1. Loading our gene expression results...")
    our_results = pd.read_csv('10x_analysis/gene_expression_check.csv')
    print(f"   - Loaded {len(our_results)} genes from our analysis")
    
    # Load paper's tables
    print("\n2. Loading paper's supplemental tables...")
    paper_s10 = load_and_clean_paper_table(
        'paper_data/sheets/Table_S10.csv',
        'Donor 2018-135 AD'
    )
    print(f"   - Table S10 (AD donor): {len(paper_s10)} genes")
    
    paper_s11 = load_and_clean_paper_table(
        'paper_data/sheets/Table_S11.csv',
        'Donor 2019-010 MCI'
    )
    print(f"   - Table S11 (MCI donor): {len(paper_s11)} genes")
    
    # Combine paper tables
    paper_combined = pd.concat([paper_s10, paper_s11], ignore_index=True)
    print(f"   - Combined: {len(paper_combined)} total entries")
    
    # Perform LEFT JOIN (all our genes, matched with paper where available)
    print("\n3. Performing LEFT JOIN (our genes + paper data)...")
    merged = our_results.merge(
        paper_combined,
        left_on='gene',
        right_on='gene',
        how='left',
        suffixes=('_ours', '_paper')
    )
    
    # Identify matches and misses
    has_paper_data = merged['p-value'].notna()
    n_matched = has_paper_data.sum()
    n_missing = (~has_paper_data).sum()
    
    print(f"   - Matched: {n_matched}/{len(our_results)} genes ({100*n_matched/len(our_results):.1f}%)")
    print(f"   - Missing from paper: {n_missing}/{len(our_results)} genes")
    
    # Save full merged results
    output_file = '10x_analysis/comparison_with_paper.csv'
    merged.to_csv(output_file, index=False)
    print(f"\n4. Saved full comparison to: {output_file}")
    
    # Create summary table
    print("\n" + "="*80)
    print("SUMMARY TABLE")
    print("="*80)
    
    # Reorder columns for readability
    summary_cols = [
        'gene',
        'fold_change', 'pvalue', 'direction',
        'pct_AD', 'pct_MCI',
        'avg_logFC', 'p-value', 'cluster', 'paper_source'
    ]
    
    # Only include columns that exist
    summary_cols = [col for col in summary_cols if col in merged.columns]
    summary = merged[summary_cols].copy()
    
    # Rename for clarity
    summary.rename(columns={
        'fold_change': 'Our_FoldChange',
        'pvalue': 'Our_Pvalue',
        'direction': 'Our_Direction',
        'pct_AD': 'Our_PctAD',
        'pct_MCI': 'Our_PctMCI',
        'avg_logFC': 'Paper_LogFC',
        'p-value': 'Paper_Pvalue',
        'cluster': 'Paper_Cluster',
        'paper_source': 'Paper_Source'
    }, inplace=True)
    
    # Sort by our p-value
    summary = summary.sort_values('Our_Pvalue')
    
    print("\n" + summary.to_string(index=False))
    
    # Save summary
    summary_file = '10x_analysis/comparison_summary.csv'
    summary.to_csv(summary_file, index=False)
    print(f"\n5. Saved summary to: {summary_file}")
    
    # Print detailed analysis
    print("\n" + "="*80)
    print("DETAILED ANALYSIS")
    print("="*80)
    
    # Genes in our results that are also in paper
    matched_genes = merged[has_paper_data]
    if len(matched_genes) > 0:
        print(f"\n✓ GENES FOUND IN PAPER ({len(matched_genes)}):")
        for _, row in matched_genes.iterrows():
            print(f"\n  {row['gene']}:")
            print(f"    Our analysis:  {row['fold_change']:.2f}× FC, p={row['pvalue']:.2e} ({row['direction']})")
            print(f"    Paper:         logFC={row['avg_logFC']:.2f}, p={row['p-value']:.2e}, cluster={row['cluster']}")
            print(f"    Source:        {row['paper_source']}")
    
    # Genes NOT in paper's cluster-enriched lists
    missing_genes = merged[~has_paper_data]
    if len(missing_genes) > 0:
        print(f"\n✗ GENES NOT IN PAPER'S CLUSTER-ENRICHED LISTS ({len(missing_genes)}):")
        for _, row in missing_genes.iterrows():
            print(f"  - {row['gene']:10s}: {row['fold_change']:.2f}× FC, p={row['pvalue']:.2e} ({row['direction']})")
        
        print("\n  Note: These genes may still be in the full dataset but were not")
        print("  identified as cluster-enriched markers by the paper's analysis.")
    
    print("\n" + "="*80)
    print("INTERPRETATION")
    print("="*80)
    print("\nThe paper identified cluster-enriched genes (markers that define specific")
    print("cell clusters), while our analysis identified differentially expressed genes")
    print("between conditions (AD vs MCI).")
    print("\nGenes appearing in BOTH analyses are particularly robust findings:")
    print("- They define specific cell populations (paper's cluster analysis)")
    print("- They are differentially expressed between AD/MCI (our DE analysis)")
    print("\nGenes only in our analysis may be:")
    print("- Expressed broadly across clusters but differ between conditions")
    print("- Not meeting the paper's threshold for cluster enrichment")
    print("- Novel findings worth investigating further")
    print("="*80)

if __name__ == '__main__':
    main()
