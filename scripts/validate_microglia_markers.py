#!/usr/bin/env python3
"""
Validate pipeline by checking expression of key microglia marker genes.

These genes should be:
1. Highly expressed (microglia identity markers)
2. Not differentially expressed between AD and Control (per paper's main finding)

Key markers:
- CX3CR1: Microglia-specific chemokine receptor
- P2RY12: Microglia purinergic receptor
- TMEM119: Microglia marker
- ITGAM (CD11B): Myeloid marker
- PTPRC (CD45): Pan-leukocyte marker
- TREM2: AD risk gene, microglial receptor
- APOE: AD risk gene
- AIF1 (IBA1): Pan-macrophage/microglia marker
- CSF1R: Microglia survival receptor
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Set style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['figure.figsize'] = (14, 10)
plt.rcParams['font.size'] = 10

# Key microglia markers to validate
MICROGLIA_MARKERS = {
    'Identity': ['CX3CR1', 'P2RY12', 'TMEM119', 'AIF1', 'CSF1R'],
    'Myeloid': ['ITGAM', 'PTPRC', 'CD68', 'FCGR3A'],
    'AD Risk': ['TREM2', 'APOE', 'CD33', 'MS4A6A'],
    'Activation': ['IL1B', 'TNF', 'IL6', 'CD86', 'CD163'],
    'Homeostatic': ['P2RY13', 'SELPLG', 'GPR34', 'OLFML3']
}

ALL_MARKERS = [gene for genes in MICROGLIA_MARKERS.values() for gene in genes]

def load_pseudobulk_counts():
    """Load the pseudobulk count matrix"""
    # Use normalized counts which has features as rows (correct orientation)
    count_file = Path('/home/jacobc/hcaTE/pseudobulk/pseudobulk_normalized_counts.csv')
    if not count_file.exists():
        print(f"ERROR: Count matrix not found at {count_file}")
        sys.exit(1)
    
    counts = pd.read_csv(count_file, index_col=0)
    print(f"Loaded normalized count matrix: {counts.shape[0]} features × {counts.shape[1]} samples")
    print(f"Note: Using DESeq2 normalized counts for visualization")
    return counts

def load_sample_info(counts):
    """Load sample condition information, filtered to samples in count matrix"""
    info_file = Path('/home/jacobc/hcaTE/pseudobulk/pseudobulk_sample_info.csv')
    info = pd.read_csv(info_file)
    
    # Filter to only samples present in counts
    info = info[info['Sample'].isin(counts.columns)]
    
    print(f"Loaded sample info: {len(info)} samples (filtered to match count matrix)")
    print(f"Conditions: {info['Condition'].value_counts().to_dict()}")
    return info

def load_de_results():
    """Load differential expression results"""
    de_file = Path('/home/jacobc/hcaTE/pseudobulk/pseudobulk_diffexp_results_AD_vs_Control.csv')
    if not de_file.exists():
        print("WARNING: DE results not found")
        return None
    
    de = pd.read_csv(de_file)
    print(f"Loaded DE results: {len(de)} features tested")
    return de

def check_marker_presence(counts, markers):
    """Check which markers are present in the count matrix"""
    present = []
    missing = []
    
    for marker in markers:
        if marker in counts.index:
            present.append(marker)
        else:
            missing.append(marker)
    
    print(f"\n{'='*80}")
    print(f"MARKER GENE PRESENCE CHECK")
    print(f"{'='*80}")
    print(f"Present: {len(present)}/{len(markers)} ({len(present)/len(markers)*100:.1f}%)")
    print(f"Missing: {len(missing)} genes")
    
    if missing:
        print(f"\nMissing genes: {', '.join(missing)}")
    
    return present, missing

def get_marker_expression(counts, sample_info, present_markers):
    """Extract expression data for present markers"""
    marker_counts = counts.loc[present_markers].copy()
    
    # Counts are already normalized by DESeq2, just log2 transform
    # DESeq2 normalization accounts for library size differences
    marker_log2 = np.log2(marker_counts + 1)
    
    # Add condition information
    marker_data = marker_log2.T
    marker_data['Condition'] = marker_data.index.map(
        sample_info.set_index('Sample')['Condition']
    )
    
    return marker_counts, marker_log2, marker_data

def calculate_marker_stats(marker_counts, sample_info):
    """Calculate summary statistics for each marker"""
    stats = []
    
    for gene in marker_counts.index:
        gene_counts = marker_counts.loc[gene]
        
        # Overall stats
        mean_count = gene_counts.mean()
        detection_rate = (gene_counts > 0).sum() / len(gene_counts) * 100
        
        # By condition
        for condition in sample_info['Condition'].unique():
            samples = sample_info[sample_info['Condition'] == condition]['Sample']
            cond_counts = gene_counts[samples]
            
            stats.append({
                'Gene': gene,
                'Condition': condition,
                'Mean_Count': cond_counts.mean(),
                'Median_Count': cond_counts.median(),
                'Detection_Rate': (cond_counts > 0).sum() / len(cond_counts) * 100,
                'N_Samples': len(cond_counts)
            })
    
    stats_df = pd.DataFrame(stats)
    return stats_df

def check_de_status(de_results, present_markers):
    """Check differential expression status of markers"""
    if de_results is None:
        return None
    
    marker_de = de_results[de_results['feature'].isin(present_markers)].copy()
    
    print(f"\n{'='*80}")
    print(f"DIFFERENTIAL EXPRESSION STATUS")
    print(f"{'='*80}")
    
    if len(marker_de) == 0:
        print("No markers found in DE results")
        return None
    
    # Categorize by significance
    marker_de['Significant'] = marker_de['padj'] < 0.05
    
    sig_up = marker_de[(marker_de['Significant']) & (marker_de['log2FoldChange'] > 0)]
    sig_down = marker_de[(marker_de['Significant']) & (marker_de['log2FoldChange'] < 0)]
    not_sig = marker_de[~marker_de['Significant']]
    
    print(f"Significantly UP in AD: {len(sig_up)} genes")
    if len(sig_up) > 0:
        for _, row in sig_up.iterrows():
            print(f"  {row['feature']}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")
    
    print(f"\nSignificantly DOWN in AD: {len(sig_down)} genes")
    if len(sig_down) > 0:
        for _, row in sig_down.iterrows():
            print(f"  {row['feature']}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.2e}")
    
    print(f"\nNOT significantly different: {len(not_sig)} genes")
    if len(not_sig) > 0:
        print("Genes:")
        for _, row in not_sig.iterrows():
            print(f"  {row['feature']}: log2FC={row['log2FoldChange']:.2f}, padj={row['padj']:.3f}")
    
    return marker_de

def plot_marker_expression(marker_data, present_markers, output_dir):
    """Create comprehensive visualization of marker expression"""
    
    # Prepare data for plotting
    plot_data = marker_data.melt(
        id_vars=['Condition'],
        value_vars=present_markers,
        var_name='Gene',
        value_name='log2NormCount'
    )
    
    # Create figure with subplots
    n_genes = len(present_markers)
    ncols = 5
    nrows = int(np.ceil(n_genes / ncols))
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(20, 4*nrows))
    axes = axes.flatten() if n_genes > 1 else [axes]
    
    # Define condition order (only include conditions present in data)
    all_conditions = ['CTR', 'CTR+', 'MCI', 'AD']
    present_conditions = [c for c in all_conditions if c in marker_data['Condition'].unique()]
    condition_order = present_conditions
    condition_colors = {'CTR': '#2ecc71', 'CTR+': '#f39c12', 'MCI': '#e74c3c', 'AD': '#c0392b'}
    
    for idx, gene in enumerate(present_markers):
        ax = axes[idx]
        gene_data = plot_data[plot_data['Gene'] == gene]
        
        # Prepare data for boxplot
        box_data = [gene_data[gene_data['Condition'] == cond]['log2NormCount'].values 
                    for cond in condition_order]
        
        # Create boxplot
        bp = ax.boxplot(box_data, labels=condition_order, patch_artist=True)
        
        # Color boxes
        for patch, cond in zip(bp['boxes'], condition_order):
            patch.set_facecolor(condition_colors.get(cond, '#95a5a6'))
            patch.set_alpha(0.7)
        
        # Add individual points
        for i, cond in enumerate(condition_order):
            cond_data = gene_data[gene_data['Condition'] == cond]['log2NormCount'].values
            x = np.random.normal(i+1, 0.04, size=len(cond_data))
            ax.plot(x, cond_data, 'ko', alpha=0.5, markersize=4)
        
        ax.set_title(gene, fontsize=12, fontweight='bold')
        ax.set_xlabel('')
        ax.set_ylabel('log2(Normalized Count + 1)' if idx % ncols == 0 else '')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3)
    
    # Hide extra subplots
    for idx in range(n_genes, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    output_file = output_dir / 'marker_expression_boxplots.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved boxplot: {output_file}")
    plt.close()
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(12, max(8, n_genes * 0.4)))
    
    # Prepare matrix for heatmap
    heatmap_data = marker_data.pivot_table(
        index='Condition',
        values=present_markers,
        aggfunc='mean'
    )
    # Reorder to match condition_order (only present conditions)
    heatmap_data = heatmap_data.reindex(condition_order)
    
    # Create heatmap using imshow
    im = ax.imshow(heatmap_data.T, aspect='auto', cmap='RdYlBu_r')
    
    # Set ticks and labels
    ax.set_xticks(np.arange(len(condition_order)))
    ax.set_yticks(np.arange(len(present_markers)))
    ax.set_xticklabels(condition_order)
    ax.set_yticklabels(present_markers)
    
    # Rotate the tick labels and set their alignment
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center")
    
    # Add text annotations
    for i in range(len(present_markers)):
        for j in range(len(condition_order)):
            text = ax.text(j, i, f'{heatmap_data.T.iloc[i, j]:.1f}',
                          ha="center", va="center", color="black", fontsize=8)
    
    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Mean log2(Normalized Count + 1)', rotation=270, labelpad=20)
    
    ax.set_title('Mean Marker Expression by Condition', fontsize=14, fontweight='bold')
    ax.set_xlabel('Condition', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)
    
    plt.tight_layout()
    output_file = output_dir / 'marker_expression_heatmap.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved heatmap: {output_file}")
    plt.close()

def create_summary_report(stats_df, marker_de, present_markers, missing_markers, output_dir):
    """Create a text summary report"""
    report_file = output_dir / 'VALIDATION_REPORT.txt'
    
    with open(report_file, 'w') as f:
        f.write('='*80 + '\n')
        f.write('MICROGLIA MARKER GENE VALIDATION REPORT\n')
        f.write('='*80 + '\n\n')
        
        f.write('OBJECTIVE:\n')
        f.write('Validate that our pipeline correctly identifies microglia and that gene\n')
        f.write('expression matches the paper\'s main finding: "Gene expression profiles did\n')
        f.write('not differ between AD donors and non-demented elderly."\n\n')
        
        f.write('='*80 + '\n')
        f.write('1. MARKER PRESENCE\n')
        f.write('='*80 + '\n')
        f.write(f'Present: {len(present_markers)} genes\n')
        f.write(f'Missing: {len(missing_markers)} genes\n\n')
        
        if missing_markers:
            f.write('Missing genes:\n')
            for gene in missing_markers:
                f.write(f'  - {gene}\n')
        
        f.write('\n' + '='*80 + '\n')
        f.write('2. EXPRESSION LEVELS (Mean counts by condition)\n')
        f.write('='*80 + '\n\n')
        
        for category, genes in MICROGLIA_MARKERS.items():
            present_in_cat = [g for g in genes if g in present_markers]
            if not present_in_cat:
                continue
                
            f.write(f'\n{category} Markers:\n')
            f.write('-'*80 + '\n')
            
            for gene in present_in_cat:
                gene_stats = stats_df[stats_df['Gene'] == gene]
                f.write(f'\n{gene}:\n')
                for _, row in gene_stats.iterrows():
                    f.write(f"  {row['Condition']}: mean={row['Mean_Count']:.1f}, "
                           f"detection={row['Detection_Rate']:.1f}%\n")
        
        f.write('\n' + '='*80 + '\n')
        f.write('3. DIFFERENTIAL EXPRESSION STATUS\n')
        f.write('='*80 + '\n\n')
        
        if marker_de is not None and len(marker_de) > 0:
            f.write('According to the paper, microglia markers should NOT be significantly\n')
            f.write('different between AD and Control.\n\n')
            
            sig_markers = marker_de[marker_de['padj'] < 0.05]
            if len(sig_markers) > 0:
                f.write(f'⚠️  WARNING: {len(sig_markers)} markers are significantly DE:\n')
                for _, row in sig_markers.iterrows():
                    direction = 'UP' if row['log2FoldChange'] > 0 else 'DOWN'
                    f.write(f"  {row['feature']}: {direction} (log2FC={row['log2FoldChange']:.2f}, "
                           f"padj={row['padj']:.2e})\n")
            else:
                f.write('✓ Good: No markers are significantly DE (matches paper)\n')
            
            not_sig = marker_de[marker_de['padj'] >= 0.05]
            if len(not_sig) > 0:
                f.write(f'\n✓ {len(not_sig)} markers are NOT significantly DE (expected):\n')
                for _, row in not_sig.iterrows():
                    f.write(f"  {row['feature']}: log2FC={row['log2FoldChange']:.2f}, "
                           f"padj={row['padj']:.3f}\n")
        else:
            f.write('No DE results available for markers\n')
        
        f.write('\n' + '='*80 + '\n')
        f.write('4. INTERPRETATION\n')
        f.write('='*80 + '\n\n')
        
        f.write('Pipeline Validation:\n')
        detection_rate = len(present_markers) / len(ALL_MARKERS) * 100
        if detection_rate >= 70:
            f.write(f'✓ PASS: {detection_rate:.1f}% of markers detected\n')
        else:
            f.write(f'⚠️  CONCERN: Only {detection_rate:.1f}% of markers detected\n')
        
        f.write('\nConsistency with Paper:\n')
        if marker_de is not None:
            sig_pct = (marker_de['padj'] < 0.05).sum() / len(marker_de) * 100
            if sig_pct < 20:
                f.write(f'✓ PASS: Only {sig_pct:.1f}% of markers are DE (matches paper\'s finding)\n')
            else:
                f.write(f'⚠️  CONCERN: {sig_pct:.1f}% of markers are DE (paper reports minimal DE)\n')
        
        f.write('\n' + '='*80 + '\n')
        f.write('FILES GENERATED\n')
        f.write('='*80 + '\n')
        f.write('- marker_expression_boxplots.png: Individual gene expression by condition\n')
        f.write('- marker_expression_heatmap.png: Mean expression heatmap\n')
        f.write('- marker_statistics.csv: Detailed statistics per gene per condition\n')
        f.write('- marker_de_results.csv: DE results for markers (if available)\n')
        f.write('- VALIDATION_REPORT.txt: This report\n')
    
    print(f"\nSaved report: {report_file}")

def main():
    print("="*80)
    print("MICROGLIA MARKER GENE VALIDATION")
    print("="*80)
    
    # Create output directory
    output_dir = Path('/home/jacobc/hcaTE/validation')
    output_dir.mkdir(exist_ok=True)
    
    # Load data
    print("\n[1/6] Loading data...")
    counts = load_pseudobulk_counts()
    sample_info = load_sample_info(counts)  # Pass counts to filter samples
    de_results = load_de_results()
    
    # Check marker presence
    print("\n[2/6] Checking marker presence...")
    present_markers, missing_markers = check_marker_presence(counts, ALL_MARKERS)
    
    if len(present_markers) == 0:
        print("\nERROR: No markers found in count matrix!")
        print("This suggests a problem with the pipeline or gene annotation.")
        sys.exit(1)
    
    # Extract expression data
    print("\n[3/6] Extracting marker expression...")
    marker_counts, marker_log2cpm, marker_data = get_marker_expression(
        counts, sample_info, present_markers
    )
    
    # Calculate statistics
    print("\n[4/6] Calculating statistics...")
    stats_df = calculate_marker_stats(marker_counts, sample_info)
    stats_df.to_csv(output_dir / 'marker_statistics.csv', index=False)
    print(f"Saved: {output_dir / 'marker_statistics.csv'}")
    
    # Check DE status
    print("\n[5/6] Checking differential expression status...")
    marker_de = check_de_status(de_results, present_markers)
    if marker_de is not None:
        marker_de.to_csv(output_dir / 'marker_de_results.csv', index=False)
        print(f"Saved: {output_dir / 'marker_de_results.csv'}")
    
    # Create plots
    print("\n[6/6] Creating visualizations...")
    plot_marker_expression(marker_data, present_markers, output_dir)
    
    # Create summary report
    print("\n[7/6] Creating summary report...")
    create_summary_report(stats_df, marker_de, present_markers, missing_markers, output_dir)
    
    print("\n" + "="*80)
    print("VALIDATION COMPLETE")
    print("="*80)
    print(f"\nResults saved to: {output_dir}")
    print("\nKey files:")
    print("  - VALIDATION_REPORT.txt: Summary of findings")
    print("  - marker_expression_boxplots.png: Expression plots")
    print("  - marker_expression_heatmap.png: Heatmap visualization")
    print("\nNext steps:")
    print("  1. Review VALIDATION_REPORT.txt to check pipeline quality")
    print("  2. If markers are present and not DE, pipeline is working correctly")
    print("  3. If markers are missing or highly DE, investigate pipeline issues")

if __name__ == '__main__':
    main()
