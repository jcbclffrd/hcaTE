#!/usr/bin/env python3
"""
Visualize SIGLEC1 expression across samples and conditions.

This script creates simple visualizations of SIGLEC1 expression patterns
to help understand the extremely low counts observed in pseudobulk analysis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def create_siglec1_plots(output_dir):
    """Create visualization plots for SIGLEC1 expression."""
    output_dir = Path(output_dir)
    
    # Load SIGLEC1 expression data
    siglec1_df = pd.read_csv(output_dir / 'siglec1_expression_by_sample.csv')
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('SIGLEC1 Expression Analysis (Pseudobulk Level)', fontsize=16, fontweight='bold')
    
    # Plot 1: Bar plot of SIGLEC1 counts by sample
    ax1 = axes[0, 0]
    samples_with_counts = siglec1_df[siglec1_df['SIGLEC1_count'] > 0]
    if len(samples_with_counts) > 0:
        samples_with_counts.plot(x='Sample', y='SIGLEC1_count', kind='bar', ax=ax1, 
                                  color=['#2ecc71' if c=='CTR+' else '#e74c3c' for c in samples_with_counts['Condition']])
        ax1.set_title('Samples with SIGLEC1 Expression (n=2)', fontweight='bold')
        ax1.set_xlabel('Sample ID')
        ax1.set_ylabel('SIGLEC1 Count')
        ax1.legend(['CTR+', 'MCI'], title='Condition')
        ax1.grid(axis='y', alpha=0.3)
    else:
        ax1.text(0.5, 0.5, 'No samples with SIGLEC1 > 0', 
                ha='center', va='center', transform=ax1.transAxes)
        ax1.set_title('Samples with SIGLEC1 Expression')
    
    # Plot 2: Distribution by condition
    ax2 = axes[0, 1]
    condition_summary = siglec1_df.groupby('Condition')['SIGLEC1_count'].sum()
    colors = {'AD': '#e74c3c', 'CTR': '#3498db', 'CTR+': '#2ecc71', 'MCI': '#f39c12'}
    condition_summary.plot(kind='bar', ax=ax2, color=[colors[c] for c in condition_summary.index])
    ax2.set_title('Total SIGLEC1 Counts by Condition', fontweight='bold')
    ax2.set_xlabel('Condition')
    ax2.set_ylabel('Total SIGLEC1 Count')
    ax2.set_xticklabels(condition_summary.index, rotation=0)
    ax2.grid(axis='y', alpha=0.3)
    
    # Plot 3: Proportion of samples with expression
    ax3 = axes[1, 0]
    condition_counts = siglec1_df.groupby('Condition').size()
    condition_with_expr = siglec1_df[siglec1_df['SIGLEC1_count'] > 0].groupby('Condition').size()
    proportion_df = pd.DataFrame({
        'Total Samples': condition_counts,
        'With SIGLEC1': condition_with_expr
    }).fillna(0)
    proportion_df['Without SIGLEC1'] = proportion_df['Total Samples'] - proportion_df['With SIGLEC1']
    
    proportion_df[['With SIGLEC1', 'Without SIGLEC1']].plot(kind='bar', stacked=True, ax=ax3, 
                                                              color=['#2ecc71', '#ecf0f1'])
    ax3.set_title('Samples with/without SIGLEC1 Expression', fontweight='bold')
    ax3.set_xlabel('Condition')
    ax3.set_ylabel('Number of Samples')
    ax3.set_xticklabels(proportion_df.index, rotation=0)
    ax3.legend(['With SIGLEC1', 'Without SIGLEC1'])
    ax3.grid(axis='y', alpha=0.3)
    
    # Plot 4: Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis('off')
    summary_stats = siglec1_df.groupby('Condition').agg({
        'SIGLEC1_count': ['sum', 'mean', 'count']
    }).round(3)
    summary_stats.columns = ['Total Counts', 'Mean', '# Samples']
    
    # Create table
    table_data = []
    table_data.append(['Condition', 'Total Counts', 'Mean', '# Samples'])
    for idx, row in summary_stats.iterrows():
        table_data.append([idx, int(row['Total Counts']), f"{row['Mean']:.3f}", int(row['# Samples'])])
    
    table = ax4.table(cellText=table_data, cellLoc='center', loc='center', 
                      colWidths=[0.25, 0.25, 0.25, 0.25])
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2)
    
    # Style header row
    for i in range(4):
        table[(0, i)].set_facecolor('#34495e')
        table[(0, i)].set_text_props(weight='bold', color='white')
    
    # Alternate row colors
    for i in range(1, len(table_data)):
        for j in range(4):
            if i % 2 == 0:
                table[(i, j)].set_facecolor('#ecf0f1')
    
    ax4.set_title('Summary Statistics by Condition', fontweight='bold', pad=20)
    
    plt.tight_layout()
    
    # Save plot
    output_path = output_dir / 'siglec1_expression_plots.png'
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nVisualization saved to: {output_path}")
    
    plt.close()
    
    # Additional: Create a simple text report
    create_text_report(siglec1_df, output_dir)


def create_text_report(siglec1_df, output_dir):
    """Create a detailed text report of SIGLEC1 findings."""
    report_path = output_dir / 'SIGLEC1_REPORT.txt'
    
    with open(report_path, 'w') as f:
        f.write("="*80 + "\n")
        f.write("SIGLEC1 EXPRESSION ANALYSIS REPORT\n")
        f.write("Pseudobulk RNA-seq from bc-SmartSeq2 Data\n")
        f.write("="*80 + "\n\n")
        
        f.write("OVERVIEW\n")
        f.write("-"*80 + "\n")
        f.write(f"Total samples analyzed: {len(siglec1_df)}\n")
        f.write(f"Total SIGLEC1 counts: {siglec1_df['SIGLEC1_count'].sum()}\n")
        f.write(f"Samples with SIGLEC1 > 0: {(siglec1_df['SIGLEC1_count'] > 0).sum()}\n")
        f.write(f"Percentage with expression: {(siglec1_df['SIGLEC1_count'] > 0).sum() / len(siglec1_df) * 100:.1f}%\n\n")
        
        f.write("CONDITION BREAKDOWN\n")
        f.write("-"*80 + "\n")
        condition_summary = siglec1_df.groupby('Condition').agg({
            'SIGLEC1_count': ['sum', 'mean', 'count']
        })
        condition_summary.columns = ['Total_Counts', 'Mean', 'N_Samples']
        f.write(condition_summary.to_string() + "\n\n")
        
        f.write("SAMPLES WITH SIGLEC1 EXPRESSION\n")
        f.write("-"*80 + "\n")
        nonzero = siglec1_df[siglec1_df['SIGLEC1_count'] > 0]
        if len(nonzero) > 0:
            for idx, row in nonzero.iterrows():
                f.write(f"  {row['Sample']}: {int(row['SIGLEC1_count'])} counts ({row['Condition']})\n")
        else:
            f.write("  No samples with SIGLEC1 expression\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("INTERPRETATION\n")
        f.write("="*80 + "\n\n")
        
        f.write("The extremely low SIGLEC1 counts observed in this pseudobulk analysis are\n")
        f.write("expected due to:\n\n")
        f.write("1. CELL-TYPE DILUTION: Each sample aggregates ~80-100 cells from multiple cell\n")
        f.write("   types (microglia, neurons, astrocytes, oligodendrocytes, etc.)\n\n")
        f.write("2. MICROGLIA SPECIFICITY: SIGLEC1 is primarily expressed in microglia, which\n")
        f.write("   typically represent 5-15% of total brain cells\n\n")
        f.write("3. PSEUDOBULK AVERAGING: Even if SIGLEC1 is highly expressed in microglia,\n")
        f.write("   the signal gets averaged across all cell types in the sample\n\n")
        f.write("RECOMMENDATION:\n")
        f.write("-"*80 + "\n")
        f.write("To properly analyze SIGLEC1 expression as shown in the paper's Figure 2:\n\n")
        f.write("1. Perform SINGLE-CELL level analysis using scTE output files\n")
        f.write("2. SUBSET to microglia cells using marker genes (CX3CR1, P2RY12, TMEM119)\n")
        f.write("3. Analyze SIGLEC1 expression within MICROGLIA by BRAIN REGION (LPS vs GFS)\n")
        f.write("4. This will eliminate the cell-type dilution effect and match the paper's\n")
        f.write("   analytical approach\n\n")
        
        f.write("="*80 + "\n")
    
    print(f"Text report saved to: {report_path}")


if __name__ == '__main__':
    import sys
    
    output_dir = '/home/jacobc/hcaTE/regional_analysis'
    
    print("="*80)
    print("SIGLEC1 EXPRESSION VISUALIZATION")
    print("="*80)
    
    try:
        create_siglec1_plots(output_dir)
        print("\n" + "="*80)
        print("Visualization complete!")
        print("="*80)
    except Exception as e:
        print(f"Error creating visualizations: {e}")
        print("\nThis is likely due to missing matplotlib/seaborn packages.")
        print("To install: pip install matplotlib seaborn")
        sys.exit(1)
