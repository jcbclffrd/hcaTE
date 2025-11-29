#!/usr/bin/env python3
"""
Aggregate single-cell scTE counts to pseudo-bulk sample-level counts

Since each bc-Smart-seq2 sample comes from a different individual,
we treat each sample as a biological replicate. We aggregate cells
within each sample to create pseudo-bulk samples.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys

# Configuration
SCTE_DIR = Path("/home/jacobc/hcaTE/scTE_output")
METADATA_FILE = Path("/home/jacobc/HumanBam2scTE/data/metadata.csv")
OUTPUT_DIR = Path("/home/jacobc/hcaTE/pseudobulk")
OUTPUT_DIR.mkdir(exist_ok=True)

def load_metadata():
    """Load sample metadata with conditions"""
    metadata = pd.read_csv(METADATA_FILE)
    print(f"Loaded metadata for {len(metadata)} samples")
    print(f"\nCondition distribution:")
    print(metadata['donor_group'].value_counts())
    return metadata

def load_scte_matrix(sample_id):
    """Load scTE output CSV for a sample"""
    csv_file = SCTE_DIR / sample_id / f"{sample_id}.csv"
    
    if not csv_file.exists():
        print(f"Warning: {csv_file} not found")
        return None
    
    # Load CSV: rows = cells, columns = genes/TEs
    df = pd.read_csv(csv_file, index_col=0)
    print(f"  {sample_id}: {df.shape[0]} cells × {df.shape[1]} features")
    
    return df

def aggregate_to_pseudobulk(metadata):
    """Aggregate cells within each sample to create pseudo-bulk samples"""
    
    # Dictionary to store pseudo-bulk counts for each sample
    pseudobulk_data = {}
    sample_conditions = {}
    
    print("\n" + "=" * 80)
    print("Aggregating single-cell data to pseudo-bulk (sample-level)...")
    print("=" * 80 + "\n")
    
    # Process each sample
    for idx, row in metadata.iterrows():
        sample_id = row['sample_id']
        condition = row['donor_group']
        
        print(f"Processing {sample_id} ({condition})")
        
        # Load single-cell matrix for this sample
        sc_matrix = load_scte_matrix(sample_id)
        if sc_matrix is None:
            continue
        
        # Sum across all cells in this sample (collapse cells → sample)
        sample_totals = sc_matrix.sum(axis=0)
        
        # Store as pseudo-bulk for this sample
        pseudobulk_data[sample_id] = sample_totals
        sample_conditions[sample_id] = condition
    
    print("\n" + "=" * 80)
    print(f"Aggregated {len(pseudobulk_data)} samples to pseudo-bulk")
    print("=" * 80)
    
    # Convert to DataFrame: rows = samples, columns = genes/TEs
    pseudobulk_df = pd.DataFrame(pseudobulk_data).T
    
    # Add condition column
    pseudobulk_df['Condition'] = pseudobulk_df.index.map(sample_conditions)
    
    # Reorder: Condition first, then all features
    cols = ['Condition'] + [c for c in pseudobulk_df.columns if c != 'Condition']
    pseudobulk_df = pseudobulk_df[cols]
    
    print("\nPseudo-bulk matrix shape:", pseudobulk_df.shape)
    print("\nSamples by condition:")
    print(pseudobulk_df['Condition'].value_counts())
    
    return pseudobulk_df

def save_pseudobulk_matrix(pseudobulk_df):
    """Save pseudo-bulk matrix in DESeq2-compatible format"""
    
    # Save as CSV (same format as expression_matrix2.csv)
    output_file = OUTPUT_DIR / "pseudobulk_expression_matrix.csv"
    pseudobulk_df.to_csv(output_file)
    print(f"\n✓ Saved: {output_file}")
    
    # Create sample info file for DESeq2
    sample_info = pd.DataFrame({
        'Sample': pseudobulk_df.index,
        'Condition': pseudobulk_df['Condition']
    })
    sample_info_file = OUTPUT_DIR / "pseudobulk_sample_info.csv"
    sample_info.to_csv(sample_info_file, index=False)
    print(f"✓ Saved: {sample_info_file}")
    
    # Summary statistics
    summary_file = OUTPUT_DIR / "pseudobulk_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PSEUDO-BULK AGGREGATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total samples: {len(pseudobulk_df)}\n")
        f.write(f"Total features: {len(pseudobulk_df.columns) - 1}\n\n")
        
        f.write("Samples by condition:\n")
        for cond, count in pseudobulk_df['Condition'].value_counts().items():
            f.write(f"  {cond}: {count} samples\n")
        
        f.write("\nTotal counts per sample:\n")
        total_counts = pseudobulk_df.drop('Condition', axis=1).sum(axis=1)
        f.write(f"  Min: {total_counts.min():,.0f}\n")
        f.write(f"  Max: {total_counts.max():,.0f}\n")
        f.write(f"  Mean: {total_counts.mean():,.0f}\n")
        f.write(f"  Median: {total_counts.median():,.0f}\n")
        
        f.write("\nDetailed counts by condition:\n")
        for cond in sorted(pseudobulk_df['Condition'].unique()):
            cond_counts = total_counts[pseudobulk_df['Condition'] == cond]
            f.write(f"\n  {cond} ({len(cond_counts)} samples):\n")
            f.write(f"    Min: {cond_counts.min():,.0f}\n")
            f.write(f"    Max: {cond_counts.max():,.0f}\n")
            f.write(f"    Mean: {cond_counts.mean():,.0f}\n")
            f.write(f"    Median: {cond_counts.median():,.0f}\n")
    
    print(f"✓ Saved: {summary_file}")

def main():
    print("\n" + "=" * 80)
    print("PSEUDO-BULK AGGREGATION")
    print("=" * 80)
    
    # Load metadata
    metadata = load_metadata()
    
    # Aggregate single-cell to pseudo-bulk
    pseudobulk_df = aggregate_to_pseudobulk(metadata)
    
    # Save outputs
    save_pseudobulk_matrix(pseudobulk_df)
    
    print("\n" + "=" * 80)
    print("✓ Pseudo-bulk aggregation complete!")
    print("=" * 80)
    print("\nOutput files:")
    print(f"  - {OUTPUT_DIR}/pseudobulk_expression_matrix.csv")
    print(f"  - {OUTPUT_DIR}/pseudobulk_sample_info.csv")
    print(f"  - {OUTPUT_DIR}/pseudobulk_summary.txt")
    print("\nNext steps:")
    print("  1. Review the summary file")
    print("  2. Run DESeq2 using the pseudo-bulk matrix")
    print("  3. Compare AD vs CTR/CTR+ conditions")

if __name__ == "__main__":
    main()
