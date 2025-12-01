#!/usr/bin/env python3
"""
Load and combine single-cell scTE data from all samples
Step 1 of single-cell clustering analysis

scTE output format:
- Each sample file contains ~80 cells (rows)
- Each cell identified by its 10bp barcode (e.g., AAATACACTC)
- ~76,394 features (columns): genes + TEs
- 158 total samples = ~12,640 total cells

Goal: Combine all 158 samples into one AnnData object for Scanpy
"""

import pandas as pd
import numpy as np
import scanpy as sc
import os
from pathlib import Path
import glob

print("=" * 80)
print("Single-Cell Data Loading - Step 1")
print("=" * 80)

# Paths
SCTE_OUTPUT_DIR = "/home/jacobc/hcaTE/scTE_output"
SAMPLE_INFO_PATH = "/home/jacobc/hcaTE/pseudobulk/pseudobulk_sample_info.csv"
OUTPUT_DIR = "/home/jacobc/hcaTE/clustering"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load sample metadata
print("\n1. Loading sample metadata...")
sample_info = pd.read_csv(SAMPLE_INFO_PATH)
print(f"   Loaded {len(sample_info)} samples with metadata")
print(f"   Conditions: {sample_info['condition'].value_counts().to_dict()}")

# Find all scTE output files
print("\n2. Finding scTE output files...")
csv_files = glob.glob(f"{SCTE_OUTPUT_DIR}/SRR*/SRR*.csv")
print(f"   Found {len(csv_files)} scTE output files")

# Load and combine all samples
print("\n3. Loading single-cell data from all samples...")
all_dfs = []
cell_metadata = []

for i, csv_file in enumerate(csv_files, 1):
    # Extract sample ID from path
    sample_id = Path(csv_file).stem
    
    if i % 20 == 0:
        print(f"   Processing sample {i}/{len(csv_files)}: {sample_id}")
    
    # Load sample
    df = pd.read_csv(csv_file, index_col=0)
    
    # Add sample ID to cell barcodes to make them unique
    df.index = [f"{sample_id}_{barcode}" for barcode in df.index]
    
    # Store dataframe
    all_dfs.append(df)
    
    # Get condition for this sample
    sample_condition = sample_info[sample_info['sample_id'] == sample_id]['condition'].values
    condition = sample_condition[0] if len(sample_condition) > 0 else "Unknown"
    
    # Create metadata for each cell in this sample
    for barcode in df.index:
        cell_metadata.append({
            'cell_id': barcode,
            'sample_id': sample_id,
            'condition': condition,
            'barcode': barcode.split('_', 1)[1]  # Original barcode without sample ID
        })

print(f"\n   Loaded {len(all_dfs)} samples")
print(f"   Total cells before combining: {sum(len(df) for df in all_dfs)}")

# Combine all samples into one matrix
print("\n4. Combining all samples into single matrix...")
combined_matrix = pd.concat(all_dfs, axis=0)
print(f"   Combined matrix shape: {combined_matrix.shape}")
print(f"   Total cells: {combined_matrix.shape[0]}")
print(f"   Total features: {combined_matrix.shape[1]}")

# Create metadata DataFrame
print("\n5. Creating cell metadata...")
cell_meta_df = pd.DataFrame(cell_metadata)
cell_meta_df = cell_meta_df.set_index('cell_id')
print(f"   Metadata shape: {cell_meta_df.shape}")

# Verify indices match
assert all(combined_matrix.index == cell_meta_df.index), "Cell IDs don't match!"

# Create AnnData object
print("\n6. Creating AnnData object...")
adata = sc.AnnData(
    X=combined_matrix.values,
    obs=cell_meta_df,
    var=pd.DataFrame(index=combined_matrix.columns)
)

# Add feature annotations (separate genes from TEs)
print("\n7. Annotating features (genes vs TEs)...")
adata.var['feature_name'] = adata.var.index
adata.var['feature_type'] = ['TE' if '#' in name else 'Gene' for name in adata.var.index]

# For TEs, extract class and family
def parse_te_name(name):
    if '#' not in name:
        return {'te_name': None, 'te_class': None, 'te_family': None}
    parts = name.split('#')
    return {
        'te_name': parts[0] if len(parts) > 0 else None,
        'te_class': parts[1] if len(parts) > 1 else None,
        'te_family': parts[2] if len(parts) > 2 else None
    }

te_info = adata.var.index.map(parse_te_name).to_list()
adata.var['te_name'] = [x['te_name'] for x in te_info]
adata.var['te_class'] = [x['te_class'] for x in te_info]
adata.var['te_family'] = [x['te_family'] for x in te_info]

print(f"   Genes: {(adata.var['feature_type'] == 'Gene').sum()}")
print(f"   TEs: {(adata.var['feature_type'] == 'TE').sum()}")

# Add basic QC metrics
print("\n8. Computing basic QC metrics...")
adata.obs['n_counts'] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
adata.obs['pct_te'] = (adata[:, adata.var['feature_type'] == 'TE'].X.sum(axis=1) / 
                         adata.obs['n_counts']) * 100

print(f"\n   QC Summary:")
print(f"   Mean counts per cell: {adata.obs['n_counts'].mean():.0f}")
print(f"   Median counts per cell: {adata.obs['n_counts'].median():.0f}")
print(f"   Mean genes per cell: {adata.obs['n_genes'].mean():.0f}")
print(f"   Median genes per cell: {adata.obs['n_genes'].median():.0f}")
print(f"   Mean % TE expression: {adata.obs['pct_te'].mean():.2f}%")

# Summary by condition
print(f"\n   Cells by condition:")
for condition, count in adata.obs['condition'].value_counts().items():
    print(f"   {condition}: {count} cells")

# Save AnnData object
print("\n9. Saving AnnData object...")
output_file = f"{OUTPUT_DIR}/raw_adata.h5ad"
adata.write_h5ad(output_file)
print(f"   Saved to: {output_file}")

# Save cell metadata separately
cell_meta_file = f"{OUTPUT_DIR}/cell_metadata.csv"
adata.obs.to_csv(cell_meta_file)
print(f"   Cell metadata saved to: {cell_meta_file}")

# Save feature metadata
feature_meta_file = f"{OUTPUT_DIR}/feature_metadata.csv"
adata.var.to_csv(feature_meta_file)
print(f"   Feature metadata saved to: {feature_meta_file}")

print("\n" + "=" * 80)
print("✓ Step 1 Complete: Single-cell data loaded and combined!")
print("=" * 80)
print(f"\nNext steps:")
print(f"  1. Review QC metrics in {cell_meta_file}")
print(f"  2. Run quality control filtering (Step 2)")
print(f"  3. Proceed to normalization and clustering")
