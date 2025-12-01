#!/usr/bin/env python3
"""
Quick test: Check data structure of scTE output
"""

import pandas as pd
import glob

# Get first 3 samples
csv_files = sorted(glob.glob("/home/jacobc/hcaTE/scTE_output/SRR*/SRR*.csv"))[:3]

print("Testing scTE output structure with first 3 samples:\n")

for csv_file in csv_files:
    sample_id = csv_file.split('/')[-1].replace('.csv', '')
    df = pd.read_csv(csv_file, index_col=0)
    
    print(f"{sample_id}:")
    print(f"  Shape: {df.shape} (cells × features)")
    print(f"  Barcodes (first 3): {df.index[:3].tolist()}")
    print(f"  Features (first 5): {df.columns[:5].tolist()}")
    print(f"  Total counts per cell (first 3): {df.sum(axis=1).head(3).tolist()}")
    print()

print(f"Total samples available: {len(glob.glob('/home/jacobc/hcaTE/scTE_output/SRR*/SRR*.csv'))}")
