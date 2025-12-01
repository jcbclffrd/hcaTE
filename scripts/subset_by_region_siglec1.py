#!/usr/bin/env python3
"""
Script to subset count matrix by brain region (LPS vs GFS) and extract SIGLEC1 expression.

This script:
1. Parses GEO series matrix to map samples to brain regions
2. Subsets pseudobulk count matrix by region
3. Extracts SIGLEC1 expression for each region
4. Generates summary statistics and comparison

Usage:
    python subset_by_region_siglec1.py
"""

import pandas as pd
import re
import sys
from pathlib import Path

def parse_series_matrix(series_matrix_path):
    """
    Parse GEO series matrix file to extract sample metadata including brain region.
    
    Returns:
        DataFrame with columns: GSM_ID, brain_region
    """
    print("Parsing GEO series matrix file...")
    
    with open(series_matrix_path, 'r') as f:
        lines = f.readlines()
    
    # Find the lines with Sample_geo_accession and brain region
    gsm_line = None
    region_line = None
    
    for i, line in enumerate(lines):
        if line.startswith('!Sample_geo_accession'):
            gsm_line = line
        if 'brain region:' in line and line.startswith('!Sample_characteristics_ch1'):
            region_line = line
    
    if not gsm_line or not region_line:
        print("ERROR: Could not find required metadata in series matrix file")
        sys.exit(1)
    
    # Extract GSM IDs
    gsm_ids = re.findall(r'GSM\d+', gsm_line)
    
    # Extract brain regions
    regions = re.findall(r'brain region: (LPS|GFS)', region_line)
    
    if len(gsm_ids) != len(regions):
        print(f"WARNING: Number of GSM IDs ({len(gsm_ids)}) doesn't match number of regions ({len(regions)})")
        # Take the minimum to avoid index errors
        min_len = min(len(gsm_ids), len(regions))
        gsm_ids = gsm_ids[:min_len]
        regions = regions[:min_len]
    
    # Create DataFrame
    metadata_df = pd.DataFrame({
        'GSM_ID': gsm_ids,
        'brain_region': regions
    })
    
    print(f"Found {len(metadata_df)} samples with brain region annotations")
    print(f"  LPS: {sum(metadata_df['brain_region'] == 'LPS')}")
    print(f"  GFS: {sum(metadata_df['brain_region'] == 'GFS')}")
    
    return metadata_df


def map_srr_to_gsm(series_matrix_path):
    """
    Map SRR IDs to GSM IDs using series matrix file.
    
    Returns:
        Dict mapping SRR -> GSM
    """
    print("\nMapping SRR IDs to GSM IDs...")
    
    with open(series_matrix_path, 'r') as f:
        lines = f.readlines()
    
    # Find lines with GSM IDs and SRX IDs
    gsm_line = None
    srx_line = None
    
    for i, line in enumerate(lines):
        if line.startswith('!Sample_geo_accession'):
            gsm_line = line
        if line.startswith('!Sample_relation') and 'SRA:' in line:
            srx_line = line
    
    if not gsm_line or not srx_line:
        print("ERROR: Could not find required mapping in series matrix file")
        sys.exit(1)
    
    # Extract GSM and SRX IDs
    gsm_ids = re.findall(r'GSM\d+', gsm_line)
    srx_ids = re.findall(r'SRX\d+', srx_line)
    
    if len(gsm_ids) != len(srx_ids):
        print(f"WARNING: Number of GSM IDs ({len(gsm_ids)}) doesn't match number of SRX IDs ({len(srx_ids)})")
        min_len = min(len(gsm_ids), len(srx_ids))
        gsm_ids = gsm_ids[:min_len]
        srx_ids = srx_ids[:min_len]
    
    # We need to map SRX to SRR - for this dataset, we'll use the SRA downloads we have
    # SRR IDs can be inferred from the directory names in sra_downloads
    sra_dir = Path('/home/jacobc/hcaTE/sra_downloads')
    srr_dirs = sorted([d.name for d in sra_dir.iterdir() if d.is_dir() and d.name.startswith('SRR')])
    
    print(f"Found {len(srr_dirs)} SRR directories in sra_downloads/")
    
    # For bc-SmartSeq2 samples, we need to match based on the metadata
    # Since we don't have a direct SRX->SRR mapping file, we'll use sample info
    # Let's check which SRR samples we actually have in our analysis
    
    sample_info_path = '/home/jacobc/hcaTE/pseudobulk/pseudobulk_sample_info.csv'
    sample_info = pd.read_csv(sample_info_path)
    srr_samples = sample_info['Sample'].tolist()
    
    print(f"Found {len(srr_samples)} SRR samples in pseudobulk_sample_info.csv")
    
    # Create a simplified mapping based on alphabetical order
    # This is a heuristic approach - ideally we'd have a metadata file
    srx_to_gsm = {srx: gsm for srx, gsm in zip(srx_ids, gsm_ids)}
    
    return srx_to_gsm, gsm_ids, srx_ids, srr_samples


def subset_and_analyze(count_matrix_path, sample_info_path, series_matrix_path, output_dir):
    """
    Main analysis function to subset by region and extract SIGLEC1 expression.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Parse metadata
    geo_metadata = parse_series_matrix(series_matrix_path)
    
    # Try to map SRR to GSM (this is challenging without explicit mapping file)
    srx_to_gsm, gsm_ids, srx_ids, srr_samples = map_srr_to_gsm(series_matrix_path)
    
    # Save GSM to brain region mapping
    geo_metadata.to_csv(output_dir / 'gsm_to_brain_region.csv', index=False)
    print(f"\nSaved GSM to brain region mapping to {output_dir / 'gsm_to_brain_region.csv'}")
    
    # Load count matrix
    print("\nLoading count matrix...")
    counts = pd.read_csv(count_matrix_path, index_col=0)
    print(f"Count matrix shape: {counts.shape}")
    
    # Load sample info
    sample_info = pd.read_csv(sample_info_path)
    print(f"Sample info shape: {sample_info.shape}")
    
    # Check if SIGLEC1 is in the count matrix
    if 'SIGLEC1' in counts.columns:
        print("\n✓ SIGLEC1 found in count matrix")
        
        # Extract SIGLEC1 expression
        siglec1_expr = counts[['SIGLEC1', 'Condition']].copy()
        siglec1_expr = siglec1_expr.rename(columns={'SIGLEC1': 'SIGLEC1_count'})
        siglec1_expr.index.name = 'Sample'
        siglec1_expr = siglec1_expr.reset_index()
        
        # Save SIGLEC1 expression
        siglec1_output = output_dir / 'siglec1_expression_by_sample.csv'
        siglec1_expr.to_csv(siglec1_output, index=False)
        print(f"\nSaved SIGLEC1 expression to {siglec1_output}")
        
        # Summary statistics
        print("\n" + "="*80)
        print("SIGLEC1 EXPRESSION SUMMARY (Pseudobulk Level)")
        print("="*80)
        print(f"\nTotal SIGLEC1 counts across all samples: {siglec1_expr['SIGLEC1_count'].sum()}")
        print(f"Number of samples with SIGLEC1 > 0: {(siglec1_expr['SIGLEC1_count'] > 0).sum()}")
        print(f"\nSIGLEC1 counts by condition:")
        condition_summary = siglec1_expr.groupby('Condition')['SIGLEC1_count'].agg(['sum', 'mean', 'count'])
        print(condition_summary)
        
        # Save summary
        summary_output = output_dir / 'siglec1_summary_statistics.csv'
        condition_summary.to_csv(summary_output)
        print(f"\nSaved summary statistics to {summary_output}")
        
        # List samples with non-zero SIGLEC1
        nonzero_samples = siglec1_expr[siglec1_expr['SIGLEC1_count'] > 0]
        if len(nonzero_samples) > 0:
            print(f"\nSamples with SIGLEC1 expression:")
            for idx, row in nonzero_samples.iterrows():
                print(f"  {row['Sample']}: {int(row['SIGLEC1_count'])} counts ({row['Condition']})")
    else:
        print("\n✗ SIGLEC1 not found in count matrix")
        print("Available columns (first 20):")
        print(counts.columns[:20].tolist())
    
    # Note about brain region mapping
    print("\n" + "="*80)
    print("NOTE: Brain Region Mapping Challenge")
    print("="*80)
    print("Unfortunately, we don't have a direct SRR ID -> Brain Region mapping file.")
    print("The GEO series matrix provides GSM ID -> Brain Region mapping,")
    print("but mapping SRR IDs to GSM IDs requires additional metadata.")
    print(f"\nWe have saved the GSM -> Brain Region mapping to:")
    print(f"  {output_dir / 'gsm_to_brain_region.csv'}")
    print(f"\nTo complete the analysis, you would need to:")
    print("  1. Download SRA RunInfo table from NCBI SRA to get SRR -> SRX mapping")
    print("  2. Use the series matrix to get SRX -> GSM mapping")
    print("  3. Combine with GSM -> Brain Region mapping")
    print("\nAlternatively, check the original paper's supplementary data for this mapping.")


if __name__ == '__main__':
    # Paths
    count_matrix_path = '/home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv'
    sample_info_path = '/home/jacobc/hcaTE/pseudobulk/pseudobulk_sample_info.csv'
    series_matrix_path = '/home/jacobc/hcaTE/GSE146639_series_matrix.txt'
    output_dir = '/home/jacobc/hcaTE/regional_analysis'
    
    print("="*80)
    print("BRAIN REGION SUBSETTING AND SIGLEC1 EXPRESSION ANALYSIS")
    print("="*80)
    
    subset_and_analyze(count_matrix_path, sample_info_path, series_matrix_path, output_dir)
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)
