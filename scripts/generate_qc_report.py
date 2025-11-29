#!/usr/bin/env python3
"""
Generate comprehensive QC report for bc-Smart-seq2 alignment and scTE quantification
"""

import os
import csv
import glob
import re
from pathlib import Path
from collections import defaultdict

# Paths
BASE_DIR = Path("/home/jacobc/hcaTE")
ALIGNED_DIR = BASE_DIR / "aligned_bams"
SCTE_DIR = BASE_DIR / "scTE_output"
OUTPUT_DIR = BASE_DIR / "qc"
OUTPUT_DIR.mkdir(exist_ok=True)

def parse_star_log(log_file):
    """Parse STAR Log.final.out file for alignment metrics"""
    metrics = {}
    try:
        with open(log_file) as f:
            for line in f:
                line = line.strip()
                if '|' in line:
                    key, value = line.split('|')
                    key = key.strip()
                    value = value.strip()
                    metrics[key] = value
    except FileNotFoundError:
        return None
    return metrics

def parse_scte_output(csv_file):
    """Parse scTE CSV output for cell and feature counts"""
    metrics = {}
    try:
        with open(csv_file) as f:
            reader = csv.reader(f)
            header = next(reader)
            
            cells = 0
            genes_features = len(header) - 1  # -1 for barcode column
            total_counts = 0
            non_zero = 0
            counts_per_cell = []
            features_per_cell = []
            
            for row in reader:
                cells += 1
                cell_counts = sum(float(x) if x else 0 for x in row[1:])
                cell_features = sum(1 for x in row[1:] if x and float(x) > 0)
                total_counts += cell_counts
                counts_per_cell.append(cell_counts)
                features_per_cell.append(cell_features)
                non_zero += cell_features
            
            if cells > 0:
                metrics['cells'] = cells
                metrics['total_features'] = genes_features
                metrics['total_counts'] = int(total_counts)
                metrics['avg_counts_per_cell'] = int(total_counts / cells)
                metrics['avg_features_per_cell'] = int(non_zero / cells)
                metrics['min_counts_per_cell'] = int(min(counts_per_cell)) if counts_per_cell else 0
                metrics['max_counts_per_cell'] = int(max(counts_per_cell)) if counts_per_cell else 0
                metrics['median_features_per_cell'] = int(sorted(features_per_cell)[cells//2]) if features_per_cell else 0
                
    except (FileNotFoundError, StopIteration, ValueError) as e:
        print(f"Error parsing {csv_file}: {e}")
        return None
    return metrics

def main():
    print("Generating comprehensive QC report...")
    print("=" * 80)
    
    # Find all samples
    sample_dirs = sorted(glob.glob(str(ALIGNED_DIR / "SRR*")))
    samples = [Path(d).name for d in sample_dirs]
    
    print(f"Found {len(samples)} samples")
    
    # Collect all metrics
    all_metrics = []
    
    for sample in samples:
        print(f"Processing {sample}...", end=' ')
        
        # Parse STAR alignment log
        star_log = ALIGNED_DIR / sample / "Log.final.out"
        star_metrics = parse_star_log(star_log)
        
        # Parse scTE output
        scte_csv = SCTE_DIR / sample / f"{sample}.csv"
        scte_metrics = parse_scte_output(scte_csv)
        
        if star_metrics and scte_metrics:
            # Extract key STAR metrics
            total_reads = star_metrics.get("Number of input reads", "0")
            unique_pct = star_metrics.get("Uniquely mapped reads %", "0%")
            multi_pct = star_metrics.get("% of reads mapped to multiple loci", "0%")
            unmapped_pct = star_metrics.get("% of reads unmapped: too short", "0%")
            
            all_metrics.append({
                'sample': sample,
                'total_reads': int(total_reads.replace(',', '')) if total_reads != "0" else 0,
                'unique_pct': float(unique_pct.rstrip('%')),
                'multi_pct': float(multi_pct.rstrip('%')),
                'unmapped_pct': float(unmapped_pct.rstrip('%')),
                **scte_metrics
            })
            print("✓")
        else:
            print("✗ (missing data)")
    
    print("\n" + "=" * 80)
    print(f"Successfully processed {len(all_metrics)} / {len(samples)} samples")
    
    # Write detailed TSV report
    output_file = OUTPUT_DIR / "alignment_scTE_qc_report.tsv"
    with open(output_file, 'w') as f:
        # Header
        f.write("sample\ttotal_reads\tunique_pct\tmulti_pct\tunmapped_pct\t")
        f.write("cells\ttotal_features\ttotal_counts\tavg_counts_per_cell\t")
        f.write("avg_features_per_cell\tmedian_features_per_cell\t")
        f.write("min_counts_per_cell\tmax_counts_per_cell\n")
        
        # Data
        for m in all_metrics:
            f.write(f"{m['sample']}\t{m['total_reads']}\t{m['unique_pct']:.2f}\t")
            f.write(f"{m['multi_pct']:.2f}\t{m['unmapped_pct']:.2f}\t")
            f.write(f"{m['cells']}\t{m['total_features']}\t{m['total_counts']}\t")
            f.write(f"{m['avg_counts_per_cell']}\t{m['avg_features_per_cell']}\t")
            f.write(f"{m['median_features_per_cell']}\t{m['min_counts_per_cell']}\t")
            f.write(f"{m['max_counts_per_cell']}\n")
    
    print(f"\nDetailed report written to: {output_file}")
    
    # Calculate summary statistics
    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    
    # Alignment stats
    total_reads = sum(m['total_reads'] for m in all_metrics)
    avg_unique = sum(m['unique_pct'] for m in all_metrics) / len(all_metrics)
    avg_multi = sum(m['multi_pct'] for m in all_metrics) / len(all_metrics)
    avg_unmapped = sum(m['unmapped_pct'] for m in all_metrics) / len(all_metrics)
    
    print("\nALIGNMENT METRICS:")
    print(f"  Total samples:          {len(all_metrics)}")
    print(f"  Total reads:            {total_reads:,}")
    print(f"  Avg unique mapped:      {avg_unique:.2f}%")
    print(f"  Avg multi-mapped:       {avg_multi:.2f}%")
    print(f"  Avg unmapped:           {avg_unmapped:.2f}%")
    print(f"  Total mapped rate:      {avg_unique + avg_multi:.2f}%")
    
    # scTE stats
    total_cells = sum(m['cells'] for m in all_metrics)
    avg_cells = total_cells / len(all_metrics)
    total_umis = sum(m['total_counts'] for m in all_metrics)
    avg_umis_per_cell = sum(m['avg_counts_per_cell'] for m in all_metrics) / len(all_metrics)
    avg_features_per_cell = sum(m['avg_features_per_cell'] for m in all_metrics) / len(all_metrics)
    median_features_per_cell = sum(m['median_features_per_cell'] for m in all_metrics) / len(all_metrics)
    
    print("\nscTE QUANTIFICATION METRICS:")
    print(f"  Total cells detected:   {total_cells:,}")
    print(f"  Avg cells/sample:       {avg_cells:.1f}")
    print(f"  Total UMI counts:       {total_umis:,}")
    print(f"  Avg UMIs/cell:          {avg_umis_per_cell:,.0f}")
    print(f"  Avg features/cell:      {avg_features_per_cell:.0f}")
    print(f"  Median features/cell:   {median_features_per_cell:.0f}")
    
    # Cell count distribution
    cell_counts = [m['cells'] for m in all_metrics]
    print(f"\nCELL COUNT DISTRIBUTION:")
    print(f"  Min cells/sample:       {min(cell_counts)}")
    print(f"  Max cells/sample:       {max(cell_counts)}")
    print(f"  Expected (84):          {sum(1 for c in cell_counts if c == 84)} samples")
    print(f"  Below expected (<84):   {sum(1 for c in cell_counts if c < 84)} samples")
    
    # Write summary report
    summary_file = OUTPUT_DIR / "qc_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("BC-SMART-SEQ2 PIPELINE QC REPORT\n")
        f.write("=" * 80 + "\n\n")
        
        f.write("ALIGNMENT METRICS:\n")
        f.write(f"  Total samples:          {len(all_metrics)}\n")
        f.write(f"  Total reads:            {total_reads:,}\n")
        f.write(f"  Avg unique mapped:      {avg_unique:.2f}%\n")
        f.write(f"  Avg multi-mapped:       {avg_multi:.2f}%\n")
        f.write(f"  Avg unmapped:           {avg_unmapped:.2f}%\n")
        f.write(f"  Total mapped rate:      {avg_unique + avg_multi:.2f}%\n\n")
        
        f.write("scTE QUANTIFICATION METRICS:\n")
        f.write(f"  Total cells detected:   {total_cells:,}\n")
        f.write(f"  Avg cells/sample:       {avg_cells:.1f}\n")
        f.write(f"  Total UMI counts:       {total_umis:,}\n")
        f.write(f"  Avg UMIs/cell:          {avg_umis_per_cell:,.0f}\n")
        f.write(f"  Avg features/cell:      {avg_features_per_cell:.0f}\n")
        f.write(f"  Median features/cell:   {median_features_per_cell:.0f}\n\n")
        
        f.write("CELL COUNT DISTRIBUTION:\n")
        f.write(f"  Min cells/sample:       {min(cell_counts)}\n")
        f.write(f"  Max cells/sample:       {max(cell_counts)}\n")
        f.write(f"  Expected (84):          {sum(1 for c in cell_counts if c == 84)} samples\n")
        f.write(f"  Below expected (<84):   {sum(1 for c in cell_counts if c < 84)} samples\n")
    
    print(f"\nSummary report written to: {summary_file}")
    print("\n" + "=" * 80)
    print("QC report generation complete!")

if __name__ == "__main__":
    main()
