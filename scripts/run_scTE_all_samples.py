#!/usr/bin/env python3
"""
Run scTE on all aligned BAM files to quantify gene and TE expression
per cell.

scTE will:
1. Use CB tags to identify cells
2. Use UB tags for UMI deduplication
3. Quantify both genes and TEs
4. Output count matrices for downstream analysis
"""

import os
import subprocess
import glob
from pathlib import Path
import time
import gc
import shutil

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths
BASE_DIR = "/home/jacobc/hcaTE"
BAM_DIR = f"{BASE_DIR}/aligned_bams"
OUTPUT_DIR = f"{BASE_DIR}/scTE_output"
LOG_DIR = f"{BASE_DIR}/qc/scTE_logs"

# scTE resources (from HumanBam2scTE project)
SCTE_BIN = "/home/jacobc/HumanBam2scTE/.venv/bin/scTE"
SCTE_INDEX = "/home/jacobc/HumanBam2scTE/annotations/hg38.exclusive.exclusive.idx"

# scTE parameters
CPU_CORES = 4  # Reduced from 20 - scTE loads index per thread (very memory intensive!)
CB_TAG = "CB"  # Cell barcode tag
UB_TAG = "UB"  # UMI barcode tag
MIN_GENES = 200  # Minimum genes per cell for filtering
EXPECT_CELLS = 100  # Expected cells per sample (~84 per pool)

# ============================================================================
# FUNCTIONS
# ============================================================================

def setup_directories():
    """Create output directories."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Log directory: {LOG_DIR}")


def find_bam_files():
    """Find all processed BAM files."""
    bam_files = []
    for sample_dir in sorted(glob.glob(f"{BAM_DIR}/SRR*")):
        sample_name = Path(sample_dir).name
        bam_file = f"{sample_dir}/{sample_name}_Aligned.sortedByCoord.out.bam"
        
        if os.path.exists(bam_file) and os.path.getsize(bam_file) > 1000000:
            bam_files.append((sample_name, bam_file))
    
    return bam_files


def run_scte(sample_name, bam_file):
    """
    Run scTE on a single BAM file.
    
    scTE will:
    - Extract cell barcodes from CB tag
    - Use UB tag for UMI deduplication
    - Count reads mapping to genes and TEs
    - Output per-cell count matrices
    """
    sample_output_dir = f"{OUTPUT_DIR}/{sample_name}"
    os.makedirs(sample_output_dir, exist_ok=True)
    
    log_file = f"{LOG_DIR}/{sample_name}_scTE.log"
    
    print(f"  Running scTE...")
    
    # scTE needs to run from within the output directory
    # and use just the sample name as prefix
    cmd = [
        SCTE_BIN,
        "-i", bam_file,
        "-o", sample_name,  # Just the sample name, not full path
        "-x", SCTE_INDEX,
        "-CB", CB_TAG,
        "-UMI", UB_TAG,
        "--min_genes", str(MIN_GENES),
        "--expect-cells", str(EXPECT_CELLS),
        "-p", str(CPU_CORES)
    ]
    
    with open(log_file, 'w') as log:
        result = subprocess.run(
            cmd,
            cwd=sample_output_dir,  # Run from output directory
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )
    
    if result.returncode != 0:
        print(f"  ✗ scTE failed (see log: {log_file})")
        return False
    
    # Check for output file (scTE creates a single CSV file)
    output_csv = f"{sample_output_dir}/{sample_name}.csv"
    
    if os.path.exists(output_csv):
        # Count detected cells
        with open(log_file, 'r') as f:
            log_content = f.read()
            if "Detect" in log_content:
                # Extract cell count from log
                import re
                match = re.search(r'Detect (\d+) cells', log_content)
                if match:
                    cell_count = match.group(1)
                    print(f"  ✓ Detected {cell_count} cells")
        
        print(f"  ✓ scTE completed successfully")
        return True
    else:
        print(f"  ✗ scTE output file missing: {output_csv}")
        return False


def process_sample(sample_name, bam_file):
    """Process a single sample with scTE."""
    print(f"\n{'='*70}")
    print(f"Processing: {sample_name}")
    print(f"{'='*70}")
    
    start_time = time.time()
    
    success = run_scte(sample_name, bam_file)
    
    # Force garbage collection and cleanup temp files after each sample
    gc.collect()
    
    # Clean up scTE temp directory if it exists
    temp_dir = f"{OUTPUT_DIR}/{sample_name}/{sample_name}_scTEtmp"
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            print(f"  Cleaned up temp directory")
        except:
            pass
    
    elapsed = time.time() - start_time
    print(f"  Time: {elapsed/60:.1f} minutes")
    
    return success


def get_completed_samples():
    """Check which samples already have scTE output."""
    completed = []
    if os.path.exists(OUTPUT_DIR):
        for sample_dir in glob.glob(f"{OUTPUT_DIR}/SRR*"):
            sample_name = Path(sample_dir).name
            output_csv = f"{sample_dir}/{sample_name}.csv"
            if os.path.exists(output_csv) and os.path.getsize(output_csv) > 1000:
                completed.append(sample_name)
    return set(completed)


def main():
    """Main pipeline."""
    print("\n" + "="*70)
    print("scTE Quantification Pipeline")
    print("Gene and TE expression per cell")
    print("="*70)
    print(f"⚠️  Using {CPU_CORES} threads (scTE is memory-intensive)")
    
    # Setup
    setup_directories()
    
    # Find BAM files
    print("\nFinding BAM files...")
    bam_files = find_bam_files()
    print(f"Found {len(bam_files)} BAM files to process")
    
    if len(bam_files) == 0:
        print("ERROR: No BAM files found!")
        return
    
    # Check for already completed samples
    completed = get_completed_samples()
    if completed:
        print(f"✓ Found {len(completed)} already completed samples")
        bam_files = [(name, bam) for name, bam in bam_files if name not in completed]
        print(f"→ {len(bam_files)} samples remaining to process")
    
    if len(bam_files) == 0:
        print("All samples already completed!")
        return
    
    # Check scTE installation
    if not os.path.exists(SCTE_BIN):
        print(f"ERROR: scTE not found at {SCTE_BIN}")
        return
    
    if not os.path.exists(SCTE_INDEX):
        print(f"ERROR: scTE index not found at {SCTE_INDEX}")
        return
    
    print(f"✓ scTE binary: {SCTE_BIN}")
    print(f"✓ scTE index: {SCTE_INDEX}")
    print(f"✓ Expected cells per sample: {EXPECT_CELLS}")
    print(f"✓ Threads per sample: {CPU_CORES}")
    
    # Ask for confirmation
    response = input(f"\nProcess {len(bam_files)} samples? (y/n): ")
    if response.lower() != 'y':
        print("Aborted.")
        return
    
    # Process all samples
    total_start = time.time()
    successful = 0
    failed = 0
    
    for i, (sample_name, bam_file) in enumerate(bam_files, 1):
        print(f"\n[{i}/{len(bam_files)}] {sample_name}")
        print(f"Progress: {successful} successful, {failed} failed")
        
        try:
            if process_sample(sample_name, bam_file):
                successful += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ✗ Error: {e}")
            failed += 1
        
        # Print estimated time remaining
        if i > 0:
            elapsed = time.time() - total_start
            avg_time = elapsed / i
            remaining = (len(bam_files) - i) * avg_time
            print(f"  ETA: {remaining/3600:.1f} hours remaining")
    
    total_time = time.time() - total_start
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total samples processed: {len(bam_files)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Already completed: {len(completed)}")
    print(f"Total time: {total_time/3600:.1f} hours")
    if len(bam_files) > 0:
        print(f"Average time: {total_time/len(bam_files)/60:.1f} minutes per sample")
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print("="*70)


if __name__ == "__main__":
    main()
