#!/usr/bin/env python3
"""
Run scTE on 10x Genomics aligned BAM files to quantify gene and TE expression
per cell.

Key differences from Smart-seq2:
- Uses CR (corrected cell barcode) instead of CB
- Uses UR (corrected UMI) instead of UB
- Much higher cell counts (~3,000 per donor)
- Lower reads per cell but more cells
"""

import os
import subprocess
import glob
from pathlib import Path
import time
import shutil

# ============================================================================
# CONFIGURATION
# ============================================================================

# Paths
BASE_DIR = "/home/jacobc/hcaTE"
BAM_DIR = f"{BASE_DIR}/10x_aligned_starsolo"
OUTPUT_DIR = f"{BASE_DIR}/10x_scTE_output"
LOG_DIR = f"{BASE_DIR}/qc/scTE_10x_logs"

# scTE resources (from HumanBam2scTE project)
SCTE_BIN = "/home/jacobc/HumanBam2scTE/.venv/bin/scTE"
SCTE_INDEX = "/home/jacobc/HumanBam2scTE/annotations/hg38.exclusive.exclusive.idx"

# scTE parameters for 10x data
CPU_CORES = 8  # Can use more for 10x since we have fewer samples
CB_TAG = "CR"  # 10x uses CR (corrected cell barcode)
UB_TAG = "UR"  # 10x uses UR (corrected UMI)
MIN_GENES = 200  # Minimum genes per cell for filtering
EXPECT_CELLS = 3000  # Expected cells per donor (~3,000 from STARsolo)

# Donors to process
DONORS = [
    "Donor2019-010",  # MCI, 3,052 cells
    "Donor2018-135"   # AD, 2,929 cells
]

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
    """Find 10x BAM files."""
    bam_files = []
    for donor in DONORS:
        bam_file = f"{BAM_DIR}/{donor}/Aligned.sortedByCoord.out.bam"
        
        if os.path.exists(bam_file):
            size_gb = os.path.getsize(bam_file) / 1e9
            print(f"  Found {donor}: {size_gb:.1f} GB")
            bam_files.append((donor, bam_file))
        else:
            print(f"  ✗ Missing: {donor}")
    
    return bam_files


def run_scte(donor_name, bam_file):
    """
    Run scTE on a single 10x BAM file.
    
    scTE will:
    - Extract cell barcodes from CR tag (corrected barcode)
    - Use UR tag for UMI deduplication (corrected UMI)
    - Count reads mapping to genes and TEs
    - Output per-cell count matrices
    """
    donor_output_dir = f"{OUTPUT_DIR}/{donor_name}"
    os.makedirs(donor_output_dir, exist_ok=True)
    
    log_file = f"{LOG_DIR}/{donor_name}_scTE.log"
    
    print(f"  Running scTE...")
    print(f"    Cell barcode tag: {CB_TAG}")
    print(f"    UMI tag: {UB_TAG}")
    print(f"    Expected cells: {EXPECT_CELLS}")
    print(f"    Min genes/cell: {MIN_GENES}")
    print(f"    Threads: {CPU_CORES}")
    
    # scTE needs to run from within the output directory
    # and use just the donor name as prefix
    cmd = [
        SCTE_BIN,
        "-i", bam_file,
        "-o", donor_name,  # Just the donor name, not full path
        "-x", SCTE_INDEX,
        "-CB", CB_TAG,      # CR for 10x (corrected barcode)
        "-UMI", UB_TAG,     # UR for 10x (corrected UMI)
        "--min_genes", str(MIN_GENES),
        "--expect-cells", str(EXPECT_CELLS),
        "-p", str(CPU_CORES)
    ]
    
    print(f"  Command: {' '.join(cmd)}")
    
    with open(log_file, 'w') as log:
        result = subprocess.run(
            cmd,
            cwd=donor_output_dir,  # Run from output directory
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True
        )
    
    if result.returncode != 0:
        print(f"  ✗ scTE failed (see log: {log_file})")
        # Print last 20 lines of log for debugging
        print("\n  Last 20 lines of log:")
        with open(log_file, 'r') as f:
            lines = f.readlines()
            for line in lines[-20:]:
                print(f"    {line.rstrip()}")
        return False
    
    # Check for output file (scTE creates a single CSV file)
    output_csv = f"{donor_output_dir}/{donor_name}.csv"
    
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
        
        # Check output file size
        size_mb = os.path.getsize(output_csv) / 1e6
        print(f"  ✓ Output file: {size_mb:.1f} MB")
        print(f"  ✓ scTE completed successfully")
        return True
    else:
        print(f"  ✗ scTE output file missing: {output_csv}")
        return False


def process_donor(donor_name, bam_file):
    """Process a single donor with scTE."""
    print(f"\n{'='*70}")
    print(f"Processing: {donor_name}")
    print(f"{'='*70}")
    
    start_time = time.time()
    
    success = run_scte(donor_name, bam_file)
    
    # Clean up scTE temp directory if it exists
    temp_dir = f"{OUTPUT_DIR}/{donor_name}/{donor_name}_scTEtmp"
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
            print(f"  Cleaned up temp directory")
        except:
            pass
    
    elapsed = time.time() - start_time
    print(f"  Time: {elapsed/60:.1f} minutes")
    
    return success


def get_completed_donors():
    """Check which donors already have scTE output."""
    completed = []
    if os.path.exists(OUTPUT_DIR):
        for donor in DONORS:
            output_csv = f"{OUTPUT_DIR}/{donor}/{donor}.csv"
            if os.path.exists(output_csv) and os.path.getsize(output_csv) > 1000:
                completed.append(donor)
    return set(completed)


def main():
    """Main pipeline."""
    print("\n" + "="*70)
    print("scTE Quantification Pipeline - 10x Genomics Data")
    print("Gene and TE expression per cell")
    print("="*70)
    print(f"Technology: 10x Genomics Chromium v2")
    print(f"Cell barcode tag: {CB_TAG} (corrected)")
    print(f"UMI tag: {UB_TAG} (corrected)")
    print(f"Threads: {CPU_CORES}")
    
    # Setup
    setup_directories()
    
    # Find BAM files
    print("\nFinding 10x BAM files...")
    bam_files = find_bam_files()
    
    if len(bam_files) == 0:
        print("ERROR: No BAM files found!")
        return
    
    # Check for already completed donors
    completed = get_completed_donors()
    if completed:
        print(f"\n✓ Found {len(completed)} already completed donors: {', '.join(completed)}")
        bam_files = [(name, bam) for name, bam in bam_files if name not in completed]
        print(f"→ {len(bam_files)} donors remaining to process")
    
    if len(bam_files) == 0:
        print("\nAll donors already completed!")
        return
    
    # Check scTE installation
    if not os.path.exists(SCTE_BIN):
        print(f"\nERROR: scTE not found at {SCTE_BIN}")
        return
    
    if not os.path.exists(SCTE_INDEX):
        print(f"ERROR: scTE index not found at {SCTE_INDEX}")
        return
    
    print(f"\n✓ scTE binary: {SCTE_BIN}")
    print(f"✓ scTE index: {SCTE_INDEX}")
    
    # STARsolo cell counts
    print("\nExpected cells from STARsolo:")
    print("  Donor2019-010 (MCI): 3,052 cells")
    print("  Donor2018-135 (AD):  2,929 cells")
    
    # Ask for confirmation
    print(f"\nWill process {len(bam_files)} donors:")
    for donor, _ in bam_files:
        print(f"  - {donor}")
    
    response = input(f"\nProceed? (y/n): ")
    if response.lower() != 'y':
        print("Aborted.")
        return
    
    # Process all donors
    total_start = time.time()
    successful = 0
    failed = 0
    
    for i, (donor_name, bam_file) in enumerate(bam_files, 1):
        print(f"\n[{i}/{len(bam_files)}]")
        
        try:
            if process_donor(donor_name, bam_file):
                successful += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    total_time = time.time() - total_start
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Total donors processed: {len(bam_files)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Already completed: {len(completed)}")
    print(f"Total time: {total_time/60:.1f} minutes ({total_time/3600:.1f} hours)")
    if len(bam_files) > 0:
        print(f"Average time: {total_time/len(bam_files)/60:.1f} minutes per donor")
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"Log directory: {LOG_DIR}")
    print("="*70)
    
    if successful > 0:
        print("\nNext steps:")
        print("1. Load count matrices into Scanpy/Seurat")
        print("2. QC filtering and normalization")
        print("3. Clustering and cell type identification")
        print("4. Compare with Smart-seq2 results")
        print("5. TE expression analysis")


if __name__ == "__main__":
    main()
