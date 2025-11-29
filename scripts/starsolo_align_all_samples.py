#!/usr/bin/env python3
"""
STARsolo Alignment for Smart-seq2 Single-Cell Data
Aligns 160 samples with 19bp barcodes + 63bp cDNA reads

Key features:
- Sequential processing (one sample at a time using all 20 cores)
- Cell barcode (CB) and UMI (UB) tagging for scTE compatibility
- Progress tracking and error handling
- BAM validation (checks for CB tags)
"""

import os
import subprocess
import glob
from pathlib import Path
import time
import shutil
import sys

# ============================================================================
# CONFIGURATION
# ============================================================================

RAM_GB = 120
CPU_CORES = 20
SRA_DIR = "/home/jacobc/hcaTE/sra_downloads"
OUTPUT_DIR = "/home/jacobc/hcaTE/aligned_bams"
GENOME_DIR = "/home/jacobc/hcaTE/star_index"
LOG_DIR = "/home/jacobc/hcaTE/qc/alignment_logs"

# Smart-seq2 barcode structure (from FASTQ analysis)
# 19bp reads: likely 16bp cell barcode + 3bp UMI/adapter
CELL_BARCODE_LENGTH = 16
UMI_LENGTH = 3

# STAR parameters for TE analysis (from HumanBam2scTE)
MULTIMAPPER_MAX = 100
MISMATCH_MAX = 999
MISMATCH_RATIO = 0.04
ANCHOR_MULTIMAP_MAX = 200

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def check_system():
    """Check system resources and prerequisites."""
    print("\nChecking system prerequisites...")
    
    # Check STAR installation
    try:
        result = subprocess.run(["STAR", "--version"], 
                              capture_output=True, text=True)
        print(f"  ✓ STAR version: {result.stdout.strip()}")
    except FileNotFoundError:
        print("  ✗ ERROR: STAR not found. Install with: sudo apt install rna-star")
        sys.exit(1)
    
    # Check samtools installation
    try:
        result = subprocess.run(["samtools", "--version"], 
                              capture_output=True, text=True)
        version = result.stdout.split('\n')[0]
        print(f"  ✓ samtools: {version}")
    except FileNotFoundError:
        print("  ✗ ERROR: samtools not found. Install with: sudo apt install samtools")
        sys.exit(1)
    
    # Check genome index
    if not os.path.exists(f"{GENOME_DIR}/genomeParameters.txt"):
        print(f"  ✗ ERROR: STAR genome index not found at {GENOME_DIR}")
        print("     Please build the index first!")
        sys.exit(1)
    print(f"  ✓ STAR genome index found")
    
    # Check disk space
    stat = shutil.disk_usage("/home/jacobc/hcaTE")
    free_gb = stat.free / (1024**3)
    print(f"  ✓ Available disk space: {free_gb:.1f} GB")
    
    # Estimate required space (5-10GB per sample)
    samples = get_available_samples(SRA_DIR)
    required_gb = len(samples) * 10
    if free_gb < required_gb:
        print(f"  ⚠ WARNING: May need ~{required_gb:.1f} GB, only {free_gb:.1f} GB available")
        response = input("    Continue anyway? (y/n): ")
        if response.lower() != 'y':
            sys.exit(1)
    
    # Create output directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    print(f"  ✓ Output directory: {OUTPUT_DIR}")
    print(f"  ✓ Log directory: {LOG_DIR}")


def get_available_samples(sra_dir):
    """Get list of samples with valid FASTQ pairs."""
    sample_dirs = sorted(glob.glob(f"{sra_dir}/SRR*"))
    valid_samples = []
    
    for sample_dir in sample_dirs:
        sample_name = Path(sample_dir).name
        r1 = f"{sample_dir}/{sample_name}_1.fastq.gz"
        r2 = f"{sample_dir}/{sample_name}_2.fastq.gz"
        
        if os.path.exists(r1) and os.path.exists(r2):
            # Check file sizes (basic validation)
            r1_size = os.path.getsize(r1)
            r2_size = os.path.getsize(r2)
            if r1_size > 1000000 and r2_size > 1000000:  # At least 1MB each
                valid_samples.append(sample_name)
            else:
                print(f"  ⚠ Skipping {sample_name}: files too small")
    
    return valid_samples


def align_sample_starsolo(sample_name):
    """
    Align single Smart-seq2 sample with STARsolo.
    
    Returns:
        bool: True if alignment successful, False otherwise
    """
    # Input files
    fastq_dir = f"{SRA_DIR}/{sample_name}"
    read1_barcodes = f"{fastq_dir}/{sample_name}_1.fastq.gz"  # 19bp barcodes
    read2_cDNA = f"{fastq_dir}/{sample_name}_2.fastq.gz"      # 63bp cDNA
    
    # Output directory
    output_dir = f"{OUTPUT_DIR}/{sample_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # Log file
    log_file = f"{LOG_DIR}/{sample_name}_alignment.log"
    
    # STARsolo command
    # NOTE: For Smart-seq2, we're using a simplified approach
    # Each sample is essentially one cell, so we'll use the sample name as the barcode
    cmd = [
        "STAR",
        "--runThreadN", str(CPU_CORES),
        "--genomeDir", GENOME_DIR,
        
        # CRITICAL: R2 (cDNA) first, R1 (barcodes) second!
        "--readFilesIn", read2_cDNA, read1_barcodes,
        "--readFilesCommand", "zcat",
        
        # Single-cell configuration - SmartSeq mode
        "--soloType", "SmartSeq",
        "--soloUMIlen", str(UMI_LENGTH),
        "--soloCBlen", str(CELL_BARCODE_LENGTH),
        "--soloStrand", "Unstranded",
        
        # No whitelist for Smart-seq2 (each sample is one cell)
        "--soloCBwhitelist", "None",
        
        # Solo features - count genes
        "--soloFeatures", "Gene",
        
        # Output settings
        "--outFileNamePrefix", f"{output_dir}/",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outBAMsortingThreadN", str(min(6, CPU_CORES)),
        
        # CRITICAL: Cell barcode and UMI tags for scTE!
        "--outSAMattributes", "NH", "HI", "AS", "nM", "CB", "UB",
        
        # Memory optimization
        "--limitBAMsortRAM", str(int(RAM_GB * 0.7 * 1e9)),
        
        # TE-friendly alignment parameters (from HumanBam2scTE)
        "--outFilterMultimapNmax", str(MULTIMAPPER_MAX),
        "--outFilterMismatchNmax", str(MISMATCH_MAX),
        "--outFilterMismatchNoverLmax", str(MISMATCH_RATIO),
        "--winAnchorMultimapNmax", str(ANCHOR_MULTIMAP_MAX),
        "--outMultimapperOrder", "Random",
        
        # Performance and cleanup
        "--outSAMunmapped", "None",  # Don't output unmapped reads
        "--outReadsUnmapped", "None",
    ]
    
    print(f"  Running STARsolo...")
    
    # Run STAR with logging
    with open(log_file, 'w') as log:
        result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    
    if result.returncode != 0:
        print(f"  ✗ STAR failed (exit code {result.returncode})")
        print(f"    See log: {log_file}")
        return False
    
    # Check for output BAM
    bam_file = f"{output_dir}/Aligned.sortedByCoord.out.bam"
    if not os.path.exists(bam_file):
        print(f"  ✗ BAM file not created")
        return False
    
    # Rename BAM for clarity
    renamed_bam = f"{output_dir}/{sample_name}_Aligned.sortedByCoord.out.bam"
    if os.path.exists(renamed_bam):
        os.remove(renamed_bam)
    os.rename(bam_file, renamed_bam)
    
    # Index BAM
    print(f"  Indexing BAM...")
    try:
        subprocess.run(["samtools", "index", renamed_bam], 
                      check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"  ✗ BAM indexing failed")
        return False
    
    print(f"  ✓ Alignment complete")
    return True


def validate_bam(bam_file, sample_name):
    """
    Validate BAM file has CB and UB tags.
    
    Returns:
        dict: Validation results with alignment stats
    """
    if not os.path.exists(bam_file):
        return {"valid": False, "error": "BAM file not found"}
    
    results = {"valid": True}
    
    # Check for CB tags
    cmd = f"samtools view {bam_file} | head -1000 | grep -o 'CB:Z:' | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    cb_count = int(result.stdout.strip())
    
    if cb_count > 0:
        results["cb_tags"] = f"✓ Present ({cb_count}/1000 reads)"
    else:
        results["cb_tags"] = "✗ Missing"
        results["valid"] = False
    
    # Check for UB tags
    cmd = f"samtools view {bam_file} | head -1000 | grep -o 'UB:Z:' | wc -l"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    ub_count = int(result.stdout.strip())
    
    if ub_count > 0:
        results["ub_tags"] = f"✓ Present ({ub_count}/1000 reads)"
    else:
        results["ub_tags"] = "✗ Missing"
        results["valid"] = False
    
    # Get alignment stats from STAR log
    log_final = f"{OUTPUT_DIR}/{sample_name}/Log.final.out"
    if os.path.exists(log_final):
        with open(log_final, 'r') as f:
            for line in f:
                if "Uniquely mapped reads %" in line:
                    results["unique_pct"] = line.split('|')[1].strip()
                elif "% of reads mapped to multiple loci" in line:
                    results["multi_pct"] = line.split('|')[1].strip()
                elif "Number of input reads" in line:
                    results["total_reads"] = line.split('|')[1].strip()
    
    return results


def print_summary(results_list):
    """Print summary of all alignments."""
    successful = sum(1 for r in results_list if r["success"])
    failed = len(results_list) - successful
    
    print("\n" + "=" * 70)
    print("ALIGNMENT SUMMARY")
    print("=" * 70)
    print(f"Total samples: {len(results_list)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    
    if failed > 0:
        print("\nFailed samples:")
        for r in results_list:
            if not r["success"]:
                print(f"  - {r['sample']}")
    
    # Calculate average alignment rates
    unique_rates = []
    for r in results_list:
        if r["success"] and "unique_pct" in r["validation"]:
            try:
                rate = float(r["validation"]["unique_pct"].strip('%'))
                unique_rates.append(rate)
            except:
                pass
    
    if unique_rates:
        avg_unique = sum(unique_rates) / len(unique_rates)
        print(f"\nAverage unique mapping rate: {avg_unique:.2f}%")
        print(f"Range: {min(unique_rates):.2f}% - {max(unique_rates):.2f}%")


def main():
    """Main alignment pipeline."""
    print("=" * 70)
    print("STARsolo Alignment for Smart-seq2 Single-Cell Data")
    print("=" * 70)
    print(f"System: {RAM_GB}GB RAM, {CPU_CORES} CPU cores")
    print(f"Processing: Sequential (one sample at a time)")
    
    # Check system
    check_system()
    
    # Get samples
    samples = get_available_samples(SRA_DIR)
    print(f"\nFound {len(samples)} samples to align")
    
    # Confirm before starting
    response = input(f"\nBegin alignment of {len(samples)} samples? (y/n): ")
    if response.lower() != 'y':
        print("Aborted.")
        sys.exit(0)
    
    # Align sequentially
    results = []
    start_time = time.time()
    
    for i, sample_name in enumerate(samples, 1):
        print(f"\n{'=' * 70}")
        print(f"[{i}/{len(samples)}] {sample_name}")
        print('=' * 70)
        
        sample_start = time.time()
        success = False
        
        try:
            # Align sample
            if align_sample_starsolo(sample_name):
                # Validate BAM
                bam_file = f"{OUTPUT_DIR}/{sample_name}/{sample_name}_Aligned.sortedByCoord.out.bam"
                validation = validate_bam(bam_file, sample_name)
                
                if validation["valid"]:
                    print(f"\n  Validation:")
                    print(f"    CB tags: {validation.get('cb_tags', 'Unknown')}")
                    print(f"    UB tags: {validation.get('ub_tags', 'Unknown')}")
                    if "unique_pct" in validation:
                        print(f"    Unique mapping: {validation['unique_pct']}")
                    if "multi_pct" in validation:
                        print(f"    Multi-mapping: {validation['multi_pct']}")
                    success = True
                else:
                    print(f"\n  ✗ Validation failed")
                
                results.append({
                    "sample": sample_name,
                    "success": success,
                    "validation": validation,
                    "time": time.time() - sample_start
                })
            else:
                results.append({
                    "sample": sample_name,
                    "success": False,
                    "error": "Alignment failed",
                    "time": time.time() - sample_start
                })
        
        except Exception as e:
            print(f"  ✗ ERROR: {e}")
            results.append({
                "sample": sample_name,
                "success": False,
                "error": str(e),
                "time": time.time() - sample_start
            })
        
        # Progress update every 10 samples
        if i % 10 == 0 or i == len(samples):
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            remaining = avg_time * (len(samples) - i)
            successful_so_far = sum(1 for r in results if r["success"])
            
            print(f"\n{'─' * 70}")
            print(f"Progress: {i}/{len(samples)} ({i*100//len(samples)}%)")
            print(f"Successful: {successful_so_far}/{i}")
            print(f"Average time per sample: {avg_time/60:.1f} minutes")
            if i < len(samples):
                print(f"Estimated time remaining: {remaining/60:.1f} minutes ({remaining/3600:.1f} hours)")
            print('─' * 70)
    
    # Final summary
    total_time = time.time() - start_time
    print_summary(results)
    print(f"\nTotal runtime: {total_time/3600:.2f} hours")
    print("=" * 70)
    
    # Save results to file
    results_file = f"{LOG_DIR}/alignment_results.txt"
    with open(results_file, 'w') as f:
        f.write("Sample\tSuccess\tUnique%\tMulti%\tCB_Tags\tUB_Tags\n")
        for r in results:
            f.write(f"{r['sample']}\t")
            f.write(f"{r['success']}\t")
            if "validation" in r:
                f.write(f"{r['validation'].get('unique_pct', 'N/A')}\t")
                f.write(f"{r['validation'].get('multi_pct', 'N/A')}\t")
                f.write(f"{r['validation'].get('cb_tags', 'N/A')}\t")
                f.write(f"{r['validation'].get('ub_tags', 'N/A')}\n")
            else:
                f.write("N/A\tN/A\tN/A\tN/A\n")
    
    print(f"\nResults saved to: {results_file}")


if __name__ == "__main__":
    main()
