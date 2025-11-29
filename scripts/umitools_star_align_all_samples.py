#!/usr/bin/env python3
"""
bc-Smart-seq2 Alignment Pipeline using UMI-tools + STAR
For 160 pooled samples (~84 cells per sample)

Workflow:
1. UMI-tools extract: Extract 7bp UMI + 10bp cell barcode from R1, add to R2 read names
2. STAR align: Map extracted R2 reads to genome with TE annotations
3. Add CB/UB tags: Convert read name info to BAM tags for scTE compatibility
4. Validate: Check CB/UB tags present

This follows the original study's barcode extraction but uses STAR instead of HISAT2
for better TE mapping, and outputs CB/UB-tagged BAMs for scTE analysis.
"""

import os
import subprocess
import glob
from pathlib import Path
import time
import shutil
import sys
import re

# ============================================================================
# CONFIGURATION
# ============================================================================

RAM_GB = 120
CPU_CORES = 20
SRA_DIR = "/home/jacobc/hcaTE/sra_downloads"
OUTPUT_DIR = "/home/jacobc/hcaTE/aligned_bams"
GENOME_DIR = "/home/jacobc/hcaTE/star_index"
LOG_DIR = "/home/jacobc/hcaTE/qc/alignment_logs"
TEMP_DIR = "/home/jacobc/hcaTE/temp"

# Barcode structure (from GSE146639 methods)
UMI_LEN = 7  # First 7bp of R1
CB_LEN = 10  # Next 10bp of R1 (positions 8-17)
BARCODE_PATTERN = "NNNNNNNCCCCCCCCCC"  # 7N + 10C = UMI + CellBarcode (extra 2bp stay on read)

# Cell barcode whitelist (84 barcodes per pool)
WHITELIST_FILE = "/home/jacobc/hcaTE/annotations/cell_barcode_whitelist.txt"

# STAR parameters for TE analysis
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
    
    # Check STAR
    try:
        result = subprocess.run(["STAR", "--version"], 
                              capture_output=True, text=True)
        print(f"  ✓ STAR version: {result.stdout.strip()}")
    except FileNotFoundError:
        print("  ✗ ERROR: STAR not found")
        sys.exit(1)
    
    # Check samtools
    try:
        result = subprocess.run(["samtools", "--version"], 
                              capture_output=True, text=True, errors='replace')
        version = result.stdout.split('\n')[0]
        print(f"  ✓ samtools: {version}")
    except FileNotFoundError:
        print("  ✗ ERROR: samtools not found")
        sys.exit(1)
    
    # Check umi_tools
    try:
        result = subprocess.run(["/home/jacobc/hcaTE/.venv/bin/umi_tools", "--version"], 
                              capture_output=True, text=True, errors='replace')
        print(f"  ✓ umi_tools: {result.stderr.strip()}")
    except FileNotFoundError:
        print("  ✗ ERROR: umi_tools not found")
        sys.exit(1)
    
    # Check genome index
    if not os.path.exists(f"{GENOME_DIR}/genomeParameters.txt"):
        print(f"  ✗ ERROR: STAR genome index not found at {GENOME_DIR}")
        sys.exit(1)
    print(f"  ✓ STAR genome index found")
    
    # Check whitelist
    if not os.path.exists(WHITELIST_FILE):
        print(f"  ✗ ERROR: Cell barcode whitelist not found at {WHITELIST_FILE}")
        sys.exit(1)
    with open(WHITELIST_FILE) as f:
        num_barcodes = len(f.readlines())
    print(f"  ✓ Cell barcode whitelist: {num_barcodes} barcodes")
    
    # Check disk space
    stat = shutil.disk_usage("/home/jacobc/hcaTE")
    free_gb = stat.free / (1024**3)
    print(f"  ✓ Available disk space: {free_gb:.1f} GB")
    
    # Create directories
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(LOG_DIR, exist_ok=True)
    os.makedirs(TEMP_DIR, exist_ok=True)
    print(f"  ✓ Output directory: {OUTPUT_DIR}")
    print(f"  ✓ Temp directory: {TEMP_DIR}")


def get_available_samples(sra_dir):
    """Get list of samples with valid FASTQ pairs."""
    sample_dirs = sorted(glob.glob(f"{sra_dir}/SRR*"))
    valid_samples = []
    
    for sample_dir in sample_dirs:
        sample_name = Path(sample_dir).name
        r1 = f"{sample_dir}/{sample_name}_1.fastq.gz"
        r2 = f"{sample_dir}/{sample_name}_2.fastq.gz"
        
        if os.path.exists(r1) and os.path.exists(r2):
            # Check file sizes
            r1_size = os.path.getsize(r1)
            r2_size = os.path.getsize(r2)
            if r1_size > 1000000 and r2_size > 1000000:
                valid_samples.append(sample_name)
    
    return valid_samples


def extract_barcodes_umis(sample_name, temp_dir):
    """
    Extract cell barcodes and UMIs from R1 using umi_tools.
    Pattern: 7bp UMI + 10bp Cell Barcode
    """
    r1_file = f"{SRA_DIR}/{sample_name}/{sample_name}_1.fastq.gz"
    r2_file = f"{SRA_DIR}/{sample_name}/{sample_name}_2.fastq.gz"
    
    r2_extracted = f"{temp_dir}/{sample_name}_R2_extracted.fastq.gz"
    extract_log = f"{LOG_DIR}/{sample_name}_umi_extract.log"
    
    print(f"  Extracting barcodes and UMIs...")
    
    cmd = [
        "/home/jacobc/hcaTE/.venv/bin/umi_tools", "extract",
        "--bc-pattern", BARCODE_PATTERN,
        "--stdin", r1_file,
        "--read2-in", r2_file,
        "--read2-out", r2_extracted,
        "--whitelist", WHITELIST_FILE,
        "--extract-method", "string",
        "-L", extract_log
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0 or not os.path.exists(r2_extracted):
        print(f"  ✗ UMI extraction failed")
        print(f"    See log: {extract_log}")
        return None
    
    return r2_extracted


def align_with_star(sample_name, extracted_r2, output_dir):
    """
    Align extracted R2 reads with STAR.
    R2 read names now contain: @READID_CELLBARCODE_UMI
    """
    print(f"  Aligning with STAR...")
    
    sample_output_dir = f"{output_dir}/{sample_name}"
    os.makedirs(sample_output_dir, exist_ok=True)
    
    cmd = [
        "STAR",
        "--runThreadN", str(CPU_CORES),
        "--genomeDir", GENOME_DIR,
        "--readFilesIn", extracted_r2,
        "--readFilesCommand", "zcat",
        
        # Output settings
        "--outFileNamePrefix", f"{sample_output_dir}/",
        "--outSAMtype", "BAM", "Unsorted",  # Unsorted for now (will sort after tagging)
        "--outBAMsortingThreadN", "6",
        
        # Memory - use most of available RAM for faster processing
        "--limitBAMsortRAM", str(int(RAM_GB * 0.85 * 1e9)),  # 102GB for BAM sorting
        "--genomeLoad", "NoSharedMemory",  # Each sample loads own genome (safer)
        
        # TE-friendly alignment
        "--outFilterMultimapNmax", str(MULTIMAPPER_MAX),
        "--outFilterMismatchNmax", str(MISMATCH_MAX),
        "--outFilterMismatchNoverLmax", str(MISMATCH_RATIO),
        "--winAnchorMultimapNmax", str(ANCHOR_MULTIMAP_MAX),
        "--outMultimapperOrder", "Random",
        
        # Output all tags
        "--outSAMattributes", "NH", "HI", "AS", "nM",
        "--outSAMunmapped", "None",
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    bam_file = f"{sample_output_dir}/Aligned.out.bam"
    if result.returncode != 0 or not os.path.exists(bam_file):
        print(f"  ✗ STAR alignment failed")
        return None
    
    return bam_file


def add_cb_ub_tags(bam_file, sample_name):
    """
    Add CB (cell barcode) and UB (UMI) tags to BAM from read names.
    Read names from umi_tools have format: @READID_CELLBARCODE_UMI
    """
    print(f"  Adding CB and UB tags...")
    
    tagged_bam = bam_file.replace("Aligned.out.bam", "Aligned.tagged.bam")
    
    # Use samtools view to read, Python to add tags, samtools view to write
    cmd_read = f"samtools view -h {bam_file}"
    cmd_write = f"samtools view -b -o {tagged_bam}"
    
    # Process BAM file
    process_read = subprocess.Popen(cmd_read.split(), stdout=subprocess.PIPE, text=True)
    process_write = subprocess.Popen(cmd_write.split(), stdin=subprocess.PIPE, text=True)
    
    for line in process_read.stdout:
        if line.startswith('@'):
            # Header line - pass through
            process_write.stdin.write(line)
        else:
            # Alignment line - add tags
            fields = line.strip().split('\t')
            read_name = fields[0]
            
            # Parse read name: READID_CELLBARCODE_UMI
            parts = read_name.split('_')
            if len(parts) >= 3:
                cell_barcode = parts[-2]
                umi = parts[-1]
                
                # Add CB and UB tags
                fields.append(f"CB:Z:{cell_barcode}")
                fields.append(f"UB:Z:{umi}")
            
            process_write.stdin.write('\t'.join(fields) + '\n')
    
    process_read.stdout.close()
    process_write.stdin.close()
    process_read.wait()
    process_write.wait()
    
    if not os.path.exists(tagged_bam):
        print(f"  ✗ Failed to add tags")
        return None
    
    return tagged_bam


def sort_and_index_bam(bam_file, sample_name):
    """Sort and index the tagged BAM file."""
    print(f"  Sorting and indexing BAM...")
    
    sorted_bam = bam_file.replace(".bam", ".sorted.bam")
    
    # Sort - use more memory per thread for faster sorting
    # Total RAM usage: 20 threads × 5G = 100GB
    cmd_sort = [
        "samtools", "sort",
        "-@", str(CPU_CORES),
        "-m", "5G",  # 5GB per thread, 20 threads = 100GB total
        "-o", sorted_bam,
        bam_file
    ]
    
    result = subprocess.run(cmd_sort, capture_output=True)
    if result.returncode != 0:
        print(f"  ✗ Sorting failed")
        return None
    
    # Index
    cmd_index = ["samtools", "index", "-@", str(min(6, CPU_CORES)), sorted_bam]
    result = subprocess.run(cmd_index, capture_output=True)
    if result.returncode != 0:
        print(f"  ✗ Indexing failed")
        return None
    
    # Rename to final name
    final_bam = sorted_bam.replace("Aligned.tagged.sorted.bam", 
                                   f"{sample_name}_Aligned.sortedByCoord.out.bam")
    os.rename(sorted_bam, final_bam)
    os.rename(sorted_bam + ".bai", final_bam + ".bai")
    
    # Clean up intermediate files
    if os.path.exists(bam_file):
        os.remove(bam_file)
    
    return final_bam


def validate_bam(bam_file, sample_name):
    """Validate BAM has CB and UB tags and get stats."""
    if not os.path.exists(bam_file):
        return {"valid": False, "error": "BAM not found"}
    
    results = {"valid": True}
    
    # Check for CB tags
    cmd = f"samtools view {bam_file} | head -1000 | grep -c 'CB:Z:' || true"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    cb_count = int(result.stdout.strip()) if result.stdout.strip() else 0
    
    if cb_count > 0:
        results["cb_tags"] = f"✓ Present ({cb_count}/1000)"
    else:
        results["cb_tags"] = "✗ Missing"
        results["valid"] = False
    
    # Check for UB tags
    cmd = f"samtools view {bam_file} | head -1000 | grep -c 'UB:Z:' || true"
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    ub_count = int(result.stdout.strip()) if result.stdout.strip() else 0
    
    if ub_count > 0:
        results["ub_tags"] = f"✓ Present ({ub_count}/1000)"
    else:
        results["ub_tags"] = "✗ Missing"
        results["valid"] = False
    
    # Get alignment stats
    log_file = f"{OUTPUT_DIR}/{sample_name}/Log.final.out"
    if os.path.exists(log_file):
        with open(log_file) as f:
            for line in f:
                if "Uniquely mapped reads %" in line:
                    results["unique_pct"] = line.split('|')[1].strip()
                elif "% of reads mapped to multiple loci" in line:
                    results["multi_pct"] = line.split('|')[1].strip()
                elif "Number of input reads" in line:
                    results["total_reads"] = line.split('|')[1].strip()
    
    return results


def process_sample(sample_name, temp_dir):
    """Process a single sample through the complete pipeline."""
    sample_start = time.time()
    
    try:
        # Step 1: Extract barcodes and UMIs
        extracted_r2 = extract_barcodes_umis(sample_name, temp_dir)
        if not extracted_r2:
            return {"success": False, "error": "Barcode extraction failed"}
        
        # Step 2: Align with STAR
        bam_file = align_with_star(sample_name, extracted_r2, OUTPUT_DIR)
        if not bam_file:
            return {"success": False, "error": "STAR alignment failed"}
        
        # Step 3: Add CB/UB tags
        tagged_bam = add_cb_ub_tags(bam_file, sample_name)
        if not tagged_bam:
            return {"success": False, "error": "Tag addition failed"}
        
        # Step 4: Sort and index
        final_bam = sort_and_index_bam(tagged_bam, sample_name)
        if not final_bam:
            return {"success": False, "error": "Sorting failed"}
        
        # Step 5: Validate
        validation = validate_bam(final_bam, sample_name)
        
        # Cleanup temp files
        if os.path.exists(extracted_r2):
            os.remove(extracted_r2)
        
        return {
            "success": validation["valid"],
            "validation": validation,
            "time": time.time() - sample_start,
            "bam_file": final_bam
        }
    
    except Exception as e:
        return {
            "success": False,
            "error": str(e),
            "time": time.time() - sample_start
        }


def print_summary(results_list):
    """Print summary of all alignments."""
    successful = sum(1 for r in results_list if r["success"])
    failed = len(results_list) - successful
    
    print("\n" + "=" * 70)
    print("ALIGNMENT SUMMARY")
    print("=" * 70)
    print(f"Total samples (pools): {len(results_list)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    print(f"Expected cells: ~{successful * 84} (84 cells per pool)")
    
    if failed > 0:
        print("\nFailed samples:")
        for r in results_list:
            if not r["success"]:
                error = r.get("error", "Unknown error")
                print(f"  - {r['sample']}: {error}")
    
    # Average alignment rates
    unique_rates = []
    for r in results_list:
        if r["success"] and "validation" in r and "unique_pct" in r["validation"]:
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
    """Main pipeline."""
    print("=" * 70)
    print("bc-Smart-seq2 Alignment Pipeline")
    print("UMI-tools + STAR + CB/UB tagging for scTE")
    print("=" * 70)
    print(f"System: {RAM_GB}GB RAM, {CPU_CORES} cores")
    print(f"Barcode structure: {UMI_LEN}bp UMI + {CB_LEN}bp Cell Barcode")
    print(f"Expected: ~84 cells per sample/pool")
    
    # Check system
    check_system()
    
    # Get samples
    samples = get_available_samples(SRA_DIR)
    print(f"\nFound {len(samples)} samples (pools) to align")
    print(f"Expected total cells: ~{len(samples) * 84}")
    
    # Confirm
    response = input(f"\nBegin alignment? (y/n): ")
    if response.lower() != 'y':
        print("Aborted.")
        sys.exit(0)
    
    # Process samples
    results = []
    start_time = time.time()
    
    for i, sample_name in enumerate(samples, 1):
        print(f"\n{'=' * 70}")
        print(f"[{i}/{len(samples)}] {sample_name}")
        print('=' * 70)
        
        result = process_sample(sample_name, TEMP_DIR)
        result["sample"] = sample_name
        results.append(result)
        
        if result["success"]:
            val = result["validation"]
            print(f"\n  Validation:")
            print(f"    CB tags: {val.get('cb_tags', 'Unknown')}")
            print(f"    UB tags: {val.get('ub_tags', 'Unknown')}")
            if "unique_pct" in val:
                print(f"    Unique mapping: {val['unique_pct']}")
            if "multi_pct" in val:
                print(f"    Multi-mapping: {val['multi_pct']}")
        else:
            print(f"\n  ✗ Failed: {result.get('error', 'Unknown error')}")
        
        # Progress update
        if i % 10 == 0 or i == len(samples):
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            remaining = avg_time * (len(samples) - i)
            successful_so_far = sum(1 for r in results if r["success"])
            
            print(f"\n{'─' * 70}")
            print(f"Progress: {i}/{len(samples)} ({i*100//len(samples)}%)")
            print(f"Successful: {successful_so_far}/{i}")
            print(f"Average time: {avg_time/60:.1f} min/sample")
            if i < len(samples):
                print(f"ETA: {remaining/60:.1f} min ({remaining/3600:.1f} hrs)")
            print('─' * 70)
    
    # Final summary
    total_time = time.time() - start_time
    print_summary(results)
    print(f"\nTotal runtime: {total_time/3600:.2f} hours")
    print("=" * 70)
    
    # Save results
    results_file = f"{LOG_DIR}/alignment_results.txt"
    with open(results_file, 'w') as f:
        f.write("Sample\tSuccess\tUnique%\tMulti%\tCB_Tags\tUB_Tags\tBAM_File\n")
        for r in results:
            f.write(f"{r['sample']}\t{r['success']}\t")
            if "validation" in r:
                v = r["validation"]
                f.write(f"{v.get('unique_pct', 'N/A')}\t")
                f.write(f"{v.get('multi_pct', 'N/A')}\t")
                f.write(f"{v.get('cb_tags', 'N/A')}\t")
                f.write(f"{v.get('ub_tags', 'N/A')}\t")
                f.write(f"{r.get('bam_file', 'N/A')}\n")
            else:
                f.write("N/A\tN/A\tN/A\tN/A\tN/A\n")
    
    print(f"\nResults saved to: {results_file}")
    print("\nReady for scTE analysis!")


if __name__ == "__main__":
    main()
