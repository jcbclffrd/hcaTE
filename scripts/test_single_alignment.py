#!/usr/bin/env python3
"""
Test STARsolo alignment on a single sample
Use this to verify the pipeline before running all 160 samples
"""

import os
import subprocess
import sys

# Configuration
SAMPLE_NAME = "SRR11271993"  # First sample for testing
RAM_GB = 120
CPU_CORES = 20
SRA_DIR = "/home/jacobc/hcaTE/sra_downloads"
OUTPUT_DIR = "/home/jacobc/hcaTE/aligned_bams"
GENOME_DIR = "/home/jacobc/hcaTE/star_index"

def test_single_sample():
    """Test alignment on a single sample."""
    
    print("=" * 70)
    print("TEST: Single Sample STARsolo Alignment")
    print("=" * 70)
    
    # Check genome index
    if not os.path.exists(f"{GENOME_DIR}/genomeParameters.txt"):
        print(f"ERROR: STAR genome index not found at {GENOME_DIR}")
        print("Please wait for the index to finish building.")
        sys.exit(1)
    
    print(f"Testing sample: {SAMPLE_NAME}")
    
    # Input files
    fastq_dir = f"{SRA_DIR}/{SAMPLE_NAME}"
    read1_barcodes = f"{fastq_dir}/{SAMPLE_NAME}_1.fastq.gz"
    read2_cDNA = f"{fastq_dir}/{SAMPLE_NAME}_2.fastq.gz"
    
    # Check inputs exist
    if not os.path.exists(read1_barcodes):
        print(f"ERROR: R1 file not found: {read1_barcodes}")
        sys.exit(1)
    if not os.path.exists(read2_cDNA):
        print(f"ERROR: R2 file not found: {read2_cDNA}")
        sys.exit(1)
    
    print("✓ Input files found")
    
    # Output directory
    output_dir = f"{OUTPUT_DIR}/{SAMPLE_NAME}_TEST"
    os.makedirs(output_dir, exist_ok=True)
    
    # STARsolo command
    cmd = [
        "STAR",
        "--runThreadN", str(CPU_CORES),
        "--genomeDir", GENOME_DIR,
        "--readFilesIn", read2_cDNA, read1_barcodes,
        "--readFilesCommand", "zcat",
        "--soloType", "SmartSeq",
        "--soloUMIlen", "3",
        "--soloCBlen", "16",
        "--soloStrand", "Unstranded",
        "--soloCBwhitelist", "None",
        "--soloFeatures", "Gene",
        "--outFileNamePrefix", f"{output_dir}/",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outBAMsortingThreadN", "6",
        "--outSAMattributes", "NH", "HI", "AS", "nM", "CB", "UB",
        "--limitBAMsortRAM", str(int(RAM_GB * 0.7 * 1e9)),
        "--outFilterMultimapNmax", "100",
        "--outFilterMismatchNmax", "999",
        "--outFilterMismatchNoverLmax", "0.04",
        "--winAnchorMultimapNmax", "200",
        "--outMultimapperOrder", "Random",
        "--outSAMunmapped", "None",
    ]
    
    print("\nRunning STAR...")
    print("Command:", " ".join(cmd[:10]), "...")
    
    # Run STAR
    result = subprocess.run(cmd)
    
    if result.returncode == 0:
        print("\n✓ STAR completed successfully")
        
        # Check output
        bam_file = f"{output_dir}/Aligned.sortedByCoord.out.bam"
        if os.path.exists(bam_file):
            print(f"✓ BAM file created: {bam_file}")
            
            # Check for CB tags
            print("\nChecking for CB tags...")
            cmd = f"samtools view {bam_file} | head -1000 | grep -o 'CB:Z:' | wc -l"
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            cb_count = int(result.stdout.strip())
            
            if cb_count > 0:
                print(f"✓ CB tags found in {cb_count}/1000 reads")
            else:
                print("✗ WARNING: No CB tags found!")
            
            # Check alignment stats
            log_file = f"{output_dir}/Log.final.out"
            if os.path.exists(log_file):
                print("\nAlignment statistics:")
                with open(log_file, 'r') as f:
                    for line in f:
                        if "Number of input reads" in line:
                            print(f"  {line.strip()}")
                        elif "Uniquely mapped reads %" in line:
                            print(f"  {line.strip()}")
                        elif "% of reads mapped to multiple loci" in line:
                            print(f"  {line.strip()}")
            
            print("\n" + "=" * 70)
            print("TEST SUCCESSFUL!")
            print("You can now run the full pipeline on all 160 samples.")
            print("=" * 70)
        else:
            print("✗ ERROR: BAM file not created")
    else:
        print(f"\n✗ STAR failed with exit code {result.returncode}")
        print(f"Check log at: {output_dir}/Log.out")

if __name__ == "__main__":
    test_single_sample()
