# Reference: Parallel Alignment Script Template

## Previous Bulk RNA-seq Alignment (HumanBam2scTE)

This project previously had a parallel alignment script for bulk RNA-seq. **DO NOT USE THIS DIRECTLY** because it's for bulk RNA-seq, not single-cell!

However, the **parallel processing structure** is useful. See `/home/jacobc/HumanBam2scTE/scripts/align_downloaded_samples.py` for reference.

## Key Features to Adapt

### 1. Sequential Processing (Not Parallel)
The script processes samples **one at a time** using full system resources:
- **RAM**: 120GB available (use ~80% for STAR)
- **CPU**: 20 cores (use all for each sample)
- **Strategy**: Process samples sequentially, not in parallel

This is **better than parallel** because:
- STAR is memory-intensive (needs 30-40GB per run)
- Single-threaded STAR with 20 cores is faster than 2 parallel runs with 10 cores each
- Simpler error handling and logging

### 2. Progress Tracking
```python
for i, sample_name in enumerate(samples, 1):
    print(f"[{i}/{len(samples)}] Aligning {sample_name}...")
    
    # ... alignment code ...
    
    # Progress update every 10 samples
    if i % 10 == 0:
        elapsed = time.time() - start_time
        avg_time = elapsed / i
        remaining = avg_time * (len(samples) - i)
        print(f"  Progress: {i}/{len(samples)} ({i*100//len(samples)}%)")
        print(f"  Estimated time remaining: {remaining/60:.1f} minutes\n")
```

### 3. Sample Discovery
```python
def get_available_samples(sra_dir):
    """Get list of samples with valid FASTQ pairs."""
    sample_dirs = sorted(glob.glob(f"{sra_dir}/SRR*"))
    valid_samples = []
    for sample_dir in sample_dirs:
        r1 = f"{sample_dir}/{Path(sample_dir).name}_1.fastq.gz"
        r2 = f"{sample_dir}/{Path(sample_dir).name}_2.fastq.gz"
        if os.path.exists(r1) and os.path.exists(r2):
            valid_samples.append(Path(sample_dir).name)
    return valid_samples
```

### 4. Error Handling
```python
successful = 0
failed = 0
failed_samples = []

for sample in samples:
    try:
        if align_sample(sample):
            successful += 1
        else:
            failed += 1
            failed_samples.append(sample)
    except Exception as e:
        print(f"ERROR: {e}")
        failed += 1
        failed_samples.append(sample)

# Report failures at end
if failed_samples:
    print(f"\nFailed samples:")
    for sample in failed_samples:
        print(f"  - {sample}")
```

## Adaptations for STARsolo

### Critical Differences

1. **Command structure**: Use `--soloType SmartSeq` instead of standard STAR
2. **Input order**: R2 (cDNA) first, then R1 (barcodes): `--readFilesIn R2 R1`
3. **Solo parameters**: Add barcode/UMI configuration
4. **Output tags**: Must include CB (cell barcode) and UB (UMI) tags

### STARsolo Command Template

```python
def align_sample_starsolo(sample_name, genome_dir, output_base, cpu_cores=20, ram_gb=120):
    """Align single-cell sample with STARsolo."""
    
    # Input files
    fastq_dir = f"/home/jacobc/hcaTE/sra_downloads/{sample_name}"
    read1_barcodes = f"{fastq_dir}/{sample_name}_1.fastq.gz"  # 19bp barcodes
    read2_cDNA = f"{fastq_dir}/{sample_name}_2.fastq.gz"      # 63bp cDNA
    
    # Output directory
    output_dir = f"{output_base}/{sample_name}"
    os.makedirs(output_dir, exist_ok=True)
    
    # STARsolo command
    cmd = [
        "STAR",
        "--runThreadN", str(cpu_cores),
        "--genomeDir", genome_dir,
        
        # IMPORTANT: R2 (cDNA) first, R1 (barcodes) second!
        "--readFilesIn", read2_cDNA, read1_barcodes,
        "--readFilesCommand", "zcat",
        
        # Single-cell configuration
        "--soloType", "SmartSeq",
        "--soloUMIlen", "10",  # Adjust based on barcode analysis
        "--soloCBlen", "16",   # Adjust based on barcode analysis
        "--soloStrand", "Unstranded",
        
        # Barcode whitelist (if you have one)
        # "--soloCBwhitelist", "/path/to/whitelist.txt",  # or "None"
        
        # Output settings
        "--outFileNamePrefix", f"{output_dir}/{sample_name}_",
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outBAMsortingThreadN", str(min(6, cpu_cores)),
        
        # Cell barcode and UMI tags (CRITICAL for scTE!)
        "--outSAMattributes", "NH", "HI", "AS", "nM", "CB", "UB",
        
        # Memory optimization
        "--limitBAMsortRAM", str(int(ram_gb * 0.8 * 1e9)),
        
        # TE-friendly alignment (from previous project)
        "--outFilterMultimapNmax", "100",
        "--outFilterMismatchNmax", "999",
        "--outFilterMismatchNoverLmax", "0.04",
        "--winAnchorMultimapNmax", "200",
        "--outMultimapperOrder", "Random",
        
        # Performance
        "--outSAMunmapped", "None",
    ]
    
    # Run STAR
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        # Index BAM
        bam_file = f"{output_dir}/{sample_name}_Aligned.sortedByCoord.out.bam"
        if os.path.exists(bam_file):
            subprocess.run(["samtools", "index", bam_file], check=True)
            return True
    
    return False
```

## Validation After Alignment

Add a validation function to check BAM files have proper tags:

```python
def validate_bam(bam_file):
    """Check if BAM has CB and UB tags."""
    cmd = ["samtools", "view", bam_file, "|", "head", "-1000", "|", "grep", "-o", "CB:Z:", "|", "wc", "-l"]
    result = subprocess.run(" ".join(cmd), shell=True, capture_output=True, text=True)
    
    cb_count = int(result.stdout.strip())
    if cb_count > 0:
        print(f"  ✓ BAM has CB tags ({cb_count}/1000 reads checked)")
        return True
    else:
        print(f"  ✗ WARNING: No CB tags found in BAM")
        return False
```

## System Resources (from HumanBam2scTE)

- **RAM**: 120GB total
- **CPU**: 20 cores
- **Storage**: Check available space first!

```python
# At start of script
RAM_GB = 120
CPU_CORES = 20

# Check available disk space
import shutil
stat = shutil.disk_usage("/home/jacobc/hcaTE")
free_gb = stat.free / (1024**3)
print(f"Available disk space: {free_gb:.1f} GB")

# Estimate: ~5-10GB per sample for BAM + temp files
required_gb = len(samples) * 10
if free_gb < required_gb:
    print(f"WARNING: May need {required_gb:.1f} GB, only {free_gb:.1f} GB available")
```

## Recommended Script Structure

```python
#!/usr/bin/env python3
"""
STARsolo Alignment for Smart-seq2 Single-Cell Data
160 samples with 19bp barcodes + 63bp cDNA
"""

import os
import subprocess
import glob
from pathlib import Path
import time
import shutil

# Configuration
RAM_GB = 120
CPU_CORES = 20
SRA_DIR = "/home/jacobc/hcaTE/sra_downloads"
OUTPUT_DIR = "/home/jacobc/hcaTE/aligned_bams"
GENOME_DIR = "/home/jacobc/hcaTE/genome/star_index"

def check_system():
    """Check system resources and prerequisites."""
    # Check STAR
    # Check samtools
    # Check disk space
    # Check genome index
    pass

def get_available_samples(sra_dir):
    """Get list of Smart-seq2 samples."""
    pass

def align_sample_starsolo(sample_name):
    """Align single sample with STARsolo."""
    pass

def validate_bam(bam_file):
    """Validate BAM has CB and UB tags."""
    pass

def main():
    print("=" * 70)
    print("STARsolo Alignment for Smart-seq2 Single-Cell Data")
    print(f"System: {RAM_GB}GB RAM, {CPU_CORES} CPU cores")
    print("=" * 70)
    
    # Check system
    check_system()
    
    # Get samples
    samples = get_available_samples(SRA_DIR)
    print(f"Found {len(samples)} samples to align")
    
    # Align sequentially
    successful = 0
    failed = 0
    start_time = time.time()
    
    for i, sample in enumerate(samples, 1):
        print(f"\n[{i}/{len(samples)}] Aligning {sample}...")
        
        if align_sample_starsolo(sample):
            bam_file = f"{OUTPUT_DIR}/{sample}/{sample}_Aligned.sortedByCoord.out.bam"
            if validate_bam(bam_file):
                successful += 1
            else:
                failed += 1
        else:
            failed += 1
        
        # Progress update
        if i % 10 == 0:
            elapsed = time.time() - start_time
            avg_time = elapsed / i
            remaining = avg_time * (len(samples) - i)
            print(f"Progress: {i}/{len(samples)} - ETA: {remaining/60:.1f} min")
    
    # Summary
    print(f"\n{'=' * 70}")
    print(f"Complete! {successful}/{len(samples)} successful")
    print(f"Total time: {(time.time() - start_time)/3600:.2f} hours")

if __name__ == "__main__":
    main()
```

---

## Tell the Agent

**"Yes, please create an alignment script now. Use the structure from `/home/jacobc/HumanBam2scTE/scripts/align_downloaded_samples.py` as a reference for parallel processing and progress tracking, but adapt it for STARsolo single-cell alignment. See `/home/jacobc/hcaTE/REFERENCE_parallel_alignment.md` for guidance on what to adapt."**
