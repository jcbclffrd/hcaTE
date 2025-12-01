#!/bin/bash

# Script: download_10x_genomics_fastq.sh
# Purpose: Download 10x Genomics FASTQ files for Alsema et al. (2020) study
# Date: December 1, 2024
# Branch: feature/10x-genomics-analysis

set -e  # Exit on error

echo "========================================"
echo "10x Genomics FASTQ Download Script"
echo "GSE146639 - Alsema et al. (2020)"
echo "========================================"
echo ""

# Configuration
THREADS=20
OUTPUT_DIR="10x_downloads/raw"

# Create output directories
mkdir -p ${OUTPUT_DIR}/Donor2019-010
mkdir -p ${OUTPUT_DIR}/Donor2018-135

echo "Output directory: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}"
echo ""

# Check if fasterq-dump is installed
if ! command -v fasterq-dump &> /dev/null; then
    echo "ERROR: fasterq-dump not found!"
    echo "Please install SRA Toolkit:"
    echo "  conda install -c bioconda sra-tools"
    exit 1
fi

# Check if pigz is installed
if ! command -v pigz &> /dev/null; then
    echo "WARNING: pigz not found, using gzip instead (slower)"
    COMPRESS_CMD="gzip"
else
    COMPRESS_CMD="pigz --processes ${THREADS}"
fi

echo "========================================"
echo "Downloading Donor 2019-010 (MCI)"
echo "GSM4403285 - 4 runs"
echo "========================================"
echo ""

# Donor 2019-010 (MCI) - 4 runs
DONOR1_RUNS=("SRR11272192" "SRR11272193" "SRR11272194" "SRR11272195")

for SRR in "${DONOR1_RUNS[@]}"; do
    echo "----------------------------------------"
    echo "Downloading: $SRR"
    echo "----------------------------------------"
    
    fasterq-dump $SRR \
        --outdir ${OUTPUT_DIR}/Donor2019-010/ \
        --split-files \
        --threads ${THREADS} \
        --progress
    
    echo "Compressing $SRR..."
    ${COMPRESS_CMD} ${OUTPUT_DIR}/Donor2019-010/${SRR}*.fastq
    
    echo "✓ Completed $SRR"
    echo ""
done

echo "========================================"
echo "Downloading Donor 2018-135 (AD)"
echo "GSM4403286 - 8 runs"
echo "========================================"
echo ""

# Donor 2018-135 (AD) - 8 runs
DONOR2_RUNS=("SRR11272196" "SRR11272197" "SRR11272198" "SRR11272199" 
             "SRR11272200" "SRR11272201" "SRR11272202" "SRR11272203")

for SRR in "${DONOR2_RUNS[@]}"; do
    echo "----------------------------------------"
    echo "Downloading: $SRR"
    echo "----------------------------------------"
    
    fasterq-dump $SRR \
        --outdir ${OUTPUT_DIR}/Donor2018-135/ \
        --split-files \
        --threads ${THREADS} \
        --progress
    
    echo "Compressing $SRR..."
    ${COMPRESS_CMD} ${OUTPUT_DIR}/Donor2018-135/${SRR}*.fastq
    
    echo "✓ Completed $SRR"
    echo ""
done

echo "========================================"
echo "Download Summary"
echo "========================================"
echo ""
echo "Donor 2019-010 (MCI) files:"
ls -lh ${OUTPUT_DIR}/Donor2019-010/*.fastq.gz 2>/dev/null || echo "No files found"
echo ""
echo "Donor 2018-135 (AD) files:"
ls -lh ${OUTPUT_DIR}/Donor2018-135/*.fastq.gz 2>/dev/null || echo "No files found"
echo ""

echo "========================================"
echo "Next Steps"
echo "========================================"
echo ""
echo "1. Verify file integrity:"
echo "   zcat ${OUTPUT_DIR}/Donor2019-010/SRR11272192_1.fastq.gz | head -n 4"
echo ""
echo "2. Merge runs for each donor:"
echo "   bash scripts/merge_10x_genomics_runs.sh"
echo ""
echo "3. Run Cell Ranger or STARsolo for alignment"
echo ""
echo "✓ Download complete!"
