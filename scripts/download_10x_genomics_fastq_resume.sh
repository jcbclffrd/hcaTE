#!/bin/bash

# Script: download_10x_genomics_fastq_resume.sh
# Purpose: Download 10x Genomics FASTQ files with resume capability
# Date: December 1, 2024
# Branch: feature/10x-genomics-analysis

set -e  # Exit on error

echo "========================================"
echo "10x Genomics FASTQ Download (Resume-capable)"
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

# Function to download a single SRR with retry
download_srr() {
    local SRR=$1
    local DONOR_DIR=$2
    local MAX_RETRIES=3
    local RETRY=0
    
    # Check if already downloaded
    if [ -f "${DONOR_DIR}/${SRR}_1.fastq.gz" ] && [ -f "${DONOR_DIR}/${SRR}_2.fastq.gz" ]; then
        echo "✓ ${SRR} already downloaded, skipping"
        return 0
    fi
    
    echo "----------------------------------------"
    echo "Downloading: $SRR"
    echo "----------------------------------------"
    
    while [ $RETRY -lt $MAX_RETRIES ]; do
        if fasterq-dump $SRR \
            --outdir ${DONOR_DIR}/ \
            --split-files \
            --threads ${THREADS} \
            --progress; then
            
            echo "Compressing $SRR..."
            gzip ${DONOR_DIR}/${SRR}*.fastq
            
            echo "✓ Completed $SRR"
            echo ""
            return 0
        else
            RETRY=$((RETRY + 1))
            echo "⚠ Download failed (attempt $RETRY/$MAX_RETRIES)"
            if [ $RETRY -lt $MAX_RETRIES ]; then
                WAIT_TIME=$((RETRY * 30))
                echo "Waiting ${WAIT_TIME}s before retry..."
                sleep $WAIT_TIME
            fi
        fi
    done
    
    echo "✗ Failed to download $SRR after $MAX_RETRIES attempts"
    echo "Continuing with remaining files..."
    echo ""
    return 1
}

echo "========================================"
echo "Downloading Donor 2019-010 (MCI)"
echo "GSM4403285 - 4 runs"
echo "========================================"
echo ""

# Donor 2019-010 (MCI) - 4 runs
DONOR1_RUNS=("SRR11272192" "SRR11272193" "SRR11272194" "SRR11272195")
DONOR1_DIR="${OUTPUT_DIR}/Donor2019-010"

for SRR in "${DONOR1_RUNS[@]}"; do
    download_srr $SRR $DONOR1_DIR
done

echo "========================================"
echo "Downloading Donor 2018-135 (AD)"
echo "GSM4403286 - 8 runs"
echo "========================================"
echo ""

# Donor 2018-135 (AD) - 8 runs
DONOR2_RUNS=("SRR11272196" "SRR11272197" "SRR11272198" "SRR11272199" 
             "SRR11272200" "SRR11272201" "SRR11272202" "SRR11272203")
DONOR2_DIR="${OUTPUT_DIR}/Donor2018-135"

for SRR in "${DONOR2_RUNS[@]}"; do
    download_srr $SRR $DONOR2_DIR
done

echo "========================================"
echo "Download Summary"
echo "========================================"
echo ""
echo "Donor 2019-010 (MCI) files:"
ls -lh ${DONOR1_DIR}/*.fastq.gz 2>/dev/null | wc -l | xargs echo "Files downloaded:"
du -sh ${DONOR1_DIR}
echo ""
echo "Donor 2018-135 (AD) files:"
ls -lh ${DONOR2_DIR}/*.fastq.gz 2>/dev/null | wc -l | xargs echo "Files downloaded:"
du -sh ${DONOR2_DIR}
echo ""

echo "========================================"
echo "Next Steps"
echo "========================================"
echo ""
echo "1. Verify file integrity:"
echo "   zcat ${DONOR1_DIR}/SRR11272192_1.fastq.gz | head -n 4"
echo ""
echo "2. Merge runs for each donor:"
echo "   bash scripts/merge_10x_genomics_runs.sh"
echo ""
echo "3. Run Cell Ranger or STARsolo for alignment"
echo ""
echo "✓ Download script complete!"
