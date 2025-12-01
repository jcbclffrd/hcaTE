#!/bin/bash

# Script: merge_10x_genomics_runs.sh
# Purpose: Merge multiple sequencing runs for each 10x Genomics donor
# Date: December 1, 2024
# Branch: feature/10x-genomics-analysis

set -e  # Exit on error

echo "========================================"
echo "10x Genomics Run Merger"
echo "GSE146639 - Alsema et al. (2020)"
echo "========================================"
echo ""

INPUT_DIR="10x_downloads/raw"
OUTPUT_DIR="10x_downloads/merged"

# Create output directory
mkdir -p ${OUTPUT_DIR}

echo "Input directory: ${INPUT_DIR}"
echo "Output directory: ${OUTPUT_DIR}"
echo ""

# Check if input files exist
if [ ! -d "${INPUT_DIR}/Donor2019-010" ] || [ ! -d "${INPUT_DIR}/Donor2018-135" ]; then
    echo "ERROR: Input directories not found!"
    echo "Please run download_10x_genomics_fastq.sh first"
    exit 1
fi

echo "========================================"
echo "Merging Donor 2019-010 (MCI)"
echo "GSM4403285 - 4 runs"
echo "========================================"
echo ""

cd ${INPUT_DIR}/Donor2019-010/

# Count files
R1_COUNT=$(ls SRR*_1.fastq.gz 2>/dev/null | wc -l)
R2_COUNT=$(ls SRR*_2.fastq.gz 2>/dev/null | wc -l)

echo "Found ${R1_COUNT} R1 files and ${R2_COUNT} R2 files"
echo ""

if [ ${R1_COUNT} -ne 4 ] || [ ${R2_COUNT} -ne 4 ]; then
    echo "WARNING: Expected 4 R1 and 4 R2 files, found ${R1_COUNT} R1 and ${R2_COUNT} R2"
    echo "Continuing with available files..."
fi

echo "Merging R1 (barcodes + UMI)..."
cat SRR11272192_1.fastq.gz SRR11272193_1.fastq.gz \
    SRR11272194_1.fastq.gz SRR11272195_1.fastq.gz \
    > ../../merged/Donor2019-010_S1_L001_R1_001.fastq.gz

echo "Merging R2 (cDNA)..."
cat SRR11272192_2.fastq.gz SRR11272193_2.fastq.gz \
    SRR11272194_2.fastq.gz SRR11272195_2.fastq.gz \
    > ../../merged/Donor2019-010_S1_L001_R2_001.fastq.gz

echo "✓ Donor 2019-010 merged successfully"
echo ""

cd ../../

echo "========================================"
echo "Merging Donor 2018-135 (AD)"
echo "GSM4403286 - 8 runs"
echo "========================================"
echo ""

cd ${INPUT_DIR}/Donor2018-135/

# Count files
R1_COUNT=$(ls SRR*_1.fastq.gz 2>/dev/null | wc -l)
R2_COUNT=$(ls SRR*_2.fastq.gz 2>/dev/null | wc -l)

echo "Found ${R1_COUNT} R1 files and ${R2_COUNT} R2 files"
echo ""

if [ ${R1_COUNT} -ne 8 ] || [ ${R2_COUNT} -ne 8 ]; then
    echo "WARNING: Expected 8 R1 and 8 R2 files, found ${R1_COUNT} R1 and ${R2_COUNT} R2"
    echo "Continuing with available files..."
fi

echo "Merging R1 (barcodes + UMI)..."
cat SRR11272196_1.fastq.gz SRR11272197_1.fastq.gz \
    SRR11272198_1.fastq.gz SRR11272199_1.fastq.gz \
    SRR11272200_1.fastq.gz SRR11272201_1.fastq.gz \
    SRR11272202_1.fastq.gz SRR11272203_1.fastq.gz \
    > ../../merged/Donor2018-135_S1_L001_R1_001.fastq.gz

echo "Merging R2 (cDNA)..."
cat SRR11272196_2.fastq.gz SRR11272197_2.fastq.gz \
    SRR11272198_2.fastq.gz SRR11272199_2.fastq.gz \
    SRR11272200_2.fastq.gz SRR11272201_2.fastq.gz \
    SRR11272202_2.fastq.gz SRR11272203_2.fastq.gz \
    > ../../merged/Donor2018-135_S1_L001_R2_001.fastq.gz

echo "✓ Donor 2018-135 merged successfully"
echo ""

cd ../../

echo "========================================"
echo "Merge Summary"
echo "========================================"
echo ""

echo "Merged files in ${OUTPUT_DIR}:"
ls -lh ${OUTPUT_DIR}/*.fastq.gz

echo ""
echo "Verifying merged files..."
echo ""

# Verify Donor 2019-010
echo "Donor 2019-010 (MCI):"
R1_READS=$(zcat ${OUTPUT_DIR}/Donor2019-010_S1_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}')
R2_READS=$(zcat ${OUTPUT_DIR}/Donor2019-010_S1_L001_R2_001.fastq.gz | wc -l | awk '{print $1/4}')
echo "  R1 reads: ${R1_READS}"
echo "  R2 reads: ${R2_READS}"

if [ "${R1_READS}" == "${R2_READS}" ]; then
    echo "  ✓ Read counts match!"
else
    echo "  ✗ WARNING: Read counts don't match!"
fi

echo ""

# Verify Donor 2018-135
echo "Donor 2018-135 (AD):"
R1_READS=$(zcat ${OUTPUT_DIR}/Donor2018-135_S1_L001_R1_001.fastq.gz | wc -l | awk '{print $1/4}')
R2_READS=$(zcat ${OUTPUT_DIR}/Donor2018-135_S1_L001_R2_001.fastq.gz | wc -l | awk '{print $1/4}')
echo "  R1 reads: ${R1_READS}"
echo "  R2 reads: ${R2_READS}"

if [ "${R1_READS}" == "${R2_READS}" ]; then
    echo "  ✓ Read counts match!"
else
    echo "  ✗ WARNING: Read counts don't match!"
fi

echo ""
echo "========================================"
echo "Next Steps"
echo "========================================"
echo ""
echo "1. Verify read structure:"
echo "   zcat ${OUTPUT_DIR}/Donor2019-010_S1_L001_R1_001.fastq.gz | head -n 8"
echo ""
echo "2. Run FastQC for quality control:"
echo "   fastqc ${OUTPUT_DIR}/*.fastq.gz"
echo ""
echo "3. Align with Cell Ranger or STARsolo:"
echo "   See INSTRUCTIONS_10x_genomics_download.md for details"
echo ""
echo "✓ Merge complete!"
