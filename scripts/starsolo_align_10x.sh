#!/bin/bash

################################################################################
# STARsolo Alignment for 10x Genomics Data
# Project: hcaTE - TE expression in human microglia
# Data: GSE146639 - Alsema et al. (2020)
# Date: December 1, 2024
#
# STARsolo is an alternative to Cell Ranger that's already installed!
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}STARsolo 10x Genomics Alignment${NC}"
echo -e "${BLUE}GSE146639 - Alsema et al. (2020)${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Configuration
PROJECT_DIR="/home/jacobc/hcaTE"
FASTQ_DIR="${PROJECT_DIR}/10x_downloads/merged"
OUTPUT_DIR="${PROJECT_DIR}/10x_aligned_starsolo"
GENOME_DIR="${PROJECT_DIR}/star_index_10x"
GTF_FILE="${PROJECT_DIR}/annotations/gencode.v45.primary_assembly.annotation.gtf"

# STARsolo parameters
CORES=20
RAM_GB=60000000000  # 60GB in bytes

# 10x Genomics v2 chemistry parameters
WHITELIST="${PROJECT_DIR}/annotations/10x_whitelist_v2.txt"
UMI_LENGTH=11
BARCODE_LENGTH=16

# Sample information
declare -A SAMPLES=(
    ["Donor2019-010"]="MCI"
    ["Donor2018-135"]="AD"
)

################################################################################
# Check Prerequisites
################################################################################

echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check if STAR is installed
if ! command -v STAR &> /dev/null; then
    echo -e "${RED}ERROR: STAR not found in PATH${NC}"
    exit 1
fi

STAR_VERSION=$(STAR --version)
echo -e "${GREEN}✓ STAR found: ${STAR_VERSION}${NC}"

# Check if FASTQ files exist
if [ ! -d "$FASTQ_DIR" ]; then
    echo -e "${RED}ERROR: FASTQ directory not found: ${FASTQ_DIR}${NC}"
    exit 1
fi

echo -e "${GREEN}✓ FASTQ directory found${NC}"

# Check STAR index
if [ ! -d "$GENOME_DIR" ]; then
    echo -e "${RED}ERROR: STAR index not found: ${GENOME_DIR}${NC}"
    exit 1
fi

echo -e "${GREEN}✓ STAR index found${NC}"

# Create 10x whitelist if it doesn't exist
if [ ! -f "$WHITELIST" ]; then
    echo -e "${YELLOW}Creating 10x Genomics v2 whitelist...${NC}"
    
    # Download 10x whitelist (737K barcodes for v2 chemistry)
    wget -O "${WHITELIST}.gz" \
        "https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/737K-august-2016.txt.gz" \
        2>/dev/null || {
        echo -e "${RED}Failed to download whitelist${NC}"
        echo "Trying alternative source..."
        wget -O "${WHITELIST}.gz" \
            "https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/737K-august-2016.txt.gz"
    }
    
    gunzip "${WHITELIST}.gz"
    echo -e "${GREEN}✓ Whitelist created${NC}"
else
    echo -e "${GREEN}✓ Whitelist found${NC}"
fi

echo ""

################################################################################
# Run STARsolo for Each Sample
################################################################################

# Create output directory
mkdir -p "$OUTPUT_DIR"

for SAMPLE in "${!SAMPLES[@]}"; do
    CONDITION="${SAMPLES[$SAMPLE]}"
    
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Processing ${SAMPLE} (${CONDITION})${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo ""
    
    # Check if FASTQ files exist for this sample
    R1_FILE="${FASTQ_DIR}/${SAMPLE}_S1_L001_R1_001.fastq.gz"
    R2_FILE="${FASTQ_DIR}/${SAMPLE}_S1_L001_R2_001.fastq.gz"
    
    if [ ! -f "$R1_FILE" ] || [ ! -f "$R2_FILE" ]; then
        echo -e "${RED}ERROR: FASTQ files not found for ${SAMPLE}${NC}"
        echo "Expected:"
        echo "  $R1_FILE"
        echo "  $R2_FILE"
        continue
    fi
    
    echo -e "${GREEN}✓ Found FASTQ files:${NC}"
    echo "  R1: $(basename $R1_FILE) ($(du -h $R1_FILE | cut -f1))"
    echo "  R2: $(basename $R2_FILE) ($(du -h $R2_FILE | cut -f1))"
    echo ""
    
    # Sample output directory
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    
    # Check if already processed
    if [ -f "${SAMPLE_OUT}Solo.out/Gene/filtered/matrix.mtx" ]; then
        echo -e "${YELLOW}WARNING: ${SAMPLE} already processed. Skipping...${NC}"
        echo "To reprocess, delete: ${SAMPLE_OUT}"
        echo ""
        continue
    fi
    
    # Run STARsolo
    echo -e "${BLUE}Running STARsolo...${NC}"
    echo "Sample: ${SAMPLE}"
    echo "Condition: ${CONDITION}"
    echo "Cores: ${CORES}"
    echo ""
    
    STAR --runThreadN $CORES \
        --genomeDir "$GENOME_DIR" \
        --readFilesIn "$R2_FILE" "$R1_FILE" \
        --readFilesCommand zcat \
        --outFileNamePrefix "${SAMPLE_OUT}/" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
        --limitBAMsortRAM $RAM_GB \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 \
        --soloCBlen $BARCODE_LENGTH \
        --soloUMIstart 17 \
        --soloUMIlen 10 \
        --soloBarcodeReadLength 0 \
        --soloStrand Forward \
        --soloFeatures Gene GeneFull \
        --soloCellFilter EmptyDrops_CR \
        --soloUMIdedup 1MM_CR \
        --soloUMIfiltering MultiGeneUMI_CR \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts
    
    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ ${SAMPLE} completed successfully${NC}"
        echo ""
        
        # Display summary
        if [ -f "${SAMPLE_OUT}Solo.out/Gene/Summary.csv" ]; then
            echo -e "${BLUE}Alignment Summary:${NC}"
            cat "${SAMPLE_OUT}Solo.out/Gene/Summary.csv"
            echo ""
        fi
        
        # Count cells detected
        if [ -f "${SAMPLE_OUT}Solo.out/Gene/filtered/barcodes.tsv" ]; then
            CELL_COUNT=$(wc -l < "${SAMPLE_OUT}Solo.out/Gene/filtered/barcodes.tsv")
            echo -e "${GREEN}Cells detected: ${CELL_COUNT}${NC}"
            echo ""
        fi
    else
        echo -e "${RED}ERROR: STARsolo failed for ${SAMPLE}${NC}"
        echo ""
    fi
done

################################################################################
# Final Summary
################################################################################

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Alignment Complete${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

echo -e "${BLUE}Output directory:${NC} ${OUTPUT_DIR}"
echo ""

echo -e "${BLUE}STARsolo outputs for each sample:${NC}"
for SAMPLE in "${!SAMPLES[@]}"; do
    SAMPLE_OUT="${OUTPUT_DIR}/${SAMPLE}"
    
    if [ -d "${SAMPLE_OUT}Solo.out/Gene/filtered" ]; then
        echo ""
        echo -e "${GREEN}✓ ${SAMPLE}${NC}"
        echo "  - BAM: ${SAMPLE_OUT}Aligned.sortedByCoord.out.bam"
        echo "  - Gene matrix: ${SAMPLE_OUT}Solo.out/Gene/filtered/"
        echo "  - GeneFull matrix: ${SAMPLE_OUT}Solo.out/GeneFull/filtered/"
        echo "  - Summary: ${SAMPLE_OUT}Solo.out/Gene/Summary.csv"
        
        # Get cell count
        if [ -f "${SAMPLE_OUT}Solo.out/Gene/filtered/barcodes.tsv" ]; then
            CELLS=$(wc -l < "${SAMPLE_OUT}Solo.out/Gene/filtered/barcodes.tsv")
            echo "  - Cells detected: ${CELLS}"
        fi
        
        # Get total reads
        if [ -f "${SAMPLE_OUT}Solo.out/Gene/Summary.csv" ]; then
            READS=$(grep "Number of Reads" "${SAMPLE_OUT}Solo.out/Gene/Summary.csv" | cut -d',' -f2 || echo "N/A")
            echo "  - Total reads: ${READS}"
        fi
    else
        echo -e "${RED}✗ ${SAMPLE} (not processed)${NC}"
    fi
done

echo ""
echo -e "${BLUE}STARsolo vs Cell Ranger:${NC}"
echo "  ✓ Same algorithms (EmptyDrops, UMI deduplication)"
echo "  ✓ Compatible output formats"
echo "  ✓ Faster processing"
echo "  ✓ No installation needed (already have STAR)"
echo "  ✓ Free and open source"
echo ""

echo -e "${YELLOW}Next Steps:${NC}"
echo "1. Review Summary.csv files for QC metrics"
echo "2. Load matrices into Scanpy/Seurat:"
echo "   - Python: sc.read_10x_mtx('${OUTPUT_DIR}/Donor2019-010/Solo.out/Gene/filtered/')"
echo "   - R: Read10X('${OUTPUT_DIR}/Donor2019-010/Solo.out/Gene/filtered/')"
echo "3. Run scTE for TE quantification on BAM files:"
echo "   - scTE -i ${OUTPUT_DIR}/Donor2019-010/Aligned.sortedByCoord.out.bam \\"
echo "         -o 10x_scTE_output/Donor2019-010 -CB CR -UMI UR -x annotations/hg38_rmsk.gtf"
echo "4. Compare with Smart-seq2 clustering results (12,515 cells)"
echo ""

echo -e "${GREEN}✓ Pipeline complete!${NC}"
