#!/bin/bash

################################################################################
# Cell Ranger Alignment for 10x Genomics Data
# Project: hcaTE - TE expression in human microglia
# Data: GSE146639 - Alsema et al. (2020)
# Date: December 1, 2024
################################################################################

set -euo pipefail

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}Cell Ranger 10x Genomics Alignment${NC}"
echo -e "${BLUE}GSE146639 - Alsema et al. (2020)${NC}"
echo -e "${BLUE}========================================${NC}"
echo ""

# Configuration
PROJECT_DIR="/home/jacobc/hcaTE"
FASTQ_DIR="${PROJECT_DIR}/10x_downloads/merged"
OUTPUT_DIR="${PROJECT_DIR}/10x_aligned"
TRANSCRIPTOME="${PROJECT_DIR}/cellranger_reference/hg38_gencode_v45_TEs"

# Cell Ranger parameters
CORES=20
MEM_GB=64
EXPECT_CELLS=3000  # Estimate based on 10x data

# Sample information
declare -A SAMPLES=(
    ["Donor2019-010"]="MCI"
    ["Donor2018-135"]="AD"
)

################################################################################
# Check Prerequisites
################################################################################

echo -e "${YELLOW}Checking prerequisites...${NC}"

# Check if Cell Ranger is installed
if ! command -v cellranger &> /dev/null; then
    echo -e "${RED}ERROR: Cell Ranger not found in PATH${NC}"
    echo ""
    echo -e "${YELLOW}Cell Ranger Installation Instructions:${NC}"
    echo ""
    echo "1. Download Cell Ranger from 10x Genomics:"
    echo "   https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest"
    echo ""
    echo "2. Extract and add to PATH:"
    echo "   tar -xzvf cellranger-8.0.1.tar.gz"
    echo "   export PATH=/path/to/cellranger-8.0.1:\$PATH"
    echo ""
    echo "3. Or use conda:"
    echo "   conda install -c bioconda cellranger"
    echo ""
    echo "4. Verify installation:"
    echo "   cellranger --version"
    echo ""
    exit 1
fi

CELLRANGER_VERSION=$(cellranger --version | grep -oP 'cellranger-\K[0-9.]+' || echo "unknown")
echo -e "${GREEN}✓ Cell Ranger found: v${CELLRANGER_VERSION}${NC}"

# Check if FASTQ files exist
if [ ! -d "$FASTQ_DIR" ]; then
    echo -e "${RED}ERROR: FASTQ directory not found: ${FASTQ_DIR}${NC}"
    echo "Please run download and merge scripts first"
    exit 1
fi

echo -e "${GREEN}✓ FASTQ directory found${NC}"

# List available FASTQ files
echo ""
echo -e "${BLUE}Available FASTQ files:${NC}"
ls -lh "$FASTQ_DIR"/*.fastq.gz 2>/dev/null || {
    echo -e "${RED}ERROR: No FASTQ files found in ${FASTQ_DIR}${NC}"
    exit 1
}

echo ""

################################################################################
# Build Reference Genome (if needed)
################################################################################

if [ ! -d "$TRANSCRIPTOME" ]; then
    echo -e "${YELLOW}Reference transcriptome not found. Building...${NC}"
    echo ""
    echo -e "${BLUE}Building Cell Ranger reference with GENCODE v45 + RepeatMasker TEs${NC}"
    
    # Create reference directory
    mkdir -p "${PROJECT_DIR}/cellranger_reference"
    
    # Check if genome and GTF exist
    GENOME_FASTA="${PROJECT_DIR}/genome/hg38.fa"
    GENES_GTF="${PROJECT_DIR}/annotations/gencode.v45.primary_assembly.annotation.gtf"
    
    if [ ! -f "$GENOME_FASTA" ]; then
        echo -e "${RED}ERROR: Genome FASTA not found: ${GENOME_FASTA}${NC}"
        exit 1
    fi
    
    if [ ! -f "$GENES_GTF" ]; then
        echo -e "${RED}ERROR: GTF annotation not found: ${GENES_GTF}${NC}"
        exit 1
    fi
    
    echo "Building Cell Ranger reference..."
    echo "This may take 1-2 hours..."
    
    cellranger mkref \
        --genome=hg38_gencode_v45_TEs \
        --fasta="$GENOME_FASTA" \
        --genes="$GENES_GTF" \
        --nthreads=$CORES \
        --memgb=$MEM_GB
    
    # Move to reference directory
    mv hg38_gencode_v45_TEs "${PROJECT_DIR}/cellranger_reference/"
    
    echo -e "${GREEN}✓ Reference genome built successfully${NC}"
else
    echo -e "${GREEN}✓ Reference transcriptome found: ${TRANSCRIPTOME}${NC}"
fi

echo ""

################################################################################
# Run Cell Ranger for Each Sample
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
    
    # Check if already processed
    if [ -d "${OUTPUT_DIR}/${SAMPLE}/outs" ]; then
        echo -e "${YELLOW}WARNING: ${SAMPLE} already processed. Skipping...${NC}"
        echo "To reprocess, delete: ${OUTPUT_DIR}/${SAMPLE}"
        echo ""
        continue
    fi
    
    # Run Cell Ranger count
    echo -e "${BLUE}Running Cell Ranger count...${NC}"
    echo "Sample: ${SAMPLE}"
    echo "Condition: ${CONDITION}"
    echo "Expected cells: ${EXPECT_CELLS}"
    echo "Cores: ${CORES}"
    echo "Memory: ${MEM_GB} GB"
    echo ""
    
    cd "$OUTPUT_DIR"
    
    cellranger count \
        --id="${SAMPLE}" \
        --transcriptome="$TRANSCRIPTOME" \
        --fastqs="$FASTQ_DIR" \
        --sample="${SAMPLE}" \
        --expect-cells=$EXPECT_CELLS \
        --localcores=$CORES \
        --localmem=$MEM_GB \
        --chemistry=auto
    
    if [ $? -eq 0 ]; then
        echo ""
        echo -e "${GREEN}✓ ${SAMPLE} completed successfully${NC}"
        echo ""
        
        # Display summary
        if [ -f "${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv" ]; then
            echo -e "${BLUE}Alignment Summary:${NC}"
            cat "${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv" | head -n 20
            echo ""
        fi
    else
        echo -e "${RED}ERROR: Cell Ranger failed for ${SAMPLE}${NC}"
        echo ""
    fi
    
    cd "$PROJECT_DIR"
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

echo -e "${BLUE}Cell Ranger outputs for each sample:${NC}"
for SAMPLE in "${!SAMPLES[@]}"; do
    if [ -d "${OUTPUT_DIR}/${SAMPLE}/outs" ]; then
        echo ""
        echo -e "${GREEN}✓ ${SAMPLE}${NC}"
        echo "  - BAM: ${OUTPUT_DIR}/${SAMPLE}/outs/possorted_genome_bam.bam"
        echo "  - Matrix: ${OUTPUT_DIR}/${SAMPLE}/outs/filtered_feature_bc_matrix/"
        echo "  - Metrics: ${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv"
        echo "  - Web Summary: ${OUTPUT_DIR}/${SAMPLE}/outs/web_summary.html"
        
        # Get cell count
        if [ -f "${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv" ]; then
            CELLS=$(grep "Estimated Number of Cells" "${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' || echo "N/A")
            READS=$(grep "Number of Reads" "${OUTPUT_DIR}/${SAMPLE}/outs/metrics_summary.csv" | cut -d',' -f2 | tr -d '"' || echo "N/A")
            echo "  - Cells detected: ${CELLS}"
            echo "  - Total reads: ${READS}"
        fi
    else
        echo -e "${RED}✗ ${SAMPLE} (not processed)${NC}"
    fi
done

echo ""
echo -e "${YELLOW}Next Steps:${NC}"
echo "1. Review web summaries: ${OUTPUT_DIR}/*/outs/web_summary.html"
echo "2. Run scTE for TE quantification on BAM files"
echo "3. Load count matrices into Scanpy/Seurat for analysis"
echo "4. Compare with Smart-seq2 clustering results (12,515 cells)"
echo ""

echo -e "${GREEN}✓ Pipeline complete!${NC}"
