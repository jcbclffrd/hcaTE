# hcaTE: Single-Cell RNA-seq Analysis of Transposable Elements

> **🤖 This pipeline was built by GitHub Copilot coding agent**

## Project Overview
This project analyzes **Smart-seq2 single-cell RNA-seq data** from human microglia in Alzheimer's Disease to quantify transposable element (TE) expression at the single-cell level.

**CRITICAL**: This is **barcoded single-cell data**, NOT bulk RNA-seq!

## Data Source
- **BioProject**: PRJNA611563
- **GEO Series**: GSE146639
- **SRA Accessions**: SRR11271993 - SRR11272311 (160 samples)
- **Title**: Single-cell RNA Sequencing of human microglia from post mortem Alzheimer's Disease CNS tissue
- **Library Type**: Smart-seq2 with barcoding (bc-Smart-seq2)
- **Sample Count**: 160 single-cell samples
- **Publication**: Olah et al. (2020) Nat Neurosci, PMID: 33192286

## Data Structure

### Read Layout (CRITICAL!)
Each sample has **two FASTQ files** with **different purposes**:

- **`SRR*_1.fastq.gz`**: **19bp cell barcodes** (16bp cell barcode + UMI)
  - Fixed length: 19bp
  - Purpose: Identify individual cells
  - Quality: Q30+ (excellent)
  
- **`SRR*_2.fastq.gz`**: **63bp cDNA sequences**
  - Fixed length: 63bp
  - Purpose: Map to genome/transcriptome
  - Quality: Q30-36 (excellent)

### Quality Control Summary
- **No trimming needed**: Excellent quality scores (Q30+) throughout
- **Minimal adapter content**: <2.2% contamination
- **Fixed-length reads**: Trimming would lose data
- **FastQC validated**: All 160 samples passed QC

## Directory Structure

```
hcaTE/
├── sra_downloads/          # 160 Smart-seq2 samples (FASTQ files)
│   ├── SRR11271993/
│   │   ├── SRR11271993_1.fastq.gz  (19bp barcodes)
│   │   └── SRR11271993_2.fastq.gz  (63bp cDNA)
│   ├── SRR11271994/
│   └── ...
├── aligned_bams/           # TODO: Create - output BAM files with cell barcodes
├── annotations/            # TODO: Copy - genome annotations
├── genome/                 # TODO: Create - reference genome
├── scripts/                # TODO: Create - analysis scripts
├── results/                # TODO: Create - scTE output, matrices
└── README.md               # This file
```

## Task: Set Up Cell Ranger / STARsolo for Alignment

### Objective
Create a pipeline to align the Smart-seq2 barcoded single-cell data using **Cell Ranger** or **STARsolo** (preferred for custom TE analysis).

### Requirements

#### 1. Choose Alignment Tool

**Option A: STARsolo (RECOMMENDED)**
- More flexible for custom TE analysis
- Can output per-cell BAMs tagged with cell barcodes
- Works well with scTE downstream
- Handles barcodes in separate file

**Option B: Cell Ranger**
- Standard industry tool
- May require barcode format conversion
- Produces standard outputs

#### 2. Handle Barcode Structure

The 19bp read structure needs to be determined:
- **Cell Barcode**: Likely 16bp (positions 1-16)
- **UMI**: Likely 10bp (positions 17-26 or similar)
- **Whitelist**: May need to extract valid cell barcodes from data

**Action Items:**
1. Analyze barcode structure (check first few samples)
2. Determine if a whitelist exists or needs to be generated
3. Configure barcode extraction parameters

#### 3. Genome Reference Setup

**Genome**: hg38 (GRCh38)
**Critical additions for TE analysis:**
- **Primary genome**: Standard hg38 from GENCODE
- **TE annotations**: RepeatMasker annotations (rmsk)
- **Gene annotations**: GENCODE v45 GTF

**Note**: Existing annotations may be available in `/home/jacobc/HumanBam2scTE/annotations/`:
- `gencode.v45.primary_assembly.annotation.gtf`
- `hg38_rmsk.bed`
- Consider copying these to `hcaTE/annotations/`

#### 4. Alignment Strategy

For **STARsolo**:
```bash
STAR --soloType SmartSeq \
     --readFilesIn SRR*_2.fastq.gz SRR*_1.fastq.gz \
     --soloUMIlen 10 \
     --soloCBlen 16 \
     --soloCBwhitelist [whitelist.txt or None] \
     --outSAMattributes NH HI AS nM CB UB \
     --outSAMtype BAM SortedByCoordinate \
     [additional parameters]
```

**Key outputs needed:**
- Cell-barcode-tagged BAM files (CB tag)
- Cell × gene count matrices
- Quality metrics per cell

#### 5. Expected Workflow

```
Step 1: Install STARsolo/Cell Ranger
Step 2: Set up genome index with TE annotations
Step 3: Determine barcode structure and create/identify whitelist
Step 4: Create alignment script for batch processing (160 samples)
Step 5: Run alignment on all samples
Step 6: Generate per-cell BAM files with CB tags
Step 7: QC: Check alignment rates, cells detected, UMI counts
```

#### 6. Downstream Integration

**Important**: The aligned BAMs will be used with **scTE** (single-cell TE quantification):
- scTE requires cell-barcode-tagged BAMs
- scTE mode: `scTE -i input.bam -o output -CB -UMI ...`
- Goal: Generate cell × TE expression matrices

#### 7. Validation Checklist

After alignment, verify:
- [ ] BAM files contain CB (cell barcode) tags
- [ ] BAM files contain UB (UMI) tags
- [ ] Multiple cells detected per sample
- [ ] Reasonable alignment rates (>50% for cDNA reads)
- [ ] Per-cell UMI counts are reasonable
- [ ] Output directory structure is organized

## Previous Errors (DO NOT REPEAT!)

❌ **Do NOT treat this as bulk RNA-seq**
❌ **Do NOT map barcode reads (R1) to genome**
❌ **Do NOT concatenate all cells together**
❌ **Do NOT use standard STAR without solo mode**

The previous pipeline (`HumanBam2scTE`) incorrectly treated this single-cell barcoded data as bulk RNA-seq, resulting in invalid differential expression results.

## Data Location

**Smart-seq2 FASTQ files**: `/home/jacobc/hcaTE/sra_downloads/`
- 160 samples
- Each with `_1.fastq.gz` (barcodes) and `_2.fastq.gz` (cDNA)

**Sample list with types**: `/home/jacobc/HumanBam2scTE/sample_types.txt`

**Original FastQC reports**: `/home/jacobc/HumanBam2scTE/fastqc_reports/`

## Questions for the Agent

1. **Barcode whitelist**: Does one exist in GEO metadata, or should we generate from data?
2. **Barcode structure**: Confirm the 16bp + UMI layout
3. **STARsolo vs Cell Ranger**: Which is preferred for TE analysis?
4. **Genome index**: Build new or use existing? Need TE integration?
5. **Batch processing**: Parallel processing strategy for 160 samples?
6. **Output format**: Cell-level BAMs or sample-level BAMs with CB tags?

## Differential Expression Analysis

### Get Sample Metadata
```bash
# Download GEO metadata for sample grouping (AD, CTR, CTR+, MCI)
wget -O GSE146639_metadata.txt "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146639&targ=self&form=text&view=full"
wget -O GSE146639_series_matrix.txt.gz "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146639/matrix/GSE146639_series_matrix.txt.gz"
gunzip GSE146639_series_matrix.txt.gz
```

### Aggregate Single-Cell to Pseudo-bulk
```bash
# Aggregate cells within each sample to create pseudo-bulk (treats samples as biological replicates)
python3 scripts/aggregate_to_pseudobulk.py
```

### Run DESeq2
```bash
# Compare AD vs Control using DESeq2
Rscript scripts/pseudobulk_diffexp_analysis.R
```

**Output**: Significant TEs in `pseudobulk/pseudobulk_diffexp_results_TEs_significant.csv`

## Contact & Context

This is part of a larger analysis to understand TE expression differences in Alzheimer's Disease microglia at single-cell resolution. The downstream tool (scTE) will quantify TEs per cell from the aligned BAMs.

---
**Date Created**: November 27, 2025
**Created By**: Analysis migration from HumanBam2scTE (bulk) → hcaTE (single-cell)
