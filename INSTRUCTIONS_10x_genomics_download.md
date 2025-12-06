# Instructions: Downloading 10x Genomics Data from NCBI

**Date Created**: December 1, 2024  
**Purpose**: Download 10x Genomics single-cell RNA-seq data for TE analysis  
**Branch**: feature/10x-genomics-analysis

---

## Overview

This document provides instructions for downloading 10x Genomics single-cell RNA-seq data from NCBI SRA for analysis of transposable element (TE) expression in human microglia.

**Key Differences from Smart-seq2:**
- 10x Genomics uses droplet-based technology (thousands of cells per sample)
- Different barcode structure (16bp cell barcode + 10-12bp UMI)
- Uses Cell Ranger or STARsolo for alignment
- Lower reads per cell but higher cell numbers

---

## Quick Start: Download Same Study 10x Genomics Data

**For the Alsema et al. (2020) 10x Genomics samples from GSE146639:**

### Automated Scripts (Recommended) 🚀

```bash
# Step 1: Download all FASTQ files (12 runs total: 4 + 8)
bash scripts/download_10x_genomics_fastq.sh

# Step 2: Merge runs for each donor
bash scripts/merge_10x_genomics_runs.sh

# Step 3: Verify merged files
zcat 10x_downloads/merged/Donor2019-010_S1_L001_R1_001.fastq.gz | head -n 8
```

**Expected Runtime**: 1-3 hours depending on network speed  
**Expected Size**: ~20-30 GB compressed

---

### Manual Download Instructions

### Option A: Download Pre-Processed Matrices (Fastest) ⚡

```bash
# Create download directory
mkdir -p 10x_downloads/

# Download Donor 2019-010 (MCI)
cd 10x_downloads/
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403285/suppl/GSM4403285_10X_Donor2019-010_barcodes.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403285/suppl/GSM4403285_10X_Donor2019-010_genes.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403285/suppl/GSM4403285_10X_Donor2019-010_matrix.mtx.gz"

# Download Donor 2018-135 (AD)
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403286/suppl/GSM4403286_10X_Donor2018-135_barcodes.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403286/suppl/GSM4403286_10X_Donor2018-135_genes.tsv.gz"
wget "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4403nnn/GSM4403286/suppl/GSM4403286_10X_Donor2018-135_matrix.mtx.gz"

# Decompress
gunzip *.gz

cd ..
```

**Matrix Format**: Standard 10x Genomics Cell Ranger output format
- `barcodes.tsv`: Cell barcodes (one per cell)
- `genes.tsv`: Gene IDs and names
- `matrix.mtx`: Sparse count matrix (Market Matrix format)

**Note**: These matrices contain **genes only**, NOT TEs. You'll need to download raw FASTQ files to quantify TEs (see Option B below).

### Option B: Download Raw FASTQ Files for TE Quantification (Required for TE Analysis) ⭐

The processed matrices above don't include TE counts. To analyze TEs, you need raw FASTQ files:

**SRR Accessions Found:**

**GSM4403285 (Donor 2019-010, MCI):**
- SRR11272192
- SRR11272193
- SRR11272194
- SRR11272195

**GSM4403286 (Donor 2018-135, AD):**
- SRR11272196
- SRR11272197
- SRR11272198
- SRR11272199
- SRR11272200
- SRR11272201
- SRR11272202
- SRR11272203

**Download Script:**

```bash
# Create download directory
mkdir -p 10x_downloads/raw/

# Download Donor 2019-010 (MCI) - 4 runs
for SRR in SRR11272192 SRR11272193 SRR11272194 SRR11272195; do
    echo "Downloading $SRR..."
    fasterq-dump $SRR \
      --outdir 10x_downloads/raw/Donor2019-010/ \
      --split-files \
      --threads 20 \
      --progress
    
    # Compress immediately to save space
    pigz --processes 20 10x_downloads/raw/Donor2019-010/${SRR}*.fastq
    
    echo "Completed $SRR"
done

# Download Donor 2018-135 (AD) - 8 runs
for SRR in SRR11272196 SRR11272197 SRR11272198 SRR11272199 SRR11272200 SRR11272201 SRR11272202 SRR11272203; do
    echo "Downloading $SRR..."
    fasterq-dump $SRR \
      --outdir 10x_downloads/raw/Donor2018-135/ \
      --split-files \
      --threads 20 \
      --progress
    
    # Compress immediately to save space
    pigz --processes 20 10x_downloads/raw/Donor2018-135/${SRR}*.fastq
    
    echo "Completed $SRR"
done
```

**Note**: Each donor was sequenced across **multiple runs** (4 runs for MCI donor, 8 runs for AD donor). These will need to be **merged** after download before alignment.

**Expected Files for 10x Genomics v2/v3:**
- `SRR_ID_1.fastq.gz`: R1 - Cell barcode (16bp) + UMI (10-12bp)
- `SRR_ID_2.fastq.gz`: R2 - cDNA sequence (90-98bp)
- `SRR_ID_3.fastq.gz`: I1 - Sample index (optional)

### Merging Multiple Runs

After downloading all runs for each donor, merge them by read type:

```bash
# Merge Donor 2019-010 (MCI)
cd 10x_downloads/raw/Donor2019-010/

# Merge R1 (barcodes + UMI)
cat SRR11272192_1.fastq.gz SRR11272193_1.fastq.gz \
    SRR11272194_1.fastq.gz SRR11272195_1.fastq.gz \
    > Donor2019-010_S1_L001_R1_001.fastq.gz

# Merge R2 (cDNA)
cat SRR11272192_2.fastq.gz SRR11272193_2.fastq.gz \
    SRR11272194_2.fastq.gz SRR11272195_2.fastq.gz \
    > Donor2019-010_S1_L001_R2_001.fastq.gz

cd ../../..

# Merge Donor 2018-135 (AD)
cd 10x_downloads/raw/Donor2018-135/

# Merge R1 (barcodes + UMI)
cat SRR11272196_1.fastq.gz SRR11272197_1.fastq.gz \
    SRR11272198_1.fastq.gz SRR11272199_1.fastq.gz \
    SRR11272200_1.fastq.gz SRR11272201_1.fastq.gz \
    SRR11272202_1.fastq.gz SRR11272203_1.fastq.gz \
    > Donor2018-135_S1_L001_R1_001.fastq.gz

# Merge R2 (cDNA)
cat SRR11272196_2.fastq.gz SRR11272197_2.fastq.gz \
    SRR11272198_2.fastq.gz SRR11272199_2.fastq.gz \
    SRR11272200_2.fastq.gz SRR11272201_2.fastq.gz \
    SRR11272202_2.fastq.gz SRR11272203_2.fastq.gz \
    > Donor2018-135_S1_L001_R2_001.fastq.gz

cd ../../..
```

**Merged File Naming**: Following 10x Genomics Cell Ranger format:
- `{SampleName}_S{SampleNumber}_L{Lane}_R{Read}_001.fastq.gz`
- Example: `Donor2019-010_S1_L001_R1_001.fastq.gz`

---

## Finding 10x Genomics Datasets

### Option 1: Same Study (Alsema et al. 2020) ⭐ **RECOMMENDED**

**Great news!** The same study (GSE146639) has **2 donors with 10x Genomics data**:

### 10x Genomics Samples Available:

1. **GSM4403285** - 10X_Donor2019-010
   - **Donor Group**: MCI (Mild Cognitive Impairment)
   - **Age**: 76, Female
   - **Brain Region**: LPS (superior parietal lobe)
   - **SRA**: SRX7878785
   - **Processed Data**: ~15,000-20,000 cells (estimated from matrix size)

2. **GSM4403286** - 10X_Donor2018-135
   - **Donor Group**: AD (Alzheimer's Disease)
   - **Age**: 81, Female
   - **Brain Region**: LPS (superior parietal lobe)
   - **SRA**: SRX7878786
   - **Processed Data**: ~10,000-15,000 cells (estimated from matrix size)

**Key Advantages of Using This Data:**
- ✅ Same study, same tissue source, same isolation protocol
- ✅ Direct comparison with Smart-seq2 results (12,515 cells analyzed)
- ✅ Same brain bank (NetherlandsBrainBank), same brain region (LPS)
- ✅ Pre-processed matrices available (already Cell Ranger output)
- ✅ Raw FASTQ files available in SRA for custom TE analysis

### Option 2: Alternative Human Microglia 10x Datasets

Search NCBI SRA for human microglia 10x Genomics datasets:

**Recommended Search Terms:**
- "human microglia 10x genomics"
- "human brain single cell 10x"
- "Alzheimer disease microglia 10x"
- "human microglia scRNA-seq"

**Example High-Quality Datasets:**

1. **Mathys et al. (2019)** - Single-cell RNA-seq of Alzheimer's brain
   - GEO: GSE138852
   - Technology: 10x Genomics 3' v2
   - ~80,000 cells from AD and control brains
   - Includes microglia subset

2. **Olah et al. (2020)** - Human microglia heterogeneity
   - GEO: GSE140511  
   - Technology: 10x Genomics
   - Purified human microglia from multiple brain regions

3. **Sankowski et al. (2019)** - Human microglia from surgical samples
   - GEO: GSE130119
   - Technology: 10x Genomics
   - Fresh human microglia

---

## Example: Downloading from SRA (Generic Workflow)

Once you've identified a 10x Genomics dataset, follow these steps:

### Step 1: Get SRA Accession Numbers

```bash
# Navigate to GEO dataset page
# Example: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138852

# Download SRA Run Table
wget -O SraRunTable.txt "https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA<PROJECT_ID>&o=acc_s%3Aa"

# Or use SRA Run Selector:
# Visit: https://www.ncbi.nlm.nih.gov/Traces/study/
# Search by: GEO accession or BioProject
# Download: "Metadata" and "Accession List"
```

### Step 2: Install SRA Toolkit

```bash
# If not already installed
conda install -c bioconda sra-tools

# Or download from NCBI:
# https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit

# Configure SRA toolkit
vdb-config --interactive
# Set cache location and download preferences
```

### Step 3: Download 10x Genomics FASTQ Files

**IMPORTANT**: 10x Genomics data has different file structure than Smart-seq2!

10x Genomics produces **3 FASTQ files per sample**:
- **R1**: Cell barcode + UMI (26-28bp)
  - 16bp cell barcode
  - 10-12bp UMI
- **R2**: cDNA sequence (90-98bp)
- **I1**: Sample index (optional, for demultiplexing)

```bash
# Create download directory
mkdir -p 10x_downloads/

# Download using fasterq-dump (faster than fastq-dump)
# Replace SRR_ACCESSION with actual accession numbers

# Download single sample
fasterq-dump SRR_ACCESSION \
  --outdir 10x_downloads/ \
  --split-files \
  --threads 20 \
  --progress

# This creates:
# SRR_ACCESSION_1.fastq (R1: barcodes + UMI)
# SRR_ACCESSION_2.fastq (R2: cDNA reads)
# SRR_ACCESSION_3.fastq (I1: sample index, if present)

# Compress files to save space
pigz --processes 20 10x_downloads/SRR_ACCESSION_*.fastq
```

### Step 4: Batch Download Multiple Samples

```bash
# Create list of SRA accessions
cat > sra_accessions_10x.txt << 'EOF'
SRR12345678
SRR12345679
SRR12345680
EOF

# Download all samples
while read SRR; do
    echo "Downloading $SRR..."
    fasterq-dump $SRR \
      --outdir 10x_downloads/${SRR}/ \
      --split-files \
      --threads 20 \
      --progress
    
    # Compress immediately
    pigz --processes 20 10x_downloads/${SRR}/*.fastq
    
    echo "Completed $SRR"
done < sra_accessions_10x.txt
```

### Step 5: Verify File Structure

```bash
# Check first sample
ls -lh 10x_downloads/SRR12345678/

# Expected output:
# SRR12345678_1.fastq.gz  (R1: barcodes + UMI)
# SRR12345678_2.fastq.gz  (R2: cDNA)
# SRR12345678_3.fastq.gz  (I1: index, optional)

# Verify read counts match
zcat 10x_downloads/SRR12345678/SRR12345678_1.fastq.gz | wc -l
zcat 10x_downloads/SRR12345678/SRR12345678_2.fastq.gz | wc -l
# Should be equal (divide by 4 for number of reads)

# Check read structure
zcat 10x_downloads/SRR12345678/SRR12345678_1.fastq.gz | head -n 2
# R1 should be ~26-28bp (16bp barcode + 10-12bp UMI)

zcat 10x_downloads/SRR12345678/SRR12345678_2.fastq.gz | head -n 2
# R2 should be ~90-98bp (cDNA sequence)
```

---

## Dataset Selection Criteria

When choosing a 10x Genomics dataset, consider:

1. **Cell Type**: Microglia-specific or contains microglia subset
2. **Sample Size**: At least 3-5 biological replicates per condition
3. **Sequencing Depth**: >20,000 reads per cell recommended
4. **Chemistry Version**: 10x v2 or v3 (v3 preferred for sensitivity)
5. **Data Availability**: Raw FASTQ files available in SRA
6. **Metadata**: Sample annotations (disease state, age, sex, brain region)
7. **Publication**: Peer-reviewed study with clear methods

---

## Recommended Dataset: Mathys et al. (2019)

**Why this dataset is excellent for TE analysis:**

- **BioProject**: PRJNA531535
- **GEO**: GSE138852
- **Publication**: Nature (2019), PMID: 31042697
- **Technology**: 10x Genomics Chromium Single Cell 3' v2
- **Sample Size**: 48 individuals (24 AD, 24 controls)
- **Cells**: ~80,000 single-nucleus RNA-seq profiles
- **Cell Types**: Includes microglia, astrocytes, neurons, oligodendrocytes
- **Brain Region**: Prefrontal cortex
- **Coverage**: ~1.5-2.5k UMIs per cell
- **Data Quality**: High-quality, well-annotated

**Microglia Subset:**
- ~5,000-10,000 microglia cells (subset by CX3CR1, P2RY12, TMEM119 expression)
- Can directly compare with Smart-seq2 results from GSE146639

### Download Command for Mathys et al. Dataset

```bash
# Create download directory
mkdir -p 10x_downloads/mathys_2019/

# Example SRA accessions from Mathys et al. (verify from SRA Run Selector)
# Visit: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA531535
# Download run table to get all SRR accessions

# Download single sample (example)
fasterq-dump SRR8718201 \
  --outdir 10x_downloads/mathys_2019/SRR8718201/ \
  --split-files \
  --threads 20 \
  --progress

pigz --processes 20 10x_downloads/mathys_2019/SRR8718201/*.fastq
```

---

## Expected File Sizes

**10x Genomics files are typically larger than Smart-seq2:**

- **Per sample**: 5-20 GB compressed (depending on cell number and depth)
- **R1 file**: ~10-20% of total size (shorter reads)
- **R2 file**: ~80-90% of total size (longer cDNA reads)
- **Total dataset**: 100-500 GB for ~48 samples

**Storage recommendations:**
- Allocate at least 500 GB for raw FASTQ files
- Additional 500 GB for aligned BAM files
- Additional 100 GB for analysis outputs

---

## Next Steps After Download

Once downloaded, proceed with:

1. **FastQC**: Quality control of FASTQ files
2. **Cell Ranger or STARsolo**: Alignment and quantification
3. **scTE**: TE quantification at single-cell level
4. **Scanpy/Seurat**: Cell clustering and analysis
5. **Compare with Smart-seq2 results** from GSE146639

See separate analysis pipeline documentation for 10x Genomics alignment and TE quantification.

---

## Troubleshooting

### Issue: Download is very slow
```bash
# Use parallel download with prefetch
prefetch --max-size 100G --progress SRR_ACCESSION

# Then convert to FASTQ
fasterq-dump SRR_ACCESSION --split-files
```

### Issue: Running out of disk space
```bash
# Download and process one sample at a time
# Delete intermediate files immediately after processing

# Or use streaming mode (no intermediate storage)
prefetch SRR_ACCESSION | fasterq-dump - --split-files
```

### Issue: Files are corrupted
```bash
# Validate downloaded files
vdb-validate SRR_ACCESSION

# Re-download if validation fails
rm -rf ~/ncbi/public/sra/SRR_ACCESSION.sra
fasterq-dump SRR_ACCESSION --split-files
```

---

## References

- 10x Genomics Support: https://support.10xgenomics.com/
- SRA Toolkit: https://github.com/ncbi/sra-tools
- Mathys et al. (2019): https://www.nature.com/articles/s41586-019-1195-2
- Cell Ranger Documentation: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger

---

**Questions or Issues?**
- Check SRA Run Selector for accurate accession numbers
- Verify file structure matches 10x Genomics specifications
- Contact: GitHub issues on jcbclffrd/hcaTE repository
