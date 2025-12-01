# Cell Ranger Setup & Alignment Guide

**Project**: hcaTE - 10x Genomics Analysis  
**Branch**: feature/10x-genomics-analysis  
**Date**: December 1, 2024

---

## Overview

This guide covers setting up Cell Ranger and running alignment for the 10x Genomics data from GSE146639.

**Data Status**: ✅ Downloaded and merged (11 GB, 2 donors)

---

## Step 1: Install Cell Ranger

### Option A: Download from 10x Genomics (Recommended)

```bash
# 1. Download Cell Ranger from 10x Genomics website
# Visit: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
# Sign up for free account if needed
# Download latest version (currently 8.0.1)

# 2. Extract
cd /home/jacobc/software/  # or your preferred location
tar -xzvf cellranger-8.0.1.tar.gz

# 3. Add to PATH
export PATH=/home/jacobc/software/cellranger-8.0.1:$PATH

# 4. Add to .bashrc for persistence
echo 'export PATH=/home/jacobc/software/cellranger-8.0.1:$PATH' >> ~/.bashrc
source ~/.bashrc

# 5. Verify installation
cellranger --version
# Should show: cellranger cellranger-8.0.1
```

### Option B: Conda Installation (Alternative)

```bash
# Create conda environment
conda create -n cellranger -c bioconda cellranger

# Activate environment
conda activate cellranger

# Verify
cellranger --version
```

**Note**: The conda version may be older than the official release.

---

## Step 2: Build Reference Genome

Cell Ranger needs a custom reference that includes both genes and TEs.

### Quick Start (Automated)

The alignment script (`scripts/cellranger_align_10x.sh`) will automatically build the reference if it doesn't exist.

### Manual Build (Optional)

```bash
cd /home/jacobc/hcaTE

# Build Cell Ranger reference with GENCODE v45 + TEs
cellranger mkref \
    --genome=hg38_gencode_v45_TEs \
    --fasta=genome/hg38.fa \
    --genes=annotations/gencode.v45.primary_assembly.annotation.gtf \
    --nthreads=20 \
    --memgb=64

# Move to reference directory
mkdir -p cellranger_reference
mv hg38_gencode_v45_TEs cellranger_reference/
```

**Expected Runtime**: 1-2 hours  
**Output Size**: ~30-40 GB

---

## Step 3: Run Cell Ranger Alignment

### Automated Pipeline (Recommended)

```bash
cd /home/jacobc/hcaTE

# Run alignment for both donors
bash scripts/cellranger_align_10x.sh
```

The script will:
1. Check for Cell Ranger installation
2. Verify merged FASTQ files exist
3. Build reference genome (if needed)
4. Run `cellranger count` for each donor
5. Generate summary statistics

### Manual Alignment (Alternative)

```bash
cd /home/jacobc/hcaTE/10x_aligned

# Donor 2019-010 (MCI)
cellranger count \
    --id=Donor2019-010 \
    --transcriptome=/home/jacobc/hcaTE/cellranger_reference/hg38_gencode_v45_TEs \
    --fastqs=/home/jacobc/hcaTE/10x_downloads/merged \
    --sample=Donor2019-010 \
    --expect-cells=3000 \
    --localcores=20 \
    --localmem=64 \
    --chemistry=auto

# Donor 2018-135 (AD)
cellranger count \
    --id=Donor2018-135 \
    --transcriptome=/home/jacobc/hcaTE/cellranger_reference/hg38_gencode_v45_TEs \
    --fastqs=/home/jacobc/hcaTE/10x_downloads/merged \
    --sample=Donor2018-135 \
    --expect-cells=3000 \
    --localcores=20 \
    --localmem=64 \
    --chemistry=auto
```

---

## Step 4: Review Results

### Cell Ranger Outputs

For each sample, Cell Ranger produces:

```
10x_aligned/Donor2019-010/outs/
├── web_summary.html              # ⭐ Interactive QC report
├── metrics_summary.csv           # Key metrics
├── possorted_genome_bam.bam      # Aligned reads (for scTE)
├── possorted_genome_bam.bam.bai  # BAM index
├── filtered_feature_bc_matrix/   # Gene expression matrix
│   ├── barcodes.tsv.gz           # Cell barcodes
│   ├── features.tsv.gz           # Gene names
│   └── matrix.mtx.gz             # Count matrix
├── raw_feature_bc_matrix/        # Unfiltered matrix (all barcodes)
├── molecule_info.h5              # Per-molecule information
└── cloupe.cloupe                 # Loupe Browser file
```

### View Web Summary

```bash
# Open in browser
firefox /home/jacobc/hcaTE/10x_aligned/Donor2019-010/outs/web_summary.html
firefox /home/jacobc/hcaTE/10x_aligned/Donor2018-135/outs/web_summary.html
```

The web summary shows:
- Number of cells detected
- Median genes per cell
- Sequencing saturation
- Reads per cell
- Mapping statistics
- t-SNE/UMAP plots

### Key Metrics to Check

```bash
# View metrics for both samples
cat 10x_aligned/Donor2019-010/outs/metrics_summary.csv
cat 10x_aligned/Donor2018-135/outs/metrics_summary.csv
```

**Expected metrics for good quality data:**
- Estimated cells: 2,000-5,000 per donor
- Median genes per cell: >1,000
- Sequencing saturation: >60%
- Reads mapped to genome: >80%
- Reads mapped to transcriptome: >50%

---

## Step 5: Next Steps - TE Quantification

After Cell Ranger alignment, we need to quantify TEs using the BAM files.

### Option A: scTE (Recommended for single-cell)

```bash
# Run scTE on Cell Ranger BAM files
scTE -i 10x_aligned/Donor2019-010/outs/possorted_genome_bam.bam \
     -o 10x_scTE_output/Donor2019-010 \
     -x annotations/hg38_rmsk.gtf \
     -CB CR \
     -UMI UR \
     -p 20

scTE -i 10x_aligned/Donor2018-135/outs/possorted_genome_bam.bam \
     -o 10x_scTE_output/Donor2018-135 \
     -x annotations/hg38_rmsk.gtf \
     -CB CR \
     -UMI UR \
     -p 20
```

### Option B: TEtranscripts (For bulk/pseudobulk)

If aggregating to pseudobulk:

```bash
TEcount \
    --BAM 10x_aligned/Donor2019-010/outs/possorted_genome_bam.bam \
    --GTF annotations/gencode.v45.primary_assembly.annotation.gtf \
    --TE annotations/hg38_rmsk.gtf \
    --format BAM \
    --mode multi
```

---

## Step 6: Load into Scanpy/Seurat

### Python (Scanpy)

```python
import scanpy as sc

# Load Cell Ranger output
adata_mci = sc.read_10x_mtx('10x_aligned/Donor2019-010/outs/filtered_feature_bc_matrix/')
adata_ad = sc.read_10x_mtx('10x_aligned/Donor2018-135/outs/filtered_feature_bc_matrix/')

# Add sample metadata
adata_mci.obs['condition'] = 'MCI'
adata_mci.obs['donor'] = 'Donor2019-010'

adata_ad.obs['condition'] = 'AD'
adata_ad.obs['donor'] = 'Donor2018-135'

# Merge
adata = adata_mci.concatenate(adata_ad)

# Load TE counts from scTE
# (Load scTE output CSV and merge with adata)
```

### R (Seurat)

```R
library(Seurat)

# Load Cell Ranger output
mci <- Read10X(data.dir = "10x_aligned/Donor2019-010/outs/filtered_feature_bc_matrix/")
ad <- Read10X(data.dir = "10x_aligned/Donor2018-135/outs/filtered_feature_bc_matrix/")

# Create Seurat objects
mci_seurat <- CreateSeuratObject(counts = mci, project = "MCI")
ad_seurat <- CreateSeuratObject(counts = ad, project = "AD")

# Merge
combined <- merge(mci_seurat, ad_seurat)
```

---

## Troubleshooting

### Issue: Cell Ranger not found

```bash
# Check PATH
echo $PATH | grep cellranger

# Add to PATH if missing
export PATH=/path/to/cellranger-8.0.1:$PATH

# Or use full path in script
/path/to/cellranger-8.0.1/cellranger count ...
```

### Issue: Out of memory

```bash
# Reduce memory allocation
cellranger count \
    --localmem=32 \
    ...

# Or run one sample at a time
```

### Issue: Low cell detection

If Cell Ranger detects fewer cells than expected:

```bash
# Force specific cell number
cellranger count \
    --force-cells=5000 \
    ...

# Or adjust expected cells
cellranger count \
    --expect-cells=5000 \
    ...
```

### Issue: Chemistry detection failed

```bash
# Manually specify chemistry
cellranger count \
    --chemistry=SC3Pv2 \
    # or SC3Pv3, SC5P-PE, SC5P-R2, etc.
    ...
```

---

## Expected Runtime & Resources

**Per sample:**
- **Runtime**: 2-4 hours (depends on cell number and read depth)
- **CPU**: 20 cores (can adjust with `--localcores`)
- **Memory**: 64 GB (can adjust with `--localmem`)
- **Storage**: ~50 GB per sample (BAM + matrices)

**Total for 2 donors:**
- **Runtime**: 4-8 hours
- **Storage**: ~100 GB

---

## Comparison with Smart-seq2

After completing alignment, you can compare:

| Feature | Smart-seq2 (bc-Smart-seq2) | 10x Genomics |
|---------|----------------------------|---------------|
| Samples | 160 pooled samples | 2 donors |
| Cells | ~12,515 cells | ~4,000-10,000 cells (estimated) |
| Reads/cell | ~500k-1M reads | ~20k-50k reads |
| Genes detected | ~4,000-6,000 genes/cell | ~1,000-2,000 genes/cell |
| TE detection | High (full-length) | Moderate (3' bias) |
| Cost | High (more reads) | Low (less reads) |
| Throughput | Low (fewer cells) | High (more cells) |

---

## References

- Cell Ranger Documentation: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger
- scTE Documentation: https://github.com/jphe/scTE
- Scanpy Tutorial: https://scanpy-tutorials.readthedocs.io/
- Seurat Tutorial: https://satijalab.org/seurat/

---

## Summary Checklist

- [ ] Cell Ranger installed and in PATH
- [ ] Reference genome built with TEs
- [ ] Both donors aligned successfully
- [ ] Web summaries reviewed
- [ ] Cell counts look reasonable
- [ ] BAM files ready for scTE
- [ ] Count matrices loaded into analysis framework

**Next**: Run scTE for TE quantification and compare with Smart-seq2 results!
