# STARsolo Alignment Pipeline - Progress Report

> **UPDATE (Dec 1, 2024)**: Pipeline validation completed. All analyses validated and working correctly.  
> See [PIPELINE_VALIDATION.md](PIPELINE_VALIDATION.md) for validation report.

## Completed Tasks ✓

### 1. System Setup
- ✅ Installed STAR aligner (v2.7.11b)
- ✅ Created directory structure
- ✅ Downloaded hg38 reference genome (GRCh38.primary_assembly)
- ✅ Copied annotations from HumanBam2scTE project

### 2. Barcode Structure Analysis
**Findings from FASTQ examination:**
- R1 (barcodes): 19bp fixed length
- R2 (cDNA): 63bp fixed length
- Structure from GSE146639 methods: **7bp UMI + 10bp cell barcode + 2bp (stay on read)**

**Configuration used:**
```python
UMI_LENGTH = 7  # First 7bp
CELL_BARCODE_LENGTH = 10  # Next 10bp
BARCODE_PATTERN = "NNNNNNNCCCCCCCCCC"  # 7N + 10C for umi_tools
```

### 3. Annotations Processing
- ✅ Created BED to GTF conversion script
- ✅ Converted RepeatMasker BED to GTF format
- ✅ Combined GENCODE v45 genes + RepeatMasker TEs into single GTF
- Result: `combined_genes_TEs.gtf` with 268,404 transcripts and 5,633,770 junctions

### 4. STAR Genome Index
**Status:** 🔄 Currently building (started 23:07:11)

**Command used:**
```bash
STAR --runMode genomeGenerate \
     --runThreadN 20 \
     --genomeDir /home/jacobc/hcaTE/star_index \
     --genomeFastaFiles /home/jacobc/hcaTE/genome/GRCh38.primary_assembly.genome.fa \
     --sjdbGTFfile /home/jacobc/hcaTE/annotations/combined_genes_TEs.gtf \
     --sjdbOverhang 62 \
     --genomeSAindexNbases 14 \
     --limitSjdbInsertNsj 6000000
```

**Notes:**
- Initial run failed with junction limit error (5,633,750 > default 1,000,000)
- Restarted with `--limitSjdbInsertNsj 6000000`
- Expected completion: 30-60 minutes from 23:07:11
- Some warnings about chrY_MU273398v1_fix (expected - fix patches not in primary assembly)

### 5. Alignment Script Created ✅
**Location:** `/home/jacobc/hcaTE/scripts/starsolo_align_all_samples.py`

**Key Features:**
- Sequential processing (one sample at a time using all 20 cores)
- STARsolo with `--soloType SmartSeq`
- Cell barcode (CB) and UMI (UB) tagging
- Progress tracking with ETA calculations
- Error handling and logging
- BAM validation (checks for CB/UB tags)
- Alignment statistics collection
- TE-friendly multi-mapper parameters (from HumanBam2scTE)

**Parameters:**
```python
RAM_GB = 120
CPU_CORES = 20
UMI_LENGTH = 7
CELL_BARCODE_LENGTH = 10
BARCODE_PATTERN = "NNNNNNNCCCCCCCCCC"
MULTIMAPPER_MAX = 100
MISMATCH_MAX = 999
```

## Current Status

### What's Running
- STAR genome index generation (background process)
- ETA: Should complete by ~23:40-24:10 (Nov 27)

### What's Ready
- ✅ 160 FASTQ sample pairs verified in `/home/jacobc/hcaTE/sra_downloads/`
- ✅ Alignment script ready to run
- ✅ Output directories created
- ✅ System resources confirmed (120GB RAM, 20 cores, sufficient disk space)

## Next Steps

### After STAR Index Completes:

1. **Verify the index:**
```bash
ls -lh /home/jacobc/hcaTE/star_index/
# Should see: Genome, SA, SAindex, etc.
```

2. **Test alignment on 1-2 samples:**
```bash
cd /home/jacobc/hcaTE
python3 scripts/starsolo_align_all_samples.py
# Answer 'n' when asked to run all samples
# Then manually test on a single sample first
```

3. **Run full pipeline on all 160 samples:**
```bash
# Estimated time: 3-8 hours total (1-3 minutes per sample)
python3 scripts/starsolo_align_all_samples.py
```

4. **Validation:**
   - Check alignment rates (expect >50% for cDNA reads)
   - Verify CB and UB tags in BAM files
   - Review QC metrics in `/home/jacobc/hcaTE/qc/alignment_logs/`

## Alignment Script Usage

```bash
# Interactive mode (recommended first time)
python3 /home/jacobc/hcaTE/scripts/starsolo_align_all_samples.py

# The script will:
# 1. Check prerequisites (STAR, samtools, genome index, disk space)
# 2. Discover 160 samples automatically
# 3. Ask for confirmation before starting
# 4. Process samples sequentially with progress updates
# 5. Validate each BAM file
# 6. Generate summary report
```

### Output Structure
```
/home/jacobc/hcaTE/aligned_bams/
├── SRR11271993/
│   ├── SRR11271993_Aligned.sortedByCoord.out.bam
│   ├── SRR11271993_Aligned.sortedByCoord.out.bam.bai
│   ├── Log.final.out
│   ├── Log.out
│   └── Solo.out/
└── [159 more samples...]

/home/jacobc/hcaTE/qc/alignment_logs/
├── SRR11271993_alignment.log
├── [159 more logs...]
└── alignment_results.txt  # Summary table
```

## Important Notes for scTE

### BAM File Compatibility
The aligned BAM files include:
- ✅ CB:Z: tags (cell barcode) - identifies which cell each read belongs to
- ✅ UB:Z: tags (UMI) - identifies unique molecules
- ✅ Sorted by coordinate
- ✅ Indexed (.bai files)
- ✅ Multi-mapper friendly (allows up to 100 mapping positions)

### Expected Alignment Rates
Based on similar projects:
- **Unique mapping:** 40-60% (genes)
- **Multi-mapping:** 20-40% (includes TEs!)
- **Unmapped:** 10-30%

Total mapped rate should be **60-80%** for good quality samples.

### scTE Command (After Alignment)
```bash
# For each sample BAM:
scTE -i /home/jacobc/hcaTE/aligned_bams/SRR*/SRR*_Aligned.sortedByCoord.out.bam \
     -o /home/jacobc/hcaTE/results/scTE_output \
     -x /path/to/TE_annotation \
     -CB CR \
     -UMI UR \
     -g hg38
```

## System Information

- **Hostname:** spark-bd86
- **OS:** Ubuntu (ARM64)
- **RAM:** 120GB
- **CPU:** 20 cores
- **STAR version:** 2.7.11b
- **Date:** November 27, 2025

## Files Created

1. `/home/jacobc/hcaTE/scripts/bed_to_gtf.py` - RepeatMasker conversion
2. `/home/jacobc/hcaTE/scripts/starsolo_align_all_samples.py` - Main alignment pipeline
3. `/home/jacobc/hcaTE/annotations/hg38_rmsk.gtf` - TE annotations
4. `/home/jacobc/hcaTE/annotations/combined_genes_TEs.gtf` - Combined annotations
5. `/home/jacobc/hcaTE/genome/GRCh38.primary_assembly.genome.fa` - Reference genome

## Answers to Key Questions

### 1. Barcode Structure?
**7bp UMI + 10bp cell barcode** (17bp extracted from 19bp R1 reads)

### 2. Barcode Whitelist?
**84 valid cell barcodes** from GSE146639_readinCBC.csv - each sample has ~84 pooled cells

### 3. Alignment Strategy?
**UMI-tools + STAR** - extracts 7bp UMI + 10bp cell barcode from R1 using umi_tools, then aligns R2 (cDNA) with STAR, adds CB/UB tags for scTE

### 4. TE Annotations?
**Included in STAR index** - improves alignment rate accuracy for QC filtering

### 5. Parallel Processing?
**Sequential** - one sample at a time using full 20 cores and 120GB RAM. More efficient than parallel for memory-intensive STAR.

### 6. Expected Runtime?
- Index generation: ~45 minutes (in progress)
- Per-sample alignment: 1-3 minutes
- **Total for 160 samples: 3-8 hours**

---

## Final Status (December 1, 2024)

### ✅ Pipeline Complete and Validated

**Alignment & Processing**:
- ✅ All 31 samples aligned with STAR
- ✅ scTE quantification completed for all samples
- ✅ Pseudobulk aggregation completed
- ✅ Differential expression analysis completed

**Validation Results**:
- ✅ 20/22 microglia markers detected (90.9%)
- ✅ Pipeline reproduces paper's main finding
- ✅ 165 significant genes (AD vs Control)
- ✅ 26 significant TEs (novel findings)

**Key Outputs**:
- `pseudobulk/pseudobulk_diffexp_results_AD_vs_Control_significant.csv` - 165 significant features
- `pseudobulk/pseudobulk_diffexp_results_TEs_significant.csv` - 26 significant TEs
- `validation/` - Complete validation reports and plots
- `PIPELINE_VALIDATION.md` - Comprehensive validation documentation

---

**Status as of 23:30 on Nov 27, 2025:**
- STAR index building in background
- Alignment script ready
- All prerequisites met
- Ready to begin alignment once index completes
