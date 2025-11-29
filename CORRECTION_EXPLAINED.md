# CRITICAL CORRECTION - bc-Smart-seq2 Data Structure

## What We Initially Thought (WRONG)

❌ Each FASTQ pair = 1 cell  
❌ R1 = 19bp: 16bp cell barcode + 3bp UMI  
❌ Use STARsolo SmartSeq mode  
❌ No whitelist needed  

## What's Actually True (CORRECT)

✅ **Each FASTQ pair = ~84 cells pooled together**  
✅ **R1 = 19bp: 7bp UMI + 10bp cell barcode + 2bp other**  
✅ **Use UMI-tools + STAR (hybrid approach)**  
✅ **Need 84-barcode whitelist from GEO**  

## The Real bc-Smart-seq2 Protocol

From the study methods (PMID: 33192286):

1. **Single cells FACS-sorted into 384-well plates**
   - Each well gets a unique barcoded oligo-dT primer
   - 84 cells are pooled into one library

2. **Barcode Structure in R1:**
   ```
   Position 1-7:   UMI (7bp) - identifies unique mRNA molecules
   Position 8-17:  Cell Barcode (10bp) - identifies which of 84 cells
   Position 18-19: Likely technical/adapter sequence
   ```

3. **Sequencing:**
   - R1: 19bp (barcode read)
   - R2: 63bp (cDNA read - this gets mapped to genome)

4. **Cell Counts:**
   - 160 samples in our data
   - ~84 cells per sample
   - **Total: ~13,440 cells** (not 160!)

## Why This Matters for scTE

scTE needs BAM files with:
- **CB tags** to identify which cell each read belongs to
- **UB tags** for UMI deduplication

With the corrected understanding:
- We need to demultiplex the 84 cells within each sample
- Extract both cell barcodes AND UMIs
- Add both as tags to the BAM

## The Corrected Pipeline

### Step 1: UMI-tools extract
```bash
umi_tools extract \
  --bc-pattern=NNNNNNNCCCCCCCCCC \  # 7N (UMI) + 10C (Cell Barcode)
  --whitelist=cell_barcode_whitelist.txt \  # 84 valid barcodes
  --error-correct-cell \  # Fix sequencing errors in barcodes
  --filter-cell-barcode \  # Keep only whitelisted barcodes
  -I R1.fastq.gz \
  --read2-in R2.fastq.gz \
  --read2-out R2_extracted.fastq.gz
```

**Output:** R2_extracted.fastq.gz with read names like:
```
@READID_CELLBARCODE_UMI
```

### Step 2: STAR align
```bash
STAR --readFilesIn R2_extracted.fastq.gz \
     --genomeDir star_index \
     --outSAMtype BAM Unsorted \
     [TE-friendly parameters...]
```

**Output:** BAM with read names containing barcode+UMI

### Step 3: Add CB/UB tags
Parse read names and add to BAM as proper tags:
```
CB:Z:CELLBARCODE
UB:Z:UMI
```

### Step 4: Sort, index, validate

**Output:** BAM files ready for scTE!

## Files We Downloaded from GEO

1. **`GSE146639_readinCBC.csv`** - Cell barcode whitelist
   - 84 barcodes, each 10bp
   - Example: CTGAATATTT, GACAGAGGGC, etc.

2. **`GSE146639_UMItools_GEO.py`** - Their processing script
   - Shows exact barcode pattern: `NNNNNNNCCCCCCCCCC`
   - Uses UMI-tools for extraction
   - Uses HISAT2 for alignment (we'll use STAR instead)
   - Generates count matrices (we need BAMs instead)

## Key Differences from Original Pipeline

| Original Study | Our Pipeline |
|---------------|-------------|
| HISAT2 align | ✅ STAR align (better for TEs) |
| featureCounts | ✅ Not needed (scTE does this) |
| UMI-tools count → matrix | ✅ Keep as BAM with CB/UB tags |
| Count matrices output | ✅ Tagged BAMs for scTE |

## Why We Needed This Correction

Our initial STARsolo script would have:
- ❌ Treated each pool as 1 cell (wrong!)
- ❌ Missed 83 cells per pool
- ❌ Lost ~13,000 cells from the analysis
- ❌ Given completely wrong results

The corrected pipeline will:
- ✅ Properly demultiplex all ~84 cells per pool
- ✅ Extract and tag both UMIs and cell barcodes
- ✅ Produce scTE-compatible BAMs
- ✅ Enable proper single-cell TE analysis

## Expected Output

After running the corrected pipeline:

```
/home/jacobc/hcaTE/aligned_bams/
├── SRR11271993/
│   ├── SRR11271993_Aligned.sortedByCoord.out.bam  ← ~84 cells in this BAM
│   ├── SRR11271993_Aligned.sortedByCoord.out.bam.bai
│   └── Log.final.out
├── SRR11271994/  ← Another ~84 cells
└── [158 more samples...]

Total: 160 pools × 84 cells/pool = ~13,440 cells
Each BAM has CB and UB tags for scTE analysis
```

## Running scTE After Alignment

```bash
# scTE can handle multiple cells in one BAM via CB tags
scTE -i /home/jacobc/hcaTE/aligned_bams/*/SRR*_Aligned.sortedByCoord.out.bam \
     -o /home/jacobc/hcaTE/results/scTE_output \
     -x hg38_rmsk.gtf \
     -CB CB \
     -UMI UB \
     -g hg38
```

scTE will:
1. Read CB tags to identify individual cells
2. Use UB tags for deduplication
3. Quantify TEs per cell
4. Output: **cell × TE expression matrix** for ~13,440 cells

---

**Date:** November 27, 2025  
**Correction made:** After examining GEO supplementary files and study methods  
**New script:** `umitools_star_align_all_samples.py`  
**Status:** Ready to run once STAR index completes
