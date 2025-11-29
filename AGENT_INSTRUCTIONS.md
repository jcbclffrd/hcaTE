# Agent Instructions: Set Up Single-Cell RNA-seq Alignment Pipeline

> **🤖 This pipeline was built by GitHub Copilot coding agent**

## Your Mission
Set up a **STARsolo** or **Cell Ranger** pipeline to align 160 Smart-seq2 single-cell RNA-seq samples with barcodes. The data is in `/home/jacobc/hcaTE/sra_downloads/` and needs proper single-cell alignment with cell barcode tagging.

## Critical Context (READ THIS FIRST!)

### Data Structure
- **160 samples** in `/home/jacobc/hcaTE/sra_downloads/`
- Each sample has **TWO FASTQ files**:
  - `SRR*_1.fastq.gz`: **19bp cell barcodes** (16bp barcode + UMI, DO NOT map to genome)
  - `SRR*_2.fastq.gz`: **63bp cDNA reads** (map these to genome)
- **Library type**: Smart-seq2 with barcoding (bc-Smart-seq2)
- **Quality**: Excellent (Q30+), no trimming needed

### What You Need to Do

1. **Install STARsolo** (preferred) or Cell Ranger
   - STARsolo is better for custom TE analysis
   - Can work with flexible barcode formats

2. **Set up genome reference**
   - Use hg38 (GRCh38)
   - Copy annotations from `/home/jacobc/HumanBam2scTE/annotations/`:
     - `gencode.v45.primary_assembly.annotation.gtf`
     - `hg38_rmsk.bed` (for TE analysis)
   - Build STAR index with TE annotations included

3. **Determine barcode structure**
   - Examine the 19bp reads to understand layout
   - Likely: 16bp cell barcode + variable UMI
   - Check if whitelist exists or generate from data

4. **Create alignment script**
   - Process all 160 samples
   - Use STARsolo with `--soloType SmartSeq`
   - Extract barcodes from R1 files
   - Map R2 files to genome
   - Output: BAM files with CB (cell barcode) and UB (UMI) tags

5. **Run alignment on all samples**
   - Parallel processing recommended (160 samples)
   - Output directory: `/home/jacobc/hcaTE/aligned_bams/`
   - Keep organized by sample

6. **Validate outputs**
   - Check BAM files have CB and UB tags
   - Verify alignment rates (>50% expected)
   - Count cells detected per sample
   - Generate QC report

## Expected STARsolo Command Structure

```bash
STAR --runThreadN 8 \
     --genomeDir /path/to/star_index \
     --readFilesIn SRR*_2.fastq.gz SRR*_1.fastq.gz \
     --readFilesCommand zcat \
     --soloType SmartSeq \
     --soloUMIlen [determine from data] \
     --soloCBlen 16 \
     --soloCBwhitelist [path or None] \
     --soloStrand Unstranded \
     --outSAMattributes NH HI AS nM CB UB \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix /home/jacobc/hcaTE/aligned_bams/SRR*/
```

## Important Notes

- **DO NOT** map the 19bp barcode reads (R1) to the genome
- **DO NOT** treat this as bulk RNA-seq (it's single-cell!)
- The barcodes need to be extracted and tagged, not mapped
- Downstream tool (scTE) requires CB-tagged BAMs

## Deliverables

When you're done, provide:
1. ✅ STAR index built with TE annotations
2. ✅ Alignment script for batch processing
3. ✅ Aligned BAM files with CB/UB tags (in `/home/jacobc/hcaTE/aligned_bams/`)
4. ✅ QC report (alignment rates, cells per sample, etc.)
5. ✅ Documentation of barcode structure determined

## Questions to Answer

As you work, figure out:
- What is the exact barcode structure in the 19bp reads?
- Is there a barcode whitelist, or should you generate one?
- What are the alignment rates for the 63bp cDNA reads?
- How many cells are detected per sample?
- Are the BAM files compatible with scTE (CB and UB tags present)?

## Read the Full Context

See `/home/jacobc/hcaTE/README.md` for complete project details, data quality info, and why this approach is necessary (previous pipeline treated it as bulk RNA-seq, which was wrong).

---

**Start by**: Reading the README.md, then examine a few sample FASTQ files to understand the barcode structure, then set up the STAR index.
