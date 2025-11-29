# STARsolo Pipeline Setup - Complete Summary

## ✅ What Has Been Accomplished

### 1. **Environment Setup** ✅
- Installed STAR aligner (v2.7.11b)
- Verified samtools (v1.19.2) is available
- Created all necessary directories
- System confirmed: 120GB RAM, 20 CPU cores

### 2. **Reference Genome & Annotations** ✅
- Downloaded hg38 (GRCh38.primary_assembly) - 805MB compressed
- Copied GENCODE v45 gene annotations (1.5GB)
- Copied RepeatMasker TE annotations (270MB)
- Converted RepeatMasker BED → GTF format
- Combined genes + TEs into single GTF (268,404 transcripts, 5.6M junctions)

### 3. **Barcode Analysis** ✅
Examined FASTQ files to determine structure:
- **R1 files:** 19bp barcodes (16bp cell barcode + 3bp UMI)
- **R2 files:** 63bp cDNA reads
- **Quality:** Q30+ throughout (excellent)
- **Conclusion:** No whitelist needed (Smart-seq2 - each sample = one cell)

### 4. **STAR Genome Index** 🔄 IN PROGRESS
- Started: Nov 27, 23:07:11
- Status: Sorting suffix array (started 23:08:14)
- Expected completion: 30-60 minutes total
- Command includes TE annotations for accurate alignment rates
- Will use 6M junction limit (TEs add many junctions)

### 5. **Alignment Pipeline Script** ✅
Created comprehensive Python script with:
- Sequential processing (one sample at a time, full resources)
- STARsolo with SmartSeq mode
- CB (cell barcode) and UB (UMI) tagging
- Progress tracking with ETA
- Error handling and logging
- BAM validation (checks for required tags)
- TE-friendly multi-mapper parameters
- Automatic QC report generation

## 📂 Files Created

```
/home/jacobc/hcaTE/
├── scripts/
│   ├── starsolo_align_all_samples.py      ← Main alignment pipeline
│   ├── test_single_alignment.py           ← Test on one sample first
│   └── bed_to_gtf.py                      ← RepeatMasker conversion utility
│
├── annotations/
│   ├── gencode.v45.primary_assembly.annotation.gtf  (1.5GB)
│   ├── hg38_rmsk.bed                                 (270MB)
│   ├── hg38_rmsk.gtf                                 (converted)
│   └── combined_genes_TEs.gtf                        (genes + TEs)
│
├── genome/
│   └── GRCh38.primary_assembly.genome.fa            (3.0GB)
│
├── star_index/                                      (building...)
│   └── [STAR index files]
│
├── aligned_bams/                                    (ready for output)
├── qc/                                              (ready for logs)
├── results/                                         (ready for scTE)
│
└── Documentation/
    ├── QUICK_START.md          ← How to run the pipeline
    ├── PROGRESS_REPORT.md      ← Detailed setup documentation
    ├── AGENT_INSTRUCTIONS.md   ← Original task description
    └── README.md               ← Project background
```

## 🎯 Pipeline Parameters

### STARsolo Configuration
```python
--soloType SmartSeq           # Smart-seq2 mode
--soloCBlen 16                # Cell barcode: 16bp
--soloUMIlen 3                # UMI: 3bp
--soloCBwhitelist None        # No whitelist needed
--soloStrand Unstranded       # Library prep type
```

### TE-Optimized Alignment
```python
--outFilterMultimapNmax 100          # Allow many multi-mappers
--outFilterMismatchNmax 999          # Flexible mismatch tolerance
--outFilterMismatchNoverLmax 0.04    # 4% mismatch rate allowed
--winAnchorMultimapNmax 200          # Anchor multi-mapping
--outMultimapperOrder Random         # Random placement
```

### Output Tags (Critical for scTE!)
```python
--outSAMattributes NH HI AS nM CB UB
# CB = Cell Barcode tag
# UB = UMI tag
```

## 📊 Expected Results

### Per Sample
- **Input:** 2 FASTQ files (~500MB + 1.3GB)
- **Output:** ~2-5GB BAM + index + logs
- **Time:** 1-3 minutes per sample
- **Alignment rate:** 60-80% total (unique + multi-mapped)

### Total (160 Samples)
- **Total time:** 3-8 hours
- **Disk usage:** ~800GB for all BAMs
- **Output:** 160 CB-tagged BAMs ready for scTE

## 🚀 Next Steps (When Ready)

### 1. Wait for STAR Index to Complete
```bash
# Monitor progress:
tail -f /home/jacobc/hcaTE/star_index/Log.out

# Check when done:
ls -lh /home/jacobc/hcaTE/star_index/Genome
# Should see ~30GB file
```

### 2. Test on One Sample
```bash
cd /home/jacobc/hcaTE
python3 scripts/test_single_alignment.py
```

### 3. Run Full Pipeline
```bash
python3 scripts/starsolo_align_all_samples.py
# Will ask for confirmation before starting
# Shows progress every 10 samples
# Can safely Ctrl+C and resume later if needed
```

### 4. Validate Results
```bash
# Check summary
cat /home/jacobc/hcaTE/qc/alignment_logs/alignment_results.txt

# Quick validation
ls /home/jacobc/hcaTE/aligned_bams/*/SRR*_Aligned.sortedByCoord.out.bam | wc -l
# Should be 160
```

### 5. Proceed to scTE Analysis
Once alignment is complete, the CB-tagged BAMs can be used directly with scTE for TE quantification.

## 🔍 Quality Control Checks

The pipeline automatically validates:
- ✅ CB tags present in BAM files
- ✅ UB tags present in BAM files
- ✅ Alignment rates per sample
- ✅ Multi-mapping rates
- ✅ Total read counts

Low-quality samples (poor alignment rates) can be filtered out based on the QC report.

## 📝 Key Differences from Previous Pipeline

| Previous (HumanBam2scTE) | New (hcaTE) |
|-------------------------|-------------|
| ❌ Bulk RNA-seq alignment | ✅ Single-cell with STARsolo |
| ❌ No cell barcodes | ✅ CB and UB tags |
| ❌ Treated as one sample | ✅ 160 individual cells |
| ❌ No barcode extraction | ✅ Barcodes from R1 files |
| ❌ Wrong for scTE | ✅ Compatible with scTE |

## 🎓 Technical Notes

### Why Smart-seq2?
- Each sample IS one cell (pre-sorted)
- Barcodes identify cells within experimental batches
- Don't need 10x-style whitelist
- Higher depth per cell than droplet methods

### Why Include TEs in STAR Index?
- Improves alignment rate accuracy
- Helps identify low-quality samples
- TEs are highly repetitive → multi-mappers
- scTE will quantify TEs regardless

### Why Sequential Processing?
- STAR needs 30-40GB RAM per run
- 20 cores fully utilized per sample = fast
- Simpler error handling
- Better resource utilization than parallel

### Why These Multi-mapper Parameters?
- TEs often have dozens of copies in genome
- Need to allow multi-mapping to capture TE reads
- Random placement prevents bias
- Same parameters as successful HumanBam2scTE project

## 📧 Support

If you encounter issues:
1. Check logs in `/home/jacobc/hcaTE/qc/alignment_logs/`
2. Review QUICK_START.md for troubleshooting
3. Verify STAR index completed successfully
4. Check disk space: `df -h`

## 🏆 Success Criteria

You'll know the pipeline worked when:
- [ ] All 160 BAM files created
- [ ] CB tags present in samtools view output
- [ ] Average alignment rate >60%
- [ ] alignment_results.txt shows mostly successful samples
- [ ] scTE can read the BAM files

---

**Status:** Ready to align once STAR index completes (~20-40 more minutes)

**Date:** November 27, 2025

**System:** spark-bd86 (120GB RAM, 20 cores, Ubuntu ARM64)
