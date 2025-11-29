# Quick Start Guide - bc-Smart-seq2 Alignment Pipeline

## ⚠️ IMPORTANT UPDATE - Corrected Pipeline

After examining the GEO data and methods, we discovered:
- Each FASTQ = ~84 cells pooled together (NOT 1 cell per FASTQ!)
- R1 contains: **7bp UMI + 10bp cell barcode** (positions 1-7 and 8-17)
- Need to use **UMI-tools** for proper barcode extraction with error correction
- Then **STAR** for alignment
- Then add **CB/UB tags** to BAM for scTE

## Current Status

🔄 **STAR genome index is building** (started Nov 27 at 23:08:14)
- Expected completion: ~30-60 minutes from start
- Monitor progress: `tail -f /home/jacobc/hcaTE/star_index/Log.out`

## Once Index Completes

### Step 1: Verify Index is Complete

```bash
# Check for the Genome file (largest file, ~30GB)
ls -lh /home/jacobc/hcaTE/star_index/Genome

# If it exists, index is complete!
```

### Step 2: Run the NEW Hybrid Pipeline (All 160 Samples)

```bash
cd /home/jacobc/hcaTE
python3 scripts/umitools_star_align_all_samples.py
```

**This uses the corrected approach:**
1. UMI-tools extract: Extracts 7bp UMI + 10bp cell barcode from R1
2. STAR align: Maps R2 (cDNA) to genome
3. Add CB/UB tags: Converts read name info to BAM tags
4. Validate: Ensures CB/UB tags present for scTE

The script will:
1. Check prerequisites automatically
2. Find all 160 samples
3. Ask for confirmation
4. Process samples sequentially (one at a time)
5. Show progress with ETA
6. Validate each BAM file
7. Generate summary report

**Estimated time:** 3-8 hours total (1-3 minutes per sample)

**Progress updates:** Every 10 samples

## Monitoring

### Check Alignment Progress
```bash
# Count completed samples
ls -d /home/jacobc/hcaTE/aligned_bams/SRR* | wc -l

# View real-time progress (if running in background)
tail -f /home/jacobc/hcaTE/qc/alignment_logs/alignment_results.txt
```

### Check Disk Space
```bash
df -h /home/jacobc/hcaTE
```

### Check a Sample's Alignment Stats
```bash
cat /home/jacobc/hcaTE/aligned_bams/SRR11271993/Log.final.out
```

## Troubleshooting

### If STAR Index Failed
```bash
# Check the log for errors
tail -50 /home/jacobc/hcaTE/star_index/Log.out

# If needed, rebuild:
rm -rf /home/jacobc/hcaTE/star_index/*
# Then re-run the STAR genome generate command from PROGRESS_REPORT.md
```

### If Alignment Fails on Some Samples
- Check logs in `/home/jacobc/hcaTE/qc/alignment_logs/`
- Failed samples are listed at the end of the run
- Can re-run just specific samples by modifying the script

### If BAMs Don't Have CB Tags
- Check STARsolo parameters in the script
- Verify barcode files (_1.fastq.gz) exist
- Review STAR log for warnings about solo parameters

## Output Structure

```
/home/jacobc/hcaTE/
├── aligned_bams/
│   ├── SRR11271993/
│   │   ├── SRR11271993_Aligned.sortedByCoord.out.bam   ← Main output
│   │   ├── SRR11271993_Aligned.sortedByCoord.out.bam.bai
│   │   ├── Log.final.out                               ← Alignment stats
│   │   └── Solo.out/                                   ← STARsolo outputs
│   └── [159 more samples...]
│
└── qc/
    └── alignment_logs/
        ├── SRR11271993_alignment.log                   ← Detailed logs
        ├── [159 more logs...]
        └── alignment_results.txt                       ← Summary table
```

## Next Steps After Alignment

### 1. Validate BAM Files
```bash
# Quick check: count BAMs with CB tags
for bam in /home/jacobc/hcaTE/aligned_bams/*/SRR*_Aligned.sortedByCoord.out.bam; do
    samtools view "$bam" | head -100 | grep -q "CB:Z:" && echo "✓ $bam" || echo "✗ $bam"
done | grep "✗" | wc -l
# Should be 0 (no missing CB tags)
```

### 2. Generate QC Report
```bash
# Compile alignment rates
grep "Uniquely mapped reads %" /home/jacobc/hcaTE/aligned_bams/*/Log.final.out > qc/unique_mapping_rates.txt
```

### 3. Run scTE for TE Quantification
```bash
# See scTE documentation for exact command
# The BAM files are now ready with CB and UB tags
```

## Files Reference

| File | Purpose |
|------|---------|
| `scripts/starsolo_align_all_samples.py` | Main pipeline - align all 160 samples |
| `scripts/test_single_alignment.py` | Test script - align one sample |
| `scripts/bed_to_gtf.py` | Convert RepeatMasker BED to GTF |
| `PROGRESS_REPORT.md` | Detailed progress and configuration |
| `AGENT_INSTRUCTIONS.md` | Original instructions |
| `README.md` | Project background |

## System Resources

- **RAM:** 120GB (using ~80-90GB per sample)
- **CPU:** 20 cores (all used per sample)
- **Disk:** Check with `df -h`
- **Processing:** Sequential (NOT parallel)

## Contact Info

- Project: hcaTE (Human Cell Atlas - Transposable Elements)
- Data: 160 Smart-seq2 samples from PRJNA611563
- Goal: Single-cell TE expression analysis in Alzheimer's Disease microglia

---

**Questions Answered:**

✅ **Barcode structure:** 16bp cell barcode + 3bp UMI (19bp total)  
✅ **Whitelist:** None needed (Smart-seq2 - each sample is one cell)  
✅ **Alignment tool:** STARsolo with SmartSeq mode  
✅ **TE annotations:** Included in index for better QC  
✅ **Output format:** CB-tagged BAMs compatible with scTE  

---

**Last Updated:** Nov 27, 2025 at 23:35
**Next Action:** Wait for STAR index to complete, then run test_single_alignment.py
