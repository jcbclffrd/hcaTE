# Brain Region Subsetting and SIGLEC1 Expression Analysis

## Overview

This analysis was performed in response to the paper's Figure 2, which shows SIGLEC1 expression differences between two brain regions:
- **LPS** (Lateral Parietal/Superior Parietal Lobe)
- **GFS** (Gyrus Frontalis Superior / Superior Frontal Gyrus)

## Key Findings

### SIGLEC1 Expression Summary (Pseudobulk Level)

**Total counts**: 5 SIGLEC1 counts across all 31 samples  
**Samples with expression**: 2 out of 31 samples

| Condition | Total Counts | Mean | # Samples |
|-----------|--------------|------|-----------|
| AD        | 0            | 0.00 | 12        |
| CTR       | 0            | 0.00 | 9         |
| CTR+      | 1            | 0.17 | 6         |
| MCI       | 4            | 1.00 | 4         |

**Samples with SIGLEC1 expression:**
1. **SRR11272123** (CTR+): 1 count
2. **SRR11272139** (MCI): 4 counts

## Why SIGLEC1 Counts Are So Low

SIGLEC1 is a **microglia-specific marker** that is highly diluted in our pseudobulk analysis:

1. **Cell-type heterogeneity**: Each bc-SmartSeq2 sample contains ~100-150 cells from multiple cell types (microglia, neurons, astrocytes, etc.)
2. **Pseudobulk aggregation**: We aggregated ~82 cells per sample on average (12,942 total cells / 157 samples)
3. **Dilution effect**: Even if SIGLEC1 is expressed in 10-20% of microglia within a sample, it gets averaged across all cell types

This is why the paper's analysis focused on:
- **Single-cell resolution** analysis
- **Microglia-specific** subsets
- **Brain region-specific** comparisons (LPS vs GFS)

## Brain Region Mapping Challenge

Unfortunately, we encountered a mapping challenge:

### What We Have:
1. **SRR IDs**: Our sample identifiers (SRR11272119, SRR11272123, etc.)
2. **GSM to Brain Region mapping**: From GEO series matrix (`gsm_to_brain_region.csv`)
   - 177 samples from **LPS** (lateral parietal/superior parietal lobe)
   - 12 samples from **GFS** (superior frontal gyrus)

### What We Need:
- **SRR to GSM mapping**: Direct link between our processed samples (SRR IDs) and the GEO sample accessions (GSM IDs)

### Why This Is Challenging:
- The GEO series matrix provides GSM → SRX (SRA experiment) mappings
- We need SRX → SRR (SRA run) mappings, which require the SRA RunInfo table
- The NCBI SRA download link was not accessible during this analysis

## Solutions to Complete Regional Analysis

### Option 1: Download SRA RunInfo (Recommended)
```bash
# Visit NCBI SRA and download RunInfo table for SRP252065
# https://www.ncbi.nlm.nih.gov/sra?term=SRP252065

# Then use the mapping:
# SRR → SRX (from RunInfo)
# SRX → GSM (from series matrix)
# GSM → Brain Region (from our gsm_to_brain_region.csv)
```

### Option 2: Check Paper's Supplementary Data
The original paper (Alsema et al. 2020) likely includes a supplementary table with:
- Sample IDs
- Brain regions (LPS vs GFS)
- Condition (AD, CTR, CTR+, MCI)
- Donor information

### Option 3: Single-Cell Level Analysis
Since SIGLEC1 has such low pseudobulk counts, a better approach would be:
1. **Extract single-cell data** from scTE output files (`scTE_output/SRR*/SRR*.csv`)
2. **Subset to microglia cells** (if cell-type annotations available)
3. **Analyze SIGLEC1 at single-cell resolution** within brain regions

This would match the paper's analytical approach more closely.

## Files Generated

1. **`siglec1_expression_by_sample.csv`**: SIGLEC1 counts for all 31 samples
2. **`siglec1_summary_statistics.csv`**: Summary stats by condition
3. **`gsm_to_brain_region.csv`**: GSM accession to brain region mapping (189 samples)
   - 177 LPS samples
   - 12 GFS samples

## Next Steps

To reproduce the paper's regional analysis, you would need to:

1. **Obtain SRR → Brain Region mapping** using one of the options above

2. **Consider single-cell analysis** instead of pseudobulk:
   ```python
   # Pseudocode for single-cell SIGLEC1 analysis
   for srr_id in samples:
       sc_data = load_scTE_output(f"scTE_output/{srr_id}/{srr_id}.csv")
       brain_region = get_brain_region(srr_id)  # LPS or GFS
       siglec1_cells = sc_data['SIGLEC1'] > 0
       
       # Analyze SIGLEC1+ cells by region and condition
   ```

3. **Cell-type annotation** (if not already done):
   - Use marker genes to identify microglia (CX3CR1, P2RY12, TMEM119, etc.)
   - Subset to microglia before analyzing SIGLEC1
   - This would dramatically reduce dilution effect

## References

- **Paper**: Alsema AM et al. (2020) "Profiling Microglia From Alzheimer's Disease Donors and Non-demented Elderly in Acute Human Postmortem Cortical Tissue" Frontiers in Molecular Neuroscience. PMID: 33192286
- **GEO Accession**: GSE146639
- **SRA Project**: SRP252065
- **Analysis Date**: November 30, 2024
