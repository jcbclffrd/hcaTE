# Pseudo-bulk Aggregation for Differential Expression Analysis

> **🤖 Built by GitHub Copilot coding agent**

## ⚠️ CRITICAL: Why Pseudo-bulk is Required

**You CANNOT run differential expression directly on single-cell data!**

### The Problem:
- **Biological replicates = Donors**, not individual cells
- Cells from the same donor are **technical replicates** (not independent)
- DESeq2/edgeR require independent biological replicates
- Running DE on individual cells violates statistical assumptions → inflated false positives

### The Solution: Pseudo-bulk Aggregation
Sum counts from all cells belonging to each donor to create donor-level pseudo-bulk samples.

---

## Overview of Your Data

### Current State (Single-Cell):
```
157 samples × 12,942 cells × 76,393 features (genes + TEs)
├── Each sample: 60-84 cells (avg 82.4 cells/sample)
├── scTE output: /home/jacobc/hcaTE/scTE_output/SRR*/SRR*.csv
└── Format: rows = cells (barcodes), columns = genes/TEs
```

### Target State (Pseudo-bulk):
```
N donors × 76,393 features
├── Each row = 1 donor (aggregated from all their cells)
├── Format: Same as expression_matrix2.csv from HumanBam2scTE
└── Compatible with existing DESeq2 script
```

---

## Step 1: Map Samples to Donors

### Find Sample Metadata

You need to determine which samples belong to which donors. Check:

```bash
# Look for sample metadata
cat /home/jacobc/HumanBam2scTE/sample_types.txt

# Or check GEO for sample information
# GSE146639: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146639
```

**What you need:**
- Sample → Donor mapping
- Sample → Condition mapping (AD vs Control)

Expected format:
```
Sample_ID    Donor_ID    Condition
SRR11271993  Donor1      AD
SRR11271994  Donor1      AD
SRR11271995  Donor2      Control
...
```

---

## Step 2: Aggregate Single-Cell Counts to Pseudo-bulk

### Create Aggregation Script

```python
#!/usr/bin/env python3
"""
Aggregate single-cell scTE counts to pseudo-bulk donor-level counts
"""

import pandas as pd
import numpy as np
import glob
from pathlib import Path

# Configuration
SCTE_DIR = Path("/home/jacobc/hcaTE/scTE_output")
OUTPUT_DIR = Path("/home/jacobc/hcaTE/pseudobulk")
OUTPUT_DIR.mkdir(exist_ok=True)

# Sample metadata (YOU NEED TO FILL THIS IN!)
# Format: {sample_id: (donor_id, condition)}
SAMPLE_METADATA = {
    'SRR11271993': ('Donor1', 'AD'),
    'SRR11271994': ('Donor1', 'AD'),
    'SRR11271995': ('Donor2', 'Control'),
    # ... ADD ALL 157 SAMPLES HERE
}

def load_scte_matrix(sample_id):
    """Load scTE output CSV for a sample"""
    csv_file = SCTE_DIR / sample_id / f"{sample_id}.csv"
    
    if not csv_file.exists():
        print(f"Warning: {csv_file} not found")
        return None
    
    # Load CSV: rows = cells, columns = genes/TEs
    df = pd.read_csv(csv_file, index_col=0)
    print(f"Loaded {sample_id}: {df.shape[0]} cells × {df.shape[1]} features")
    
    return df

def aggregate_to_pseudobulk():
    """Aggregate all single-cell samples to pseudo-bulk by donor"""
    
    # Dictionary to accumulate counts per donor
    donor_counts = {}
    donor_conditions = {}
    
    print("=" * 80)
    print("Aggregating single-cell data to pseudo-bulk...")
    print("=" * 80)
    
    for sample_id, (donor_id, condition) in SAMPLE_METADATA.items():
        print(f"\nProcessing {sample_id} → {donor_id} ({condition})")
        
        # Load single-cell matrix for this sample
        sc_matrix = load_scte_matrix(sample_id)
        if sc_matrix is None:
            continue
        
        # Sum across all cells in this sample (axis=0 for columns)
        sample_totals = sc_matrix.sum(axis=0)
        
        # Accumulate into donor
        if donor_id not in donor_counts:
            donor_counts[donor_id] = sample_totals
            donor_conditions[donor_id] = condition
        else:
            donor_counts[donor_id] += sample_totals
    
    print("\n" + "=" * 80)
    print(f"Aggregated to {len(donor_counts)} donors")
    print("=" * 80)
    
    # Convert to DataFrame: rows = donors, columns = genes/TEs
    pseudobulk_df = pd.DataFrame(donor_counts).T
    
    # Add condition column
    pseudobulk_df['Condition'] = pseudobulk_df.index.map(donor_conditions)
    
    # Reorder: Condition first, then all features
    cols = ['Condition'] + [c for c in pseudobulk_df.columns if c != 'Condition']
    pseudobulk_df = pseudobulk_df[cols]
    
    print("\nPseudo-bulk matrix shape:", pseudobulk_df.shape)
    print("\nDonors by condition:")
    print(pseudobulk_df['Condition'].value_counts())
    
    return pseudobulk_df

def save_pseudobulk_matrix(pseudobulk_df):
    """Save pseudo-bulk matrix in DESeq2-compatible format"""
    
    # Save as CSV (same format as expression_matrix2.csv)
    output_file = OUTPUT_DIR / "pseudobulk_expression_matrix.csv"
    pseudobulk_df.to_csv(output_file)
    print(f"\n✓ Saved: {output_file}")
    
    # Create sample info file for DESeq2
    sample_info = pd.DataFrame({
        'Donor': pseudobulk_df.index,
        'Condition': pseudobulk_df['Condition']
    })
    sample_info_file = OUTPUT_DIR / "pseudobulk_sample_info.csv"
    sample_info.to_csv(sample_info_file, index=False)
    print(f"✓ Saved: {sample_info_file}")
    
    # Summary statistics
    summary_file = OUTPUT_DIR / "pseudobulk_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("PSEUDO-BULK AGGREGATION SUMMARY\n")
        f.write("=" * 80 + "\n\n")
        
        f.write(f"Total donors: {len(pseudobulk_df)}\n")
        f.write(f"Total features: {len(pseudobulk_df.columns) - 1}\n\n")
        
        f.write("Donors by condition:\n")
        for cond, count in pseudobulk_df['Condition'].value_counts().items():
            f.write(f"  {cond}: {count} donors\n")
        
        f.write("\nTotal counts per donor:\n")
        total_counts = pseudobulk_df.drop('Condition', axis=1).sum(axis=1)
        f.write(f"  Min: {total_counts.min():,.0f}\n")
        f.write(f"  Max: {total_counts.max():,.0f}\n")
        f.write(f"  Mean: {total_counts.mean():,.0f}\n")
        f.write(f"  Median: {total_counts.median():,.0f}\n")
    
    print(f"✓ Saved: {summary_file}")

def main():
    print("\nStarting pseudo-bulk aggregation...")
    
    # Aggregate single-cell to pseudo-bulk
    pseudobulk_df = aggregate_to_pseudobulk()
    
    # Save outputs
    save_pseudobulk_matrix(pseudobulk_df)
    
    print("\n" + "=" * 80)
    print("✓ Pseudo-bulk aggregation complete!")
    print("=" * 80)
    print("\nNext steps:")
    print("1. Review: /home/jacobc/hcaTE/pseudobulk/pseudobulk_summary.txt")
    print("2. Run DESeq2 on: /home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv")
    print("3. Use the existing diffexp_analysis.R script from HumanBam2scTE")

if __name__ == "__main__":
    main()
```

**Save as:** `/home/jacobc/hcaTE/scripts/aggregate_to_pseudobulk.py`

---

## Step 3: Run the Aggregation

```bash
cd /home/jacobc/hcaTE

# Make script executable
chmod +x scripts/aggregate_to_pseudobulk.py

# Run aggregation
python3 scripts/aggregate_to_pseudobulk.py
```

**Expected output:**
```
pseudobulk/
├── pseudobulk_expression_matrix.csv    # Main matrix for DESeq2
├── pseudobulk_sample_info.csv          # Sample metadata
└── pseudobulk_summary.txt               # QC summary
```

---

## Step 4: Run Differential Expression with DESeq2

### Option A: Use Existing Script (Recommended)

Copy and modify the DESeq2 script from HumanBam2scTE:

```bash
# Copy existing script
cp /home/jacobc/HumanBam2scTE/diffexp_analysis.R \
   /home/jacobc/hcaTE/scripts/pseudobulk_diffexp_analysis.R

# Edit to point to new pseudo-bulk matrix
# Change input file from expression_matrix2.csv to:
# /home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv
```

### Option B: Create New DESeq2 Script

```R
#!/usr/bin/env Rscript

library(DESeq2)
library(dplyr)
library(ggplot2)

# Load pseudo-bulk matrix
counts <- read.csv("/home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv", 
                   row.names = 1)

# Separate condition column
conditions <- counts$Condition
counts <- counts %>% select(-Condition)

# Transpose: DESeq2 expects genes as rows, samples as columns
counts <- t(counts)

# Create sample metadata
colData <- data.frame(
  Donor = colnames(counts),
  Condition = conditions,
  row.names = colnames(counts)
)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = colData,
  design = ~ Condition
)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("Condition", "AD", "Control"))

# Save results
write.csv(as.data.frame(res), 
          "/home/jacobc/hcaTE/pseudobulk/DESeq2_results.csv")

# Summary
summary(res, alpha = 0.05)
```

---

## Step 5: Quality Control Checks

### Check 1: Donor Distribution
```bash
# Should see balanced number of samples per donor
cat /home/jacobc/hcaTE/pseudobulk/pseudobulk_summary.txt
```

### Check 2: Library Sizes
```R
# In R, check if donors have similar total counts
counts <- read.csv("/home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv", 
                   row.names = 1)
library_sizes <- rowSums(counts[, -1])  # Exclude Condition column
summary(library_sizes)
```

**Expected:** Similar library sizes across donors (within ~2-3x range)

### Check 3: Number of Replicates
```bash
# DESeq2 requires ≥3 replicates per condition
# Check your design
```

---

## Common Issues and Solutions

### Issue 1: Unbalanced Replicates
**Problem:** 3 AD donors, 15 Control donors  
**Solution:** DESeq2 can handle unbalanced, but power may be limited for smaller group

### Issue 2: Low Read Counts After Aggregation
**Problem:** Some donors have very few total reads  
**Solution:** Filter out donors with <500,000 total counts

### Issue 3: Missing Sample Metadata
**Problem:** Don't know which samples belong to which donors  
**Solution:** Download metadata from GEO GSE146639

### Issue 4: Batch Effects
**Problem:** Donors sequenced in different batches  
**Solution:** Add batch to DESeq2 design: `design = ~ Batch + Condition`

---

## Validation

Before running DESeq2, validate your pseudo-bulk matrix:

```python
import pandas as pd

# Load pseudo-bulk
pb = pd.read_csv("/home/jacobc/hcaTE/pseudobulk/pseudobulk_expression_matrix.csv", 
                 index_col=0)

# Check dimensions
print(f"Donors: {len(pb)}")
print(f"Features: {len(pb.columns) - 1}")  # -1 for Condition column

# Check conditions
print("\nCondition distribution:")
print(pb['Condition'].value_counts())

# Check if counts are reasonable
total_counts = pb.drop('Condition', axis=1).sum(axis=1)
print(f"\nTotal counts per donor:")
print(f"  Min: {total_counts.min():,.0f}")
print(f"  Max: {total_counts.max():,.0f}")
print(f"  Mean: {total_counts.mean():,.0f}")

# Should be in millions (summing ~80 cells per sample, multiple samples per donor)
```

---

## Expected Results

### Pseudo-bulk Matrix Format:
```
             Condition    Gene1    Gene2    TE1      TE2      ...
Donor1       AD           12500    8900     450      120      ...
Donor2       AD           15200    9500     380      95       ...
Donor3       Control      13800    9200     420      110      ...
...
```

### File Sizes:
- `pseudobulk_expression_matrix.csv`: ~150-200 MB
- Similar to `expression_matrix2.csv` from HumanBam2scTE

---

## Next Steps After Aggregation

1. ✅ **Validate pseudo-bulk matrix** (check donor counts, conditions)
2. ✅ **Run DESeq2** on pseudo-bulk matrix
3. ✅ **Compare results** to original HumanBam2scTE (should be similar but correct)
4. ✅ **Analyze TE families** separately if needed
5. ✅ **Generate volcano plots** and heatmaps

---

## Key Differences from Original Analysis

| Aspect | HumanBam2scTE (WRONG) | hcaTE (CORRECT) |
|--------|----------------------|-----------------|
| **Input data** | Treated as bulk RNA-seq | Properly handled as single-cell |
| **Alignment** | Standard STAR | UMI-tools + STAR with CB/UB tags |
| **Quantification** | FeatureCounts | scTE with cell barcodes |
| **Replicates** | Unknown (treated samples as replicates?) | Donors = biological replicates |
| **Aggregation** | None (bulk approach) | Cell → Sample → Donor pseudo-bulk |
| **Statistical power** | Inflated (pseudo-replication) | Correct (donor-level) |

---

## Questions?

If you need help:
1. Check sample metadata in GEO GSE146639
2. Look at `sample_types.txt` from HumanBam2scTE
3. Verify donor assignments with publication methods

---

**Date Created:** November 29, 2025  
**GitHub Repository:** https://github.com/jcbclffrd/hcaTE
