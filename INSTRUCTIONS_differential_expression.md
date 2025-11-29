# Instructions: Differential Expression Analysis on Single-Cell Data

## Context

You've successfully completed the scTE pipeline and have **single-cell count matrices** with:
- **Rows**: Cell barcodes (individual cells)
- **Columns**: Features (TEs and genes)
- **Data**: Very sparse (typical for single-cell data)

## CRITICAL: Single-Cell vs Bulk Analysis

The differential expression script in `/home/jacobc/HumanBam2scTE/scripts/diffexp_analysis.R` was designed for **BULK RNA-seq data**, where:
- Rows = samples (donors)
- Columns = features (TEs/genes)

**For single-cell data, you need to adapt the approach!**

## Two Approaches for Single-Cell Differential Expression

### Approach 1: Pseudo-bulk Aggregation (RECOMMENDED)

**Aggregate cells by donor** to create pseudo-bulk samples, then use the existing DESeq2 script.

#### Steps:

1. **Aggregate cells by donor/sample**
   ```R
   # For each donor, sum counts across all cells from that donor
   # This creates one "pseudo-bulk" sample per donor
   
   library(data.table)
   library(Matrix)
   
   # Load scTE output (sparse matrix format)
   # Assuming you have: barcodes × features matrix
   
   # Extract donor ID from cell barcodes (format: SRR*_BARCODE)
   cell_metadata <- data.table(
     barcode = rownames(count_matrix),
     sample_id = sub("_.*", "", rownames(count_matrix))  # Extract SRR ID
   )
   
   # Aggregate: sum counts for each donor
   pseudo_bulk <- sapply(unique(cell_metadata$sample_id), function(donor) {
     cells_from_donor <- cell_metadata[sample_id == donor, barcode]
     colSums(count_matrix[cells_from_donor, , drop = FALSE])
   })
   
   # Now pseudo_bulk is: donors × features (compatible with DESeq2 script!)
   ```

2. **Save pseudo-bulk matrix**
   ```R
   # Convert to data.table with sample column
   pseudo_bulk_dt <- as.data.table(t(pseudo_bulk))
   pseudo_bulk_dt[, sample := rownames(pseudo_bulk)]
   
   # Save in same format as expression_matrix2.csv
   fwrite(pseudo_bulk_dt, "expression_matrix_pseudobulk.csv")
   ```

3. **Run the existing DESeq2 script**
   ```bash
   cd /home/jacobc/hcaTE
   
   # Copy the DESeq2 script
   cp /home/jacobc/HumanBam2scTE/scripts/diffexp_analysis.R ./scripts/
   
   # Create metadata file (if not already done)
   # metadata.csv needs: sample_id, donor_group (AD/CTR/etc), brain_region
   
   # Run analysis on pseudo-bulk data
   Rscript scripts/diffexp_analysis.R
   ```

**Advantages of pseudo-bulk:**
- Uses established DESeq2 workflow
- Accounts for donor-level variation
- Reduces sparsity by aggregating cells
- Matches biological replicates (donors, not cells)
- The existing script will work with minimal modification

### Approach 2: True Single-Cell Analysis (MORE COMPLEX)

Use single-cell-specific tools like **Seurat**, **MAST**, or **scVI** that handle sparsity and cell-level variation.

**This requires different R packages and a completely different script!**

## Recommended Workflow: Pseudo-bulk + DESeq2

### Step 1: Create Pseudo-bulk Aggregation Script

Create `/home/jacobc/hcaTE/scripts/aggregate_to_pseudobulk.R`:

```R
#!/usr/bin/env Rscript
# Aggregate single-cell counts to pseudo-bulk by donor

library(data.table)
library(Matrix)

cat("Creating pseudo-bulk aggregation from single-cell data...\n\n")

# ===== Load scTE output =====
cat("Loading scTE count matrix...\n")

# TODO: Adjust path to your scTE output
# scTE outputs are typically in scte_results/
# Format might be: combined_matrix.csv or individual sample matrices

# If you have individual sample matrices:
scte_dir <- "scte_results"
sample_files <- list.files(scte_dir, pattern = "*.csv", full.names = TRUE)

# Load and combine all samples
all_cells_list <- lapply(sample_files, function(f) {
  dt <- fread(f)
  # Assuming first column is cell barcode, rest are features
  # Add sample prefix to barcodes
  sample_id <- sub("\\.csv$", "", basename(f))
  dt[, cell_barcode := paste0(sample_id, "_", V1)]
  dt[, V1 := NULL]
  return(dt)
})

# Combine all cells
combined_matrix <- rbindlist(all_cells_list, fill = TRUE)

cat("  Total cells:", nrow(combined_matrix), "\n")
cat("  Total features:", ncol(combined_matrix) - 1, "\n\n")

# ===== Aggregate by donor =====
cat("Aggregating cells by donor (pseudo-bulk)...\n")

# Extract sample/donor ID from cell barcode
combined_matrix[, sample_id := sub("_.*", "", cell_barcode)]

# Aggregate: sum counts for each donor
pseudo_bulk <- combined_matrix[, lapply(.SD, sum), by = sample_id, .SDcols = -"cell_barcode"]

cat("  Pseudo-bulk samples:", nrow(pseudo_bulk), "\n\n")

# ===== Format for DESeq2 script =====
# Rename sample_id column to "sample" to match expected format
setnames(pseudo_bulk, "sample_id", "sample")

# Save
fwrite(pseudo_bulk, "expression_matrix_pseudobulk.csv")
cat("Pseudo-bulk matrix saved to: expression_matrix_pseudobulk.csv\n")
cat("  Dimensions:", nrow(pseudo_bulk), "samples x", ncol(pseudo_bulk)-1, "features\n\n")

cat("Next step: Run DESeq2 differential expression analysis\n")
cat("  Rscript scripts/diffexp_analysis.R\n")
```

### Step 2: Prepare Metadata

Create or copy metadata file with donor information:

```bash
cd /home/jacobc/hcaTE

# If metadata exists in old project:
cp /home/jacobc/HumanBam2scTE/metadata.csv ./

# Or create new metadata.csv with columns:
# sample_id, donor_group, brain_region
# Example:
# sample_id,donor_group,brain_region
# SRR11272119,AD,hippocampus
# SRR11272120,CTR,hippocampus
# ...
```

### Step 3: Adapt DESeq2 Script

Copy and modify the script to use pseudo-bulk data:

```bash
cd /home/jacobc/hcaTE
mkdir -p scripts

# Copy the DESeq2 script
cp /home/jacobc/HumanBam2scTE/scripts/diffexp_analysis.R scripts/

# Modify line 27 to use new file:
# Change: expr_matrix <- fread("expression_matrix2.csv")
# To:     expr_matrix <- fread("expression_matrix_pseudobulk.csv")
```

### Step 4: Run the Analysis

```bash
cd /home/jacobc/hcaTE

# Step 1: Aggregate to pseudo-bulk
Rscript scripts/aggregate_to_pseudobulk.R

# Step 2: Run DESeq2 differential expression
Rscript scripts/diffexp_analysis.R
```

## Expected Outputs

After running, you should get:
- `diffexp_results_AD_vs_nonAD.csv` - Full differential expression results
- `diffexp_results_AD_vs_nonAD_significant.csv` - Significant TEs (padj < 0.05)
- `volcano_plot_AD_vs_nonAD.png` - Volcano plot visualization
- `ma_plot_AD_vs_nonAD.png` - MA plot visualization
- `diffexp_session_info.txt` - R session info

## Quality Control Checks

Before running DESeq2, verify:

1. **Cell counts per donor**: Are there enough cells per sample?
   ```R
   table(combined_matrix$sample_id)
   # Should have >100 cells per donor ideally
   ```

2. **Sparsity after aggregation**: Pseudo-bulk should be less sparse
   ```R
   # Before aggregation (single-cell)
   sum(combined_matrix[, -"cell_barcode"] == 0) / prod(dim(combined_matrix[, -"cell_barcode"]))
   
   # After aggregation (pseudo-bulk)
   sum(pseudo_bulk[, -"sample"] == 0) / prod(dim(pseudo_bulk[, -"sample"]))
   ```

3. **Library sizes**: Check total counts per pseudo-bulk sample
   ```R
   pseudo_bulk[, total_counts := rowSums(.SD), .SDcols = -"sample"]
   summary(pseudo_bulk$total_counts)
   # Should be in millions for good pseudo-bulk
   ```

## Troubleshooting

**Issue: Not enough cells per donor**
- Consider filtering out low-quality cells first
- May need to group by donor + condition instead of individual donors

**Issue: Still too sparse after aggregation**
- Filter low-count features before DESeq2
- The script already does this (line ~118): `keep <- rowSums(counts(dds) >= 10) >= 3`

**Issue: Batch effects between samples**
- Consider adding sequencing batch as covariate in design formula
- May need to run batch correction (ComBat-seq or similar)

## Alternative: True Single-Cell Analysis

If pseudo-bulk isn't appropriate for your question, consider:
- **Seurat + FindMarkers**: Popular single-cell DE tool
- **MAST**: Accounts for bimodality in single-cell data
- **scVI**: Deep learning approach for single-cell
- **Wilcoxon rank-sum test**: Non-parametric, good for sparse data

Let me know if you need help setting up true single-cell DE analysis!

---

## Summary

**Recommended approach:**
1. Run `aggregate_to_pseudobulk.R` to create pseudo-bulk samples
2. Use existing `diffexp_analysis.R` script on aggregated data
3. Interpret results at donor level (not cell level)

This gives you statistical power by treating donors as biological replicates (correct!) rather than cells (pseudoreplication).
