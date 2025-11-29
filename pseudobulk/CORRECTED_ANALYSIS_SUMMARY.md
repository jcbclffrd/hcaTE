# Corrected Differential Expression Analysis Summary

## Issue Identified and Fixed

### The Problem
The original DESeq2 analysis script incorrectly classified features as genes vs TEs:
```r
# INCORRECT (old version):
results_dt[, feature_type := ifelse(grepl("^ENS", feature), "Gene", "TE")]
```

This assumed gene names would start with "ENS" (Ensembl IDs), but scTE outputs **gene symbols** (e.g., HSPA6, JUN, COA1), not Ensembl IDs. As a result:
- **All genes were misclassified as TEs**
- The "155 significant TEs" were actually mostly genes
- True TE differential expression was masked

### The Fix
```r
# CORRECT (new version):
# TEs in scTE output have format: NAME#CLASS#FAMILY (e.g., MER61B#LTR#ERV1)
# Genes are gene symbols without "#" (e.g., HSPA6, JUN, A1BG)
results_dt[, feature_type := ifelse(grepl("#", feature), "TE", "Gene")]
```

## Corrected Results

### Overall Statistics
- **Total features analyzed:** 15,661 (after filtering low counts)
  - 13,648 genes
  - 2,013 TEs
- **Significant features (padj < 0.05):** 165
  - 139 significant genes (84.2%)
  - 26 significant TEs (15.8%)
- **Direction of change:**
  - 135 upregulated in AD (81.8%)
  - 30 downregulated in AD (18.2%)

### Significant Transposable Elements (n=26)

#### Upregulated in AD (n=9)
1. **HERVFH21-int#LTR#ERV1** - log2FC: 3.81, padj: 8.3e-04
2. **LTR47A#LTR#ERVL** - log2FC: 3.71, padj: 0.038
3. **HERVK11-int#LTR#ERVK** - log2FC: 3.80, padj: 0.034
4. **MER61B#LTR#ERV1** - log2FC: 3.49, padj: 8.1e-04
5. **(CTAT)n#Simple_repeat** - log2FC: 3.39, padj: 0.0016
6. **MER11A#LTR#ERVK** - log2FC: 2.85, padj: 0.0073
7. **LTR41#LTR#ERVL** - log2FC: 2.01, padj: 0.040
8. **MER57A1#LTR#ERV1** - log2FC: 1.79, padj: 0.039
9. **L4_B_Mam#LINE#RTE-X** - log2FC: 1.14, padj: 0.027

#### Downregulated in AD (n=17)
1. **(AATGGAATGG)n#Simple_repeat** - log2FC: -4.77, padj: 0.015
2. **(GGAAT)n#Simple_repeat** - log2FC: -3.50, padj: 0.046
3. **MER4A1#LTR#ERV1** - log2FC: -3.16, padj: 0.0037
4. **(GAATG)n#Satellite** - log2FC: -2.71, padj: 0.050
5. **(ATGGA)n#Simple_repeat** - log2FC: -2.45, padj: 0.034
6. **(ATTCC)n#Simple_repeat** - log2FC: -2.32, padj: 0.013
7. **(AACA)n#Simple_repeat** - log2FC: -2.19, padj: 0.015
8. **(CATTC)n#Satellite** - log2FC: -2.03, padj: 0.046
9. **MLT1G3#LTR#ERVL-MaLR** - log2FC: -1.73, padj: 0.027
10. **ERVL-int#LTR#ERVL** - log2FC: -1.69, padj: 0.016
11. **(TTTTAAT)n#Simple_repeat** - log2FC: -1.39, padj: 0.029
12. **(TGAA)n#Simple_repeat** - log2FC: -1.19, padj: 0.050
13. **(TTTTA)n#Simple_repeat** - log2FC: -1.13, padj: 0.035
14. **(A)n#Simple_repeat** - log2FC: -0.80, padj: 0.016
15. **(T)n#Simple_repeat** - log2FC: -0.84, padj: 0.0058
16. **L1MA9#LINE#L1** - log2FC: -0.73, padj: 0.050
17. **AluSc8#SINE#Alu** - log2FC: -0.69, padj: 0.034

### TE Class Distribution in Significant Results

**LTR retrotransposons:** 9 TEs
- ERV1: 4 (MER61B, MER4A1, MER57A1, HERVFH21-int)
- ERVL: 4 (LTR47A, LTR41, ERVL-int, MLT1G3)
- ERVK: 2 (MER11A, HERVK11-int)

**Simple repeats:** 11 TEs
- Various microsatellites and low-complexity sequences
- Mostly downregulated in AD

**LINE elements:** 2 TEs
- L4_B_Mam (RTE-X family) - upregulated
- L1MA9 (L1 family) - downregulated

**SINE elements:** 1 TE
- AluSc8 (Alu family) - downregulated

**Satellites:** 2 TEs
- Both downregulated in AD

## Top Significant Genes (Not TEs!)

The following genes showed the strongest differential expression (these were previously misclassified as TEs):

### Top 5 Upregulated Genes in AD:
1. **JUN** - log2FC: 6.40, padj: 8.3e-06 (transcription factor, stress response)
2. **ENSG00000289750** - log2FC: 6.04, padj: 1.7e-04 (novel gene)
3. **HSPA6** - log2FC: 5.81, padj: 4.7e-06 (heat shock protein)
4. **COA1** - log2FC: 5.80, padj: 1.5e-05 (cytochrome c oxidase assembly)
5. **PCYOX1** - log2FC: 5.61, padj: 1.1e-03 (prenylcysteine oxidase)

### Top 5 Downregulated Genes in AD:
1. **MIR181A1HG** - log2FC: -5.18, padj: 0.016 (microRNA host gene)
2. **SYNE2** - log2FC: -3.64, padj: 0.016 (nuclear envelope protein)
3. **ZDHHC14** - log2FC: -2.80, padj: 0.014 (palmitoyl transferase)
4. **SLC6A1** - log2FC: -1.65, padj: 4.1e-04 (GABA transporter)

## Biological Interpretation

### Genes:
- Strong upregulation of stress response genes (HSPA6, JUN)
- Mitochondrial dysfunction markers (COA1)
- Changes in neurotransmitter transport (SLC6A1)

### TEs:
- **HERV activation:** Multiple HERV/LTR elements upregulated (HERVFH21, HERVK11)
- **Simple repeats:** Predominantly downregulated, possibly reflecting global chromatin changes
- **LINE/SINE elements:** Mixed patterns (L4_B upregulated, L1MA9/Alu downregulated)

## Files Generated

- `pseudobulk_diffexp_results_AD_vs_Control.csv` - All 15,661 features
- `pseudobulk_diffexp_results_TEs_significant.csv` - 26 significant TEs (CORRECTED)
- `pseudobulk_diffexp_results_AD_vs_Control_significant.csv` - All 165 significant features
- `volcano_plot_pseudobulk_AD_vs_Control.png` - Updated with correct classification
- `volcano_plot_TEs_only.png` - TE-specific plot with correct TEs
- `ma_plot_pseudobulk_AD_vs_Control.png` - Mean-expression plot

## Comparison: Old vs New

| Metric | Incorrect Analysis | Corrected Analysis |
|--------|-------------------|-------------------|
| "Significant TEs" | 155 | 26 |
| Actual significant genes | 0 (misclassified) | 139 |
| True TE upregulated | Unknown | 9 |
| True TE downregulated | Unknown | 17 |
| Top hit | "HSPA6 (TE)" | "JUN (Gene)" |

## Conclusion

The corrected analysis reveals:
1. **26 true TEs** show significant differential expression in AD vs Control
2. **139 genes** (not TEs) are significantly differentially expressed
3. LTR/HERV elements show upregulation in AD microglia
4. Simple repeats and some LINE/SINE elements are downregulated
5. The original finding of HSPA6 and JUN upregulation is **correct**, but these are **genes, not TEs**

This correction is critical for proper biological interpretation and downstream validation experiments.
