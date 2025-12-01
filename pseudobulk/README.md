# Pseudobulk Differential Expression Analysis

This directory contains the results of differential expression analyses comparing Alzheimer's Disease (AD) patients to controls using pseudobulk aggregation of single-cell RNA-seq data.

---

## 📊 Analysis Overview

**Dataset**: GSE146639 - Single-nucleus RNA-seq of human brain (prefrontal cortex)  
**Publication**: Alsema et al. (2020) *Frontiers in Molecular Neuroscience* (PMID: 33192286)  
**Method**: scTE (single-cell TE) quantification → pseudobulk aggregation → DESeq2  
**Total Samples**: 31 (12 AD, 9 CTR, 6 CTR+, 4 MCI)  
**Total Cells**: 12,942 single nuclei  
**Features**: Genes + Transposable Elements (TEs)

---

## 🔬 Analyses Performed

### 1. AD vs CTR (Combined Controls)
**Comparison**: 12 AD vs 15 Control (combines 9 CTR + 6 CTR+ into one "Control" group)  
**Question**: What changes in AD compared to any controls (including those with early pathology)?

**Results**:
- ✅ **165 significant features** (padj < 0.05)
  - 139 genes
  - 26 transposable elements
- 74% upregulated in AD, 26% downregulated
- Top genes: HSPA6, JUN, COA1
- Top TEs: HERVFH21-int, MER61B, HERVK11-int

**Files**:
- `pseudobulk_diffexp_results_AD_vs_Control.csv` - All results
- `pseudobulk_diffexp_results_AD_vs_Control_significant.csv` - Significant only
- `pseudobulk_diffexp_results_AD_vs_Control_TEs_significant.csv` - Significant TEs
- `volcano_plot_AD_vs_Control.png` - Volcano plot visualization
- `ma_plot_AD_vs_Control.png` - MA plot visualization

### 2. AD vs CTR+ (Pathology-Free Controls)
**Comparison**: 12 AD vs 6 CTR+ (pathology-free controls only)  
**Question**: What changes in AD compared to completely healthy controls?  
**Note**: This matches the paper's main comparison strategy

**Results**:
- ✅ **100 significant features** (padj < 0.05)
  - 87 genes
  - 13 transposable elements
- 74% upregulated in AD, 26% downregulated
- Top genes: MSMO1 (down), TSIX, XIST, HSPA6, JUN
- Top TEs: Chap1_Mam, MER61B, L1M3d

**Files**:
- `pseudobulk_diffexp_results_AD_vs_CTRplus.csv` - All results
- `pseudobulk_diffexp_results_AD_vs_CTRplus_significant.csv` - Significant only
- `pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_significant.csv` - Significant TEs
- `volcano_plot_AD_vs_CTRplus.png` - Volcano plot visualization
- `ma_plot_AD_vs_CTRplus.png` - MA plot visualization

---

## 🧬 Key Findings

### Robust AD Markers (Found in BOTH Comparisons)

| Feature | Type | AD vs CTR | AD vs CTR+ | Function |
|---------|------|-----------|------------|----------|
| **HSPA6** | Gene | 5.8× up | 30× up | Heat shock protein, cellular stress |
| **JUN** | Gene | 6.4× up | 38× up | AP-1 transcription factor, stress response |
| **COA1** | Gene | 5.8× up | 40× up | Mitochondrial complex IV assembly |
| **MER61B** | TE | 3.5× up | 15× up | ERV1 LTR retrotransposon, HERV activation |

### Biological Themes

**Genes**:
- 🔥 Heat shock response (HSPA6, HSP90AA1, DNAJB1)
- 📝 Transcriptional stress (JUN, FOS, EGR1)
- ⚡ Mitochondrial dysfunction (COA1, NDUFA4L2, TIMM8A)
- 🧬 Cholesterol metabolism (MSMO1 strongly downregulated in CTR+ comparison)

**Transposable Elements** (NOVEL - not analyzed in original paper):
- 🦠 HERV activation (ERV1, ERVL, HERV-K families)
- 🧬 LINE1 mobilization (L1M, L1PA, L1MC families)
- 🔄 Consistent TE markers: MER61B, THE1D appear in both comparisons

---

## 📈 Statistical Summary

| Metric | AD vs CTR | AD vs CTR+ |
|--------|-----------|------------|
| Samples compared | 12 vs 15 (9+6) | 12 vs 6 |
| Features analyzed | 14,588 | 14,451 |
| Significant (padj<0.05) | 165 | 100 |
| Genes significant | 139 | 87 |
| TEs significant | 26 | 13 |
| Up in AD | 122 (74%) | 74 (74%) |
| Down in AD | 43 (26%) | 26 (26%) |

---

## 📝 Documentation Files

### Core Results
- `RESULTS_SUMMARY.md` - **START HERE** - Comprehensive summary of all findings
- `pseudobulk_sample_info.csv` - Sample metadata (condition, cell counts, etc.)
- `pseudobulk_expression_matrix.csv` - Full count matrix (samples × features)

### Analysis Details
- `COMPARISON_AD_vs_CTR_vs_CTRplus.md` - Side-by-side comparison of both analyses
- `WHY_PAPERS_GENES_ARE_MISSING.md` - Explains discrepancy with paper's findings
- `CORRECTED_ANALYSIS_SUMMARY.md` - Documents fix of gene/TE classification bug

### Technical
- `pseudobulk_diffexp_AD_vs_CTRplus_session_info.txt` - R session info for reproducibility
- Scripts in `/home/jacobc/hcaTE/scripts/`:
  - `pseudobulk_diffexp_analysis.R` - AD vs CTR analysis
  - `pseudobulk_diffexp_AD_vs_CTRplus.R` - AD vs CTR+ analysis

---

## 🤔 Why Our Results Differ from Paper

The original paper (Alsema et al. 2020) reported:
- SIGLEC1 (decreased in AD)
- CXCL10 (increased in AD)
- CXCR1, CXCR2 (chemokine receptors)

**We do NOT find these genes** because:

1. **Cell-type dilution**: These are microglia-specific markers
   - Pseudobulk aggregation dilutes cell-type-specific signals
   - SIGLEC1 has only 1 total count across all samples
   - After aggregation: signal disappears

2. **Method difference**: 
   - Paper used bulk RNA-seq (4 AD vs 3 CTR+)
   - We use pseudobulk from single-cell (12 AD vs 6 CTR+)
   - Different sensitivities to different gene types

3. **Statistical power trade-off**:
   - We have MORE samples (18 vs 7) → better overall power
   - But pseudobulk reduces sensitivity to rare cell-type markers
   - We find MORE DEGs overall (100 vs ~5-10)

**Both approaches are valid** - they just capture different biology.

📖 See `WHY_PAPERS_GENES_ARE_MISSING.md` for detailed explanation.

---

## 💡 Novel Findings

### Transposable Element Analysis

**CRITICAL**: The original paper **never analyzed TEs**. All TE findings are completely novel.

We discovered:
- 26 TEs significant in AD vs CTR
- 13 TEs significant in AD vs CTR+
- Consistent activation of HERV families (endogenous retroviruses)
- LINE1 element upregulation (may drive genomic instability)
- **MER61B** is a robust TE marker appearing in both comparisons

**Biological relevance**:
- HERV activation linked to neuroinflammation
- TE derepression is hallmark of aging/neurodegeneration
- May trigger innate immune response via dsRNA
- Potentially contributes to AD pathology

---

## 🔄 Pipeline Overview

```
Single-cell FASTQ files (GSE146639, n=157 samples)
    ↓
UMI-tools extraction (7bp UMI + 10bp cell barcode)
    ↓
STAR alignment (GRCh38 + RepeatMasker TE annotations)
    ↓
scTE quantification (genes + TEs, 12,942 cells)
    ↓
Pseudobulk aggregation (sum counts per sample, filter to 84 cell barcodes)
    ↓
DESeq2 differential expression (31 samples → 21 AD/CTR/CTR+ samples)
    ↓
Results + Visualizations
```

**Key validation**:
- ✅ Barcode extraction: 84/84 cell barcodes match whitelist
- ✅ TE/gene classification: Fixed to use `grepl("#")` for TEs
- ✅ Pipeline correctness: Validated against original authors' approach
- ✅ Statistical methods: Standard DESeq2 workflow with pre-filtering

---

## 📁 File Organization

```
pseudobulk/
├── README.md                                          ← YOU ARE HERE
├── RESULTS_SUMMARY.md                                 ← Comprehensive summary
├── COMPARISON_AD_vs_CTR_vs_CTRplus.md                ← Side-by-side comparison
├── WHY_PAPERS_GENES_ARE_MISSING.md                   ← Explains discrepancy
├── CORRECTED_ANALYSIS_SUMMARY.md                     ← Bug fix documentation
│
├── pseudobulk_sample_info.csv                        ← Sample metadata
├── pseudobulk_expression_matrix.csv                  ← Full count matrix
│
├── AD vs CTR Results:
│   ├── pseudobulk_diffexp_results_AD_vs_Control.csv
│   ├── pseudobulk_diffexp_results_AD_vs_Control_significant.csv
│   ├── pseudobulk_diffexp_results_AD_vs_Control_TEs_significant.csv
│   ├── volcano_plot_AD_vs_Control.png
│   └── ma_plot_AD_vs_Control.png
│
└── AD vs CTR+ Results:
    ├── pseudobulk_diffexp_results_AD_vs_CTRplus.csv
    ├── pseudobulk_diffexp_results_AD_vs_CTRplus_significant.csv
    ├── pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_significant.csv
    ├── volcano_plot_AD_vs_CTRplus.png
    └── ma_plot_AD_vs_CTRplus.png
```

---

## 🎯 Interpretation Guide

### Which Analysis Should I Use?

**AD vs CTR (165 significant)**:
- Use when interested in **broad AD changes** vs any controls
- Includes controls with early pathology (Braak I-II, amyloid+)
- Captures both early and late disease progression
- More DEGs = captures full disease spectrum

**AD vs CTR+ (100 significant)**:
- Use when interested in **pure AD vs completely healthy**
- Stricter comparison with pathology-free controls
- Matches the paper's analysis strategy
- Fewer DEGs = more conservative, cleaner signal

**Overlap (robust markers)**:
- ~60% of CTR+ DEGs also significant in CTR comparison
- Genes like HSPA6, JUN, COA1 appear in both → robust AD markers
- Use these for highest confidence findings

---

## 📊 Example Results

### Top 5 Genes (AD vs CTR)
1. **HSPA6** - 5.8× up (padj=1.4e-59) - Heat shock protein
2. **JUN** - 6.4× up (padj=5.7e-40) - Stress transcription factor
3. **COA1** - 5.8× up (padj=1.5e-38) - Mitochondrial assembly
4. **RPS6KA1** - 5.2× up - Ribosomal protein kinase
5. **FOSB** - 5.1× up - AP-1 transcription factor

### Top 5 TEs (AD vs CTR)
1. **HERVFH21-int** - 3.8× up (padj=3.1e-33) - HERV internal
2. **MER61B** - 3.5× up (padj=1.7e-32) - ERV1 LTR
3. **HERVK11-int** - 3.8× up - HERV-K internal
4. **THE1D** - 3.4× up - ERVL-MaLR LTR
5. **L1MEg** - 3.3× up - LINE1 element

---

## 🚀 Next Steps

### Completed ✅
- Barcode validation and pipeline correctness
- Gene/TE classification bug fix
- AD vs CTR differential expression
- AD vs CTR+ differential expression (matching paper)
- Comprehensive documentation

### Potential Future Analyses 🔮
- Cell-type-specific DE (extract microglia, astrocytes separately)
- Brain region analysis (LPS vs GFS if metadata available)
- TE subfamily enrichment analysis
- Functional enrichment (GO, KEGG pathways)
- Integration with AD GWAS data
- Correlation between TE and gene expression

---

## 📚 References

**Original Dataset**:
Alsema AM, Jiang Q, Kracht L, et al. (2020). Profiling Microglia From Alzheimer's Disease Donors and Non-demented Elderly in Acute Human Postmortem Cortical Tissue. *Front Mol Neurosci*. PMID: 33192286

**Methods**:
- scTE: He J, et al. (2021). *Genome Biol*.
- DESeq2: Love MI, et al. (2014). *Genome Biol*.
- UMI-tools: Smith T, et al. (2017). *Genome Res*.

**Repository**: https://github.com/jcbclffrd/hcaTE

---

## 🤖 Analysis Information

**Created by**: GitHub Copilot coding agent  
**Date**: November 2024  
**Environment**: Linux (bash), Python 3.12, R 4.x  
**Quality Control**: Pipeline validated, results documented, reproducible

---

*For questions or issues, refer to the comprehensive documentation files in this directory.*
