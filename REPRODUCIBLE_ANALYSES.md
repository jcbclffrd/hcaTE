# Reproducible Analyses from Alsema et al. 2020

## Paper: "Profiling Microglia From Alzheimer's Disease Donors and Non-demented Elderly"
- **PMID**: 33192286
- **DOI**: 10.3389/fnmol.2020.00134
- **GEO**: GSE146639

---

## What We've Already Done ✓

### 1. **Differential Expression Analysis (Pseudobulk)**
- **Status**: COMPLETED
- **Results**: 
  - AD vs Control: **165 significant features** (padj < 0.05)
  - AD vs CTR+: **100 significant features** (padj < 0.05)  
  - Transposable Elements: **26 significant TEs** (padj < 0.05)
- **Files**: 
  - `pseudobulk/pseudobulk_diffexp_results_AD_vs_Control_significant.csv`
  - `pseudobulk/pseudobulk_diffexp_results_AD_vs_CTRplus_significant.csv`
  - `pseudobulk/pseudobulk_diffexp_results_TEs_significant.csv`

### Key Findings So Far:
- **Top TE upregulated in AD**: MER61B (LTR/ERV1, log2FC = 3.49, padj < 0.001)
- **Top TE downregulated in AD**: MER4A1 (LTR/ERV1, log2FC = -3.16, padj < 0.01)
- Most significant TEs are LTR retrotransposons and simple repeats

---

## Potential Reproducible Analyses

### 🎯 **Analysis 1: Gene Expression Comparison (Paper's Main Finding)**

**Paper's Finding**: 
> "Gene expression profiles and subcluster composition of microglia did not differ between AD donors and non-demented elderly in bulk RNA sequencing nor in single-cell sequencing."

**What We Can Do**:
1. ✅ Already completed differential expression (pseudobulk)
2. 🔄 **Compare our results to paper's claim**:
   - Check if microglia markers are differentially expressed
   - Look at specific genes mentioned in paper (CD11B, CD45, CX3CR1, etc.)
   - Validate that most genes show no significant changes

**Script to Create**:
```python
# Compare microglia marker genes between conditions
# - Extract specific marker genes from our DE results
# - Create volcano plot showing lack of major DE
# - Statistical summary matching paper's approach
```

**Difficulty**: ⭐ Easy - we already have the data

---

### 🎯 **Analysis 2: Transposable Element Expression in AD**

**Paper's Context**: 
The paper focuses on microglia but doesn't specifically highlight TEs. However, TE dysregulation in neurodegeneration is a hot topic.

**What We Found**:
- **26 significant TEs** in AD vs Control
- ERV families (endogenous retroviruses) are prominent
- Simple repeats also show differential expression

**What We Can Do**:
1. ✅ Already identified significant TEs
2. 🆕 **TE family enrichment analysis**:
   - Group TEs by family (LTR, LINE, SINE, etc.)
   - Test if specific TE families are enriched in AD
   - Create summary plots by TE class

3. 🆕 **TE expression correlation with known AD genes**:
   - Correlate TE expression with APOE, TREM2, etc.
   - Look for co-expression patterns

**Script to Create**:
```python
# TE family enrichment and visualization
# - Categorize TEs by class/family
# - Fisher's exact test for enrichment
# - Heatmap of top TEs across conditions
# - Correlation analysis with AD genes
```

**Difficulty**: ⭐⭐ Moderate - requires TE annotation parsing

---

### 🎯 **Analysis 3: Regional Differences (LPS vs GFS)**

**Paper's Design**:
- 13 samples from **LPS** (Lateral/Superior Parietal Lobe)
- 12 samples from **GFS** (Superior Frontal Gyrus)

**Our Challenge**:
- ❌ Need SRR → Brain Region mapping (see `regional_analysis/README.md`)
- Mapping blocked by lack of SRA RunInfo table

**What We Can Do**:
1. 🔄 **Download proper sample mappings**:
   - Get SRA RunInfo for SRP252065
   - Create SRR → GSM → Brain Region mapping
   
2. 🆕 **Compare LPS vs GFS expression**:
   - Differential expression between brain regions
   - Check if TEs differ by region
   - Subset by condition (AD vs CTR within each region)

**Expected Result**: Paper suggests minimal regional differences in microglia

**Difficulty**: ⭐⭐⭐ Hard - requires external data download

---

### 🎯 **Analysis 4: Single-Cell Heterogeneity (Most Exciting!)**

**Paper's Finding**:
> "We identified 7 human microglial subpopulations with heterogeneity in gene expression."

**Our Data**:
- We have **scTE output** for 160 samples
- Each sample has ~80-100 single cells
- Total: ~12,942 cells

**What We Can Do**:
1. 🆕 **Single-cell clustering analysis**:
   - Load scTE single-cell matrices
   - Perform dimensionality reduction (PCA, UMAP)
   - Cluster cells (Seurat or Scanpy)
   - Identify microglial subpopulations
   
2. 🆕 **TE expression at single-cell level**:
   - Which TEs are cell-type specific?
   - Are certain TE families expressed in specific microglia subtypes?
   - Do AD samples have different subpopulation proportions?

3. 🆕 **Validate paper's 7 subpopulations**:
   - Can we identify similar clusters?
   - Check marker genes for each subtype
   - Compare AD vs Control subpopulation composition

**Script Requirements**:
```python
# Single-cell analysis workflow
# 1. Load all scTE outputs (160 samples × ~80 cells each)
# 2. Create combined cell × gene matrix
# 3. Normalize and filter
# 4. PCA + UMAP
# 5. Clustering (Leiden/Louvain)
# 6. Marker gene identification per cluster
# 7. TE expression per cluster
# 8. Proportion analysis (AD vs CTR)
```

**Difficulty**: ⭐⭐⭐⭐ Very Hard - requires single-cell analysis framework

---

### 🎯 **Analysis 5: Specific Gene/TE Validation**

**What We Can Do**:
1. 🆕 **Reproduce specific gene expression plots**:
   - Check expression of key genes mentioned in paper
   - Create condition-specific boxplots
   - Statistical tests matching paper's methods

2. 🆕 **Focus on microglia markers**:
   - CX3CR1, P2RY12, TMEM119, CD11B (ITGAM), CD45 (PTPRC)
   - TREM2, APOE (AD risk genes)
   - Inflammatory markers (IL1B, TNF, etc.)

3. 🆕 **TE-specific plots**:
   - Show our top TEs (MER61B, HERVFH21, MER4A1)
   - Compare to known TE expression in neurodegeneration literature

**Script to Create**:
```python
# Gene/TE validation plots
# - Extract specific genes from count matrix
# - Create publication-quality boxplots
# - Add statistical annotations (p-values)
# - Compare to paper's figures
```

**Difficulty**: ⭐ Easy - straightforward plotting

---

## Recommended Next Steps (Ranked by Value)

### 🥇 **HIGHEST PRIORITY: Analysis 4 (Single-Cell Clustering)**
**Why**: 
- This is the paper's main contribution
- We have all the data needed
- Would reveal TE expression heterogeneity at single-cell level
- Most novel and publishable

**Requirements**:
- Python packages: scanpy, anndata, umap-learn
- Computational: ~16GB RAM, 1-2 hours runtime
- Files: `scTE_output/SRR*/SRR*.csv` (already have)

---

### 🥈 **MEDIUM PRIORITY: Analysis 2 (TE Family Enrichment)**
**Why**:
- Builds on existing pseudobulk results
- Quick to implement
- Provides mechanistic insights
- Could lead to novel findings

**Requirements**:
- Parse TE annotations from RepeatMasker
- Group by TE class/family
- Statistical testing (Fisher's exact)

---

### 🥉 **LOWER PRIORITY: Analysis 5 (Gene Validation)**
**Why**:
- Easy to do
- Validates our pipeline
- Less novel but good sanity check

**Requirements**:
- Extract key genes from pseudobulk matrix
- Create plots with ggplot2 or matplotlib

---

## Data We Already Have

✅ **scTE outputs**: 160 samples × ~80 cells = ~12,942 single cells  
✅ **Pseudobulk matrices**: Aggregated gene + TE counts  
✅ **Sample metadata**: Condition (AD, CTR, CTR+, MCI)  
✅ **Differential expression results**: 165 significant features (AD vs Control)  
✅ **TE annotations**: RepeatMasker bed file  
✅ **Gene annotations**: GENCODE v45 GTF  

## Data We Still Need

❌ **Brain region mapping**: SRR → LPS/GFS (for Analysis 3)  
❌ **Cell-type annotations**: Which cells are microglia? (for Analysis 4)  
❌ **Paper's supplementary tables**: For direct comparison  

---

## Conclusion

The **most impactful** analysis would be **#4: Single-Cell Clustering**, as it would:
1. Validate the paper's finding of 7 microglial subpopulations
2. Reveal TE expression patterns at single-cell resolution
3. Show which subtypes are altered in AD
4. Potentially discover novel TE-expressing cell states

This would be a publishable extension of the original paper focusing on transposable elements.

**Estimated Time**:
- Analysis 5 (Gene validation): 1-2 hours
- Analysis 2 (TE enrichment): 3-4 hours  
- Analysis 4 (Single-cell): 1-2 days (includes learning scanpy/Seurat)
- Analysis 3 (Regional): Blocked until we get mapping data

**Which would you like to pursue?**
