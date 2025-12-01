#!/usr/bin/env Rscript
# Differential Expression Analysis: AD vs CTR+ (to match paper's analysis)
# Using DESeq2 for transposable element expression analysis
# 🤖 Built by GitHub Copilot coding agent

# Install/load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE))
    BiocManager::install("DESeq2")

library(DESeq2)
library(data.table)
library(ggplot2)

cat("===== Differential Expression Analysis: AD vs CTR+ =====\n")
cat("(Matching the paper's main comparison)\n\n")

# Set working directory
setwd("/home/jacobc/hcaTE/pseudobulk")

# ========== 1. Load Pseudo-bulk Expression Matrix ==========
cat("Loading pseudo-bulk expression matrix...\n")
expr_matrix <- fread("pseudobulk_expression_matrix.csv")
cat("  Matrix dimensions:", nrow(expr_matrix), "samples x", ncol(expr_matrix)-1, "features\n")

# ========== 2. Extract Sample Info ==========
cat("\nExtracting sample information...\n")

# First column is sample IDs, Condition column is the group
sample_ids <- expr_matrix$V1
conditions <- expr_matrix$Condition

# Remove Condition column and sample ID column to get pure counts
count_data <- as.matrix(expr_matrix[, -c("V1", "Condition")])
rownames(count_data) <- sample_ids

# Transpose: DESeq2 needs features as rows, samples as columns
count_data <- t(count_data)

cat("  Count matrix dimensions:", nrow(count_data), "features x", ncol(count_data), "samples\n")

# ========== 3. Create Sample Metadata ==========
cat("\nCreating sample metadata...\n")

coldata <- data.frame(
  sample_id = sample_ids,
  condition = conditions,
  row.names = sample_ids
)

cat("  Sample distribution by condition:\n")
print(table(coldata$condition))
cat("\n")

# ========== 4. Define Comparison Groups: AD vs CTR+ ONLY ==========
cat("Setting up comparison: AD vs CTR+ (matching paper)\n")

# Filter to ONLY AD and CTR+ samples (exclude CTR and MCI)
coldata_filtered <- coldata[coldata$condition %in% c("AD", "CTR+"), ]
count_data_filtered <- count_data[, rownames(coldata_filtered)]

cat("  Final sample sizes:\n")
cat("    AD:", sum(coldata_filtered$condition == "AD"), "samples\n")
cat("    CTR+:", sum(coldata_filtered$condition == "CTR+"), "samples\n")
cat("    (Excluded CTR:", sum(coldata$condition == "CTR"), "samples)\n")
cat("    (Excluded MCI:", sum(coldata$condition == "MCI"), "samples)\n\n")

# Set reference level (CTR+ is reference, AD is treatment)
coldata_filtered$condition <- factor(coldata_filtered$condition, 
                                     levels = c("CTR+", "AD"))

# ========== 5. Create DESeq2 Dataset ==========
cat("Creating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = coldata_filtered,
  design = ~ condition
)

cat("  DESeq2 dataset created\n\n")

# ========== 6. Pre-filtering: Remove Low-Count Features ==========
cat("Pre-filtering low-count features...\n")

keep <- rowSums(counts(dds)) >= 10
dds_filtered <- dds[keep, ]

cat("  Original features:", format(nrow(dds), big.mark=","), "\n")
cat("  After filtering:", format(nrow(dds_filtered), big.mark=","), "\n")
cat("  Removed:", format(nrow(dds) - nrow(dds_filtered), big.mark=","), "low-count features\n\n")

# ========== 7. Run DESeq2 Differential Expression Analysis ==========
cat("Running DESeq2 differential expression analysis...\n")
cat("  This may take several minutes...\n\n")

dds_result <- DESeq(dds_filtered)

cat("  DESeq2 analysis complete!\n\n")

# ========== 8. Extract Results ==========
cat("Extracting differential expression results...\n")

results <- results(dds_result, 
                   contrast = c("condition", "AD", "CTR+"),
                   alpha = 0.05)

cat("  Results summary:\n")
summary(results)
cat("\n")

# ========== 9. Convert to Data Table and Add Metadata ==========
cat("Processing results...\n")

results_dt <- as.data.table(results, keep.rownames = "feature")

# CRITICAL: Correct TE/Gene classification
# TEs in scTE have format: NAME#CLASS#FAMILY (e.g., MER61B#LTR#ERV1)
# Genes are gene symbols without "#" (e.g., HSPA6, JUN, A1BG)
results_dt[, feature_type := ifelse(grepl("#", feature), "TE", "Gene")]

# Add significance categories
results_dt[, significance := ifelse(padj < 0.001, "padj < 0.001",
                                   ifelse(padj < 0.01, "padj < 0.01",
                                         ifelse(padj < 0.05, "padj < 0.05", "Not significant")))]

# Add regulation direction
results_dt[, regulation := ifelse(is.na(padj), "Not significant",
                                 ifelse(padj < 0.05 & log2FoldChange > 0, "Up in AD",
                                       ifelse(padj < 0.05 & log2FoldChange < 0, "Down in AD", 
                                             "Not significant")))]

# Add absolute log2FC for easier filtering
results_dt[, abs_log2FC := abs(log2FoldChange)]

# ========== 10. Summary Statistics ==========
cat("\n===== SUMMARY STATISTICS =====\n")
cat("Total features analyzed:", nrow(results_dt), "\n")
cat("  Genes:", sum(results_dt$feature_type == "Gene"), "\n")
cat("  TEs:", sum(results_dt$feature_type == "TE"), "\n\n")

cat("Significant features (padj < 0.05):", 
    sum(results_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("  Genes:", sum(results_dt$feature_type == "Gene" & results_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("  TEs:", sum(results_dt$feature_type == "TE" & results_dt$padj < 0.05, na.rm = TRUE), "\n\n")

cat("Regulation direction (significant only):\n")
cat("  Up in AD:", sum(results_dt$regulation == "Up in AD", na.rm = TRUE), "\n")
cat("  Down in AD:", sum(results_dt$regulation == "Down in AD", na.rm = TRUE), "\n\n")

# ========== 11. Save Results ==========
cat("Saving results...\n")

# All results
fwrite(results_dt, "pseudobulk_diffexp_results_AD_vs_CTRplus.csv")
cat("  All results: pseudobulk_diffexp_results_AD_vs_CTRplus.csv\n")

# Significant results only
results_sig <- results_dt[padj < 0.05 & !is.na(padj)]
results_sig <- results_sig[order(padj)]
fwrite(results_sig, "pseudobulk_diffexp_results_AD_vs_CTRplus_significant.csv")
cat("  Significant: pseudobulk_diffexp_results_AD_vs_CTRplus_significant.csv\n")

# TEs only (all)
results_te <- results_dt[feature_type == "TE"]
results_te <- results_te[order(padj)]
fwrite(results_te, "pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_only.csv")
cat("  TEs only: pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_only.csv\n")

# Significant TEs only
results_te_sig <- results_dt[feature_type == "TE" & padj < 0.05 & !is.na(padj)]
results_te_sig <- results_te_sig[order(padj)]
fwrite(results_te_sig, "pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_significant.csv")
cat("  Significant TEs: pseudobulk_diffexp_results_AD_vs_CTRplus_TEs_significant.csv\n\n")

# ========== 12. Generate Volcano Plot ==========
cat("Generating volcano plot...\n")

# Prepare data for plotting
plot_data <- results_dt[!is.na(padj)]
plot_data[, neg_log10_padj := -log10(padj)]

# Define significance colors
plot_data[, plot_color := ifelse(significance == "padj < 0.001", "High (p < 0.001)",
                                 ifelse(significance == "padj < 0.01", "Medium (p < 0.01)",
                                       ifelse(significance == "padj < 0.05", "Low (p < 0.05)", 
                                             "Not significant")))]

plot_data$plot_color <- factor(plot_data$plot_color,
                               levels = c("Not significant", "Low (p < 0.05)", 
                                        "Medium (p < 0.01)", "High (p < 0.001)"))

# Create volcano plot
p1 <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_padj, color = plot_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not significant" = "gray70",
                                "Low (p < 0.05)" = "skyblue",
                                "Medium (p < 0.01)" = "orange",
                                "High (p < 0.001)" = "red3")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  labs(title = "Volcano Plot: AD vs CTR+ (Matching Paper)",
       subtitle = paste0("AD (n=", sum(coldata_filtered$condition == "AD"), ") vs CTR+ (n=", 
                        sum(coldata_filtered$condition == "CTR+"), ") | ",
                        sum(results_dt$padj < 0.05, na.rm = TRUE), 
                        " significant features (", 
                        sum(results_te_sig$padj < 0.05, na.rm = TRUE),
                        " TEs)"),
       x = "Log2 Fold Change (AD vs CTR+)",
       y = "-Log10(Adjusted P-value)",
       color = "Significance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        legend.position = "right")

ggsave("volcano_plot_AD_vs_CTRplus.png", p1, width = 10, height = 7, dpi = 300)
cat("  Saved: volcano_plot_AD_vs_CTRplus.png\n\n")

# ========== 13. Generate MA Plot ==========
cat("Generating MA plot...\n")

p2 <- ggplot(plot_data, aes(x = log10(baseMean + 1), y = log2FoldChange, color = plot_color)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Not significant" = "gray70",
                                "Low (p < 0.05)" = "skyblue",
                                "Medium (p < 0.01)" = "orange",
                                "High (p < 0.001)" = "red3")) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  labs(title = "MA Plot: AD vs CTR+ (Matching Paper)",
       subtitle = paste0("AD (n=", sum(coldata_filtered$condition == "AD"), ") vs CTR+ (n=", 
                        sum(coldata_filtered$condition == "CTR+"), ")"),
       x = "Log10(Mean Expression + 1)",
       y = "Log2 Fold Change (AD vs CTR+)",
       color = "Significance") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10),
        legend.position = "right")

ggsave("ma_plot_AD_vs_CTRplus.png", p2, width = 10, height = 7, dpi = 300)
cat("  Saved: ma_plot_AD_vs_CTRplus.png\n\n")

# ========== 14. Compare with Paper's Findings ==========
cat("===== COMPARISON WITH PAPER =====\n")
cat("Paper's reported genes (AD vs CTR+):\n")
cat("  - SIGLEC1 (decreased in AD, FDR=0.002)\n")
cat("  - CXCL10 (increased in AD, FDR=0.02)\n")
cat("  - CXCR1, CXCR2 (chemokine receptors)\n\n")

# Check if we found these genes
paper_genes <- c("SIGLEC1", "CXCL10", "CXCR1", "CXCR2")
cat("Our findings for paper's reported genes:\n")
for (gene in paper_genes) {
  gene_result <- results_dt[feature == gene]
  if (nrow(gene_result) > 0) {
    cat(sprintf("  %s: log2FC=%.2f, padj=%.4f, %s\n", 
                gene, 
                gene_result$log2FoldChange, 
                ifelse(is.na(gene_result$padj), 1, gene_result$padj),
                gene_result$regulation))
  } else {
    cat(sprintf("  %s: Not found in results\n", gene))
  }
}

cat("\n===== ANALYSIS COMPLETE =====\n")
cat("Output files saved to: /home/jacobc/hcaTE/pseudobulk/\n")
cat("Files with '_AD_vs_CTRplus' suffix contain this analysis\n\n")

# Save session info
sink("pseudobulk_diffexp_AD_vs_CTRplus_session_info.txt")
sessionInfo()
sink()
