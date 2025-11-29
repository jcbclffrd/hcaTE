#!/usr/bin/env Rscript
# Differential Expression Analysis: Pseudo-bulk scRNA-seq
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

cat("===== Differential Expression Analysis: Pseudo-bulk scRNA-seq =====\n\n")

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
cat("  Total counts:", format(sum(count_data), big.mark=","), "\n")
cat("  Non-zero values:", format(sum(count_data > 0), big.mark=","), "\n\n")

# ========== 3. Create Sample Metadata ==========
cat("Creating sample metadata...\n")

coldata <- data.frame(
  sample_id = sample_ids,
  condition = conditions,
  row.names = sample_ids
)

cat("  Sample distribution by condition:\n")
print(table(coldata$condition))
cat("\n")

# ========== 4. Define Comparison Groups ==========
cat("Setting up comparison: AD vs Controls\n")

# Combine CTR and CTR+ as control group for maximum power
coldata$disease_group <- ifelse(coldata$condition == "AD", "AD", 
                               ifelse(coldata$condition %in% c("CTR", "CTR+"), "Control", "MCI"))

# Filter out MCI samples for cleaner AD vs Control comparison
coldata_filtered <- coldata[coldata$disease_group != "MCI", ]
count_data_filtered <- count_data[, rownames(coldata_filtered)]

cat("  Final sample sizes:\n")
cat("    AD:", sum(coldata_filtered$disease_group == "AD"), "samples\n")
cat("    Control:", sum(coldata_filtered$disease_group == "Control"), "samples\n")
cat("    (Excluded MCI:", sum(coldata$disease_group == "MCI"), "samples)\n\n")

# Set reference level
coldata_filtered$disease_group <- factor(coldata_filtered$disease_group, 
                                        levels = c("Control", "AD"))

# ========== 5. Create DESeq2 Dataset ==========
cat("Creating DESeq2 dataset...\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_data_filtered,
  colData = coldata_filtered,
  design = ~ disease_group
)

cat("  DESeq2 dataset created\n\n")

# ========== 6. Pre-filtering: Remove Low-Count Features ==========
cat("Pre-filtering low-count features...\n")

# Keep features with at least 10 counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds_filtered <- dds[keep, ]

cat("  Original features:", format(nrow(dds), big.mark=","), "\n")
cat("  After filtering:", format(nrow(dds_filtered), big.mark=","), "\n")
cat("  Removed:", format(nrow(dds) - nrow(dds_filtered), big.mark=","), "low-count features\n\n")

# ========== 7. Run DESeq2 Differential Expression Analysis ==========
cat("Running DESeq2 differential expression analysis...\n")
cat("  This may take several minutes for large datasets...\n\n")

dds_de <- DESeq(dds_filtered)

cat("  DESeq2 analysis complete!\n\n")

# ========== 8. Extract Results ==========
cat("Extracting differential expression results...\n")

# Get results: AD vs Control (log2FoldChange = log2(AD/Control))
results <- results(dds_de, contrast = c("disease_group", "AD", "Control"))

# Convert to data.table
results_dt <- as.data.table(results, keep.rownames = "feature")

cat("  Total features tested:", format(nrow(results_dt), big.mark=","), "\n")
cat("  Features with p-value < 0.05:", sum(results_dt$pvalue < 0.05, na.rm = TRUE), "\n")
cat("  Features with padj < 0.05:", sum(results_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("  Features with padj < 0.01:", sum(results_dt$padj < 0.01, na.rm = TRUE), "\n\n")

# ========== 9. Annotate Results ==========
cat("Annotating results...\n")

results_dt[, significance := "Not significant"]
results_dt[padj < 0.05, significance := "padj < 0.05"]
results_dt[padj < 0.01, significance := "padj < 0.01"]
results_dt[padj < 0.001, significance := "padj < 0.001"]

results_dt[, regulation := ifelse(log2FoldChange > 0, "Up in AD", "Down in AD")]
results_dt[, abs_log2FC := abs(log2FoldChange)]

# Annotate gene vs TE
results_dt[, feature_type := ifelse(grepl("^ENS", feature), "Gene", "TE")]

cat("  Feature type distribution:\n")
print(table(results_dt$feature_type))
cat("\n")

# ========== 10. Sort Results ==========
setorder(results_dt, padj, -abs_log2FC)

# ========== 11. Summary Statistics ==========
cat("===== DIFFERENTIAL EXPRESSION SUMMARY =====\n\n")

cat("Total features analyzed:", format(nrow(results_dt), big.mark=","), "\n")
cat("  - Genes:", sum(results_dt$feature_type == "Gene", na.rm = TRUE), "\n")
cat("  - TEs:", sum(results_dt$feature_type == "TE", na.rm = TRUE), "\n\n")

cat("Significant features (padj < 0.05):", sum(results_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("  - Upregulated in AD:", 
    sum(results_dt$padj < 0.05 & results_dt$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  - Downregulated in AD:", 
    sum(results_dt$padj < 0.05 & results_dt$log2FoldChange < 0, na.rm = TRUE), "\n\n")

cat("Significant TEs (padj < 0.05):", 
    sum(results_dt$padj < 0.05 & results_dt$feature_type == "TE", na.rm = TRUE), "\n")
cat("Significant genes (padj < 0.05):", 
    sum(results_dt$padj < 0.05 & results_dt$feature_type == "Gene", na.rm = TRUE), "\n\n")

# Top upregulated in AD
cat("Top 10 upregulated in AD (padj < 0.05):\n")
top_up <- results_dt[padj < 0.05 & log2FoldChange > 0][1:min(10, .N)]
print(top_up[, .(feature, feature_type, log2FoldChange, padj)])
cat("\n")

# Top downregulated in AD
cat("Top 10 downregulated in AD (padj < 0.05):\n")
top_down <- results_dt[padj < 0.05 & log2FoldChange < 0][1:min(10, .N)]
print(top_down[, .(feature, feature_type, log2FoldChange, padj)])
cat("\n")

# Top TEs specifically
cat("Top 10 upregulated TEs in AD (padj < 0.05):\n")
top_te_up <- results_dt[padj < 0.05 & log2FoldChange > 0 & feature_type == "TE"][1:min(10, .N)]
print(top_te_up[, .(feature, log2FoldChange, padj)])
cat("\n")

# ========== 12. Save Results ==========
cat("Saving results...\n")

# Save all results
fwrite(results_dt, "pseudobulk_diffexp_results_AD_vs_Control.csv")
cat("  Full results: pseudobulk_diffexp_results_AD_vs_Control.csv\n")

# Save significant only
results_sig <- results_dt[padj < 0.05]
fwrite(results_sig, "pseudobulk_diffexp_results_AD_vs_Control_significant.csv")
cat("  Significant results: pseudobulk_diffexp_results_AD_vs_Control_significant.csv\n")

# Save TE-specific results
results_te <- results_dt[feature_type == "TE"]
fwrite(results_te, "pseudobulk_diffexp_results_TEs_only.csv")
cat("  TE-only results: pseudobulk_diffexp_results_TEs_only.csv\n")

# Save significant TEs only
results_te_sig <- results_dt[padj < 0.05 & feature_type == "TE"]
fwrite(results_te_sig, "pseudobulk_diffexp_results_TEs_significant.csv")
cat("  Significant TEs: pseudobulk_diffexp_results_TEs_significant.csv\n\n")

# ========== 13. Generate Volcano Plot ==========
cat("Generating volcano plot...\n")

volcano_data <- results_dt[!is.na(padj)]
volcano_data[, neg_log10_padj := -log10(padj)]

# Color by significance and feature type
volcano_data[, plot_color := "Not significant"]
volcano_data[padj < 0.05 & feature_type == "TE", plot_color := "Significant TE"]
volcano_data[padj < 0.05 & feature_type == "Gene", plot_color := "Significant Gene"]

p_volcano <- ggplot(volcano_data, aes(x = log2FoldChange, y = neg_log10_padj, color = plot_color)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = c("Not significant" = "gray70",
                                "Significant Gene" = "blue",
                                "Significant TE" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(title = "Volcano Plot: AD vs Control (Pseudo-bulk)",
       subtitle = paste0(sum(results_dt$padj < 0.05, na.rm = TRUE), 
                        " significant features (", 
                        sum(results_te_sig$padj < 0.05, na.rm = TRUE),
                        " TEs, ",
                        sum(results_dt$padj < 0.05 & results_dt$feature_type == "Gene", na.rm = TRUE),
                        " genes)"),
       x = "log2 Fold Change (AD / Control)",
       y = "-log10(adjusted p-value)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("volcano_plot_pseudobulk_AD_vs_Control.png", p_volcano, width = 12, height = 9, dpi = 300)
cat("  Volcano plot: volcano_plot_pseudobulk_AD_vs_Control.png\n\n")

# ========== 14. Generate MA Plot ==========
cat("Generating MA plot...\n")

p_ma <- ggplot(results_dt[!is.na(padj)], aes(x = log10(baseMean + 1), 
                                              y = log2FoldChange, 
                                              color = padj < 0.05)) +
  geom_point(alpha = 0.4, size = 1.2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"),
                     labels = c("TRUE" = "padj < 0.05", "FALSE" = "Not significant")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(title = "MA Plot: AD vs Control (Pseudo-bulk)",
       subtitle = "Mean expression vs log2 fold change",
       x = "log10(Mean Expression + 1)",
       y = "log2 Fold Change (AD / Control)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("ma_plot_pseudobulk_AD_vs_Control.png", p_ma, width = 12, height = 9, dpi = 300)
cat("  MA plot: ma_plot_pseudobulk_AD_vs_Control.png\n\n")

# ========== 15. TE-specific Volcano Plot ==========
cat("Generating TE-specific volcano plot...\n")

te_volcano_data <- results_dt[feature_type == "TE" & !is.na(padj)]
te_volcano_data[, neg_log10_padj := -log10(padj)]

p_te_volcano <- ggplot(te_volcano_data, aes(x = log2FoldChange, y = neg_log10_padj, 
                                             color = padj < 0.05)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "gray60"),
                     labels = c("TRUE" = "padj < 0.05", "FALSE" = "Not significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray40") +
  labs(title = "Volcano Plot: Transposable Elements (AD vs Control)",
       subtitle = paste0(sum(te_volcano_data$padj < 0.05, na.rm = TRUE), 
                        " significant TEs out of ",
                        nrow(te_volcano_data), " analyzed"),
       x = "log2 Fold Change (AD / Control)",
       y = "-log10(adjusted p-value)",
       color = "Significance") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave("volcano_plot_TEs_only.png", p_te_volcano, width = 12, height = 9, dpi = 300)
cat("  TE volcano plot: volcano_plot_TEs_only.png\n\n")

# ========== 16. Save Normalized Counts ==========
cat("Saving normalized counts for downstream analysis...\n")

normalized_counts <- counts(dds_de, normalized = TRUE)
normalized_counts_dt <- as.data.table(normalized_counts, keep.rownames = "feature")

fwrite(normalized_counts_dt, "pseudobulk_normalized_counts.csv")
cat("  Normalized counts: pseudobulk_normalized_counts.csv\n\n")

# ========== 17. Save Session Info ==========
cat("Saving session info...\n")
sink("pseudobulk_diffexp_session_info.txt")
sessionInfo()
sink()
cat("  Session info: pseudobulk_diffexp_session_info.txt\n\n")

cat("===== ANALYSIS COMPLETE =====\n\n")

cat("Output files:\n")
cat("  Results:\n")
cat("    - pseudobulk_diffexp_results_AD_vs_Control.csv (all features)\n")
cat("    - pseudobulk_diffexp_results_AD_vs_Control_significant.csv (padj < 0.05)\n")
cat("    - pseudobulk_diffexp_results_TEs_only.csv (all TEs)\n")
cat("    - pseudobulk_diffexp_results_TEs_significant.csv (significant TEs)\n")
cat("  Plots:\n")
cat("    - volcano_plot_pseudobulk_AD_vs_Control.png\n")
cat("    - ma_plot_pseudobulk_AD_vs_Control.png\n")
cat("    - volcano_plot_TEs_only.png\n")
cat("  Data:\n")
cat("    - pseudobulk_normalized_counts.csv\n")
cat("    - pseudobulk_diffexp_session_info.txt\n\n")

cat("Summary:\n")
cat("  Comparison: AD (", sum(coldata_filtered$disease_group == "AD"), 
    " samples) vs Control (", sum(coldata_filtered$disease_group == "Control"), 
    " samples)\n", sep="")
cat("  Significant features:", sum(results_dt$padj < 0.05, na.rm = TRUE), "\n")
cat("  Significant TEs:", sum(results_te_sig$padj < 0.05, na.rm = TRUE), "\n\n")

cat("Next steps:\n")
cat("  1. Review significant TEs in pseudobulk_diffexp_results_TEs_significant.csv\n")
cat("  2. Examine volcano plots for overall patterns\n")
cat("  3. Consider TE family-level analysis\n")
cat("  4. Validate top candidates\n\n")
