#!/usr/bin/env Rscript
# Parse Table 1.xlsx from Alsema et al. 2020 paper
# Extract all sheets and search for heat shock proteins

library(readxl)
library(dplyr)
library(tidyr)

# ============================================================================
# CONFIGURATION
# ============================================================================

EXCEL_FILE <- "/home/jacobc/hcaTE/paper_data/Table 1.xlsx"
OUTPUT_DIR <- "/home/jacobc/hcaTE/paper_data"

# Genes of interest (heat shock proteins and stress response)
GENES_OF_INTEREST <- c(
  "HSPA6", "HSPA1A", "HSPA1B", "HSPA8", "HSPA9", 
  "HSP90AA1", "HSP90AB1", "HSPB1", "HSPD1", "HSPE1",
  "DNAJB1", "DNAJA1", "JUN", "FOS", "COA1"
)

# ============================================================================
# FUNCTIONS
# ============================================================================

cat("================================================================================\n")
cat("Parsing Table 1.xlsx from Alsema et al. 2020\n")
cat("================================================================================\n\n")

# Check if file exists
if (!file.exists(EXCEL_FILE)) {
  stop("ERROR: Excel file not found: ", EXCEL_FILE)
}

cat("1. Reading Excel file...\n")
cat("   File:", EXCEL_FILE, "\n")
cat("   Size:", file.info(EXCEL_FILE)$size / 1024^2, "MB\n\n")

# Get all sheet names
sheet_names <- excel_sheets(EXCEL_FILE)
cat("   Found", length(sheet_names), "sheets:\n")
for (i in seq_along(sheet_names)) {
  cat("      ", i, ". ", sheet_names[i], "\n", sep = "")
}

cat("\n2. Extracting all sheets...\n")

# Create output directory for individual sheets
sheets_dir <- file.path(OUTPUT_DIR, "sheets")
dir.create(sheets_dir, showWarnings = FALSE, recursive = TRUE)

# Store all sheets
all_sheets <- list()
sheet_summaries <- data.frame(
  sheet_name = character(),
  n_rows = integer(),
  n_cols = integer(),
  column_names = character(),
  stringsAsFactors = FALSE
)

for (sheet_name in sheet_names) {
  cat("   Reading:", sheet_name, "...")
  
  # Read sheet
  tryCatch({
    df <- read_excel(EXCEL_FILE, sheet = sheet_name)
    all_sheets[[sheet_name]] <- df
    
    # Save as CSV
    output_file <- file.path(sheets_dir, paste0(gsub("[^A-Za-z0-9_]", "_", sheet_name), ".csv"))
    write.csv(df, output_file, row.names = FALSE)
    
    cat(" ", nrow(df), "rows,", ncol(df), "cols\n")
    
    # Store summary
    sheet_summaries <- rbind(sheet_summaries, data.frame(
      sheet_name = sheet_name,
      n_rows = nrow(df),
      n_cols = ncol(df),
      column_names = paste(colnames(df), collapse = " | "),
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat(" ERROR:", e$message, "\n")
  })
}

# Save sheet summary
write.csv(sheet_summaries, file.path(OUTPUT_DIR, "sheet_summaries.csv"), row.names = FALSE)
cat("\n   ✓ All sheets saved to:", sheets_dir, "\n")
cat("   ✓ Summary saved to:", file.path(OUTPUT_DIR, "sheet_summaries.csv"), "\n")

# ============================================================================
# SEARCH FOR GENES OF INTEREST
# ============================================================================

cat("\n3. Searching for genes of interest across all sheets...\n")
cat("   Looking for:", paste(GENES_OF_INTEREST, collapse = ", "), "\n\n")

results <- list()

for (sheet_name in names(all_sheets)) {
  df <- all_sheets[[sheet_name]]
  
  # Try to find gene column (common names: Gene, gene, Gene.name, symbol, etc.)
  gene_col <- NULL
  possible_gene_cols <- c("Gene", "gene", "Gene.name", "gene_name", "Gene_name", 
                          "symbol", "Symbol", "SYMBOL", "geneName", "GeneSymbol",
                          "Gene symbol", "gene symbol")
  
  for (col in possible_gene_cols) {
    if (col %in% colnames(df)) {
      gene_col <- col
      break
    }
  }
  
  # If no obvious gene column, check first few columns for gene-like content
  if (is.null(gene_col)) {
    for (i in 1:min(3, ncol(df))) {
      col_name <- colnames(df)[i]
      # Check if first few values look like gene names
      sample_vals <- head(df[[col_name]], 10)
      if (is.character(sample_vals) || is.factor(sample_vals)) {
        # Check if they match typical gene name patterns (all caps, 2-10 chars)
        matches <- sum(grepl("^[A-Z0-9]{2,10}$", as.character(sample_vals)))
        if (matches >= 5) {  # At least half should match
          gene_col <- col_name
          break
        }
      }
    }
  }
  
  if (is.null(gene_col)) {
    cat("   [", sheet_name, "] - No gene column found\n", sep = "")
    next
  }
  
  cat("   [", sheet_name, "] - Using column '", gene_col, "':\n", sep = "")
  
  # Search for genes of interest
  found_genes <- df %>%
    filter(!!sym(gene_col) %in% GENES_OF_INTEREST)
  
  if (nrow(found_genes) > 0) {
    cat("      ✓ FOUND ", nrow(found_genes), " gene(s)!\n", sep = "")
    
    # Print details
    for (i in 1:nrow(found_genes)) {
      gene <- found_genes[[gene_col]][i]
      cat("         ", gene, "\n", sep = "")
      
      # Print key columns if they exist
      for (col in c("logFC", "log2FoldChange", "adj.P.Val", "padj", "pvalue", "P.Value", "FDR")) {
        if (col %in% colnames(found_genes)) {
          val <- found_genes[[col]][i]
          cat("            ", col, ": ", val, "\n", sep = "")
        }
      }
    }
    
    # Store results
    results[[sheet_name]] <- found_genes
    
  } else {
    cat("      - None found\n")
  }
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n4. Saving results...\n")

if (length(results) > 0) {
  # Combine all results
  all_results <- bind_rows(
    lapply(names(results), function(sheet) {
      results[[sheet]] %>%
        mutate(sheet_name = sheet, .before = 1)
    })
  )
  
  output_file <- file.path(OUTPUT_DIR, "heat_shock_proteins_in_paper.csv")
  write.csv(all_results, output_file, row.names = FALSE)
  cat("   ✓ Heat shock protein results saved to:", output_file, "\n")
  
  # Print summary
  cat("\n================================================================================\n")
  cat("SUMMARY: Heat Shock Proteins Found in Paper\n")
  cat("================================================================================\n\n")
  
  print(all_results)
  
} else {
  cat("   ⚠ WARNING: No genes of interest found in any sheet!\n")
}

cat("\n================================================================================\n")
cat("Analysis complete!\n")
cat("================================================================================\n\n")

cat("Output files:\n")
cat("  - Individual sheets: ", sheets_dir, "/\n", sep = "")
cat("  - Sheet summaries: ", file.path(OUTPUT_DIR, "sheet_summaries.csv"), "\n", sep = "")
if (length(results) > 0) {
  cat("  - HSP results: ", file.path(OUTPUT_DIR, "heat_shock_proteins_in_paper.csv"), "\n", sep = "")
}
cat("\n")
