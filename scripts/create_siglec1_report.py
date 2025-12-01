#!/usr/bin/env python3
"""
Create a text report summarizing SIGLEC1 expression
"""

import pandas as pd
from pathlib import Path

output_dir = Path('/home/jacobc/hcaTE/regional_analysis')
siglec1_df = pd.read_csv(output_dir / 'siglec1_expression_by_sample.csv')

report_path = output_dir / 'SIGLEC1_REPORT.txt'

with open(report_path, 'w') as f:
    f.write('='*80 + '\n')
    f.write('SIGLEC1 EXPRESSION ANALYSIS REPORT\n')
    f.write('Pseudobulk RNA-seq from bc-SmartSeq2 Data\n')
    f.write('='*80 + '\n\n')
    
    f.write('OVERVIEW\n')
    f.write('-'*80 + '\n')
    f.write(f'Total samples analyzed: {len(siglec1_df)}\n')
    f.write(f'Total SIGLEC1 counts: {siglec1_df["SIGLEC1_count"].sum()}\n')
    f.write(f'Samples with SIGLEC1 > 0: {(siglec1_df["SIGLEC1_count"] > 0).sum()}\n')
    f.write(f'Percentage with expression: {(siglec1_df["SIGLEC1_count"] > 0).sum() / len(siglec1_df) * 100:.1f}%\n\n')
    
    f.write('CONDITION BREAKDOWN\n')
    f.write('-'*80 + '\n')
    condition_summary = siglec1_df.groupby('Condition').agg({
        'SIGLEC1_count': ['sum', 'mean', 'count']
    })
    condition_summary.columns = ['Total_Counts', 'Mean', 'N_Samples']
    f.write(condition_summary.to_string() + '\n\n')
    
    f.write('SAMPLES WITH SIGLEC1 EXPRESSION\n')
    f.write('-'*80 + '\n')
    nonzero = siglec1_df[siglec1_df['SIGLEC1_count'] > 0]
    if len(nonzero) > 0:
        for idx, row in nonzero.iterrows():
            f.write(f'  {row["Sample"]}: {int(row["SIGLEC1_count"])} counts ({row["Condition"]})\n')
    else:
        f.write('  No samples with SIGLEC1 expression\n')
    
    f.write('\n' + '='*80 + '\n')
    f.write('INTERPRETATION\n')
    f.write('='*80 + '\n\n')
    
    f.write('The extremely low SIGLEC1 counts observed in this pseudobulk analysis are\n')
    f.write('expected due to:\n\n')
    f.write('1. CELL-TYPE DILUTION: Each sample aggregates ~80-100 cells from multiple cell\n')
    f.write('   types (microglia, neurons, astrocytes, oligodendrocytes, etc.)\n\n')
    f.write('2. MICROGLIA SPECIFICITY: SIGLEC1 is primarily expressed in microglia, which\n')
    f.write('   typically represent 5-15% of total brain cells\n\n')
    f.write('3. PSEUDOBULK AVERAGING: Even if SIGLEC1 is highly expressed in microglia,\n')
    f.write('   the signal gets averaged across all cell types in the sample\n\n')
    f.write('RECOMMENDATION:\n')
    f.write('-'*80 + '\n')
    f.write('To properly analyze SIGLEC1 expression as shown in the paper\'s Figure 2:\n\n')
    f.write('1. Perform SINGLE-CELL level analysis using scTE output files\n')
    f.write('2. SUBSET to microglia cells using marker genes (CX3CR1, P2RY12, TMEM119)\n')
    f.write('3. Analyze SIGLEC1 expression within MICROGLIA by BRAIN REGION (LPS vs GFS)\n')
    f.write('4. This will eliminate the cell-type dilution effect and match the paper\'s\n')
    f.write('   analytical approach\n\n')
    
    f.write('='*80 + '\n')

print(f'Text report created successfully at: {report_path}')
print()
print('Report contents:')
print()
with open(report_path, 'r') as f:
    print(f.read())
